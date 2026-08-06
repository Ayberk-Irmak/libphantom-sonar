// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — active sonar echo synthesis.
//
// Given a Pulse Descriptor Word from the analyser, generate the echo a target
// of a chosen apparent range, strength and radial velocity would have produced.
// This is the countermeasure half of the library, and it only works because the
// analyser ran first: everything here is parameterised by what was intercepted.
//
// Three things it does, all of them ordinary physics:
//
//   delay       an apparent range offset dd costs dt = 2 dd / c, two-way
//   Doppler     a target closing at v scales the echo by (1+v/c)/(1-v/c),
//               which is ~1 + 2v/c -- the two-way factor, not the one-way one
//   strength    target strength in dB, applied as an amplitude ratio
//
// and a "ghost swarm": several such echoes in one pulse train, to place more
// apparent targets than there are vehicles.
//
// WHAT IT DOES NOT DO. Anti-phase cancellation of a hull's acoustic
// cross-section is bounded by timing accuracy in a way that rules it out in
// practice -- 20 dB costs under two millimetres of path error. `null_gain_db()`
// at the bottom measures that rather than claiming around it.
#ifndef PHANTOM_ECHO_SYNTH_HPP
#define PHANTOM_ECHO_SYNTH_HPP

#include "phantom/eigenray.hpp"
#include "phantom/ping_analyzer.hpp"
#include "phantom/types.hpp"
#include "phantom/waveform.hpp"

#include <span>

namespace phantom {

// One synthetic return.
struct EchoSpec {
    // Apparent range offset from the true reflector, in metres. Positive puts
    // the ghost further away. Converted to delay as 2*dr/c.
    Real range_offset_m = 0;
    // Delay added directly, on top of the range offset. This is the multipath
    // route: a traced path's arrival time is a time, not a range, and forcing
    // it through 2*dr/c would assume a straight line at a nominal sound speed
    // -- exactly the assumption the ray tracer exists to remove.
    Real extra_delay_s = 0;
    // Echo level relative to the intercepted ping, in dB. Negative is quieter.
    Real target_strength_db = 0;
    // Radial velocity to impose, m/s, positive closing.
    Real radial_velocity_mps = 0;
    // Multiplies the echo duration. >1 mimics an extended target whose returns
    // from bow and stern smear the pulse; 1 leaves the length alone.
    Real length_scale = 1;
};

// Two-way delay for an apparent range offset: dt = 2 dr / c.
[[nodiscard]] constexpr Real echo_delay_s(Real range_offset_m, Real sound_speed_mps) noexcept {
    return (sound_speed_mps > 0) ? (2 * range_offset_m) / sound_speed_mps : 0;
}

// Two-way Doppler scale for a target closing at v:
//
//     alpha = (c + v) / (c - v)  ~  1 + 2 v/c
//
// The exact form is used rather than the linear approximation because v/c is
// large underwater -- at 30 m/s the two differ by 0.4%, which is 4 samples of
// delay over a 20 ms pulse and enough to matter for a matched filter.
[[nodiscard]] Real echo_doppler_scale(Real radial_velocity_mps, Real sound_speed_mps) noexcept;

// Renders one echo into `out`, starting at index 0. The caller places it in the
// output stream at `pdw.toa_s + echo_delay_s(...)`.
//
// Uses the intercepted PDW's waveform parameters, so the echo is coherent with
// the ping that provoked it -- which is exactly what makes a matched filter at
// the far end accept it.
std::size_t synthesize_echo(const PulseDescriptor& pdw,
                            const EchoSpec& echo,
                            Real sample_rate_hz,
                            Real sound_speed_mps,
                            std::span<Real> out) noexcept;

// Renders a ghost swarm: `echoes.size()` returns summed into `out`, each at its
// own delay, strength and velocity. `out` is zeroed first and must be long
// enough for the latest ghost; `required_length()` reports how long that is.
//
// Returns the number of samples written, or 0 if `out` is too short.
std::size_t synthesize_swarm(const PulseDescriptor& pdw,
                             std::span<const EchoSpec> echoes,
                             Real sample_rate_hz,
                             Real sound_speed_mps,
                             std::span<Real> out) noexcept;

// Builds one echo per eigenray, so the arrival structure comes from the traced
// paths rather than from hand-picked delays.
//
// Delays are referenced to the FIRST path, so the swarm starts at zero and the
// caller places it wherever the reply is transmitted. Levels are the two-way
// transmission loss (out along the path and back along the same one) relative
// to the strongest path, plus `target_strength_db`.
//
// This is a same-path round trip: a monostatic geometry where the echo returns
// the way the ping came. A bistatic case needs two path sets and is not what
// this does.
//
// Returns the number of echoes written.
std::size_t echoes_from_eigenrays(std::span<const Eigenray> paths,
                                  Real range_m,
                                  Real frequency_hz,
                                  Real source_speed_mps,
                                  Real target_strength_db,
                                  std::span<EchoSpec> out) noexcept;

// Samples `synthesize_swarm` needs for this PDW and echo set.
[[nodiscard]] std::size_t swarm_length(const PulseDescriptor& pdw,
                                       std::span<const EchoSpec> echoes,
                                       Real sample_rate_hz,
                                       Real sound_speed_mps) noexcept;

// ---------------------------------------------------------------------------
// Anti-phase cancellation, and what it actually costs
// ---------------------------------------------------------------------------
//
// The idea that turns up in every countermeasure specification: transmit the
// inverse of the incoming ping and cancel your own echo. It works. The question
// is what it costs, and the answer is best given in millimetres.
//
// Residual power for e(t) - e(t - tau) over a flat band of width B centred on
// f_c:
//
//     G^2 = 2 - 2 cos(2 pi f_c tau) sinc(B tau)
//
// For any realistic sonar (f_c >> B) the cosine dominates and the bandwidth
// term is negligible -- at 12 kHz centre with a 12 kHz sweep and a 2 us error,
// widening the band from 0 to 12 kHz moves the residual by 0.3 dB. So
// cancellation is set almost entirely by TIMING, not by how wideband the ping
// is. Solving the narrowband limit gives the number an engineer needs:
//
//     tau_max = asin(10^(G/20) / 2) / (pi f_c)
//
// At 12 kHz, 20 dB of cancellation needs tau < 1.3 us -- under two millimetres
// of path error at 1500 m/s. And past a quarter period the canceller ADDS
// energy: a mistimed null makes the vehicle louder, by up to 6 dB. That is the
// outcome a naive implementation ships.
//
// Two effects this does NOT model, both of which make the real figure worse:
//
//   - Hull scattering is distributed over many wavelengths, so the echo leaves
//     with a bearing-dependent phase that one point projector cannot match at
//     more than one bearing at a time.
//   - A canceller cannot emit the inverse of a sample it has not yet received.
//     tau therefore has a floor set by acquisition and processing latency,
//     which on any real receiver is orders of magnitude above 1.3 us.
//
// Treat the output as an upper bound on achievable cancellation, not a
// prediction of it.

// Residual level in dB relative to the uncancelled echo, for a residual timing
// error `timing_error_s`. Negative is cancellation; 0 dB is no effect;
// positive means the canceller has made the target louder.
[[nodiscard]] Real null_gain_db(Real bandwidth_hz,
                                Real centre_freq_hz,
                                Real timing_error_s) noexcept;

// Largest timing error that still achieves at least `target_db` (negative) of
// cancellation at `centre_freq_hz`, in the narrowband limit. Multiply by the
// sound speed to get the path-length budget. Returns 0 for a non-negative
// target or a non-positive frequency.
[[nodiscard]] Real max_timing_error_s(Real target_db, Real centre_freq_hz) noexcept;

}  // namespace phantom

#endif  // PHANTOM_ECHO_SYNTH_HPP
