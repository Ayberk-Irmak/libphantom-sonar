// SPDX-License-Identifier: Apache-2.0
#include "phantom/echo_synth.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

// Rebuilds the transmitted waveform from what the analyser reported. The PDW
// carries centre frequency and bandwidth rather than the endpoints, because
// that is what a classifier estimates; the endpoints come back out by the same
// convention PulseSpec uses.
PulseSpec spec_from_descriptor(const PulseDescriptor& pdw) noexcept {
    PulseSpec s;
    s.type = pdw.type;
    s.duration_s = pdw.duration_s;
    s.amplitude = kOne;
    const Real half = pdw.bandwidth_hz / kTwo;
    switch (pdw.type) {
        case PulseType::Cw:
            s.f_start_hz = pdw.centre_freq_hz;
            s.f_end_hz = pdw.centre_freq_hz;
            break;
        case PulseType::LfmUp:
        case PulseType::Hfm:
            s.f_start_hz = pdw.centre_freq_hz - half;
            s.f_end_hz = pdw.centre_freq_hz + half;
            break;
        case PulseType::LfmDown:
            s.f_start_hz = pdw.centre_freq_hz + half;
            s.f_end_hz = pdw.centre_freq_hz - half;
            break;
        case PulseType::Unknown:
            break;
    }
    return s;
}

inline Real sinc_pi(Real x) noexcept {
    // sin(pi x) / (pi x), continuous at 0.
    if (std::fabs(x) < static_cast<Real>(1e-9)) return kOne;
    const Real a = kPi * x;
    return std::sin(a) / a;
}

}  // namespace

Real echo_doppler_scale(Real radial_velocity_mps, Real sound_speed_mps) noexcept {
    if (!(sound_speed_mps > kZero)) return kOne;
    const Real denom = sound_speed_mps - radial_velocity_mps;
    if (!(denom > kZero)) return kOne;   // at or beyond the sound speed
    return (sound_speed_mps + radial_velocity_mps) / denom;
}

std::size_t synthesize_echo(const PulseDescriptor& pdw,
                            const EchoSpec& echo,
                            Real sample_rate_hz,
                            Real sound_speed_mps,
                            std::span<Real> out) noexcept {
    if (!(sample_rate_hz > kZero) || !(sound_speed_mps > kZero)) return 0;
    if (!(echo.length_scale > kZero)) return 0;

    PulseSpec spec = spec_from_descriptor(pdw);
    spec.duration_s = pdw.duration_s * echo.length_scale;
    if (!spec.valid()) return 0;

    // Two-way Doppler: the echo is scaled once on the way out and once on the
    // way back, so alpha = (c+v)/(c-v), not (1 + v/c).
    const Real alpha = echo_doppler_scale(echo.radial_velocity_mps, sound_speed_mps);
    const Real doppler = alpha - kOne;

    spec.amplitude = pdw.amplitude
                   * std::pow(static_cast<Real>(10), echo.target_strength_db / static_cast<Real>(20));
    if (!(spec.amplitude > kZero)) return 0;

    return render_real_doppler(spec, sample_rate_hz, doppler, out);
}

std::size_t echoes_from_eigenrays(std::span<const Eigenray> paths,
                                 Real range_m,
                                 Real frequency_hz,
                                 Real source_speed_mps,
                                 Real target_strength_db,
                                 std::span<EchoSpec> out) noexcept {
    if (paths.empty() || out.empty()) return 0;

    // Two-way loss along each path, and the earliest arrival to reference to.
    Real earliest = paths[0].travel_time_s;
    Real best_loss = kZero;
    bool have_best = false;
    for (const Eigenray& p : paths) {
        if (p.travel_time_s < earliest) earliest = p.travel_time_s;
        if (p.near_caustic) continue;   // level not trustworthy; see Eigenray
        const Real loss = kTwo * transmission_loss_db(p, range_m, frequency_hz,
                                                      source_speed_mps);
        if (!have_best || loss < best_loss) { best_loss = loss; have_best = true; }
    }
    if (!have_best) return 0;

    std::size_t written = 0;
    for (const Eigenray& p : paths) {
        if (written >= out.size()) break;
        if (p.near_caustic) continue;
        const Real loss = kTwo * transmission_loss_db(p, range_m, frequency_hz,
                                                      source_speed_mps);
        EchoSpec& e = out[written++];
        e.range_offset_m = kZero;
        // Two-way travel time along this path, relative to the first arrival.
        e.extra_delay_s = kTwo * (p.travel_time_s - earliest);
        e.target_strength_db = target_strength_db - (loss - best_loss);
        e.radial_velocity_mps = kZero;
        e.length_scale = kOne;
    }
    return written;
}

std::size_t swarm_length(const PulseDescriptor& pdw,
                         std::span<const EchoSpec> echoes,
                         Real sample_rate_hz,
                         Real sound_speed_mps) noexcept {
    if (!(sample_rate_hz > kZero) || !(sound_speed_mps > kZero) || echoes.empty()) return 0;

    std::size_t furthest = 0;
    for (const EchoSpec& e : echoes) {
        if (!(e.length_scale > kZero)) continue;
        const Real delay = echo_delay_s(e.range_offset_m, sound_speed_mps) + e.extra_delay_s;
        if (delay < kZero) continue;   // a ghost cannot arrive before the ping
        const Real alpha = echo_doppler_scale(e.radial_velocity_mps, sound_speed_mps);
        const Real dur = pdw.duration_s * e.length_scale / alpha;
        const auto end = static_cast<std::size_t>((delay + dur) * sample_rate_hz) + 2;
        if (end > furthest) furthest = end;
    }
    return furthest;
}

std::size_t synthesize_swarm(const PulseDescriptor& pdw,
                             std::span<const EchoSpec> echoes,
                             Real sample_rate_hz,
                             Real sound_speed_mps,
                             std::span<Real> out) noexcept {
    const std::size_t need = swarm_length(pdw, echoes, sample_rate_hz, sound_speed_mps);
    if (need == 0 || out.size() < need) return 0;

    for (std::size_t i = 0; i < need; ++i) out[i] = kZero;

    for (const EchoSpec& e : echoes) {
        const Real delay = echo_delay_s(e.range_offset_m, sound_speed_mps) + e.extra_delay_s;
        if (delay < kZero) continue;
        const auto offset = static_cast<std::size_t>(delay * sample_rate_hz);
        if (offset >= need) continue;

        // Each ghost is rendered into the tail of the caller's buffer and
        // summed in place, so no scratch and no allocation are needed.
        const std::size_t room = need - offset;
        const std::size_t n = synthesize_echo(pdw, e, sample_rate_hz, sound_speed_mps,
                                              out.subspan(offset, room));
        (void)n;
    }
    return need;
}

Real null_gain_db(Real bandwidth_hz, Real centre_freq_hz, Real timing_error_s) noexcept {
    if (!(centre_freq_hz > kZero)) return kZero;
    const Real b = (bandwidth_hz > kZero) ? bandwidth_hz : kZero;
    const Real tau = timing_error_s;

    // Residual r(t) = e(t) - e(t - tau); over a flat band of width B centred on
    // f_c the residual-to-original power ratio integrates to
    //     G^2 = 2 - 2 cos(2 pi f_c tau) sinc(B tau)
    const Real ratio = kTwo - kTwo * std::cos(kTwo * kPi * centre_freq_hz * tau)
                                  * sinc_pi(b * tau);
    if (!(ratio > kZero)) return static_cast<Real>(-200);   // perfect null
    const Real db = static_cast<Real>(10) * std::log10(ratio);
    return (db < static_cast<Real>(-200)) ? static_cast<Real>(-200) : db;
}

Real max_timing_error_s(Real target_db, Real centre_freq_hz) noexcept {
    if (!(centre_freq_hz > kZero) || !(target_db < kZero)) return kZero;
    // Narrowband limit: G^2 = 4 sin^2(pi f_c tau), so
    //     tau = asin(10^(G/20) / 2) / (pi f_c)
    const Real amp = std::pow(static_cast<Real>(10), target_db / static_cast<Real>(20));
    const Real x = amp / kTwo;
    if (!(x < kOne)) return kZero;
    return std::asin(x) / (kPi * centre_freq_hz);
}

}  // namespace phantom
