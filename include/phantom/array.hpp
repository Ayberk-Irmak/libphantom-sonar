// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — line arrays, beamforming and bearing.
//
// Everything before this release was single-channel: a hydrophone knows WHEN a
// ping arrived and what shape it was, but not where it came from. That is the
// largest structural gap left, and it is also where v0.6's result gets its
// teeth -- reverberation scales with the ensonified area, and a narrower beam
// is the only way to shrink the azimuthal half of it.
//
// Angles are measured from BROADSIDE: 0 is perpendicular to the array, +/-pi/2
// is endfire. That convention makes sin(theta) the natural variable, and every
// formula below is written in it.
//
// Narrowband throughout. A wideband array steers by delay rather than by phase,
// and the difference matters once the signal bandwidth approaches the inverse
// of the array's traversal time -- see `narrowband_bandwidth_limit_hz`.
#ifndef PHANTOM_ARRAY_HPP
#define PHANTOM_ARRAY_HPP

#include "phantom/fft.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom {

// A uniform line array.
struct LineArray {
    std::size_t element_count = 0;
    Real spacing_m = 0;

    [[nodiscard]] constexpr bool valid() const noexcept {
        return element_count >= 2 && spacing_m > 0;
    }
    // Distance between the first and last element.
    [[nodiscard]] constexpr Real aperture_m() const noexcept {
        return valid() ? static_cast<Real>(element_count - 1) * spacing_m : 0;
    }
};

[[nodiscard]] constexpr Real wavelength_m(Real frequency_hz, Real sound_speed_mps) noexcept {
    return (frequency_hz > 0) ? sound_speed_mps / frequency_hz : 0;
}

// Array factor magnitude, normalised to 1 at the steer angle:
//
//   psi = k d (sin(look) - sin(steer)),   k = 2 pi / lambda
//   B   = | sin(N psi / 2) / (N sin(psi / 2)) |
//
// This is the whole spatial response of a uniformly shaded line array. Its
// closed forms -- unity at the steer angle, nulls at psi = 2 pi m / N, and a
// first sidelobe approaching -13.26 dB -- are what the test suite checks
// against, because a beam pattern that is subtly wrong still looks like a beam.
[[nodiscard]] Real array_factor(const LineArray& array, Real lambda_m,
                                Real steer_rad, Real look_rad) noexcept;

// Offset from the steer angle to the first null, in radians.
//
// Exact rather than the small-angle approximation: the null is where
// sin(look) - sin(steer) = lambda / (N d), so its angular offset depends on
// where the beam is pointing. A beam steered towards endfire is far broader
// than the same beam at broadside, which the lambda/(N d cos theta) form
// captures and lambda/(N d) does not. Returns 0 when the null would require
// |sin| > 1, i.e. when the beam runs off the end of the visible region.
[[nodiscard]] Real first_null_offset_rad(const LineArray& array, Real lambda_m,
                                         Real steer_rad) noexcept;

// -3 dB beamwidth in radians. Uses the standard 0.886 lambda / (N d cos theta),
// which is the half-power width of a uniformly shaded line array.
[[nodiscard]] Real beamwidth_3db_rad(const LineArray& array, Real lambda_m,
                                     Real steer_rad) noexcept;

// Largest element spacing with no grating lobe when steering out to
// `max_steer_rad`:  d <= lambda / (1 + |sin(max_steer)|).
//
// Steering to endfire therefore demands d <= lambda/2 -- the familiar rule, and
// a hard constraint on how much aperture a given element count can buy.
[[nodiscard]] Real max_spacing_no_grating_lobes_m(Real lambda_m, Real max_steer_rad) noexcept;

// True if the array has a grating lobe when steered to `steer_rad`.
[[nodiscard]] bool has_grating_lobe(const LineArray& array, Real lambda_m,
                                    Real steer_rad) noexcept;

// Array gain against spatially white noise: 10 log10(N).
//
// Signal adds coherently across elements (amplitude x N) and noise
// incoherently (amplitude x sqrt(N)), so power SNR improves by N. This assumes
// the noise is uncorrelated between elements, which isotropic ambient roughly
// is at d >= lambda/2 and directional interference emphatically is not.
[[nodiscard]] Real array_gain_db(const LineArray& array) noexcept;

// Highest signal bandwidth for which phase steering is a good approximation to
// delay steering, at a given steer angle: the array's traversal time is
// (N-1) d sin(theta) / c, and phase steering breaks down once the signal
// decorrelates across it. Taken as 1 / (10 * traversal), i.e. a tenth of the
// inverse traversal time, which keeps the edge-element phase error small.
//
// A 32-element half-wavelength array at 10 kHz steered to 45 degrees traverses
// in 1.1 ms, so phase steering holds to about 900 Hz of bandwidth -- far less
// than the 12 kHz chirps the analyser uses. Beamform per frequency bin, or
// steer with delays.
[[nodiscard]] Real narrowband_bandwidth_limit_hz(const LineArray& array,
                                                 Real steer_rad,
                                                 Real sound_speed_mps) noexcept;

// Cramer-Rao lower bound on bearing, radians.
//
//   var(theta) >= 6 / (rho (k d cos theta)^2 N (N^2 - 1))
//
// with `rho` the per-element SNR in power. Derived exactly as the arrival-time
// bound in docs/math_spec.md 7.2 was, and it is the spatial twin of it: the
// element index takes the place of time, and Sum n'^2 = N(N^2-1)/12 about the
// array centre takes the place of the waveform's mean-square bandwidth.
//
// Two consequences worth reading off it. Bearing accuracy improves as
// N^(-3/2), not N^(-1/2), because more elements buy both more signal AND more
// aperture. And it degrades as 1/cos(theta): a target at endfire is far harder
// to locate than one at broadside, because the array's projected aperture
// shrinks.
[[nodiscard]] Real bearing_crlb_rad(const LineArray& array, Real lambda_m,
                                    Real bearing_rad, Real element_snr) noexcept;

// Fills `out` with the complex element responses to a unit plane wave from
// `bearing_rad`. Returns the element count, or 0 on bad input.
std::size_t synthesize_plane_wave(const LineArray& array, Real lambda_m,
                                  Real bearing_rad, Real amplitude,
                                  std::span<Complex> out) noexcept;

// Conventional (delay-and-sum, phase-shaded) beamforming. Fills `out_power`
// with |beam|^2 at `out_power.size()` steering angles evenly spanning
// [angle_min_rad, angle_max_rad]. Returns the count written.
std::size_t beamform_power(const LineArray& array, Real lambda_m,
                           std::span<const Complex> elements,
                           Real angle_min_rad, Real angle_max_rad,
                           std::span<Real> out_power) noexcept;

// Bearing from a beam scan: the peak cell, refined by a parabola through its
// neighbours. `angle_min_rad`/`angle_max_rad` must match the scan that produced
// `power`.
[[nodiscard]] Real estimate_bearing_rad(std::span<const Real> power,
                                        Real angle_min_rad,
                                        Real angle_max_rad) noexcept;

// Split-beam bearing: beamform the two halves of the array separately and read
// the bearing off the phase between them.
//
//   dphi = k D (sin(theta) - sin(steer)),   D = (N/2) d
//
// Unambiguous exactly out to the first null, since |dphi| = pi there -- so it
// refines a bearing already known to within the mainlobe, and cannot find one
// outside it. Much cheaper than a fine scan, and the standard technique for a
// compact array.
//
// Returns the estimated bearing in radians, or `steer_rad` if the input is
// unusable.
[[nodiscard]] Real split_beam_bearing_rad(const LineArray& array, Real lambda_m,
                                          std::span<const Complex> elements,
                                          Real steer_rad) noexcept;

}  // namespace phantom

#endif  // PHANTOM_ARRAY_HPP
