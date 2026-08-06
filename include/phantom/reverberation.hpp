// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — reverberation.
//
// v0.5 capped the surface scattering loss and said the energy "goes into a
// diffuse field the ray model does not carry". This is that field coming back.
//
// Reverberation is the sum of backscatter from every scatterer the pulse
// illuminates, and at short range it is what an active sonar actually competes
// against -- not the ambient noise the analyser assumes. Getting that wrong
// changes the answer to the only question that matters: can you detect the
// target, and does a bigger transmitter help?
//
// It does not, and this header shows why in one line. Writing the echo and the
// reverberation with the same source level and the same two-way loss:
//
//     EL = SL - 2 TL + TS
//     RL = SL - 2 TL + S_s + 10 log10(A)
//     EL - RL = TS - S_s - 10 log10(A)
//
// SL and TL cancel *exactly*. In a reverberation-limited geometry, doubling the
// transmit power raises the target and the background together and buys
// nothing. Only a shorter pulse or a narrower beam -- both of which shrink the
// ensonified area A -- improve detection. That is the whole design argument for
// pulse compression, and it is why the analyser's waveforms are chirps.
#ifndef PHANTOM_REVERBERATION_HPP
#define PHANTOM_REVERBERATION_HPP

#include "phantom/types.hpp"

#include <span>

namespace phantom {

// ---------------------------------------------------------------------------
// Scattering strengths, dB re 1 m^2 per m^2 of surface (or per m^3 of volume)
// ---------------------------------------------------------------------------

// Lambert's law for bottom backscatter:
//
//     S_s = mu + 10 log10(sin^2 theta) = mu + 20 log10(sin theta)
//
// `mu_db` is the Lambert coefficient, about -27 dB for many sediments (the
// figure usually attributed to Mackenzie). The sin^2 dependence is the whole
// content of Lambert's law: a facet scatters in proportion to its projected
// area both on transmit and on receive.
//
// Empirical and angle-only. It carries no frequency dependence, which real
// bottoms do have, and it fails near normal incidence where a smooth bottom
// gives a specular return far above Lambert.
[[nodiscard]] Real lambert_bottom_scattering_db(Real mu_db, Real grazing_angle_rad) noexcept;

// Chapman-Harris (1962) empirical fit for wind-driven surface backscatter:
//
//     beta = 158 (v f^(1/3))^(-0.58)          v in knots, f in Hz
//     S_s  = 3.3 beta log10(theta/30) - 42.4 log10(beta) + 2.6    theta in degrees
//
// The formula is written out here so a reader can check it against the paper.
// Its coefficients have NOT been verified against the original publication in
// this repository -- see docs/validation.md. What is verified is that it
// behaves: monotone in wind and in grazing angle, and in the -20 to -70 dB band
// that surface backscatter occupies.
[[nodiscard]] Real chapman_harris_surface_scattering_db(Real wind_speed_knots,
                                                        Real frequency_hz,
                                                        Real grazing_angle_rad) noexcept;

[[nodiscard]] constexpr Real mps_to_knots(Real mps) noexcept {
    return mps * static_cast<Real>(1.943844);
}

// ---------------------------------------------------------------------------
// Ensonified geometry
// ---------------------------------------------------------------------------

// The patch of boundary contributing at range r: an annulus segment of
// azimuthal width r*phi and radial depth c*tau/2.
//
// The c*tau/2 is the reason pulse compression matters here. It is the RANGE
// resolution, so a chirp's effective tau is 1/B rather than its length -- a
// 20 ms pulse with 12 kHz of bandwidth shrinks the ensonified area by a factor
// of 240, which is 24 dB straight off the reverberation.
[[nodiscard]] Real ensonified_area_m2(Real range_m, Real beamwidth_rad,
                                      Real pulse_length_s, Real sound_speed_mps) noexcept;

// The shell of water contributing at range r: r^2 * Omega * c*tau/2.
[[nodiscard]] Real ensonified_volume_m3(Real range_m, Real solid_angle_sr,
                                        Real pulse_length_s, Real sound_speed_mps) noexcept;

// ---------------------------------------------------------------------------
// The sonar equation terms
// ---------------------------------------------------------------------------

// RL = SL - 2 TL + S_s + 10 log10(A).
//
// With spherical spreading this falls as 30 log10(r) for boundary reverberation
// -- two-way loss is 40 log10(r) and the area grows as 10 log10(r). Volume
// reverberation falls as 20 log10(r), because the volume grows as r^2. Those
// two exponents are the signature of which mechanism dominates, and the test
// suite checks both.
[[nodiscard]] Real reverberation_level_db(Real source_level_db,
                                          Real transmission_loss_db,
                                          Real scattering_strength_db,
                                          Real scatterer_extent) noexcept;

// EL - RL = TS - S_s - 10 log10(A).
//
// Note what is absent: source level and transmission loss. They cancel. This is
// the quantity that decides detection in a reverberation-limited geometry, and
// no amount of transmit power changes it.
[[nodiscard]] Real echo_to_reverberation_ratio_db(Real target_strength_db,
                                                  Real scattering_strength_db,
                                                  Real scatterer_extent) noexcept;

// Range at which reverberation drops to the ambient noise level -- beyond it
// the geometry is noise-limited and a bigger transmitter helps again.
//
// Solves RL(r) = NL for boundary reverberation (30 log10 r) by bisection.
// Returns 0 if reverberation is already below the noise at `min_range_m`.
[[nodiscard]] Real reverberation_limited_range_m(Real source_level_db,
                                                 Real scattering_strength_db,
                                                 Real beamwidth_rad,
                                                 Real pulse_length_s,
                                                 Real sound_speed_mps,
                                                 Real noise_level_db,
                                                 Real min_range_m = 10,
                                                 Real max_range_m = 100000) noexcept;

// ---------------------------------------------------------------------------
// A reverberation time series, for feeding the detector
// ---------------------------------------------------------------------------

struct ReverbProfile {
    Real source_level_db = 200;
    Real scattering_strength_db = -30;   // constant with angle, for simplicity
    Real beamwidth_rad = static_cast<Real>(0.2);
    Real pulse_length_s = static_cast<Real>(0.02);
    Real sound_speed_mps = 1500;
    Real noise_level_db = 60;            // ambient floor the reverberation decays into
    Real reference_level_db = 200;       // the dB value mapped to unit amplitude
};

// Fills `out` with the reverberation AMPLITUDE envelope for each sample,
// starting at `start_time_s` after transmission. Range is c*t/2.
//
// This is the level, not a realisation: multiply by a zero-mean random series
// to get a signal. Reverberation is the sum of many independent scatterers, so
// by the central limit theorem the result is Gaussian and its envelope Rayleigh
// -- which is what CFAR's exponential-power assumption expects.
//
// SIMPLIFICATION: a true reverberation series is the scatterer field CONVOLVED
// with the transmitted pulse, so it is correlated over the pulse length. White
// noise scaled by this envelope is correlated over one sample instead. After
// matched filtering both have a correlation time of 1/B, so the detector sees
// the right statistics; before it, they differ. Do not use this to study
// anything that depends on the pre-correlation spectrum.
std::size_t reverberation_envelope(const ReverbProfile& profile,
                                   Real sample_rate_hz,
                                   Real start_time_s,
                                   std::span<Real> out) noexcept;

}  // namespace phantom

#endif  // PHANTOM_REVERBERATION_HPP
