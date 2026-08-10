// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — wideband beamforming, shading, and adaptive beams.
//
// v0.7 built an array and then showed why it could not be used: phase steering
// holds to 91 Hz of bandwidth at 45 degrees, and this library's own waveforms
// are 12 kHz chirps. This header fixes that, and joins the array to the pulse
// analyser so a Pulse Descriptor Word can finally carry a bearing.
//
// Three pieces:
//
//   WIDEBAND   steer per frequency bin. A phase shift that is wrong across a
//              band is exactly right within one bin, so transforming, steering
//              bin by bin and transforming back is exact wideband beamforming
//              rather than an approximation to it.
//
//   SHADING    a uniform array pays -13.26 dB sidelobes. A window buys tens of
//              dB of sidelobe suppression for a wider mainlobe and a little
//              gain -- the trade every array designer makes, quantified here.
//
//   ADAPTIVE   conventional beamforming cannot resolve two sources inside a
//              beamwidth however much SNR you give it: the resolution is set by
//              the aperture, full stop. MVDR can, by placing nulls rather than
//              scanning a fixed pattern.
#ifndef PHANTOM_BEAMFORMER_HPP
#define PHANTOM_BEAMFORMER_HPP

#include "phantom/array.hpp"
#include "phantom/fft.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom {

// ---------------------------------------------------------------------------
// Shading
// ---------------------------------------------------------------------------

enum class Shading : std::uint8_t { Uniform = 0, Hann, Hamming, Blackman };

[[nodiscard]] const char* shading_name(Shading s) noexcept;

// Weight for element `index` of `count`, in [0, 1].
[[nodiscard]] Real shading_weight(Shading shading, std::size_t index,
                                  std::size_t count) noexcept;

// Array factor with shading. The shaded pattern is the discrete-time Fourier
// transform of the window, so its sidelobe level IS the window's -- which is
// why the numbers below are the familiar ones and not something array-specific.
//
//   Uniform   -13.3 dB sidelobes, narrowest mainlobe
//   Hann      -31.5 dB, about 1.6x the width
//   Hamming   -42.7 dB, about 1.5x
//   Blackman  -58.1 dB, about 1.9x
[[nodiscard]] Real shaded_array_factor(const LineArray& array, Shading shading,
                                       Real lambda_m, Real steer_rad,
                                       Real look_rad) noexcept;

// Shading loss: the gain given up relative to uniform weighting, in dB.
// Equal to -10 log10( (sum w)^2 / (N sum w^2) ), always <= 0.
[[nodiscard]] Real shading_loss_db(Shading shading, std::size_t count) noexcept;

// ---------------------------------------------------------------------------
// Wideband beamforming
// ---------------------------------------------------------------------------
//
// Two stages, because the memory-versus-time trade is real and belongs to the
// caller. Transforming every element once costs N*M complex values of storage;
// re-transforming per beam costs nothing but runs N times as many FFTs per
// beam. For a 24-element array and 61 beams that is 1464 transforms rather than
// 24, so the storage is almost always the right side of the trade.

// Transforms each element's block into `spectra`, which must be N*fft_size long
// and is laid out element-major. `scratch` is fft_size long.
//
// `elements` is element-major too: element n occupies
// [n*fft_size, (n+1)*fft_size).
bool prepare_element_spectra(const FftView& fft,
                             std::size_t element_count,
                             std::span<const Real> elements,
                             std::span<Complex> spectra,
                             std::span<Complex> scratch) noexcept;

// One beam from prepared spectra, steered to `steer_rad`.
//
// Bin k carries frequency k*fs/M for k < M/2 and (k-M)*fs/M above it; the
// negative half must be steered with the negative frequency or the beam comes
// out real-valued nonsense. Element n is delayed by x_n sin(theta)/c, applied
// as exp(-j 2 pi f_k tau_n) -- an exact fractional delay, not an interpolation.
bool beamform_wideband(const FftView& fft,
                       const LineArray& array,
                       std::span<const Complex> spectra,
                       Real sample_rate_hz,
                       Real sound_speed_mps,
                       Real steer_rad,
                       Shading shading,
                       std::span<Complex> work,
                       std::span<Real> out_beam) noexcept;

// ---------------------------------------------------------------------------
// Adaptive beamforming (MVDR / Capon)
// ---------------------------------------------------------------------------

// Zeroes an N*N covariance accumulator.
void covariance_clear(std::size_t n, std::span<Complex> cov) noexcept;

// Accumulates one snapshot: R += x x^H.
bool covariance_accumulate(std::span<const Complex> snapshot,
                           std::size_t n, std::span<Complex> cov) noexcept;

// Divides the accumulator by the snapshot count.
void covariance_normalise(std::size_t n, std::size_t snapshots,
                          std::span<Complex> cov) noexcept;

// MVDR spatial spectrum: P(theta) = 1 / (a^H R^-1 a).
//
// `loading` is diagonal loading as a FRACTION of the mean diagonal power.
// It is not optional: with fewer snapshots than elements R is singular, and
// even with enough it is ill-conditioned whenever a source is strong. Loading
// trades resolution for stability, and 0.01 is a reasonable default.
//
// `work` must be at least n*n + 2*n complex: the Cholesky factor plus two
// vectors. Nothing is allocated.
//
// Returns the number of angles written, 0 on bad input or a non-positive-
// definite covariance (which means the loading was too small).
std::size_t mvdr_power(const LineArray& array, Real lambda_m,
                       std::span<const Complex> cov,
                       Real loading,
                       Real angle_min_rad, Real angle_max_rad,
                       std::span<Complex> work,
                       std::span<Real> out_power) noexcept;

// Forward-backward spatial smoothing.
//
// MVDR fails on COHERENT sources -- and the multipath this same library
// produces in v0.4 is exactly that: one arrival reaching the array by two
// paths, perfectly correlated. Two coherent arrivals make the covariance
// rank-deficient in a way loading cannot fix, and the adaptive beamformer nulls
// the target along with the interferer.
//
// The standard remedy: average the covariances of overlapping subarrays, which
// randomises the relative phase between the coherent pair and restores rank.
// Forward-backward adds the exchange-reversed conjugate, doubling the effective
// snapshot count for free on a uniform line array.
//
//   R_smooth = (1 / 2K) sum_k [ R_k + J conj(R_k) J ]
//
// The cost is aperture: an N-element array smoothed with subarrays of L
// elements has the resolution of an L-element one. Resolving P coherent sources
// needs L > P, so the trade is real and the caller must make it.
//
// `cov` is N*N, `out` is L*L. Returns false on bad sizes.
bool spatial_smooth(std::span<const Complex> cov, std::size_t n,
                    std::size_t subarray, std::span<Complex> out) noexcept;

// Angular separation below which conventional beamforming cannot resolve two
// equal sources, whatever the SNR: the Rayleigh criterion, one null-to-peak
// spacing. MVDR is not bound by it, which is the whole point of it.
[[nodiscard]] Real conventional_resolution_limit_rad(const LineArray& array,
                                                     Real lambda_m,
                                                     Real bearing_rad) noexcept;

}  // namespace phantom

#endif  // PHANTOM_BEAMFORMER_HPP
