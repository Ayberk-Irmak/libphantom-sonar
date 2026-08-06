// SPDX-License-Identifier: Apache-2.0
#include "phantom/matched_filter.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);

}  // namespace

bool matched_filter_prepare(const FftView& fft,
                            std::span<const Complex> replica,
                            std::span<Complex> spectrum_out,
                            MatchedFilter& out) noexcept {
    if (!fft.valid()) return false;
    const std::size_t m = fft.size;
    if (spectrum_out.size() != m) return false;
    if (replica.empty() || replica.size() > m) return false;

    Real energy = kZero;
    for (std::size_t i = 0; i < replica.size(); ++i) {
        spectrum_out[i] = replica[i];
        const Real re = replica[i].real();
        const Real im = replica[i].imag();
        energy += re * re + im * im;
    }
    for (std::size_t i = replica.size(); i < m; ++i) spectrum_out[i] = Complex(kZero, kZero);
    if (!(energy > kZero)) return false;

    fft_forward(fft, spectrum_out);
    // Correlation, not convolution: Y = X * conj(R).
    for (Complex& c : spectrum_out) c = std::conj(c);

    out.spectrum = std::span<const Complex>(spectrum_out.data(), spectrum_out.size());
    out.replica_length = replica.size();
    out.replica_energy = energy;
    return true;
}

bool matched_filter_from_pulse(const FftView& fft,
                               const PulseSpec& spec,
                               Real sample_rate_hz,
                               std::span<Complex> replica_scratch,
                               std::span<Complex> spectrum_out,
                               MatchedFilter& out) noexcept {
    const std::size_t n = render_analytic(spec, sample_rate_hz, replica_scratch);
    if (n == 0) return false;
    return matched_filter_prepare(fft, replica_scratch.subspan(0, n), spectrum_out, out);
}

std::size_t matched_filter_apply(const FftView& fft,
                                 const MatchedFilter& mf,
                                 std::span<const Real> block,
                                 std::span<Complex> scratch,
                                 std::span<Complex> out) noexcept {
    if (!fft.valid() || !mf.valid()) return 0;
    const std::size_t m = fft.size;
    if (block.size() != m || scratch.size() != m) return 0;
    if (mf.spectrum.size() != m) return 0;

    const std::size_t lags = matched_filter_stride(m, mf.replica_length);
    if (lags == 0 || out.size() < lags) return 0;

    for (std::size_t i = 0; i < m; ++i) scratch[i] = Complex(block[i], kZero);
    fft_forward(fft, scratch);

    for (std::size_t i = 0; i < m; ++i) {
        const Real ar = scratch[i].real(), ai = scratch[i].imag();
        const Real br = mf.spectrum[i].real(), bi = mf.spectrum[i].imag();
        scratch[i] = Complex(ar * br - ai * bi, ar * bi + ai * br);
    }
    fft_inverse(fft, scratch);

    for (std::size_t i = 0; i < lags; ++i) out[i] = scratch[i];
    return lags;
}

std::size_t correlate_direct(std::span<const Real> signal,
                             std::span<const Complex> replica,
                             std::span<Complex> out) noexcept {
    if (replica.empty() || signal.size() < replica.size()) return 0;
    const std::size_t lags = signal.size() - replica.size() + 1;
    if (out.size() < lags) return 0;

    for (std::size_t n = 0; n < lags; ++n) {
        Real re = kZero;
        Real im = kZero;
        for (std::size_t k = 0; k < replica.size(); ++k) {
            const Real x = signal[n + k];
            re += x * replica[k].real();
            im -= x * replica[k].imag();  // conj(r[k])
        }
        out[n] = Complex(re, im);
    }
    return lags;
}

Real parabolic_peak_offset(std::span<const Real> magnitude, std::size_t peak_index) noexcept {
    if (peak_index == 0 || peak_index + 1 >= magnitude.size()) return kZero;
    const Real ym = magnitude[peak_index - 1];
    const Real y0 = magnitude[peak_index];
    const Real yp = magnitude[peak_index + 1];
    const Real denom = ym - static_cast<Real>(2) * y0 + yp;
    // A non-negative denominator means the samples do not bracket a maximum
    // (flat or a minimum), so there is nothing to refine.
    if (!(denom < kZero)) return kZero;
    const Real delta = static_cast<Real>(0.5) * (ym - yp) / denom;
    return clamp(delta, static_cast<Real>(-0.5), static_cast<Real>(0.5));
}

}  // namespace phantom
