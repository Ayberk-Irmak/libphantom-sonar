// SPDX-License-Identifier: Apache-2.0
#include "phantom/beamformer.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

// Complex Cholesky, R = L L^H, in place over the lower triangle.
// Returns false if R is not positive definite, which for a covariance matrix
// means the diagonal loading was too small for the snapshot count.
bool cholesky(std::span<Complex> a, std::size_t n) noexcept {
    for (std::size_t j = 0; j < n; ++j) {
        Real d = a[j * n + j].real();
        for (std::size_t k = 0; k < j; ++k) {
            d -= std::norm(a[j * n + k]);
        }
        if (!(d > kZero)) return false;
        const Real ljj = std::sqrt(d);
        a[j * n + j] = Complex(ljj, kZero);
        for (std::size_t i = j + 1; i < n; ++i) {
            Complex s = a[i * n + j];
            for (std::size_t k = 0; k < j; ++k) {
                s -= a[i * n + k] * std::conj(a[j * n + k]);
            }
            a[i * n + j] = s / ljj;
        }
    }
    return true;
}

// Solves L L^H x = b with L lower-triangular, in place on x.
void cholesky_solve(std::span<const Complex> l, std::size_t n,
                    std::span<Complex> x) noexcept {
    // Forward substitution: L y = b
    for (std::size_t i = 0; i < n; ++i) {
        Complex s = x[i];
        for (std::size_t k = 0; k < i; ++k) s -= l[i * n + k] * x[k];
        x[i] = s / l[i * n + i].real();
    }
    // Back substitution: L^H x = y
    for (std::size_t ii = n; ii-- > 0;) {
        Complex s = x[ii];
        for (std::size_t k = ii + 1; k < n; ++k) s -= std::conj(l[k * n + ii]) * x[k];
        x[ii] = s / l[ii * n + ii].real();
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Shading
// ---------------------------------------------------------------------------

const char* shading_name(Shading s) noexcept {
    switch (s) {
        case Shading::Uniform:  return "uniform";
        case Shading::Hann:     return "Hann";
        case Shading::Hamming:  return "Hamming";
        case Shading::Blackman: return "Blackman";
    }
    return "unknown";
}

Real shading_weight(Shading shading, std::size_t index, std::size_t count) noexcept {
    if (count == 0 || index >= count) return kZero;
    if (shading == Shading::Uniform || count == 1) return kOne;

    // Symmetric windows over [0, 1]. The endpoints of Hann land on zero, which
    // is correct: those elements contribute nothing and the effective aperture
    // is slightly smaller than the physical one.
    const Real x = static_cast<Real>(index) / static_cast<Real>(count - 1);
    const Real a = kTwo * kPi * x;
    switch (shading) {
        case Shading::Uniform:
            return kOne;
        case Shading::Hann:
            return static_cast<Real>(0.5) * (kOne - std::cos(a));
        case Shading::Hamming:
            return static_cast<Real>(0.54) - static_cast<Real>(0.46) * std::cos(a);
        case Shading::Blackman:
            return static_cast<Real>(0.42) - static_cast<Real>(0.5) * std::cos(a)
                 + static_cast<Real>(0.08) * std::cos(kTwo * a);
    }
    return kOne;
}

Real shaded_array_factor(const LineArray& array, Shading shading,
                         Real lambda_m, Real steer_rad, Real look_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return kZero;
    const Real k = kTwo * kPi / lambda_m;
    const Real psi = k * array.spacing_m * (std::sin(look_rad) - std::sin(steer_rad));

    Real re = kZero;
    Real im = kZero;
    Real norm = kZero;
    for (std::size_t i = 0; i < array.element_count; ++i) {
        const Real w = shading_weight(shading, i, array.element_count);
        const Real ph = psi * static_cast<Real>(i);
        re += w * std::cos(ph);
        im += w * std::sin(ph);
        norm += w;
    }
    if (!(norm > kZero)) return kZero;
    return std::sqrt(re * re + im * im) / norm;
}

Real shading_loss_db(Shading shading, std::size_t count) noexcept {
    if (count == 0) return kZero;
    Real sum = kZero;
    Real sum_sq = kZero;
    for (std::size_t i = 0; i < count; ++i) {
        const Real w = shading_weight(shading, i, count);
        sum += w;
        sum_sq += w * w;
    }
    if (!(sum > kZero) || !(sum_sq > kZero)) return kZero;
    // Coherent gain squared over incoherent, relative to uniform.
    const Real eff = (sum * sum) / (static_cast<Real>(count) * sum_sq);
    return static_cast<Real>(10) * std::log10(eff);
}

// ---------------------------------------------------------------------------
// Wideband beamforming
// ---------------------------------------------------------------------------

bool prepare_element_spectra(const FftView& fft,
                             std::size_t element_count,
                             std::span<const Real> elements,
                             std::span<Complex> spectra,
                             std::span<Complex> scratch) noexcept {
    if (!fft.valid() || element_count == 0) return false;
    const std::size_t m = fft.size;
    if (elements.size() < element_count * m) return false;
    if (spectra.size() < element_count * m) return false;
    if (scratch.size() < m) return false;

    for (std::size_t n = 0; n < element_count; ++n) {
        for (std::size_t i = 0; i < m; ++i) {
            scratch[i] = Complex(elements[n * m + i], kZero);
        }
        fft_forward(fft, scratch.subspan(0, m));
        for (std::size_t i = 0; i < m; ++i) spectra[n * m + i] = scratch[i];
    }
    return true;
}

bool beamform_wideband(const FftView& fft,
                       const LineArray& array,
                       std::span<const Complex> spectra,
                       Real sample_rate_hz,
                       Real sound_speed_mps,
                       Real steer_rad,
                       Shading shading,
                       std::span<Complex> work,
                       std::span<Real> out_beam) noexcept {
    if (!fft.valid() || !array.valid()) return false;
    if (!(sample_rate_hz > kZero) || !(sound_speed_mps > kZero)) return false;
    const std::size_t m = fft.size;
    if (spectra.size() < array.element_count * m) return false;
    if (work.size() < m || out_beam.size() < m) return false;

    for (std::size_t i = 0; i < m; ++i) work[i] = Complex(kZero, kZero);

    const Real centre = static_cast<Real>(array.element_count - 1) / kTwo;
    const Real sin_steer = std::sin(steer_rad);
    const Real df = sample_rate_hz / static_cast<Real>(m);

    Real weight_sum = kZero;
    for (std::size_t n = 0; n < array.element_count; ++n) {
        weight_sum += shading_weight(shading, n, array.element_count);
    }
    if (!(weight_sum > kZero)) return false;

    for (std::size_t n = 0; n < array.element_count; ++n) {
        const Real w = shading_weight(shading, n, array.element_count)
                     / weight_sum;
        if (w == kZero) continue;
        const Real x = (static_cast<Real>(n) - centre) * array.spacing_m;
        const Real tau = x * sin_steer / sound_speed_mps;

        for (std::size_t bin = 0; bin < m; ++bin) {
            // Bins above M/2 carry NEGATIVE frequencies. Steering them with the
            // positive value leaves the beam Hermitian-broken and the inverse
            // transform comes back as noise that still looks like a signal.
            const Real f = (bin <= m / 2)
                         ? static_cast<Real>(bin) * df
                         : (static_cast<Real>(bin) - static_cast<Real>(m)) * df;
            const Real ph = -kTwo * kPi * f * tau;
            const Real c = std::cos(ph);
            const Real s = std::sin(ph);
            const Complex& e = spectra[n * m + bin];
            const Real er = e.real();
            const Real ei = e.imag();
            work[bin] += Complex(w * (er * c - ei * s), w * (er * s + ei * c));
        }
    }

    fft_inverse(fft, work.subspan(0, m));
    for (std::size_t i = 0; i < m; ++i) out_beam[i] = work[i].real();
    return true;
}

// ---------------------------------------------------------------------------
// MVDR
// ---------------------------------------------------------------------------

void covariance_clear(std::size_t n, std::span<Complex> cov) noexcept {
    if (cov.size() < n * n) return;
    for (std::size_t i = 0; i < n * n; ++i) cov[i] = Complex(kZero, kZero);
}

bool covariance_accumulate(std::span<const Complex> snapshot,
                           std::size_t n, std::span<Complex> cov) noexcept {
    if (n == 0 || snapshot.size() < n || cov.size() < n * n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            cov[i * n + j] += snapshot[i] * std::conj(snapshot[j]);
        }
    }
    return true;
}

void covariance_normalise(std::size_t n, std::size_t snapshots,
                          std::span<Complex> cov) noexcept {
    if (snapshots == 0 || cov.size() < n * n) return;
    const Real inv = kOne / static_cast<Real>(snapshots);
    for (std::size_t i = 0; i < n * n; ++i) cov[i] *= inv;
}

std::size_t mvdr_power(const LineArray& array, Real lambda_m,
                       std::span<const Complex> cov,
                       Real loading,
                       Real angle_min_rad, Real angle_max_rad,
                       std::span<Complex> work,
                       std::span<Real> out_power) noexcept {
    if (!array.valid() || !(lambda_m > kZero) || out_power.empty()) return 0;
    if (!(angle_max_rad > angle_min_rad)) return 0;
    const std::size_t n = array.element_count;
    if (cov.size() < n * n || work.size() < n * n + 2 * n) return 0;

    std::span<Complex> l = work.subspan(0, n * n);
    std::span<Complex> a = work.subspan(n * n, n);
    std::span<Complex> x = work.subspan(n * n + n, n);

    // Diagonal loading as a fraction of the mean diagonal power. Without it a
    // covariance built from fewer snapshots than elements is singular, and one
    // with a strong source is ill-conditioned regardless.
    Real trace = kZero;
    for (std::size_t i = 0; i < n; ++i) trace += cov[i * n + i].real();
    const Real eps = (loading > kZero) ? loading * trace / static_cast<Real>(n) : kZero;

    for (std::size_t i = 0; i < n * n; ++i) l[i] = cov[i];
    for (std::size_t i = 0; i < n; ++i) l[i * n + i] += Complex(eps, kZero);
    if (!cholesky(l, n)) return 0;

    const Real k = kTwo * kPi / lambda_m;
    const Real centre = static_cast<Real>(n - 1) / kTwo;
    const std::size_t n_angles = out_power.size();
    const Real step = (n_angles > 1)
                    ? (angle_max_rad - angle_min_rad) / static_cast<Real>(n_angles - 1)
                    : kZero;

    for (std::size_t ai = 0; ai < n_angles; ++ai) {
        const Real steer = angle_min_rad + step * static_cast<Real>(ai);
        const Real s = std::sin(steer);
        for (std::size_t i = 0; i < n; ++i) {
            const Real pos = (static_cast<Real>(i) - centre) * array.spacing_m;
            const Real ph = k * pos * s;
            a[i] = Complex(std::cos(ph), std::sin(ph));
            x[i] = a[i];
        }
        cholesky_solve(l, n, x);

        // a^H R^-1 a is real and positive for a Hermitian positive-definite R;
        // any imaginary part is rounding.
        Real denom = kZero;
        for (std::size_t i = 0; i < n; ++i) {
            denom += (std::conj(a[i]) * x[i]).real();
        }
        out_power[ai] = (denom > kZero) ? kOne / denom : kZero;
    }
    return n_angles;
}

bool spatial_smooth(std::span<const Complex> cov, std::size_t n,
                    std::size_t subarray, std::span<Complex> out) noexcept {
    if (n < 2 || subarray < 1 || subarray > n) return false;
    if (cov.size() < n * n || out.size() < subarray * subarray) return false;

    const std::size_t k = n - subarray + 1;   // number of subarrays
    for (std::size_t i = 0; i < subarray * subarray; ++i) out[i] = Complex(kZero, kZero);

    for (std::size_t s = 0; s < k; ++s) {
        for (std::size_t i = 0; i < subarray; ++i) {
            for (std::size_t j = 0; j < subarray; ++j) {
                out[i * subarray + j] += cov[(i + s) * n + (j + s)];
            }
        }
    }

    // Backward: J conj(R) J, with J the exchange matrix, so
    //     R_fb[i][j] = R_f[i][j] + conj(R_f[L-1-i][L-1-j])
    // and since R_f is Hermitian that conjugate is just R_f[L-1-j][L-1-i].
    //
    // Written that way each entry pairs with exactly one other, (i,j) with
    // (L-1-j, L-1-i), and BOTH receive the same sum -- forward-backward
    // smoothing produces a persymmetric matrix, which is the check that the
    // indices are the right way round. Handling the pair together is also what
    // makes the in-place update safe: the naive version overwrites entries it
    // has yet to read, and the result is a covariance that still looks
    // plausible and puts every bearing in the wrong place.
    for (std::size_t i = 0; i < subarray; ++i) {
        for (std::size_t j = 0; j < subarray; ++j) {
            const std::size_t idx = i * subarray + j;
            const std::size_t pidx = (subarray - 1 - j) * subarray + (subarray - 1 - i);
            if (idx > pidx) continue;
            const Complex sum = out[idx] + out[pidx];
            out[idx] = sum;
            out[pidx] = sum;
        }
    }

    const Real inv = kOne / static_cast<Real>(2 * k);
    for (std::size_t i = 0; i < subarray * subarray; ++i) out[i] *= inv;
    return true;
}

Real conventional_resolution_limit_rad(const LineArray& array, Real lambda_m,
                                       Real bearing_rad) noexcept {
    if (!array.valid() || !(lambda_m > kZero)) return kZero;
    const Real cos_th = std::fabs(std::cos(bearing_rad));
    if (!(cos_th > static_cast<Real>(1e-6))) return kPi;
    const Real n = static_cast<Real>(array.element_count);
    // One null-to-peak spacing: lambda / (N d cos theta).
    return lambda_m / (n * array.spacing_m * cos_th);
}

}  // namespace phantom
