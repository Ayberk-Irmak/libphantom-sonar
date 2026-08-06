// SPDX-License-Identifier: Apache-2.0
#include "phantom/fft.hpp"

#include <cmath>

namespace phantom {
namespace {

// Butterflies are written on real and imaginary parts by hand rather than with
// std::complex operator*. The standard requires operator* to handle infinities
// and NaNs specially, which on several compilers emits a call to a runtime
// helper (__mulsc3) instead of four multiplies -- a large, silent cost in the
// innermost loop of the whole DSP chain.
inline void butterfly(Complex& lo, Complex& hi, Real wr, Real wi) noexcept {
    const Real hr = hi.real();
    const Real hi_ = hi.imag();
    const Real vr = hr * wr - hi_ * wi;
    const Real vi = hr * wi + hi_ * wr;
    const Real lr = lo.real();
    const Real li = lo.imag();
    lo = Complex(lr + vr, li + vi);
    hi = Complex(lr - vr, li - vi);
}

void bit_reverse(std::span<Complex> a) noexcept {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            const Complex t = a[i];
            a[i] = a[j];
            a[j] = t;
        }
    }
}

// Shared kernel. `Conjugate` flips the sign of the twiddle imaginary part,
// turning the forward transform into the unscaled inverse.
//
// It is a template parameter, not a runtime flag, on purpose: as a `bool`
// argument the sign select sat inside the innermost loop, costing a branch per
// butterfly and blocking vectorisation of the whole kernel. Hoisting it to
// compile time is worth roughly 3x on an 8192-point transform.
template <bool Conjugate>
void transform_impl(const FftView& plan, std::span<Complex> data) noexcept {
    const std::size_t n = data.size();
    bit_reverse(data);

    Complex* __restrict a = data.data();
    const Complex* __restrict tw = plan.twiddles.data();

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const std::size_t stride = n / len;
        for (std::size_t k = 0; k < half; ++k) {
            // The twiddle is constant across every block at this stage, so it
            // is loaded once here rather than once per butterfly.
            const Complex w = tw[k * stride];
            const Real wr = w.real();
            const Real wi = Conjugate ? -w.imag() : w.imag();
            for (std::size_t base = 0; base < n; base += len) {
                butterfly(a[base + k], a[base + k + half], wr, wi);
            }
        }
    }
}

}  // namespace

bool fft_init(std::span<Complex> twiddles) noexcept {
    const std::size_t half = twiddles.size();
    if (half == 0) return false;
    const std::size_t n = half * 2;
    if (!is_power_of_two(n)) return false;

    // Computed from the angle each time rather than by repeated complex
    // multiplication: recurrence accumulates phase error that would show up as
    // a raised noise floor in the matched filter, which is exactly where it
    // would be hardest to attribute.
    const Real scale = static_cast<Real>(-2) * kPi / static_cast<Real>(n);
    for (std::size_t k = 0; k < half; ++k) {
        const Real angle = scale * static_cast<Real>(k);
        twiddles[k] = Complex(std::cos(angle), std::sin(angle));
    }
    return true;
}

void fft_forward(const FftView& plan, std::span<Complex> data) noexcept {
    if (!plan.valid() || data.size() != plan.size) return;
    transform_impl<false>(plan, data);
}

void fft_inverse(const FftView& plan, std::span<Complex> data) noexcept {
    if (!plan.valid() || data.size() != plan.size) return;
    transform_impl<true>(plan, data);
    const Real inv = static_cast<Real>(1) / static_cast<Real>(data.size());
    for (Complex& c : data) c *= inv;
}

}  // namespace phantom
