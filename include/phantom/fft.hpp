// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — in-place radix-2 FFT, zero allocation.
//
// Written rather than pulled in because the library promises no dependencies,
// and every stage of the ping-analysis chain needs a transform. It is a
// textbook iterative Cooley-Tukey: decimation-in-time, bit-reversal
// permutation, precomputed twiddles. No cleverness, verified against a direct
// DFT to machine precision.
//
// Buffers are caller-owned, matching the rest of the library: the twiddle table
// is filled once into storage you provide, and every transform runs in place.
#ifndef PHANTOM_FFT_HPP
#define PHANTOM_FFT_HPP

#include "phantom/types.hpp"

#include <array>
#include <complex>
#include <span>

namespace phantom {

using Complex = std::complex<Real>;

// Non-owning view of a prepared transform. `twiddles` holds n/2 entries.
struct FftView {
    std::span<const Complex> twiddles;
    std::size_t size = 0;  // transform length, a power of two

    [[nodiscard]] constexpr bool valid() const noexcept {
        return size >= 2 && (size & (size - 1)) == 0 && twiddles.size() == size / 2;
    }
};

[[nodiscard]] constexpr bool is_power_of_two(std::size_t n) noexcept {
    return n >= 2 && (n & (n - 1)) == 0;
}

// Fills `twiddles` with W_n^k = exp(-2*pi*i*k/n) for k in [0, n/2).
// Returns false when the span length is not a valid half-transform size.
bool fft_init(std::span<Complex> twiddles) noexcept;

// In-place forward transform. `data.size()` must equal `plan.size`.
// Unnormalised: forward followed by inverse reproduces the input exactly.
void fft_forward(const FftView& plan, std::span<Complex> data) noexcept;

// In-place inverse transform, scaled by 1/n.
void fft_inverse(const FftView& plan, std::span<Complex> data) noexcept;

// Owning fixed-size plan. The twiddle table costs N/2 complex values.
template <std::size_t N>
class FftPlan {
    static_assert(is_power_of_two(N), "FFT size must be a power of two >= 2");

  public:
    FftPlan() noexcept { fft_init(std::span<Complex>(twiddles_.data(), twiddles_.size())); }

    [[nodiscard]] FftView view() const noexcept {
        return FftView{std::span<const Complex>(twiddles_.data(), twiddles_.size()), N};
    }

    void forward(std::span<Complex> data) const noexcept { fft_forward(view(), data); }
    void inverse(std::span<Complex> data) const noexcept { fft_inverse(view(), data); }

    [[nodiscard]] static constexpr std::size_t size() noexcept { return N; }

  private:
    std::array<Complex, N / 2> twiddles_{};
};

}  // namespace phantom

#endif  // PHANTOM_FFT_HPP
