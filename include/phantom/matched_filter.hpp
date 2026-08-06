// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — replica correlation (matched filtering).
//
// The matched filter is the optimal detector for a known waveform in white
// Gaussian noise, and it is what turns a long low-amplitude ping into a short
// tall peak. Correlating a REAL received signal against a COMPLEX (analytic)
// replica yields the complex envelope directly, so no separate basebanding
// stage is needed:
//
//     y[n] = sum_k  x[n+k] * conj(r[k])
//
// Implemented by FFT overlap-save. A direct correlation costs L multiplies per
// output sample -- for a 100 ms chirp at 192 kHz that is 19200, or ~3.7 GMAC/s,
// which no embedded target will do. The transform reduces it to O(log M).
//
// Circular-correlation aliasing confines the valid outputs to lags [0, M-L]:
// for those, the window x[n .. n+L-1] never wraps. That is M-L+1 usable lags
// per block, and it is exactly how far the caller must advance the input.
#ifndef PHANTOM_MATCHED_FILTER_HPP
#define PHANTOM_MATCHED_FILTER_HPP

#include "phantom/fft.hpp"
#include "phantom/types.hpp"
#include "phantom/waveform.hpp"

#include <span>

namespace phantom {

// A prepared replica: its conjugated spectrum, its length and its energy.
struct MatchedFilter {
    std::span<const Complex> spectrum;  // conj(FFT(zero-padded replica)), length M
    std::size_t replica_length = 0;     // L
    Real replica_energy = 0;            // sum |r[k]|^2

    [[nodiscard]] constexpr bool valid() const noexcept {
        return replica_length >= 1 && spectrum.size() >= replica_length
            && replica_energy > 0;
    }
};

// Number of usable correlation lags produced from one block of `fft_size`
// samples by a replica of `replica_length`. Also the block advance.
[[nodiscard]] constexpr std::size_t matched_filter_stride(std::size_t fft_size,
                                                          std::size_t replica_length) noexcept {
    return (fft_size >= replica_length) ? fft_size - replica_length + 1 : 0;
}

// Prepares a replica for repeated use. `spectrum_out` must be `fft_size` long
// and is retained by the returned MatchedFilter, so it must outlive it.
bool matched_filter_prepare(const FftView& fft,
                            std::span<const Complex> replica,
                            std::span<Complex> spectrum_out,
                            MatchedFilter& out) noexcept;

// Convenience: renders `spec` and prepares it in one step. `replica_scratch`
// must hold at least the pulse length.
bool matched_filter_from_pulse(const FftView& fft,
                               const PulseSpec& spec,
                               Real sample_rate_hz,
                               std::span<Complex> replica_scratch,
                               std::span<Complex> spectrum_out,
                               MatchedFilter& out) noexcept;

// Correlates one block. `block` and `scratch` are both `fft_size` long; `out`
// receives the first `matched_filter_stride()` lags. Returns the lag count.
std::size_t matched_filter_apply(const FftView& fft,
                                 const MatchedFilter& mf,
                                 std::span<const Real> block,
                                 std::span<Complex> scratch,
                                 std::span<Complex> out) noexcept;

// Direct O(N*L) correlation. Reference implementation: the FFT path is checked
// against this in the test suite. Do not use it in a control loop.
std::size_t correlate_direct(std::span<const Real> signal,
                             std::span<const Complex> replica,
                             std::span<Complex> out) noexcept;

// Refines a peak to sub-sample resolution by fitting a parabola through the
// three magnitude samples straddling it. Returns the offset in samples from
// `peak_index`, in [-0.5, +0.5]; 0 when the peak sits on a block edge.
[[nodiscard]] Real parabolic_peak_offset(std::span<const Real> magnitude,
                                         std::size_t peak_index) noexcept;

}  // namespace phantom

#endif  // PHANTOM_MATCHED_FILTER_HPP
