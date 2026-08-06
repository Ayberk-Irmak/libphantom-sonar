// SPDX-License-Identifier: Apache-2.0
#include "phantom/ping_analyzer.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// Floor for the CFAR noise estimate. On a synthetic noise-free block the
// training cells can be exactly zero, and every cell would then clear any
// finite threshold.
constexpr Real kNoiseFloor = static_cast<Real>(1e-30);

}  // namespace

bool AnalyzerView::valid() const noexcept {
    return !templates.empty() && fft.valid() && sample_rate_hz > kZero
        && max_replica_length >= 1 && max_replica_length <= fft.size;
}

std::size_t AnalyzerView::stride() const noexcept {
    if (!valid()) return 0;
    return matched_filter_stride(fft.size, max_replica_length);
}

bool AnalyzerWorkspace::valid(std::size_t fft_size) const noexcept {
    return fft_scratch.size() >= fft_size
        && correlation.size() >= fft_size
        && best_power.size() >= fft_size
        && best_template.size() >= fft_size
        && prefix.size() >= fft_size + 1;
}

std::size_t suggested_cfar_guard(const AnalyzerView& bank) noexcept {
    return bank.valid() ? bank.max_replica_length : 0;
}

Real suggested_dead_time_s(const AnalyzerView& bank) noexcept {
    if (!bank.valid()) return kZero;
    return static_cast<Real>(bank.max_replica_length) / bank.sample_rate_hz;
}

Real cfar_alpha(std::size_t training_cells, Real pfa) noexcept {
    if (training_cells == 0 || !(pfa > kZero) || !(pfa < kOne)) return kZero;
    const Real n = static_cast<Real>(training_cells);
    return n * (std::pow(pfa, -kOne / n) - kOne);
}

std::size_t analyze_block(const AnalyzerView& bank,
                          const DetectorConfig& cfg,
                          std::span<const Real> block,
                          Real block_start_time_s,
                          const AnalyzerWorkspace& work,
                          std::span<PulseDescriptor> out) noexcept {
    if (!bank.valid()) return 0;
    const std::size_t m = bank.fft.size;
    if (block.size() != m || !work.valid(m) || out.empty()) return 0;

    const std::size_t lags = bank.stride();
    if (lags == 0) return 0;
    if (bank.templates.size() > 0xFFFFu) return 0;

    // --- 1. Run the bank, keeping the winner per lag ------------------------
    for (std::size_t i = 0; i < lags; ++i) {
        work.best_power[i] = kZero;
        work.best_template[i] = 0;
    }

    for (std::size_t t = 0; t < bank.templates.size(); ++t) {
        const MatchedFilter& mf = bank.templates[t].filter;
        const std::size_t produced = matched_filter_apply(bank.fft, mf, block,
                                                          work.fft_scratch.subspan(0, m),
                                                          work.correlation);
        if (produced < lags) continue;

        // Peaks are compared across templates as |y|^2 / E, i.e. normalised by
        // replica energy. Without that a longer replica would win on every lag
        // purely for being longer, and the bank would always report the longest
        // pulse it knows about.
        const Real inv_energy = kOne / mf.replica_energy;
        for (std::size_t i = 0; i < lags; ++i) {
            const Real re = work.correlation[i].real();
            const Real im = work.correlation[i].imag();
            const Real p = (re * re + im * im) * inv_energy;
            if (p > work.best_power[i]) {
                work.best_power[i] = p;
                work.best_template[i] = static_cast<std::uint16_t>(t);
            }
        }
    }

    // --- 2. Prefix sums, so each CFAR window costs O(1) ---------------------
    work.prefix[0] = kZero;
    for (std::size_t i = 0; i < lags; ++i) {
        work.prefix[i + 1] = work.prefix[i] + work.best_power[i];
    }

    const std::size_t guard = cfg.cfar_guard;
    const std::size_t train = (cfg.cfar_train > 0) ? cfg.cfar_train : 1;

    // Cell-averaging CFAR over the training cells either side of the guard
    // band, truncated at the block edges rather than wrapped.
    auto noise_estimate = [&](std::size_t i) noexcept -> Real {
        const std::size_t lo_hi = (i > guard) ? (i - guard) : 0;
        const std::size_t lo_lo = (lo_hi > train) ? (lo_hi - train) : 0;
        const std::size_t hi_lo = (i + guard + 1 < lags) ? (i + guard + 1) : lags;
        const std::size_t hi_hi = (hi_lo + train < lags) ? (hi_lo + train) : lags;

        const std::size_t count = (lo_hi - lo_lo) + (hi_hi - hi_lo);
        if (count == 0) return kNoiseFloor;
        const Real sum = (work.prefix[lo_hi] - work.prefix[lo_lo])
                       + (work.prefix[hi_hi] - work.prefix[hi_lo]);
        const Real mean = sum / static_cast<Real>(count);
        return (mean > kNoiseFloor) ? mean : kNoiseFloor;
    };

    // --- 3. Threshold, pick peaks, enforce dead time ------------------------
    const Real fs = bank.sample_rate_hz;
    auto dead_samples = static_cast<std::size_t>(
        (cfg.dead_time_s > kZero) ? cfg.dead_time_s * fs : kZero);
    if (dead_samples == 0) dead_samples = 1;

    std::size_t written = 0;
    std::size_t i = 0;
    while (i < lags && written < out.size()) {
        const Real noise = noise_estimate(i);
        if (!(work.best_power[i] > cfg.threshold_alpha * noise)) {
            ++i;
            continue;
        }

        // Take the strongest cell inside the dead-time window, not the first
        // one over threshold: a compressed pulse spans several cells and the
        // leading edge is not its peak.
        const std::size_t window_end = (i + dead_samples < lags) ? (i + dead_samples) : lags;
        std::size_t peak = i;
        for (std::size_t j = i; j < window_end; ++j) {
            if (work.best_power[j] > work.best_power[peak]) peak = j;
        }

        const std::uint16_t t = work.best_template[peak];
        const PulseTemplate& tpl = bank.templates[t];

        // Sub-sample refinement on magnitude, which is what the parabola
        // assumption is a good fit for near a compressed mainlobe.
        Real offset = kZero;
        if (peak > 0 && peak + 1 < lags) {
            const Real a = std::sqrt(work.best_power[peak - 1]);
            const Real b = std::sqrt(work.best_power[peak]);
            const Real c = std::sqrt(work.best_power[peak + 1]);
            const Real denom = a - static_cast<Real>(2) * b + c;
            if (denom < kZero) {
                offset = clamp(static_cast<Real>(0.5) * (a - c) / denom,
                               static_cast<Real>(-0.5), static_cast<Real>(0.5));
            }
        }

        const Real peak_noise = noise_estimate(peak);
        const Real normalised_peak = work.best_power[peak];

        PulseDescriptor& pdw = out[written++];
        pdw.toa_s = block_start_time_s
                  + (static_cast<Real>(peak) + offset) / fs;
        pdw.snr_db = static_cast<Real>(10) * std::log10(normalised_peak / peak_noise);
        // Undo the energy normalisation to recover the raw peak, then invert
        // |y_peak| = A * E / 2 for the received amplitude.
        pdw.peak_magnitude = std::sqrt(normalised_peak * tpl.filter.replica_energy);
        pdw.amplitude = static_cast<Real>(2) * pdw.peak_magnitude / tpl.filter.replica_energy;
        pdw.type = tpl.spec.type;
        pdw.centre_freq_hz = tpl.spec.centre_frequency_hz();
        pdw.bandwidth_hz = tpl.spec.bandwidth_hz();
        pdw.duration_s = tpl.spec.duration_s;
        pdw.chirp_rate_hz_s = tpl.spec.chirp_rate_hz_s();
        pdw.template_index = t;

        i = peak + dead_samples;
    }

    return written;
}

}  // namespace phantom
