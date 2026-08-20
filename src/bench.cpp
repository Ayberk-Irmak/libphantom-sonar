// SPDX-License-Identifier: Apache-2.0
#include "phantom/bench.hpp"

#include <algorithm>
#include <cmath>

namespace phantom::bench {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// Energy in a frequency band, from a real block, via the caller's transform.
Real band_energy(const FftView& fft, std::span<const Real> x,
                 Real fs, Real lo_hz, Real hi_hz,
                 std::span<Complex> scratch) noexcept {
    const std::size_t n = fft.size;
    if (scratch.size() < n) return kZero;
    Real total = kZero;
    const std::size_t hop = n / 2;
    std::size_t blocks = 0;
    for (std::size_t start = 0; start + n <= x.size(); start += hop) {
        for (std::size_t i = 0; i < n; ++i) {
            // Hann window: without it the spectral leakage from a loud
            // out-of-band transient lands squarely in the band being measured,
            // which is how a startup click passes for signal.
            const Real w = static_cast<Real>(0.5)
                         * (kOne - std::cos(static_cast<Real>(2) * kPi
                            * static_cast<Real>(i) / static_cast<Real>(n - 1)));
            scratch[i] = Complex(x[start + i] * w, kZero);
        }
        fft_forward(fft, scratch.subspan(0, n));
        const Real df = fs / static_cast<Real>(n);
        for (std::size_t k = 1; k < n / 2; ++k) {
            const Real f = static_cast<Real>(k) * df;
            if (f < lo_hz || f > hi_hz) continue;
            total += std::norm(scratch[k]);
        }
        ++blocks;
    }
    return (blocks > 0) ? total / static_cast<Real>(blocks) : kZero;
}

Real to_db(Real power) noexcept {
    constexpr Real kFloor = static_cast<Real>(1e-30);
    return static_cast<Real>(10) * std::log10((power > kFloor) ? power : kFloor);
}

}  // namespace

bool qualify_channel(const FftView& fft,
                     std::span<const Real> active,
                     std::span<const Real> silent,
                     const QualifyConfig& cfg,
                     std::span<Complex> work,
                     QualifyResult& out) noexcept {
    const std::size_t n = fft.size;
    if (!fft.valid() || work.size() < n) return false;
    if (active.size() < n || silent.size() < n) return false;
    if (!(cfg.sample_rate_hz > kZero)) return false;
    if (!(cfg.band_high_hz > cfg.band_low_hz)) return false;

    const Real ea = band_energy(fft, active, cfg.sample_rate_hz,
                                cfg.band_low_hz, cfg.band_high_hz, work);
    const Real es = band_energy(fft, silent, cfg.sample_rate_hz,
                                cfg.band_low_hz, cfg.band_high_hz, work);

    std::size_t clipped = 0;
    for (const Real v : active) {
        // 16-bit full scale in normalised units, with a sample of slack.
        if (std::fabs(v) >= static_cast<Real>(0.9999)) ++clipped;
    }

    out.active_db = to_db(ea);
    out.silent_db = to_db(es);
    out.excess_db = out.active_db - out.silent_db;
    out.clipped_fraction = static_cast<Real>(clipped) / static_cast<Real>(active.size());
    out.usable = (out.excess_db >= cfg.min_excess_db)
              && (out.clipped_fraction <= cfg.max_clipped_fraction);
    return true;
}

TwoDistanceResult two_distance_solve(Real d1, Real t1, Real d2, Real t2,
                                     Real timing_resolution_s) noexcept {
    TwoDistanceResult r;
    const Real dd = d2 - d1;
    const Real dt = t2 - t1;
    if (!(std::fabs(dt) > timing_resolution_s)) return r;   // not resolvable
    if (!(std::fabs(dd) > kZero)) return r;
    const Real c = dd / dt;
    if (!(c > kZero)) return r;                             // delays went the wrong way
    r.sound_speed_mps = c;
    r.fixed_latency_s = t1 - d1 / c;
    r.valid = true;
    return r;
}

std::size_t estimate_impulse_response(const FftView& fft,
                                      std::span<const Real> probe,
                                      std::span<const Real> recording,
                                      Real epsilon,
                                      std::span<Complex> work,
                                      std::span<Real> out) noexcept {
    const std::size_t n = fft.size;
    if (!fft.valid() || work.size() < 2 * n || out.size() < n) return 0;
    if (probe.empty() || recording.empty()) return 0;
    if (probe.size() > n || recording.size() > n) return 0;
    if (!(epsilon > kZero)) return 0;

    std::span<Complex> S = work.subspan(0, n);
    std::span<Complex> R = work.subspan(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        S[i] = Complex(i < probe.size() ? probe[i] : kZero, kZero);
        R[i] = Complex(i < recording.size() ? recording[i] : kZero, kZero);
    }
    fft_forward(fft, S);
    fft_forward(fft, R);

    // Regularisation floor from the probe's own mean power. Using an absolute
    // constant instead would make the result depend on the recording's gain,
    // which is exactly the thing a channel estimate must not do.
    Real mean_power = kZero;
    for (std::size_t i = 0; i < n; ++i) mean_power += std::norm(S[i]);
    mean_power /= static_cast<Real>(n);
    const Real floor_power = epsilon * mean_power;
    if (!(floor_power > kZero)) return 0;

    for (std::size_t i = 0; i < n; ++i) {
        const Real denom = std::norm(S[i]) + floor_power;
        R[i] = std::conj(S[i]) * R[i] / denom;
    }
    fft_inverse(fft, R);
    for (std::size_t i = 0; i < n; ++i) out[i] = R[i].real();
    return n;
}

std::size_t find_arrivals(std::span<const Real> h, Real fs, Real min_separation_s,
                          Real threshold_db, std::span<Arrival> out) noexcept {
    if (h.empty() || out.empty() || !(fs > kZero)) return 0;
    const auto guard = static_cast<std::size_t>(min_separation_s * fs);
    if (guard == 0) return 0;

    Real peak = kZero;
    for (const Real v : h) peak = std::max(peak, std::fabs(v));
    if (!(peak > kZero)) return 0;
    const Real threshold = peak * std::pow(static_cast<Real>(10),
                                           threshold_db / static_cast<Real>(20));

    // Greedy: take the strongest remaining sample, then exclude its guard band.
    // O(n * out.size()) with no scratch, which is the shape this library uses
    // everywhere a "sort by strength" would otherwise need an allocation.
    std::size_t found = 0;
    static_cast<void>(threshold);
    std::size_t taken[64];
    const std::size_t cap = (out.size() < 64) ? out.size() : std::size_t{64};
    while (found < cap) {
        Real best = kZero;
        std::size_t best_i = h.size();
        for (std::size_t i = 0; i < h.size(); ++i) {
            const Real v = std::fabs(h[i]);
            if (v <= best || v < threshold) continue;
            bool blocked = false;
            for (std::size_t k = 0; k < found; ++k) {
                const std::size_t d = (i > taken[k]) ? i - taken[k] : taken[k] - i;
                if (d < guard) { blocked = true; break; }
            }
            if (blocked) continue;
            best = v;
            best_i = i;
        }
        if (best_i >= h.size()) break;
        taken[found] = best_i;
        out[found].delay_s = static_cast<Real>(best_i) / fs;
        out[found].amplitude = best;
        out[found].relative_db = static_cast<Real>(20) * std::log10(best / peak);
        ++found;
    }
    return found;
}

}  // namespace phantom::bench
