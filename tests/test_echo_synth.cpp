// SPDX-License-Identifier: Apache-2.0
// Doppler bank and echo synthesis verification.
//
// The echo tests are closed-loop: a ping is detected, an echo is synthesised
// from the resulting descriptor, and the echo is fed back through the analyser.
// If the delay, amplitude or Doppler that comes back out does not match what
// went in, one of the two halves is wrong -- and because they were written
// against the same physics, a shared sign error shows up as a round trip that
// does not close rather than as two tests that both pass.
#include "framework.hpp"

#include "phantom/echo_synth.hpp"
#include "phantom/ping_analyzer.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr Real        kFs  = 96000;
constexpr std::size_t kFft = 8192;
constexpr Real        kC   = 1500;

AnalyzerScratch<kFft>  g_scratch;
std::array<Real, kFft> g_block;

PulseSpec make(PulseType t, Real f0, Real f1, Real dur) {
    PulseSpec s;
    s.type = t;
    s.f_start_hz = f0;
    s.f_end_hz = f1;
    s.duration_s = dur;
    return s;
}

const PulseSpec kLfm = make(PulseType::LfmUp, 8000, 20000, static_cast<Real>(0.02));
const PulseSpec kHfm = make(PulseType::Hfm,   8000, 20000, static_cast<Real>(0.02));
const PulseSpec kCw  = make(PulseType::Cw,   12000, 12000, static_cast<Real>(0.02));

DetectorConfig cfg_for(const AnalyzerView& v) {
    DetectorConfig cfg;
    cfg.cfar_guard = suggested_cfar_guard(v);
    cfg.cfar_train = 256;
    cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(1e-6));
    cfg.dead_time_s = suggested_dead_time_s(v);
    return cfg;
}

void place(std::span<const Real> sig, std::size_t at, Real gain = 1) {
    for (std::size_t i = 0; i < sig.size() && at + i < kFft; ++i) {
        g_block[at + i] += gain * sig[i];
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Doppler bank
// ---------------------------------------------------------------------------

PT_TEST(doppler_tolerance_matches_the_measured_loss) {
    // The bin-spacing formula is only useful if it predicts the real loss in
    // the regime a bank is designed for. Checked directly against the matched
    // filter: at the tolerance the formula reports, the peak must have dropped
    // by about the requested amount.
    static FftPlan<kFft> plan;
    static std::array<Complex, kFft> rep{}, spec{};
    static std::array<Complex, kFft> out{};
    static std::array<Real, kFft> tmp{};

    auto measured_loss_db = [&](const PulseSpec& s, Real delta) {
        const std::size_t l = render_analytic(s, kFs, rep);
        MatchedFilter mf;
        matched_filter_prepare(plan.view(), std::span<const Complex>(rep.data(), l), spec, mf);

        // Placed well inside the block on purpose. An HFM's Doppler response is
        // a large DELAY -- tau = (alpha-1)/(alpha k f0), which is -3.3 ms at
        // delta = 0.109 -- so a pulse near the start of the block would have its
        // peak pushed to a negative lag and clipped, and the "loss" measured
        // would be the test's framing rather than the waveform's.
        auto peak = [&](Real d) {
            g_block.fill(0);
            const std::size_t n = render_real_doppler(s, kFs, d, tmp);
            place(std::span<const Real>(tmp.data(), n), 2500);
            const std::size_t lags = matched_filter_apply(plan.view(), mf, g_block,
                                                          g_scratch.view().fft_scratch, out);
            double best = 0;
            for (std::size_t i = 0; i < lags; ++i) {
                best = std::max(best, static_cast<double>(std::abs(out[i])));
            }
            return best;
        };
        return 20.0 * std::log10(peak(delta) / peak(0));
    };

    std::printf("       %-8s %6s %14s %12s %12s\n",
                "waveform", "L (dB)", "delta", "v (m/s)", "measured dB");
    for (const auto& [name, s] : {std::pair<const char*, PulseSpec>{"LFM", kLfm},
                                  {"HFM", kHfm}, {"CW", kCw}}) {
        for (double want : {1.0, 3.0}) {
            const Real tol = doppler_tolerance(s, static_cast<Real>(want));
            const double got = -measured_loss_db(s, tol);
            std::printf("       %-8s %6.1f %14.3g %12.2f %12.2f\n",
                        name, want, static_cast<double>(tol),
                        static_cast<double>(tol) * 1500.0, got);
            // Within a factor of two of the requested loss is enough to space
            // bins sensibly; the formula is a second-order expansion, not an
            // identity.
            PT_CHECK(got > 0.4 * want);
            PT_CHECK(got < 2.5 * want);
        }
    }
}

PT_TEST(doppler_bank_sizing_reflects_the_waveform) {
    // The headline consequence: an HFM needs far fewer Doppler bins than an
    // LFM of the same time-bandwidth product, because time scaling maps the
    // HFM family onto itself. Both sweep 8-20 kHz for 20 ms here.
    // static, not automatic: a 64-template bank at this transform size is 8.4 MB
    // of spectra and would blow the stack. See the note on PulseBank.
    static PulseBank<64, kFft> bank(kFs, kC);
    const std::size_t lfm_bins = bank.doppler_bins_required(kLfm, -20, 20, 1);
    const std::size_t hfm_bins = bank.doppler_bins_required(kHfm, -20, 20, 1);
    const std::size_t cw_bins  = bank.doppler_bins_required(kCw,  -20, 20, 1);

    std::printf("       +/- 20 m/s at 1 dB straddling loss:\n");
    std::printf("         LFM %zu bins, HFM %zu bins, CW %zu bins\n",
                lfm_bins, hfm_bins, cw_bins);
    PT_CHECK(lfm_bins >= 4);
    // Two bins for the HFM is just the two endpoints: one interval covers the
    // whole 40 m/s span, where the LFM needs several and the CW needs an order
    // of magnitude more. That ratio is the design argument for HFM.
    PT_CHECK(hfm_bins == 2);
    PT_CHECK(2 * hfm_bins <= lfm_bins);
    PT_CHECK(cw_bins > 2 * lfm_bins);
}

PT_TEST(doppler_bank_recovers_radial_velocity) {
    // A moving transmitter must be both detected AND have its speed reported,
    // to within the bin spacing.
    static PulseBank<48, kFft> bank(kFs, kC);
    bank.clear();
    const std::size_t added = bank.add_doppler_bank(kLfm, -12, 12, 1);
    const double spacing = static_cast<double>(bank.doppler_bin_spacing_mps(kLfm, -12, 12, 1));
    std::printf("       LFM bank: %zu bins over +/- 12 m/s, spacing %.2f m/s\n",
                added, spacing);
    PT_CHECK(added >= 3);
    PT_CHECK(spacing > 0);

    const AnalyzerView v = bank.view();
    const DetectorConfig cfg = cfg_for(v);
    std::array<PulseDescriptor, 8> pdw{};
    static std::array<Real, kFft> tmp{};
    pt::Rng rng(606);

    double worst = 0.0;
    for (double v_true : {-10.0, -5.0, 0.0, 5.0, 10.0}) {
        g_block.fill(0);
        for (Real& x : g_block) x = static_cast<Real>(0.05 * rng.normal());
        const std::size_t n = render_real_doppler(kLfm, kFs,
                                                  static_cast<Real>(v_true / 1500.0), tmp);
        place(std::span<const Real>(tmp.data(), n), 1500);

        const std::size_t m = analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw);
        PT_CHECK(m >= 1);
        if (m == 0) continue;
        const double est = static_cast<double>(pdw[0].radial_velocity_mps);
        std::printf("       true %+6.1f m/s -> reported %+6.1f m/s%s\n",
                    v_true, est, pdw[0].doppler_resolved ? "" : "  (unresolved)");
        PT_CHECK(pdw[0].doppler_resolved);
        worst = std::max(worst, std::fabs(est - v_true));
    }
    // A bin-quantised estimate cannot do better than half the bin spacing;
    // the bound is derived from the bank rather than hardcoded, so it stays
    // honest if the tolerance formula ever changes.
    std::printf("       worst velocity error %.2f m/s (half-bin is %.2f)\n",
                worst, spacing / 2.0);
    PT_CHECK(worst <= spacing / 2.0 + 0.51);
}

PT_TEST(doppler_bank_recovers_the_peak_a_zero_doppler_bank_loses) {
    // The point of the bank. A 10 m/s target against a zero-Doppler LFM
    // template loses several dB; the matching bin gets it back.
    // 0.3 dB design, so the bins are dense enough that a 5 m/s target has a
    // non-zero bin nearer than the zero-Doppler one. At the 1 dB spacing it
    // would not, and the bank could not help -- which is the real lesson:
    // coverage is only as good as the bin density you paid for.
    static PulseBank<48, kFft> matched(kFs, kC);
    matched.clear();
    const std::size_t nbins = matched.add_doppler_bank(kLfm, -12, 12,
                                                       static_cast<Real>(0.3));
    std::printf("       fine bank: %zu bins, spacing %.2f m/s\n", nbins,
                static_cast<double>(matched.doppler_bin_spacing_mps(
                    kLfm, -12, 12, static_cast<Real>(0.3))));

    static PulseBank<4, kFft> plain(kFs, kC);
    plain.clear();
    plain.add(kLfm);

    static std::array<Real, kFft> tmp{};
    std::array<PulseDescriptor, 8> pdw{};

    auto peak_amplitude = [&](const AnalyzerView& v, double v_true) {
        g_block.fill(0);
        const std::size_t n = render_real_doppler(kLfm, kFs,
                                                  static_cast<Real>(v_true / 1500.0), tmp);
        place(std::span<const Real>(tmp.data(), n), 1500);
        const DetectorConfig cfg = cfg_for(v);
        const std::size_t m = analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw);
        return (m > 0) ? static_cast<double>(pdw[0].peak_magnitude) : 0.0;
    };

    const double ref = peak_amplitude(plain.view(), 0.0);
    std::printf("       %8s %14s %14s\n", "v (m/s)", "no bank (dB)", "with bank (dB)");
    for (double v_true : {5.0, 10.0}) {
        const double without = 20.0 * std::log10(peak_amplitude(plain.view(), v_true) / ref);
        const double with = 20.0 * std::log10(peak_amplitude(matched.view(), v_true) / ref);
        std::printf("       %8.1f %14.2f %14.2f\n", v_true, without, with);
        // The bank must recover most of whatever was lost. Demanding a fixed
        // number of dB would fail at speeds where there was less than that to
        // recover in the first place.
        const double lost = -without;
        const double recovered = with - without;
        PT_CHECK(recovered >= 0.5 * lost - 0.05);
        PT_CHECK(with > -1.5);            // and lands close to the matched peak
    }
}

// ---------------------------------------------------------------------------
// Echo synthesis
// ---------------------------------------------------------------------------

PT_TEST(echo_delay_is_two_way) {
    // dt = 2 dr / c. A ghost 150 m further out is 200 ms later, not 100 ms.
    PT_CHECK_NEAR(echo_delay_s(150, 1500), 0.2, pt::tol(1e-12, 1e-8));
    PT_CHECK_NEAR(echo_delay_s(0, 1500), 0.0, 1e-12);
    PT_CHECK_NEAR(echo_delay_s(-75, 1500), -0.1, pt::tol(1e-12, 1e-8));
    PT_CHECK_NEAR(echo_delay_s(150, 0), 0.0, 1e-12);   // no sound speed, no answer
}

PT_TEST(echo_doppler_is_two_way_and_uses_the_exact_form) {
    // alpha = (c+v)/(c-v), which is ~1 + 2v/c but not equal to it. At 30 m/s
    // the difference is 0.4% -- four samples of delay across a 20 ms pulse, so
    // enough to matter to a matched filter.
    const double c = 1500.0;
    for (double v : {5.0, 15.0, 30.0}) {
        const double exact = static_cast<double>(echo_doppler_scale(static_cast<Real>(v), 1500));
        const double linear = 1.0 + 2.0 * v / c;
        std::printf("       v = %5.1f m/s: exact %.6f, linear approx %.6f (%.3f%% apart)\n",
                    v, exact, linear, 100.0 * (exact - linear) / (exact - 1.0));
        PT_CHECK_REL(exact, (c + v) / (c - v), pt::tol(1e-12, 1e-7));
        PT_CHECK(exact > linear);   // the exact form is always the larger
    }
    PT_CHECK_NEAR(echo_doppler_scale(0, 1500), 1.0, pt::tol(1e-12, 1e-7));
    // Beyond the sound speed the expression is meaningless; it must not
    // return a negative or infinite scale.
    PT_CHECK_NEAR(echo_doppler_scale(2000, 1500), 1.0, pt::tol(1e-12, 1e-7));
}

PT_TEST(echo_round_trip_closes_on_delay_and_amplitude) {
    // Closed loop: detect a ping, synthesise a ghost from the descriptor, feed
    // it back, and check the analyser reports the range offset and target
    // strength that were asked for.
    static PulseBank<8, kFft> bank(kFs, kC);
    bank.clear();
    bank.add(kCw);
    bank.add(kLfm);
    bank.add(kHfm);
    const AnalyzerView v = bank.view();
    const DetectorConfig cfg = cfg_for(v);

    std::array<PulseDescriptor, 8> pdw{};
    static std::array<Real, kFft> tmp{};

    // --- intercept -----------------------------------------------------
    g_block.fill(0);
    const std::size_t n = render_real(kLfm, kFs, tmp);
    place(std::span<const Real>(tmp.data(), n), 400, static_cast<Real>(0.8));
    PT_CHECK(analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw) >= 1);
    const PulseDescriptor ping = pdw[0];
    PT_CHECK(ping.type == PulseType::LfmUp);
    PT_CHECK_NEAR(ping.amplitude, 0.8, 0.05);

    // --- synthesise and re-detect ---------------------------------------
    std::printf("       %10s %10s %14s %12s\n",
                "range (m)", "TS (dB)", "delay (ms)", "amp ratio");
    for (double range_m : {20.0, 45.0}) {
        for (double ts_db : {0.0, -6.0}) {
            EchoSpec e;
            e.range_offset_m = static_cast<Real>(range_m);
            e.target_strength_db = static_cast<Real>(ts_db);

            static std::array<Real, kFft> echo{};
            const std::size_t m = synthesize_echo(ping, e, kFs, kC, echo);
            PT_CHECK(m > 0);

            // A trace of noise, because a perfectly silent block gives CFAR a
            // degenerate noise estimate and no real hydrophone ever does.
            static pt::Rng rng(9090);
            for (Real& x : g_block) x = static_cast<Real>(0.01 * rng.normal());
            const auto at = static_cast<std::size_t>(
                (static_cast<double>(ping.toa_s) + 2.0 * range_m / 1500.0) * 96000.0);
            place(std::span<const Real>(echo.data(), m), at);

            const std::size_t k = analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw);
            PT_CHECK(k >= 1);
            if (k == 0) continue;

            const double delay_ms = (static_cast<double>(pdw[0].toa_s)
                                   - static_cast<double>(ping.toa_s)) * 1e3;
            const double ratio = static_cast<double>(pdw[0].amplitude)
                               / static_cast<double>(ping.amplitude);
            std::printf("       %10.1f %10.1f %14.3f %12.3f\n",
                        range_m, ts_db, delay_ms, ratio);

            PT_CHECK_NEAR(delay_ms, 2.0 * range_m / 1500.0 * 1e3, 0.05);
            PT_CHECK_NEAR(ratio, std::pow(10.0, ts_db / 20.0), 0.06);
            PT_CHECK(pdw[0].type == PulseType::LfmUp);
        }
    }
}

PT_TEST(echo_round_trip_closes_on_doppler) {
    // The same loop, but for velocity: impose a radial velocity on the ghost
    // and check a Doppler bank reads it back. The echo carries the TWO-way
    // scale, so a bank built for one-way transmitter motion must report
    // approximately 2v -- and that factor of two is exactly the kind of thing
    // this round trip exists to catch.
    static PulseBank<48, kFft> bank(kFs, kC);
    bank.clear();
    bank.add_doppler_bank(kHfm, -1, 1, 1);        // HFM: one bin covers it
    const Real design = static_cast<Real>(0.3);
    const std::size_t lfm_bins = bank.add_doppler_bank(kLfm, -30, 30, design);
    const double spacing = static_cast<double>(
        bank.doppler_bin_spacing_mps(kLfm, -30, 30, design));
    std::printf("       LFM Doppler bank: %zu bins, spacing %.2f m/s\n", lfm_bins, spacing);
    PT_CHECK(lfm_bins >= 8);

    const AnalyzerView v = bank.view();
    const DetectorConfig cfg = cfg_for(v);
    std::array<PulseDescriptor, 8> pdw{};
    static std::array<Real, kFft> tmp{}, echo{};

    g_block.fill(0);
    const std::size_t n = render_real(kLfm, kFs, tmp);
    place(std::span<const Real>(tmp.data(), n), 300);
    PT_CHECK(analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw) >= 1);
    const PulseDescriptor ping = pdw[0];

    std::printf("       %12s %16s %16s\n",
                "ghost v(m/s)", "two-way scale", "bank reads (m/s)");
    for (double v_ghost : {5.0, 10.0}) {
        EchoSpec e;
        e.radial_velocity_mps = static_cast<Real>(v_ghost);
        const std::size_t m = synthesize_echo(ping, e, kFs, kC, echo);
        PT_CHECK(m > 0);

        g_block.fill(0);
        place(std::span<const Real>(echo.data(), m), 1200);
        const std::size_t k = analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw);
        PT_CHECK(k >= 1);
        if (k == 0) continue;

        const double alpha = static_cast<double>(echo_doppler_scale(static_cast<Real>(v_ghost), 1500));
        const double reported = static_cast<double>(pdw[0].radial_velocity_mps);
        std::printf("       %12.1f %16.6f %16.1f\n", v_ghost, alpha, reported);

        // The bank's bins are labelled in one-way velocity, so a two-way echo
        // scale of alpha reads back as (alpha-1)*c, i.e. about 2*v_ghost.
        // The discriminating claim: the reading matches the TWO-way scale, not
        // the one-way velocity of the ghost. A missing factor of two here is
        // the single most likely bug in the whole echo path.
        PT_CHECK_NEAR(reported, (alpha - 1.0) * 1500.0, spacing / 2.0 + 0.51);
        PT_CHECK(std::fabs(reported - (alpha - 1.0) * 1500.0)
                 < std::fabs(reported - v_ghost));
    }
}

PT_TEST(echo_ghost_swarm_places_every_target) {
    // Several ghosts from one intercepted ping, each at its own range and
    // strength, all detected at the right places.
    static PulseBank<8, kFft> bank(kFs, kC);
    bank.clear();
    bank.add(kCw);
    bank.add(kLfm);
    bank.add(kHfm);
    const AnalyzerView v = bank.view();
    DetectorConfig cfg = cfg_for(v);
    cfg.dead_time_s = static_cast<Real>(0.004);

    std::array<PulseDescriptor, 16> pdw{};
    static std::array<Real, kFft> tmp{};

    g_block.fill(0);
    const std::size_t n = render_real(kLfm, kFs, tmp);
    place(std::span<const Real>(tmp.data(), n), 200);
    PT_CHECK(analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw) >= 1);
    const PulseDescriptor ping = pdw[0];

    // Designated initialisers on purpose. EchoSpec gained `extra_delay_s` as its
    // second member in v0.4, and every positional initialiser silently became a
    // different echo -- a -8 dB target strength was read as a -8 SECOND delay,
    // and the ghost was dropped for arriving before the ping. Naming the fields
    // makes the next insertion a no-op instead of a puzzle.
    const EchoSpec swarm[] = {
        {.range_offset_m = 10, .target_strength_db = 0},
        {.range_offset_m = 25, .target_strength_db = -4},
        {.range_offset_m = 45, .target_strength_db = -8},
    };
    const std::span<const EchoSpec> echoes(swarm, 3);

    static std::array<Real, kFft> out{};
    const std::size_t need = swarm_length(ping, echoes, kFs, kC);
    std::printf("       swarm needs %zu samples (%.1f ms)\n",
                need, static_cast<double>(need) / 96000.0 * 1e3);
    PT_CHECK(need > 0);
    PT_CHECK(need < kFft);
    const std::size_t written = synthesize_swarm(ping, echoes, kFs, kC, out);
    PT_CHECK(written == need);

    g_block.fill(0);
    place(std::span<const Real>(out.data(), written), 100);
    const std::size_t k = analyze_block(v, cfg, g_block, 0, g_scratch.view(), pdw);
    std::printf("       %zu ghosts transmitted, %zu detected\n", echoes.size(), k);
    PT_CHECK(k >= 3);

    for (std::size_t g = 0; g < 3; ++g) {
        const double want_ms = 2.0 * static_cast<double>(swarm[g].range_offset_m) / 1500.0 * 1e3;
        bool found = false;
        for (std::size_t i = 0; i < k; ++i) {
            const double got_ms = (static_cast<double>(pdw[i].toa_s) - 100.0 / 96000.0) * 1e3;
            if (std::fabs(got_ms - want_ms) < 0.1) {
                std::printf("         range %5.0f m -> %7.3f ms (want %7.3f), amp %.3f\n",
                            static_cast<double>(swarm[g].range_offset_m), got_ms, want_ms,
                            static_cast<double>(pdw[i].amplitude));
                found = true;
                break;
            }
        }
        PT_CHECK(found);
    }

    // A swarm that does not fit is a refusal, not an overrun.
    std::array<Real, 64> tiny{};
    PT_CHECK(synthesize_swarm(ping, echoes, kFs, kC, tiny) == 0);
}

PT_TEST(echo_length_scale_stretches_the_return) {
    // An extended target smears the return, because bow and stern are at
    // different ranges. Doubling the length must double the detected duration
    // -- and must NOT change the arrival time of the leading edge.
    static PulseBank<8, kFft> bank(kFs, kC);
    bank.clear();
    bank.add(kLfm);
    const AnalyzerView v = bank.view();

    std::array<PulseDescriptor, 8> pdw{};
    static std::array<Real, kFft> tmp{}, echo{};

    g_block.fill(0);
    const std::size_t n = render_real(kLfm, kFs, tmp);
    place(std::span<const Real>(tmp.data(), n), 200);
    PT_CHECK(analyze_block(v, cfg_for(v), g_block, 0, g_scratch.view(), pdw) >= 1);
    const PulseDescriptor ping = pdw[0];

    EchoSpec e;
    const std::size_t m1 = synthesize_echo(ping, e, kFs, kC, echo);
    e.length_scale = 2;
    const std::size_t m2 = synthesize_echo(ping, e, kFs, kC, echo);
    std::printf("       length scale 1 -> %zu samples, 2 -> %zu samples\n", m1, m2);
    PT_CHECK(m1 > 0);
    PT_CHECK_REL(static_cast<double>(m2), 2.0 * static_cast<double>(m1), 0.01);

    e.length_scale = 0;                 // nonsensical: refuse
    PT_CHECK(synthesize_echo(ping, e, kFs, kC, echo) == 0);
    e.length_scale = 1;
    PT_CHECK(synthesize_echo(ping, e, 0, kC, echo) == 0);      // no sample rate
    PT_CHECK(synthesize_echo(ping, e, kFs, 0, echo) == 0);     // no sound speed
}

// ---------------------------------------------------------------------------
// Anti-phase cancellation: the honest limits
// ---------------------------------------------------------------------------

PT_TEST(anti_phase_null_requires_millimetre_timing) {
    // The claim that turns up in every countermeasure specification -- transmit
    // the inverse and cancel your own echo -- is quantified here rather than
    // argued about.
    //
    // Residual power for e(t) - e(t - tau) over a flat band:
    //     G^2 = 2 - 2 cos(2 pi f_c tau) sinc(B tau)
    // At tau = 0 the null is perfect. The question is how fast it degrades.
    const Real fc = 12000;
    const Real bw = 12000;

    PT_CHECK(null_gain_db(bw, fc, 0) < -100);   // perfect timing, perfect null

    std::printf("       12 kHz centre, 12 kHz band:\n");
    std::printf("       %14s %14s %16s\n", "timing (us)", "path (mm)", "residual (dB)");
    double previous = -1000.0;
    for (double tau_us : {0.1, 0.5, 1.0, 5.0, 10.0, 20.8}) {
        const Real tau = static_cast<Real>(tau_us * 1e-6);
        const double g = static_cast<double>(null_gain_db(bw, fc, tau));
        std::printf("       %14.1f %14.2f %16.2f\n", tau_us, tau_us * 1e-6 * 1500.0 * 1e3, g);
        PT_CHECK(g > previous);        // monotone worsening over the first lobe
        previous = g;
    }

    // 20 dB of cancellation needs about a microsecond, i.e. 1.5 mm of path
    // error at 1500 m/s. That is the number that kills the idea in practice.
    Real tau_20db = 0;
    for (int i = 1; i <= 2000; ++i) {
        const Real tau = static_cast<Real>(i) * static_cast<Real>(1e-8);
        if (null_gain_db(bw, fc, tau) > static_cast<Real>(-20)) break;
        tau_20db = tau;
    }
    std::printf("       -20 dB needs tau < %.2f us = %.2f mm of path error\n",
                static_cast<double>(tau_20db) * 1e6,
                static_cast<double>(tau_20db) * 1500.0 * 1e3);
    PT_CHECK(tau_20db > 0);
    PT_CHECK(tau_20db < static_cast<Real>(3e-6));

    // Past a quarter period the "canceller" is adding energy, not removing it:
    // it makes the vehicle louder. A naive implementation ships exactly this.
    const double quarter = 1.0 / (4.0 * 12000.0);
    const double at_quarter = static_cast<double>(
        null_gain_db(bw, fc, static_cast<Real>(quarter)));
    std::printf("       at a quarter period (%.1f us): %+.2f dB -- louder, not quieter\n",
                quarter * 1e6, at_quarter);
    PT_CHECK(at_quarter > 0.0);
}

PT_TEST(anti_phase_null_is_a_timing_problem_not_a_bandwidth_one) {
    // The bandwidth term is real but negligible for any sonar with f_c >> B,
    // which is all of them. Pinned here so the header's claim is measured
    // rather than asserted.
    const Real fc = 12000;
    const Real tau = static_cast<Real>(2e-6);   // 3 mm of path error

    std::printf("       2 us (3 mm) timing error at 12 kHz centre:\n");
    std::printf("       %14s %16s\n", "bandwidth (Hz)", "residual (dB)");
    const double at_zero = static_cast<double>(null_gain_db(0, fc, tau));
    double at_wide = 0;
    for (double b : {0.0, 1000.0, 6000.0, 12000.0}) {
        at_wide = static_cast<double>(null_gain_db(static_cast<Real>(b), fc, tau));
        std::printf("       %14.0f %16.2f\n", b, at_wide);
    }
    // Widening from nothing to a full-octave sweep changes the answer by less
    // than a dB: timing dominates.
    PT_CHECK(std::fabs(at_wide - at_zero) < 1.0);

    // The closed-form timing budget, checked against the numerical model.
    std::printf("       %12s %14s %14s\n", "target (dB)", "tau max (us)", "path (mm)");
    for (double target : {-6.0, -10.0, -20.0, -30.0}) {
        const Real t = max_timing_error_s(static_cast<Real>(target), fc);
        std::printf("       %12.0f %14.3f %14.3f\n",
                    target, static_cast<double>(t) * 1e6,
                    static_cast<double>(t) * 1500.0 * 1e3);
        PT_CHECK(t > 0);
        // At exactly that error the narrowband residual must equal the target.
        PT_CHECK_NEAR(null_gain_db(0, fc, t), target, 0.05);
        // A little more error must miss it.
        PT_CHECK(static_cast<double>(null_gain_db(0, fc, t * static_cast<Real>(1.2))) > target);
    }

    // 20 dB of cancellation costs under two millimetres of path accuracy.
    const Real t20 = max_timing_error_s(static_cast<Real>(-20), fc);
    PT_CHECK(static_cast<double>(t20) * 1500.0 * 1e3 < 2.0);

    // Rejections.
    PT_CHECK(max_timing_error_s(static_cast<Real>(6), fc) == 0);   // positive target
    PT_CHECK(max_timing_error_s(static_cast<Real>(-20), 0) == 0);  // no centre frequency
    PT_CHECK(null_gain_db(12000, 0, tau) == 0);
}
