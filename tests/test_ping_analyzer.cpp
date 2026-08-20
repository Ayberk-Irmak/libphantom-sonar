// SPDX-License-Identifier: Apache-2.0
// Ping analyser verification.
//
// The detector is statistical, so it is checked statistically: measured false
// alarm rate against the CFAR design value, a classification confusion matrix,
// and -- the strongest available check -- the arrival-time estimator's variance
// against the Cramer-Rao lower bound.
//
// The CRLB test is the one worth reading. It is not a comparison against a
// recorded baseline or against another implementation; it is a comparison
// against the best variance ANY unbiased estimator could achieve for this
// waveform at this SNR. An estimator sitting near the bound is doing its job;
// one comfortably below it is a bug in the measurement.
#include "framework.hpp"

#include "phantom/ping_analyzer.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr Real        kFs  = 96000;
constexpr std::size_t kFft = 8192;

PulseBank<8, kFft>*    g_bank = nullptr;
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

// The four waveform hypotheses the bank is asked to distinguish.
const PulseSpec kCw      = make(PulseType::Cw,      12000, 12000, static_cast<Real>(0.02));
const PulseSpec kLfmUp   = make(PulseType::LfmUp,    8000, 20000, static_cast<Real>(0.02));
const PulseSpec kLfmDown = make(PulseType::LfmDown, 20000,  8000, static_cast<Real>(0.02));
const PulseSpec kHfm     = make(PulseType::Hfm,      8000, 20000, static_cast<Real>(0.02));

void build_bank() {
    static PulseBank<8, kFft> bank(kFs);
    bank.clear();
    bank.add(kCw);
    bank.add(kLfmUp);
    bank.add(kLfmDown);
    bank.add(kHfm);
    g_bank = &bank;
}

void fill_noise(pt::Rng& rng, double sigma) {
    for (Real& v : g_block) v = static_cast<Real>(sigma * rng.normal());
}

void add_pulse(const PulseSpec& spec, std::size_t at, Real amplitude) {
    std::array<Real, kFft> tmp{};
    const std::size_t n = render_real(spec, kFs, tmp);
    for (std::size_t i = 0; i < n && at + i < kFft; ++i) {
        g_block[at + i] += amplitude * tmp[i];
    }
}

DetectorConfig default_cfg() {
    DetectorConfig cfg;
    // The guard MUST clear the widest response in the bank. The CW does not
    // compress, so its correlation triangle is as wide as the pulse itself --
    // see analyzer_cfar_guard_must_clear_the_response below, which pins the
    // failure that a too-small guard produces.
    cfg.cfar_guard = 1920;
    cfg.cfar_train = 256;
    cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(1e-6));
    cfg.dead_time_s = static_cast<Real>(0.005);
    return cfg;
}

}  // namespace

PT_TEST(analyzer_bank_construction) {
    build_bank();
    PT_CHECK(g_bank->size() == 4);
    PT_CHECK(g_bank->max_replica_length() == 1920);   // 20 ms at 96 kHz
    PT_CHECK(g_bank->stride() == kFft - 1920 + 1);
    PT_CHECK(g_bank->view().valid());

    // A pulse longer than the transform cannot be prepared.
    PulseBank<2, 1024> small(kFs);
    PT_CHECK(!small.add(make(PulseType::LfmUp, 8000, 20000, static_cast<Real>(0.5))));
    PT_CHECK(small.add(make(PulseType::Cw, 12000, 12000, static_cast<Real>(0.005))));
    PT_CHECK(small.size() == 1);
    // ... and neither can an invalid spec.
    PT_CHECK(!small.add(make(PulseType::LfmUp, 20000, 8000, static_cast<Real>(0.005))));
}

PT_TEST(analyzer_finds_a_clean_pulse_at_the_right_time) {
    build_bank();
    const DetectorConfig cfg = default_cfg();
    std::array<PulseDescriptor, 8> pdw{};

    pt::Rng rng(11);
    for (std::size_t at : {std::size_t(500), std::size_t(2000), std::size_t(4000)}) {
        fill_noise(rng, 0.02);
        add_pulse(kLfmUp, at, 1);
        const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                            g_scratch.view(), pdw);
        PT_CHECK(n >= 1);
        if (n == 0) continue;
        const double expected_toa = static_cast<double>(at) / static_cast<double>(kFs);
        std::printf("       pulse at %5zu: ToA %.6f s (truth %.6f), %s, SNR %.1f dB\n",
                    at, static_cast<double>(pdw[0].toa_s), expected_toa,
                    pulse_type_name(pdw[0].type), static_cast<double>(pdw[0].snr_db));
        PT_CHECK(pdw[0].type == PulseType::LfmUp);
        PT_CHECK_NEAR(pdw[0].toa_s, expected_toa, 2.0 / static_cast<double>(kFs));
        // Amplitude estimate recovers the transmitted unit amplitude.
        PT_CHECK_NEAR(pdw[0].amplitude, 1.0, 0.1);
    }
}

PT_TEST(analyzer_classifies_the_waveform_family) {
    // Confusion matrix over the four hypotheses. CW versus the chirps is easy;
    // up versus down versus hyperbolic is the part that needs the bank to be
    // energy-normalised, or the longest replica simply always wins.
    build_bank();
    const DetectorConfig cfg = default_cfg();
    std::array<PulseDescriptor, 8> pdw{};
    pt::Rng rng(2024);

    const PulseSpec* truths[] = {&kCw, &kLfmUp, &kLfmDown, &kHfm};
    const char* names[] = {"CW", "LFM-up", "LFM-down", "HFM"};
    std::size_t correct_total = 0, trials_total = 0;

    std::printf("       %-10s %8s %8s\n", "truth", "correct", "trials");
    for (std::size_t k = 0; k < 4; ++k) {
        std::size_t correct = 0;
        const std::size_t trials = 25;
        for (std::size_t t = 0; t < trials; ++t) {
            fill_noise(rng, 0.5);
            add_pulse(*truths[k], 1500, static_cast<Real>(0.6));
            const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                                g_scratch.view(), pdw);
            if (n >= 1 && pdw[0].type == truths[k]->type) ++correct;
        }
        std::printf("       %-10s %8zu %8zu\n", names[k], correct, trials);
        correct_total += correct;
        trials_total += trials;
        PT_CHECK(correct >= 20);
    }
    std::printf("       overall %zu/%zu\n", correct_total, trials_total);
    PT_CHECK(correct_total >= trials_total - 8);
}

PT_TEST(analyzer_separates_two_pulses) {
    // Two arrivals in one block: a strong direct path and a weaker echo. This
    // is what the dead-time logic exists for -- without it a single compressed
    // pulse reports once per cell over threshold.
    build_bank();
    DetectorConfig cfg = default_cfg();
    cfg.dead_time_s = static_cast<Real>(0.003);
    std::array<PulseDescriptor, 8> pdw{};

    pt::Rng rng(77);
    fill_noise(rng, 0.05);
    add_pulse(kLfmUp, 800, 1);
    add_pulse(kLfmUp, 3400, static_cast<Real>(0.4));

    const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                        g_scratch.view(), pdw);
    std::printf("       reported %zu detections\n", n);
    for (std::size_t i = 0; i < n && i < 6; ++i) {
        std::printf("         ToA %.6f s (sample %6.0f)  amp %.3f  SNR %5.1f dB  %s\n",
                    static_cast<double>(pdw[i].toa_s),
                    static_cast<double>(pdw[i].toa_s) * static_cast<double>(kFs),
                    static_cast<double>(pdw[i].amplitude),
                    static_cast<double>(pdw[i].snr_db),
                    pulse_type_name(pdw[i].type));
    }

    // Both true arrivals must be present, correctly typed and correctly scaled.
    // The bank may also emit cross-template ghosts -- see the next test, which
    // measures them -- so this looks for the truths rather than demanding an
    // exact count.
    auto find = [&](double sample, PulseType type, double amp) {
        for (std::size_t i = 0; i < n; ++i) {
            const double s_i = static_cast<double>(pdw[i].toa_s) * static_cast<double>(kFs);
            if (std::fabs(s_i - sample) <= 2.0 && pdw[i].type == type
                && std::fabs(static_cast<double>(pdw[i].amplitude) - amp) < 0.05) {
                return true;
            }
        }
        return false;
    };
    PT_CHECK(find(800.0, PulseType::LfmUp, 1.0));
    PT_CHECK(find(3400.0, PulseType::LfmUp, 0.4));
    PT_CHECK(n >= 2);
}

PT_TEST(analyzer_reports_cross_template_ghosts) {
    // A known and unavoidable property of a matched filter bank whose templates
    // are not mutually orthogonal: a real arrival also lights up the OTHER
    // templates, at a shifted lag, and those responses are reported as extra
    // detections.
    //
    // Here the LFM-up and the HFM sweep the same 8-20 kHz band, so an LFM
    // arrival produces an HFM response about 10 dB down at a different lag.
    // Nothing is wrong: the bank faithfully reports that an HFM replica does
    // partially match an LFM signal. Removing it needs association logic across
    // detections, which is a tracker's job and is v0.3 work.
    //
    // The number is pinned here so a future change to the bank's normalisation
    // cannot quietly make it worse.
    build_bank();
    DetectorConfig cfg = default_cfg();
    cfg.dead_time_s = static_cast<Real>(0.003);
    std::array<PulseDescriptor, 8> pdw{};

    pt::Rng rng(4711);
    fill_noise(rng, 0.02);
    add_pulse(kLfmUp, 2000, 1);
    const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                        g_scratch.view(), pdw);

    // Find the true detection and the strongest ghost.
    double truth_amp = 0.0, ghost_amp = 0.0;
    std::size_t ghosts = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double s_i = static_cast<double>(pdw[i].toa_s) * static_cast<double>(kFs);
        if (std::fabs(s_i - 2000.0) <= 2.0 && pdw[i].type == PulseType::LfmUp) {
            truth_amp = static_cast<double>(pdw[i].amplitude);
        } else {
            ++ghosts;
            ghost_amp = std::max(ghost_amp, static_cast<double>(pdw[i].amplitude));
        }
    }
    PT_CHECK(truth_amp > 0.9);
    std::printf("       one LFM arrival -> %zu detections (%zu ghosts), "
                "strongest ghost %.1f dB below the truth\n",
                n, ghosts, 20.0 * std::log10(ghost_amp / truth_amp));
    if (ghosts > 0) {
        // Ghosts must stay well below the real arrival, so an amplitude gate
        // can suppress them when the dynamic range of interest allows one.
        PT_CHECK(20.0 * std::log10(ghost_amp / truth_amp) < -8.0);
    }

    // A bank of mutually dissimilar waveforms does not have the problem: with
    // the HFM removed, the LFM-down and CW responses are far too weak to cross
    // the threshold.
    static PulseBank<8, kFft> clean(kFs);
    clean.clear();
    clean.add(kCw);
    clean.add(kLfmUp);
    clean.add(kLfmDown);
    fill_noise(rng, 0.02);
    add_pulse(kLfmUp, 2000, 1);
    const std::size_t n2 = analyze_block(clean.view(), cfg, g_block, 0, g_scratch.view(), pdw);
    std::printf("       same signal, bank without the HFM -> %zu detections\n", n2);
    PT_CHECK(n2 == 1);
    PT_CHECK(pdw[0].type == PulseType::LfmUp);
}

PT_TEST(analyzer_false_alarm_rate_tracks_the_cfar_design) {
    // Noise only. The measured false alarm rate should sit in the vicinity of
    // the design Pfa. Exactness is not expected: CA-CFAR's Pfa formula assumes
    // independent cells, and a matched filter's output is correlated over its
    // mainlobe, so the realised rate is somewhat higher. What must hold is that
    // the knob works -- a tighter Pfa must produce fewer alarms.
    build_bank();
    std::array<PulseDescriptor, 64> pdw{};
    pt::Rng rng(31415);

    std::printf("       %10s %12s %14s\n", "design Pfa", "alpha", "alarms/block");
    double previous = 1e9;
    for (double pfa : {1e-3, 1e-5, 1e-7}) {
        DetectorConfig cfg = default_cfg();
        cfg.cfar_guard = 32;
        cfg.cfar_train = 256;
        cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(pfa));
        cfg.dead_time_s = static_cast<Real>(0.001);

        std::size_t alarms = 0;
        const std::size_t blocks = 30;
        for (std::size_t b = 0; b < blocks; ++b) {
            fill_noise(rng, 1.0);
            alarms += analyze_block(g_bank->view(), cfg, g_block, 0, g_scratch.view(), pdw);
        }
        const double per_block = static_cast<double>(alarms) / static_cast<double>(blocks);
        std::printf("       %10.0e %12.2f %14.2f\n", pfa,
                    static_cast<double>(cfg.threshold_alpha), per_block);
        PT_CHECK(per_block <= previous + 0.05);   // monotone in the design Pfa
        previous = per_block;
    }
    PT_CHECK(previous < 1.0);   // the tightest setting is genuinely quiet
}

PT_TEST(analyzer_toa_variance_approaches_the_cramer_rao_bound) {
    // The strongest check available: the estimator's variance against the best
    // any unbiased estimator could achieve.
    //
    // For delay estimation of a known signal in white Gaussian noise of
    // per-sample variance sigma^2, with delay measured in samples,
    //     var(tau) >= sigma^2 / F,      F = sum_n (ds/dn)^2
    // and for s[n] = A cos(phi_n) with instantaneous frequency f_n,
    //     F = (A^2/2) sum_n (2 pi f_n / fs)^2
    //
    // WHICH bound applies depends on the receiver. That distinction is the
    // whole point of this test:
    //
    //   COHERENT bound   -- frequencies measured about ZERO. Achievable only by
    //                       a receiver that tracks absolute carrier phase.
    //   ENVELOPE bound   -- frequencies measured about the CENTRE frequency,
    //                       i.e. driven by the RMS bandwidth. This is the bound
    //                       for a magnitude (envelope) detector, which is what
    //                       this analyser is.
    //
    // For an 8-20 kHz sweep the two differ by a factor of
    // f_rms / B_rms = 14.4 kHz / 3.46 kHz = 4.16 in standard deviation.
    // Comparing an envelope detector against the coherent bound would make a
    // perfectly efficient estimator look four times worse than it is.
    //
    // Carrier-coherent processing is not merely unimplemented here: it needs
    // absolute phase, which Doppler destroys, so the envelope bound is the
    // relevant one for a passive intercept receiver.
    build_bank();
    DetectorConfig cfg = default_cfg();
    cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(1e-4));

    const PulseSpec& spec = kLfmUp;
    const std::size_t at = 2000;
    const double amplitude = 1.0;
    const double sigma = 3.0;                    // input SNR = -12.6 dB

    const std::size_t l = pulse_length(spec, kFs);
    const double f_c = static_cast<double>(spec.centre_frequency_hz());

    double fisher_coherent = 0.0;
    double fisher_envelope = 0.0;
    for (std::size_t i = 0; i < l; ++i) {
        const double f = static_cast<double>(
            pulse_frequency(spec, static_cast<Real>(static_cast<double>(i) / static_cast<double>(kFs))));
        const double w0 = 2.0 * 3.14159265358979323846 * f / static_cast<double>(kFs);
        const double wc = 2.0 * 3.14159265358979323846 * (f - f_c) / static_cast<double>(kFs);
        fisher_coherent += w0 * w0;
        fisher_envelope += wc * wc;
    }
    fisher_coherent *= amplitude * amplitude / 2.0;
    fisher_envelope *= amplitude * amplitude / 2.0;

    // Cross-check the coherent Fisher information against a spectral
    // derivative of the rendered waveform: two independent routes to the same
    // number, so a slip in either is visible.
    {
        static std::array<Real, kFft> replica{};
        render_real(spec, kFs, replica);
        static FftPlan<kFft> plan;
        static std::array<Complex, kFft> deriv{};
        for (std::size_t i = 0; i < kFft; ++i) {
            deriv[i] = Complex((i < l) ? replica[i] : static_cast<Real>(0), 0);
        }
        plan.forward(deriv);
        for (std::size_t k = 0; k < kFft; ++k) {
            const double kk = (k <= kFft / 2) ? static_cast<double>(k)
                                              : static_cast<double>(k) - static_cast<double>(kFft);
            const double w = 2.0 * 3.14159265358979323846 * kk / static_cast<double>(kFft);
            deriv[k] *= Complex(0, static_cast<Real>(w));
        }
        plan.inverse(deriv);
        double spectral = 0.0;
        for (std::size_t i = 0; i < l; ++i) spectral += static_cast<double>(std::norm(deriv[i]));
        std::printf("       Fisher: analytic %.1f, spectral derivative %.1f\n",
                    fisher_coherent, spectral);
        PT_CHECK_REL(spectral, fisher_coherent, 0.02);
    }

    const double crlb_coherent_s =
        std::sqrt(sigma * sigma / fisher_coherent) / static_cast<double>(kFs);
    const double crlb_envelope_s =
        std::sqrt(sigma * sigma / fisher_envelope) / static_cast<double>(kFs);

    // --- Monte Carlo -------------------------------------------------------
    pt::Rng rng(20260806);
    std::array<PulseDescriptor, 8> pdw{};
    const std::size_t trials = 600;
    std::size_t used = 0;
    double sum = 0.0, sum_sq = 0.0;

    for (std::size_t t = 0; t < trials; ++t) {
        fill_noise(rng, sigma);
        add_pulse(spec, at, static_cast<Real>(amplitude));
        const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                            g_scratch.view(), pdw);
        if (n == 0) continue;
        const double err = static_cast<double>(pdw[0].toa_s)
                         - static_cast<double>(at) / static_cast<double>(kFs);
        // Gross outliers are detection failures, not estimator inefficiency;
        // the CRLB says nothing about them. Counted and reported, not hidden.
        if (std::fabs(err) > 20.0 / static_cast<double>(kFs)) continue;
        sum += err;
        sum_sq += err * err;
        ++used;
    }

    PT_CHECK(used > trials * 3 / 4);
    const double mean = sum / static_cast<double>(used);
    const double sd = std::sqrt(sum_sq / static_cast<double>(used) - mean * mean);

    std::printf("       input SNR %.1f dB, %zu/%zu trials used, bias %+.3f us\n",
                10.0 * std::log10(amplitude * amplitude / 2.0 / (sigma * sigma)),
                used, trials, mean * 1e6);
    std::printf("       measured sigma_toa   = %.3f us\n", sd * 1e6);
    std::printf("       coherent CRLB        = %.3f us   (ratio %.2f -- not the applicable bound)\n",
                crlb_coherent_s * 1e6, sd / crlb_coherent_s);
    std::printf("       envelope CRLB        = %.3f us   (ratio %.2f)\n",
                crlb_envelope_s * 1e6, sd / crlb_envelope_s);

    PT_CHECK(std::fabs(mean) < 1.5 * sd);            // unbiased at its own scale
    PT_CHECK(sd > 0.8 * crlb_envelope_s);            // cannot beat its own bound
    PT_CHECK(sd < 1.35 * crlb_envelope_s);           // and is essentially efficient
    // The two bounds must be far enough apart for this to be a real claim.
    PT_CHECK(crlb_envelope_s > 3.0 * crlb_coherent_s);
}

PT_TEST(analyzer_estimator_stays_efficient_across_snr) {
    // Efficiency must hold over a range, not at one lucky operating point. A
    // constant ratio to the bound across 20 dB is the signature of an estimator
    // limited by noise rather than by a systematic interpolation error.
    build_bank();
    DetectorConfig cfg = default_cfg();
    cfg.threshold_alpha = cfar_alpha(512, static_cast<Real>(1e-4));

    const PulseSpec& spec = kLfmUp;
    const std::size_t at = 2000;
    const std::size_t l = pulse_length(spec, kFs);
    const double f_c = static_cast<double>(spec.centre_frequency_hz());

    double fisher_envelope = 0.0;
    for (std::size_t i = 0; i < l; ++i) {
        const double f = static_cast<double>(
            pulse_frequency(spec, static_cast<Real>(static_cast<double>(i) / static_cast<double>(kFs))));
        const double wc = 2.0 * 3.14159265358979323846 * (f - f_c) / static_cast<double>(kFs);
        fisher_envelope += wc * wc;
    }
    fisher_envelope *= 0.5;

    pt::Rng rng(99991);
    std::array<PulseDescriptor, 8> pdw{};
    std::printf("       %8s %12s %12s %8s\n", "sigma", "measured us", "env CRLB", "ratio");
    for (double sigma : {3.0, 1.0, 0.3}) {
        std::size_t used = 0;
        double sum = 0.0, sum_sq = 0.0;
        for (std::size_t t = 0; t < 300; ++t) {
            fill_noise(rng, sigma);
            add_pulse(spec, at, 1);
            const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, 0,
                                                g_scratch.view(), pdw);
            if (n == 0) continue;
            const double err = static_cast<double>(pdw[0].toa_s)
                             - static_cast<double>(at) / static_cast<double>(kFs);
            if (std::fabs(err) > 20.0 / static_cast<double>(kFs)) continue;
            sum += err;
            sum_sq += err * err;
            ++used;
        }
        if (used < 100) { PT_CHECK(false); continue; }
        const double mean = sum / static_cast<double>(used);
        const double sd = std::sqrt(sum_sq / static_cast<double>(used) - mean * mean);
        const double crlb = std::sqrt(sigma * sigma / fisher_envelope) / static_cast<double>(kFs);
        std::printf("       %8.2f %12.4f %12.4f %8.2f\n", sigma, sd * 1e6, crlb * 1e6, sd / crlb);
        PT_CHECK(sd > 0.8 * crlb);
        PT_CHECK(sd < 1.35 * crlb);
    }
}

PT_TEST(analyzer_cfar_guard_must_clear_the_response) {
    // Pins the failure mode suggested_cfar_guard() exists to prevent. With a
    // guard narrower than the CW's correlation triangle, the pulse leaks into
    // its own training cells, lifts its own threshold, and vanishes -- silently.
    build_bank();
    std::array<PulseDescriptor, 8> pdw{};
    pt::Rng rng(4242);

    PT_CHECK(suggested_cfar_guard(g_bank->view()) == g_bank->max_replica_length());

    DetectorConfig narrow = default_cfg();
    narrow.cfar_guard = 32;                       // far smaller than the 1920-sample response
    DetectorConfig wide = default_cfg();
    wide.cfar_guard = suggested_cfar_guard(g_bank->view());

    std::size_t found_narrow = 0, found_wide = 0;
    for (std::size_t t = 0; t < 20; ++t) {
        fill_noise(rng, 0.5);
        add_pulse(kCw, 1500, static_cast<Real>(0.6));
        static std::array<Real, kFft> saved{};
        saved = g_block;
        if (analyze_block(g_bank->view(), narrow, g_block, 0, g_scratch.view(), pdw) > 0) ++found_narrow;
        g_block = saved;
        if (analyze_block(g_bank->view(), wide, g_block, 0, g_scratch.view(), pdw) > 0) ++found_wide;
    }
    std::printf("       CW detected: %zu/20 with a 32-cell guard, %zu/20 with %zu\n",
                found_narrow, found_wide, wide.cfar_guard);
    PT_CHECK(found_narrow == 0);   // the silent failure, reproduced
    PT_CHECK(found_wide >= 19);    // and fixed by the suggested guard

    // The chirps compress to ~fs/B cells, so they survive the narrow guard --
    // which is exactly why this bug hides: three of four waveforms still work.
    std::size_t chirp_narrow = 0;
    for (std::size_t t = 0; t < 20; ++t) {
        fill_noise(rng, 0.5);
        add_pulse(kLfmUp, 1500, static_cast<Real>(0.6));
        if (analyze_block(g_bank->view(), narrow, g_block, 0, g_scratch.view(), pdw) > 0) ++chirp_narrow;
    }
    std::printf("       LFM detected: %zu/20 with the same 32-cell guard\n", chirp_narrow);
    PT_CHECK(chirp_narrow >= 19);
}

PT_TEST(analyzer_rejects_bad_inputs) {
    build_bank();
    const DetectorConfig cfg = default_cfg();
    std::array<PulseDescriptor, 4> pdw{};

    // Wrong block length.
    std::array<Real, 128> short_block{};
    PT_CHECK(analyze_block(g_bank->view(), cfg, short_block, 0, g_scratch.view(), pdw) == 0);

    // Empty output span.
    PT_CHECK(analyze_block(g_bank->view(), cfg, g_block, 0, g_scratch.view(),
                           std::span<PulseDescriptor>()) == 0);

    // Empty bank.
    PulseBank<4, kFft> empty(kFs);
    PT_CHECK(!empty.view().valid());
    PT_CHECK(analyze_block(empty.view(), cfg, g_block, 0, g_scratch.view(), pdw) == 0);

    // Silence produces nothing, and does not divide by a zero noise estimate.
    for (Real& v : g_block) v = 0;
    PT_CHECK(analyze_block(g_bank->view(), cfg, g_block, 0, g_scratch.view(), pdw) == 0);

    PT_CHECK(cfar_alpha(0, static_cast<Real>(1e-6)) == 0);
    PT_CHECK(cfar_alpha(512, 0) == 0);
    PT_CHECK(cfar_alpha(512, 2) == 0);
    PT_CHECK(cfar_alpha(512, static_cast<Real>(1e-6)) > 0);
}

PT_TEST(analyzer_timestamps_are_block_relative) {
    // Streaming: the caller advances by stride() and passes the timestamp of
    // block[0]. The reported ToA must be absolute, not block-relative.
    build_bank();
    const DetectorConfig cfg = default_cfg();
    std::array<PulseDescriptor, 4> pdw{};

    pt::Rng rng(5);
    fill_noise(rng, 0.02);
    add_pulse(kHfm, 1000, 1);

    const Real t0 = static_cast<Real>(12.5);
    const std::size_t n = analyze_block(g_bank->view(), cfg, g_block, t0, g_scratch.view(), pdw);
    PT_CHECK(n >= 1);
    if (n >= 1) {
        PT_CHECK_NEAR(pdw[0].toa_s, 12.5 + 1000.0 / static_cast<double>(kFs),
                      2.0 / static_cast<double>(kFs));
        PT_CHECK(pdw[0].type == PulseType::Hfm);
        PT_CHECK_NEAR(pdw[0].bandwidth_hz, 12000.0, 1.0);
        PT_CHECK_NEAR(pdw[0].duration_s, 0.02, 1e-9);
    }
}
