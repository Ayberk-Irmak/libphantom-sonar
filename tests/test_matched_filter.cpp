// SPDX-License-Identifier: Apache-2.0
// Matched filter verification.
//
// The FFT path is checked against a direct O(N*L) correlation that shares no
// code with it, and the detection performance is checked against the closed
// forms that define what a matched filter is supposed to do: peak height equals
// half the replica energy, compressed width is 1/B, coherent gain is L/2.
//
// The Doppler tests matter more here than they would in radar. With c ~ 1500
// m/s, v/c is roughly 500x larger than for an airborne radar at the same speed,
// so waveform time-scaling is a first-order effect rather than a correction.
#include "framework.hpp"

#include "phantom/matched_filter.hpp"
#include "phantom/waveform.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr Real        kFs   = 96000;
constexpr std::size_t kFft  = 8192;

// Shared workspace. Static so the test binary does not put megabytes on the
// stack; the library itself never allocates either way.
FftPlan<kFft>                   g_plan;
std::array<Real, kFft>          g_signal;
std::array<Complex, kFft>       g_replica;
std::array<Complex, kFft>       g_spectrum;
std::array<Complex, kFft>       g_scratch;
std::array<Complex, kFft>       g_out;
std::array<Real, kFft>          g_mag;

// Places a pulse into an otherwise silent block at sample offset `at`.
std::size_t place_pulse(const PulseSpec& spec, std::size_t at, Real amplitude,
                        Real doppler = 0) {
    g_signal.fill(0);
    std::array<Real, kFft> tmp{};
    const std::size_t n = (doppler == 0)
                        ? render_real(spec, kFs, tmp)
                        : render_real_doppler(spec, kFs, doppler, tmp);
    for (std::size_t i = 0; i < n && at + i < kFft; ++i) {
        g_signal[at + i] = amplitude * tmp[i];
    }
    return n;
}

struct Peak {
    std::size_t index = 0;
    double value = 0;
};

Peak find_peak(std::span<const Complex> y) {
    Peak p;
    for (std::size_t i = 0; i < y.size(); ++i) {
        const double v = static_cast<double>(std::abs(y[i]));
        if (v > p.value) { p.value = v; p.index = i; }
    }
    return p;
}

PulseSpec lfm(Real f0, Real f1, Real dur) {
    PulseSpec s;
    s.type = (f1 > f0) ? PulseType::LfmUp : PulseType::LfmDown;
    s.f_start_hz = f0;
    s.f_end_hz = f1;
    s.duration_s = dur;
    return s;
}

}  // namespace

PT_TEST(matched_filter_fft_matches_direct_correlation) {
    // The structural test: two independent implementations of the same maths.
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.005));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    PT_CHECK(l > 0);

    place_pulse(s, 1234, 1);
    // Add clutter so the comparison is not just on a clean pulse.
    pt::Rng rng(31337);
    for (Real& v : g_signal) v += static_cast<Real>(0.3 * rng.normal());

    MatchedFilter mf;
    PT_CHECK(matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l),
                                    g_spectrum, mf));

    const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    PT_CHECK(lags == matched_filter_stride(kFft, l));

    static std::array<Complex, kFft> ref{};
    const std::size_t ref_lags = correlate_direct(g_signal, std::span<const Complex>(g_replica.data(), l), ref);
    PT_CHECK(ref_lags == lags);

    double scale = 0.0, worst = 0.0;
    for (std::size_t i = 0; i < lags; ++i) scale = std::max(scale, static_cast<double>(std::abs(ref[i])));
    for (std::size_t i = 0; i < lags; ++i) {
        worst = std::max(worst, static_cast<double>(std::abs(g_out[i] - ref[i])) / scale);
    }
    std::printf("       FFT vs direct correlation: worst relative error = %.3g\n", worst);
    PT_CHECK(worst < pt::tol(1e-12, 1e-3));
}

PT_TEST(matched_filter_peak_lands_on_the_true_delay) {
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.005));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    for (std::size_t delay : {std::size_t(0), std::size_t(1), std::size_t(777),
                              std::size_t(3000), std::size_t(6000)}) {
        place_pulse(s, delay, 1);
        const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        if (delay >= lags) continue;
        const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));
        PT_CHECK(p.index == delay);
    }
}

PT_TEST(matched_filter_peak_height_is_half_the_replica_energy) {
    // Correlating a real cos() against a complex exp() replica:
    //   y = sum a_k cos(phi_k) * a_k exp(-j phi_k)
    //     = sum a_k^2 (cos^2 phi_k - j cos phi_k sin phi_k)
    //     = sum a_k^2 ( (1 + cos 2phi_k)/2 - j sin(2phi_k)/2 )
    //     -> E/2, since the double-frequency terms average away.
    // So a peak at E/2 confirms the replica normalisation end to end.
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.01));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);
    PT_CHECK_REL(mf.replica_energy, static_cast<double>(l), pt::tol(1e-12, 1e-4));

    place_pulse(s, 500, 1);
    const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));

    std::printf("       replica energy %.1f, peak |y| = %.2f (expected E/2 = %.1f)\n",
                static_cast<double>(mf.replica_energy), p.value,
                static_cast<double>(mf.replica_energy) / 2.0);
    PT_CHECK_REL(p.value, static_cast<double>(mf.replica_energy) / 2.0, 0.01);

    // Doubling the received amplitude doubles the peak: the filter is linear.
    place_pulse(s, 500, 2);
    matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const Peak p2 = find_peak(std::span<const Complex>(g_out.data(), lags));
    PT_CHECK_REL(p2.value, 2.0 * p.value, 0.01);
}

PT_TEST(matched_filter_compresses_the_pulse_to_one_over_bandwidth) {
    // The reason chirps are transmitted at all: a 10 ms pulse resolves range to
    // 1/B, not to its own length. At B = 12 kHz that is 83 us -- a 120-fold
    // improvement over the uncompressed pulse.
    const Real bw = 12000;
    const Real dur = static_cast<Real>(0.01);
    const PulseSpec s = lfm(8000, 8000 + bw, dur);
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    place_pulse(s, 2000, 1);
    const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));

    // Walk out from the peak to the -3 dB points.
    const double half = p.value / std::sqrt(2.0);
    std::size_t lo = p.index, hi = p.index;
    while (lo > 0 && static_cast<double>(std::abs(g_out[lo])) > half) --lo;
    while (hi + 1 < lags && static_cast<double>(std::abs(g_out[hi])) > half) ++hi;

    const double width_samples = static_cast<double>(hi - lo);
    const double width_s = width_samples / static_cast<double>(kFs);
    const double expected_s = 1.0 / static_cast<double>(bw);
    std::printf("       -3 dB width %.1f samples = %.1f us  (1/B = %.1f us), "
                "pulse was %.0f us\n",
                width_samples, width_s * 1e6, expected_s * 1e6,
                static_cast<double>(dur) * 1e6);
    // The rectangular-taper mainlobe is 0.886/B between -3 dB points; allow for
    // the sample grid on top of that.
    PT_CHECK(width_s > 0.6 * expected_s);
    PT_CHECK(width_s < 1.6 * expected_s);
    PT_CHECK(width_s < static_cast<double>(dur) / 50.0);
}

PT_TEST(matched_filter_coherent_gain_is_half_the_replica_length) {
    // Output SNR / input SNR:
    //   |y_peak|^2 = (A L / 2)^2         signal
    //   E|y_noise|^2 = sigma^2 * L       noise, white; sum(cos^2 + sin^2) = L
    //   -> output SNR = A^2 L / (4 sigma^2)
    //   input SNR = (A^2/2) / sigma^2
    //   -> gain = L/2
    // The factor of two is the half of a real signal's energy that lies in the
    // negative-frequency half-plane and is discarded by an analytic replica.
    //
    // Signal and noise are measured SEPARATELY rather than from one noisy
    // realisation: the signal peak is deterministic, and averaging noise-only
    // power over many lags gives a far tighter estimate than 60 trials of a
    // noncentral chi-square would.
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.01));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    const double sigma = 1.0;
    const double amplitude = 0.05;
    const std::size_t at = 1500;

    // Signal only.
    place_pulse(s, at, static_cast<Real>(amplitude));
    const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const double signal_power = static_cast<double>(std::norm(g_out[at]));

    // Noise only.
    pt::Rng rng(24680);
    double noise_sum = 0.0;
    std::size_t noise_count = 0;
    for (std::size_t t = 0; t < 40; ++t) {
        for (Real& v : g_signal) v = static_cast<Real>(sigma * rng.normal());
        matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        for (std::size_t i = 0; i < lags; i += 5) {
            noise_sum += static_cast<double>(std::norm(g_out[i]));
            ++noise_count;
        }
    }
    const double noise_power = noise_sum / static_cast<double>(noise_count);

    const double out_snr = signal_power / noise_power;
    const double in_snr = (amplitude * amplitude / 2.0) / (sigma * sigma);
    const double gain = out_snr / in_snr;
    const double expected = static_cast<double>(l) / 2.0;

    std::printf("       L = %zu: noise power %.1f (theory sigma^2*L = %.1f)\n",
                l, noise_power, static_cast<double>(l));
    std::printf("       input SNR %.1f dB -> output SNR %.1f dB, gain %.1f (theory L/2 = %.1f)\n",
                10.0 * std::log10(in_snr), 10.0 * std::log10(out_snr), gain, expected);
    PT_CHECK_REL(noise_power, static_cast<double>(l), 0.05);
    PT_CHECK(gain > 0.95 * expected);
    PT_CHECK(gain < 1.05 * expected);
}

PT_TEST(matched_filter_detects_a_pulse_buried_in_noise) {
    // The practical consequence of that gain: a pulse 17 dB BELOW the noise is
    // still found, and found in the right cell.
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.01));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    const double amplitude = 0.3;   // input SNR = -13.5 dB, output ~ +13 dB
    const std::size_t at = 1500;
    pt::Rng rng(1357);
    std::size_t hits = 0;
    const std::size_t trials = 100;

    for (std::size_t t = 0; t < trials; ++t) {
        place_pulse(s, at, static_cast<Real>(amplitude));
        for (Real& v : g_signal) v += static_cast<Real>(rng.normal());
        const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));
        if (p.index + 1 >= at && p.index <= at + 1) ++hits;
    }
    std::printf("       input SNR %.1f dB: %zu/%zu peaks within one sample of truth\n",
                10.0 * std::log10(amplitude * amplitude / 2.0), hits, trials);
    // Not 100%: at ~13 dB output SNR a noise cell occasionally out-peaks the
    // pulse. Asserting perfection here would only mean the SNR was chosen too
    // generously to be testing anything.
    PT_CHECK(hits >= 90);
}

PT_TEST(matched_filter_cw_is_doppler_sensitive) {
    // A CW pulse is a narrow spike in frequency, so a Doppler shift walks it
    // straight out of the filter's passband. The first null is at a shift of
    // 1/T Hz, i.e. at a closing speed of c/(f0*T).
    const Real f0 = 12000;
    const Real dur = static_cast<Real>(0.05);
    PulseSpec s;
    s.type = PulseType::Cw;
    s.f_start_hz = f0;
    s.duration_s = dur;

    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    place_pulse(s, 100, 1, 0);
    std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const double ref_peak = find_peak(std::span<const Complex>(g_out.data(), lags)).value;

    // Null when the Doppler frequency shift reaches 1/T: v = c/(f0*T).
    const double v_null = 1500.0 / (static_cast<double>(f0) * static_cast<double>(dur));
    std::printf("       CW %.0f Hz x %.0f ms: predicted first null at %.2f m/s\n",
                static_cast<double>(f0), static_cast<double>(dur) * 1e3, v_null);

    place_pulse(s, 100, 1, static_cast<Real>(v_null / 1500.0));
    lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
    const double null_peak = find_peak(std::span<const Complex>(g_out.data(), lags)).value;

    const double loss_db = 20.0 * std::log10(null_peak / ref_peak);
    std::printf("       peak at that speed: %.1f dB relative to zero Doppler\n", loss_db);
    // Not a perfect null: the received pulse is time-scaled, so it is also
    // shorter than the replica, which partially fills the sinc null in.
    PT_CHECK(loss_db < -8.0);

    // The contrast is the point. An LFM over the same band and duration barely
    // notices the same speed.
    const PulseSpec chirp = lfm(f0 - 4000, f0 + 4000, dur);
    static std::array<Complex, kFft> rep2{}, spec2{};
    const std::size_t l2 = render_analytic(chirp, kFs, rep2);
    MatchedFilter mf2;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(rep2.data(), l2), spec2, mf2);

    place_pulse(chirp, 100, 1, 0);
    lags = matched_filter_apply(g_plan.view(), mf2, g_signal, g_scratch, g_out);
    const double chirp_ref = find_peak(std::span<const Complex>(g_out.data(), lags)).value;
    place_pulse(chirp, 100, 1, static_cast<Real>(v_null / 1500.0));
    lags = matched_filter_apply(g_plan.view(), mf2, g_signal, g_scratch, g_out);
    const double chirp_dop = find_peak(std::span<const Complex>(g_out.data(), lags)).value;
    const double chirp_loss = 20.0 * std::log10(chirp_dop / chirp_ref);
    std::printf("       LFM over the same band at the same speed: %.2f dB\n", chirp_loss);
    PT_CHECK(chirp_loss > -1.0);
    PT_CHECK(chirp_loss > loss_db + 6.0);
}

PT_TEST(matched_filter_lfm_range_doppler_coupling_matches_theory) {
    // An LFM keeps its peak under Doppler but BIASES the delay estimate. The
    // textbook narrowband result is dt = -f_doppler/mu with f_doppler taken at
    // the CENTRE frequency. That is the wrong formula underwater.
    //
    // Here Doppler is a time SCALING, not a frequency shift. With alpha = 1+v/c
    // the received phase is
    //     phi_rx(t) = 2*pi (alpha f0 t + mu alpha^2 t^2 / 2)
    // so the received chirp rate is mu*alpha^2 and the received duration T/alpha.
    // Minimising the residual phase against the replica over the overlap -- the
    // least-squares linear fit to t^2 on [0,T'] has slope T' -- gives
    //
    //     dt = -(v/c) * f_end / mu,     f_end = f0 + mu*T
    //
    // the END of the sweep, not its centre. For an 8-16 kHz upsweep that is a
    // 33% error that does not vanish at low speed. The cases below discriminate
    // between the two formulas, including a downsweep where they share a sign
    // but differ in magnitude.
    struct Case { Real f0, f1, dur; };
    const Case cases[] = {
        {8000, 16000, static_cast<Real>(0.04)},
        {10000, 12000, static_cast<Real>(0.04)},
        {16000, 8000, static_cast<Real>(0.04)},   // downsweep
    };

    static std::array<Complex, kFft> rep{}, spec{};
    for (const Case& c : cases) {
        const PulseSpec s = lfm(c.f0, c.f1, c.dur);
        const double mu = static_cast<double>(s.chirp_rate_hz_s());
        const double f_end = static_cast<double>(c.f1);
        const double f_c = static_cast<double>(s.centre_frequency_hz());

        const std::size_t l = render_analytic(s, kFs, rep);
        MatchedFilter mf;
        matched_filter_prepare(g_plan.view(), std::span<const Complex>(rep.data(), l), spec, mf);

        const std::size_t at = 2000;
        place_pulse(s, at, 1, 0);
        std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        const double ref_peak = find_peak(std::span<const Complex>(g_out.data(), lags)).value;
        PT_CHECK(find_peak(std::span<const Complex>(g_out.data(), lags)).index == at);

        std::printf("       %.0f-%.0f Hz, mu = %+.3g Hz/s\n",
                    static_cast<double>(c.f0), static_cast<double>(c.f1), mu);
        for (double v : {2.0, 5.0, 10.0}) {
            const double delta = v / 1500.0;
            place_pulse(s, at, 1, static_cast<Real>(delta));
            lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
            const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));

            const double measured = (static_cast<double>(p.index) - static_cast<double>(at))
                                  / static_cast<double>(kFs);
            const double wideband = -delta * f_end / mu;
            const double narrowband = -delta * f_c / mu;
            const double loss_db = 20.0 * std::log10(p.value / ref_peak);
            std::printf("       v=%4.1f m/s: %+8.1f us   wideband %+8.1f   "
                        "narrowband %+8.1f   peak %+5.2f dB\n",
                        v, measured * 1e6, wideband * 1e6, narrowband * 1e6, loss_db);

            // The wideband expression is itself a least-squares phase fit, and
            // the measured peak is quantised to the sample grid, so a few
            // samples of slack is honest. It is still far tighter than the gap
            // to the narrowband formula, which the next check makes explicit.
            PT_CHECK_NEAR(measured, wideband, 3.0 / static_cast<double>(kFs));
            PT_CHECK(loss_db > -6.0);   // the peak survives: "Doppler tolerant"
        }
    }

    // And the discriminating claim: on the widest-sweep case the narrowband
    // formula is wrong by far more than the tolerance above allows. Without
    // this, passing the test would only show the tolerance was generous.
    {
        const PulseSpec s = lfm(8000, 16000, static_cast<Real>(0.04));
        const double mu = static_cast<double>(s.chirp_rate_hz_s());
        const double delta = 10.0 / 1500.0;
        const double wideband = -delta * 16000.0 / mu;
        const double narrowband = -delta * static_cast<double>(s.centre_frequency_hz()) / mu;
        const double gap_samples = std::fabs(wideband - narrowband) * static_cast<double>(kFs);
        std::printf("       narrowband formula is off by %.1f samples on the 8-16 kHz case\n",
                    gap_samples);
        PT_CHECK(gap_samples > 6.0);
    }
}

PT_TEST(matched_filter_hfm_beats_lfm_under_doppler) {
    // The reason HFM is the underwater default. A time-scaled HFM is, to first
    // order, a DELAYED HFM -- the scaling maps the waveform family onto itself
    // -- so the matched filter peak survives essentially intact. An LFM's peak
    // degrades steadily; a CW's collapses (see the CW test above).
    const Real f0 = 8000, f1 = 20000;
    const Real dur = static_cast<Real>(0.06);   // TB = 720

    PulseSpec lfm_spec = lfm(f0, f1, dur);
    PulseSpec hfm_spec = lfm_spec;
    hfm_spec.type = PulseType::Hfm;

    auto peak_loss_db = [&](const PulseSpec& spec, double v) {
        static std::array<Complex, kFft> rep{};
        static std::array<Complex, kFft> spec_buf{};
        const std::size_t l = render_analytic(spec, kFs, rep);
        MatchedFilter mf;
        matched_filter_prepare(g_plan.view(), std::span<const Complex>(rep.data(), l), spec_buf, mf);

        place_pulse(spec, 200, 1, 0);
        std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        const double ref = find_peak(std::span<const Complex>(g_out.data(), lags)).value;

        place_pulse(spec, 200, 1, static_cast<Real>(v / 1500.0));
        lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        const double got = find_peak(std::span<const Complex>(g_out.data(), lags)).value;
        return 20.0 * std::log10(got / ref);
    };

    std::printf("       TB = %.0f, 60 ms sweep over 8-20 kHz\n",
                static_cast<double>(lfm_spec.time_bandwidth_product()));
    std::printf("       %8s %10s %10s %13s\n", "v (m/s)", "LFM (dB)", "HFM (dB)", "HFM margin");
    for (double v : {5.0, 15.0, 30.0}) {
        const double a_db = peak_loss_db(lfm_spec, v);
        const double b_db = peak_loss_db(hfm_spec, v);
        std::printf("       %8.1f %10.2f %10.2f %10.2f dB\n", v, a_db, b_db, b_db - a_db);
        PT_CHECK(b_db > -1.0);          // HFM holds its peak across the range
        PT_CHECK(b_db > a_db + 3.0);    // and beats the LFM by a clear margin
    }
    // At ordinary AUV speeds the LFM is already losing several dB.
    PT_CHECK(peak_loss_db(lfm_spec, 15.0) < -6.0);
}

PT_TEST(matched_filter_parabolic_interpolation_refines_the_peak) {
    // Sub-sample delay estimation. A pulse placed at a fractional offset must
    // be recovered to well under one sample, or the range estimate is quantised
    // to 1/fs -- 1.6 cm of two-way range at 96 kHz, but 8 m at 100 Hz.
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.005));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;
    matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l), g_spectrum, mf);

    double worst = 0.0;
    for (double frac : {0.0, 0.2, 0.35, -0.3, 0.5}) {
        // Render the pulse shifted by a fractional sample by evaluating its
        // phase at offset times.
        g_signal.fill(0);
        const std::size_t at = 1000;
        const Real dt = static_cast<Real>(1) / kFs;
        for (std::size_t i = 0; i < l; ++i) {
            const Real t = (static_cast<Real>(i) - static_cast<Real>(frac)) * dt;
            if (t < 0 || t > s.duration_s) continue;
            g_signal[at + i] = std::cos(pulse_phase(s, t));
        }

        const std::size_t lags = matched_filter_apply(g_plan.view(), mf, g_signal, g_scratch, g_out);
        for (std::size_t i = 0; i < lags; ++i) g_mag[i] = std::abs(g_out[i]);
        const Peak p = find_peak(std::span<const Complex>(g_out.data(), lags));
        const Real off = parabolic_peak_offset(std::span<const Real>(g_mag.data(), lags), p.index);

        const double estimate = static_cast<double>(p.index) + static_cast<double>(off);
        const double truth = static_cast<double>(at) + frac;
        worst = std::max(worst, std::fabs(estimate - truth));
        std::printf("       true offset %+.2f -> estimated %+.3f samples\n",
                    frac, estimate - static_cast<double>(at));
    }
    PT_CHECK(worst < 0.1);
}

PT_TEST(matched_filter_rejects_bad_buffers) {
    const PulseSpec s = lfm(8000, 20000, static_cast<Real>(0.005));
    const std::size_t l = render_analytic(s, kFs, g_replica);
    MatchedFilter mf;

    // Spectrum buffer must be exactly the transform length.
    std::array<Complex, 64> small{};
    PT_CHECK(!matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l),
                                     small, mf));
    // Replica longer than the transform cannot be zero-padded into it.
    PT_CHECK(!matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), kFft + 1),
                                     g_spectrum, mf));
    PT_CHECK(matched_filter_prepare(g_plan.view(), std::span<const Complex>(g_replica.data(), l),
                                    g_spectrum, mf));
    // Wrong block length is a refusal, not an overrun.
    std::array<Real, 128> short_block{};
    PT_CHECK(matched_filter_apply(g_plan.view(), mf, short_block, g_scratch, g_out) == 0);
    // An empty replica has no energy and cannot be a filter.
    PT_CHECK(!matched_filter_prepare(g_plan.view(), std::span<const Complex>(), g_spectrum, mf));

    PT_CHECK(matched_filter_stride(1024, 100) == 925);
    PT_CHECK(matched_filter_stride(64, 100) == 0);
}
