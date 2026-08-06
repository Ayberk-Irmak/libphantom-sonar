// SPDX-License-Identifier: Apache-2.0
// Waveform synthesis verification.
//
// The phase functions are checked by differentiating the specified phase
// numerically and comparing against the specified instantaneous frequency.
// That catches an integration constant or a factor of two, which is the classic
// chirp bug and is invisible in a plot.
#include "framework.hpp"

#include "phantom/waveform.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {
constexpr Real kFs = 96000;
constexpr double kPiD = 3.14159265358979323846;
}  // namespace

PT_TEST(waveform_cw_frequency_is_constant) {
    PulseSpec s;
    s.type = PulseType::Cw;
    s.f_start_hz = 12000;
    s.duration_s = static_cast<Real>(0.01);
    PT_CHECK(s.valid());
    PT_CHECK_NEAR(s.bandwidth_hz(), 0.0, 1e-12);
    PT_CHECK_NEAR(s.chirp_rate_hz_s(), 0.0, 1e-12);
    PT_CHECK_NEAR(s.centre_frequency_hz(), 12000.0, 1e-9);

    for (int i = 0; i <= 10; ++i) {
        const auto t = static_cast<Real>(static_cast<double>(i) * 0.001);
        PT_CHECK_REL(pulse_phase(s, t), 2.0 * kPiD * 12000.0 * static_cast<double>(t),
                     pt::tol(1e-13, 1e-5));
        PT_CHECK_NEAR(pulse_frequency(s, t), 12000.0, 1e-9);
    }
}

PT_TEST(waveform_lfm_instantaneous_frequency_is_linear) {
    PulseSpec s;
    s.type = PulseType::LfmUp;
    s.f_start_hz = 8000;
    s.f_end_hz = 20000;
    s.duration_s = static_cast<Real>(0.02);
    PT_CHECK(s.valid());
    PT_CHECK_NEAR(s.bandwidth_hz(), 12000.0, 1e-9);
    PT_CHECK_NEAR(s.centre_frequency_hz(), 14000.0, 1e-9);
    PT_CHECK_NEAR(s.chirp_rate_hz_s(), 600000.0, 1e-3);
    PT_CHECK_NEAR(s.time_bandwidth_product(), 240.0, 1e-9);

    const double h = 1e-7;
    double worst = 0.0;
    for (int i = 1; i < 20; ++i) {
        const auto t = static_cast<Real>(static_cast<double>(i) * 0.001);
        const double dphi = static_cast<double>(pulse_phase(s, static_cast<Real>(static_cast<double>(t) + h)))
                          - static_cast<double>(pulse_phase(s, static_cast<Real>(static_cast<double>(t) - h)));
        const double f_num = dphi / (2.0 * h) / (2.0 * kPiD);
        const double f_spec = static_cast<double>(pulse_frequency(s, t));
        worst = std::max(worst, std::fabs(f_num - f_spec) / f_spec);
        PT_CHECK_REL(f_spec, 8000.0 + 600000.0 * static_cast<double>(t), pt::tol(1e-9, 1e-4));
    }
    std::printf("       LFM: worst |dphase/dt - 2*pi*f| / f = %.3g\n", worst);
    // The finite difference is taken on the phase as evaluated in the working
    // precision, so in a float build the step h=1e-7 sits near the noise floor
    // of the difference itself. The bound is loosened for that, not because
    // the phase function is any less exact.
    PT_CHECK(worst < pt::tol(1e-6, 5e-2));
}

PT_TEST(waveform_hfm_has_linear_period) {
    // The defining property of an HFM: the PERIOD is linear in time, so 1/f(t)
    // is a straight line. Together with the LFM test this pins that the two are
    // genuinely different families rather than a copy-paste.
    PulseSpec s;
    s.type = PulseType::Hfm;
    s.f_start_hz = 8000;
    s.f_end_hz = 20000;
    s.duration_s = static_cast<Real>(0.02);
    PT_CHECK(s.valid());

    const double inv0 = 1.0 / 8000.0;
    const double inv1 = 1.0 / 20000.0;
    for (int i = 0; i <= 20; ++i) {
        const double frac = static_cast<double>(i) / 20.0;
        const auto t = static_cast<Real>(frac * 0.02);
        PT_CHECK_REL(1.0 / static_cast<double>(pulse_frequency(s, t)),
                     inv0 + (inv1 - inv0) * frac, pt::tol(1e-9, 1e-4));
    }
    PT_CHECK_REL(pulse_frequency(s, 0), 8000.0, pt::tol(1e-9, 1e-5));
    PT_CHECK_REL(pulse_frequency(s, s.duration_s), 20000.0, pt::tol(1e-9, 1e-5));

    const double h = 1e-7;
    double worst = 0.0;
    for (int i = 1; i < 20; ++i) {
        const auto t = static_cast<Real>(static_cast<double>(i) * 0.001);
        const double dphi = static_cast<double>(pulse_phase(s, static_cast<Real>(static_cast<double>(t) + h)))
                          - static_cast<double>(pulse_phase(s, static_cast<Real>(static_cast<double>(t) - h)));
        const double f_num = dphi / (2.0 * h) / (2.0 * kPiD);
        const double f_spec = static_cast<double>(pulse_frequency(s, t));
        worst = std::max(worst, std::fabs(f_num - f_spec) / f_spec);
    }
    std::printf("       HFM: worst |dphase/dt - 2*pi*f| / f = %.3g\n", worst);
    PT_CHECK(worst < pt::tol(1e-6, 5e-2));

    // Sweeping the same band, an HFM lingers at low frequencies, so its
    // mid-sweep frequency sits below the LFM's arithmetic mean.
    PulseSpec lfm = s;
    lfm.type = PulseType::LfmUp;
    const auto mid = static_cast<Real>(0.01);
    std::printf("       mid-sweep: HFM %.1f Hz vs LFM %.1f Hz\n",
                static_cast<double>(pulse_frequency(s, mid)),
                static_cast<double>(pulse_frequency(lfm, mid)));
    PT_CHECK(pulse_frequency(s, mid) < pulse_frequency(lfm, mid));
}

PT_TEST(waveform_hfm_degenerates_to_cw_gracefully) {
    // As the sweep vanishes the closed-form phase divides by k -> 0. The
    // implementation must fall back to the linear form rather than emit inf.
    PulseSpec s;
    s.type = PulseType::Hfm;
    s.f_start_hz = 10000;
    s.f_end_hz = static_cast<Real>(10000.0000001);
    s.duration_s = static_cast<Real>(0.01);
    for (int i = 0; i <= 10; ++i) {
        const auto t = static_cast<Real>(static_cast<double>(i) * 0.001);
        const double ph = static_cast<double>(pulse_phase(s, t));
        PT_CHECK(std::isfinite(ph));
        PT_CHECK_REL(ph, 2.0 * kPiD * 10000.0 * static_cast<double>(t), pt::tol(1e-5, 1e-3));
    }
}

PT_TEST(waveform_render_lengths_and_amplitude) {
    PulseSpec s;
    s.type = PulseType::LfmUp;
    s.f_start_hz = 5000;
    s.f_end_hz = 15000;
    s.duration_s = static_cast<Real>(0.01);   // 960 samples at 96 kHz
    s.amplitude = 2;

    std::array<Real, 2048> real_buf{};
    std::array<Complex, 2048> cplx_buf{};

    PT_CHECK(pulse_length(s, kFs) == 960);
    PT_CHECK(render_real(s, kFs, real_buf) == 960);
    PT_CHECK(render_analytic(s, kFs, cplx_buf) == 960);

    // The analytic envelope is flat at the specified amplitude; the real
    // waveform is its real part, so it never exceeds it.
    for (std::size_t i = 0; i < 960; ++i) {
        PT_CHECK_NEAR(std::abs(cplx_buf[i]), 2.0, pt::tol(1e-12, 1e-5));
        PT_CHECK(std::fabs(static_cast<double>(real_buf[i])) <= 2.0 + 1e-6);
        PT_CHECK_NEAR(real_buf[i], cplx_buf[i].real(), pt::tol(1e-12, 1e-5));
    }

    std::array<Real, 16> tiny{};
    PT_CHECK(render_real(s, kFs, tiny) == 0);   // refusal, not an overrun
}

PT_TEST(waveform_tapers_have_expected_shape) {
    const Real T = 1;
    PT_CHECK_NEAR(taper_weight(Taper::Hann, 0, T), 0.0, 1e-12);
    PT_CHECK_NEAR(taper_weight(Taper::Hann, T, T), 0.0, pt::tol(1e-12, 1e-6));
    PT_CHECK_NEAR(taper_weight(Taper::Hann, static_cast<Real>(0.5), T), 1.0, pt::tol(1e-12, 1e-6));
    PT_CHECK_NEAR(taper_weight(Taper::Rectangular, static_cast<Real>(0.5), T), 1.0, 1e-12);
    PT_CHECK_NEAR(taper_weight(Taper::Rectangular, static_cast<Real>(1.5), T), 0.0, 1e-12);
    PT_CHECK_NEAR(taper_weight(Taper::Tukey25, static_cast<Real>(0.5), T), 1.0, 1e-12);
    PT_CHECK_NEAR(taper_weight(Taper::Tukey25, 0, T), 0.0, 1e-12);
    PT_CHECK(taper_weight(Taper::Tukey25, static_cast<Real>(0.06), T) > 0);
    PT_CHECK(taper_weight(Taper::Tukey25, static_cast<Real>(0.06), T) < 1);
}

PT_TEST(waveform_doppler_scales_time_not_just_frequency) {
    // A closing target compresses the pulse: the received copy is SHORTER.
    // Modelling Doppler as a pure frequency shift -- the habit carried over
    // from radar -- would leave the duration unchanged, which this pins down.
    PulseSpec s;
    s.type = PulseType::LfmUp;
    s.f_start_hz = 8000;
    s.f_end_hz = 16000;
    s.duration_s = static_cast<Real>(0.02);

    std::array<Real, 8192> rx{};
    std::array<Real, 8192> plain{};
    const Real doppler = static_cast<Real>(10.0 / 1500.0);   // 10 m/s at 1500 m/s

    const std::size_t n0 = render_real(s, kFs, plain);
    const std::size_t n1 = render_real_doppler(s, kFs, doppler, rx);
    std::printf("       10 m/s closing: %zu -> %zu samples (%.2f%% compression)\n",
                n0, n1, 100.0 * (1.0 - static_cast<double>(n1) / static_cast<double>(n0)));
    PT_CHECK(n1 < n0);
    PT_CHECK_REL(static_cast<double>(n1),
                 static_cast<double>(n0) / (1.0 + static_cast<double>(doppler)), 1e-3);

    // Zero Doppler must reproduce the transmitted waveform exactly.
    const std::size_t n2 = render_real_doppler(s, kFs, 0, rx);
    PT_CHECK(n2 == n0);
    for (std::size_t i = 0; i < n2; ++i) PT_CHECK_NEAR(rx[i], plain[i], pt::tol(1e-12, 1e-5));
}

PT_TEST(waveform_rejects_invalid_specs) {
    PulseSpec s;
    s.type = PulseType::LfmUp;
    s.f_start_hz = 10000;
    s.f_end_hz = 5000;                   // a downsweep declared as an upsweep
    s.duration_s = static_cast<Real>(0.01);
    PT_CHECK(!s.valid());

    s.type = PulseType::Hfm;
    s.f_end_hz = 0;                      // 1/f undefined
    PT_CHECK(!s.valid());

    s.type = PulseType::Cw;
    s.duration_s = 0;
    PT_CHECK(!s.valid());

    s.duration_s = static_cast<Real>(0.01);
    PT_CHECK(s.valid());
    std::array<Real, 64> buf{};
    PT_CHECK(render_real(s, 0, buf) == 0);   // no sample rate
}
