// SPDX-License-Identifier: Apache-2.0
// Reverberation.
//
// Three closed forms carry this module. The 30 log10(r) and 20 log10(r) decay
// laws are the signature of which mechanism dominates, and the cancellation of
// source level from the echo-to-reverberation ratio is the result that decides
// how an active sonar should be designed.
#include "framework.hpp"

#include "phantom/ping_analyzer.hpp"
#include "phantom/reverberation.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {
constexpr double kPiD = 3.14159265358979323846;
}

// ---------------------------------------------------------------------------
// Scattering strengths
// ---------------------------------------------------------------------------

PT_TEST(lambert_law_is_twenty_log_sine) {
    // S_s = mu + 20 log10(sin theta). The sin^2 in intensity is the projected
    // area entering twice, once on transmit and once on receive.
    const Real mu = -27;
    std::printf("       mu = -27 dB\n       %10s %14s %14s\n",
                "grazing", "Lambert (dB)", "mu+20log10 sin");
    for (int deg : {5, 15, 30, 60, 90}) {
        const double th = static_cast<double>(deg) * kPiD / 180.0;
        const double got = static_cast<double>(
            lambert_bottom_scattering_db(mu, deg2rad(static_cast<Real>(deg))));
        const double expect = -27.0 + 20.0 * std::log10(std::sin(th));
        std::printf("       %9dd %14.3f %14.3f\n", deg, got, expect);
        PT_CHECK_NEAR(got, expect, pt::tol(1e-9, 1e-4));
    }
    // At normal incidence the sine is one, so the strength IS the coefficient.
    PT_CHECK_NEAR(lambert_bottom_scattering_db(mu, kHalfPi), -27.0, pt::tol(1e-9, 1e-4));
    // Monotone: a steeper look sees more backscatter.
    for (int deg = 2; deg < 90; ++deg) {
        PT_CHECK(lambert_bottom_scattering_db(mu, deg2rad(static_cast<Real>(deg + 1)))
               > lambert_bottom_scattering_db(mu, deg2rad(static_cast<Real>(deg))));
    }
    // Grazing to zero does not run to negative infinity.
    PT_CHECK(lambert_bottom_scattering_db(mu, 0) > -200);
    PT_CHECK(std::isfinite(static_cast<double>(lambert_bottom_scattering_db(mu, 0))));
}

PT_TEST(chapman_harris_behaves_like_surface_backscatter) {
    // The coefficients are NOT verified against the original paper -- see
    // docs/validation.md. What is verified here is that the fit behaves: the
    // right band, monotone in wind, monotone in angle.
    std::printf("       %8s %10s %10s %14s\n", "wind(kt)", "f (Hz)", "graze", "S_s (dB)");
    for (double v : {5.0, 15.0, 30.0}) {
        for (int deg : {5, 20, 45}) {
            const double s = static_cast<double>(chapman_harris_surface_scattering_db(
                static_cast<Real>(v), 2000, deg2rad(static_cast<Real>(deg))));
            std::printf("       %8.0f %10.0f %9dd %14.2f\n", v, 2000.0, deg, s);
            // Surface backscatter lives in this band; outside it something is
            // wrong with the formula or its units. The lower edge is generous
            // because the near-calm, near-grazing corner (5 knots at 5 degrees
            // gives -83 dB) is exactly where an empirical fit is least
            // trustworthy and the physical return is genuinely tiny.
            PT_CHECK(s > -95.0);
            PT_CHECK(s < -10.0);
        }
    }

    // More wind, more scattering.
    for (double v = 3; v < 40; v += 1) {
        PT_CHECK(chapman_harris_surface_scattering_db(static_cast<Real>(v + 1), 2000,
                                                      deg2rad(static_cast<Real>(20)))
               > chapman_harris_surface_scattering_db(static_cast<Real>(v), 2000,
                                                      deg2rad(static_cast<Real>(20))));
    }
    // Steeper grazing, more backscatter.
    for (int deg = 2; deg < 80; ++deg) {
        PT_CHECK(chapman_harris_surface_scattering_db(10, 2000, deg2rad(static_cast<Real>(deg + 1)))
               > chapman_harris_surface_scattering_db(10, 2000, deg2rad(static_cast<Real>(deg))));
    }

    // The angular term vanishes at 30 degrees by construction, so the value
    // there depends only on wind and frequency -- a structural check on the
    // log10(theta/30) form.
    for (double v : {8.0, 20.0}) {
        const double beta = 158.0 * std::pow(v * std::cbrt(2000.0), -0.58);
        const double expect = -42.4 * std::log10(beta) + 2.6;
        const double got = static_cast<double>(chapman_harris_surface_scattering_db(
            static_cast<Real>(v), 2000, deg2rad(static_cast<Real>(30))));
        PT_CHECK_NEAR(got, expect, pt::tol(1e-9, 1e-3));
    }

    PT_CHECK(chapman_harris_surface_scattering_db(0, 2000, deg2rad(static_cast<Real>(20))) < -100);
    PT_CHECK(chapman_harris_surface_scattering_db(10, 0, deg2rad(static_cast<Real>(20))) < -100);
    PT_CHECK(chapman_harris_surface_scattering_db(10, 2000, 0) < -100);
    PT_CHECK_NEAR(mps_to_knots(1), 1.943844, 1e-6);
}

// ---------------------------------------------------------------------------
// Geometry and the range laws
// ---------------------------------------------------------------------------

PT_TEST(ensonified_extent_scales_as_the_geometry_says) {
    // Area grows linearly with range, volume quadratically. Both scale with the
    // range resolution c*tau/2.
    const Real phi = static_cast<Real>(0.2);
    const Real tau = static_cast<Real>(0.01);
    const Real c = 1500;

    PT_CHECK_REL(ensonified_area_m2(1000, phi, tau, c),
                 1000.0 * 0.2 * (1500.0 * 0.01 / 2.0), pt::tol(1e-12, 1e-5));
    PT_CHECK_REL(ensonified_area_m2(2000, phi, tau, c) / ensonified_area_m2(1000, phi, tau, c),
                 2.0, pt::tol(1e-12, 1e-5));
    PT_CHECK_REL(ensonified_volume_m3(2000, phi, tau, c) / ensonified_volume_m3(1000, phi, tau, c),
                 4.0, pt::tol(1e-12, 1e-5));
    // Halving the pulse halves both.
    PT_CHECK_REL(ensonified_area_m2(1000, phi, tau / 2, c) / ensonified_area_m2(1000, phi, tau, c),
                 0.5, pt::tol(1e-12, 1e-5));

    PT_CHECK(ensonified_area_m2(0, phi, tau, c) == 0);
    PT_CHECK(ensonified_area_m2(1000, 0, tau, c) == 0);
    PT_CHECK(ensonified_area_m2(1000, phi, 0, c) == 0);
    PT_CHECK(ensonified_volume_m3(1000, 0, tau, c) == 0);
}

PT_TEST(boundary_reverberation_falls_as_thirty_log_range) {
    // Two-way spreading costs 40 log10(r); the ensonified area gives back
    // 10 log10(r). The 30 log10(r) that remains is the signature of boundary
    // reverberation, and it is how you tell it from the volume kind.
    const Real sl = 200, ss = -30, phi = static_cast<Real>(0.2), tau = static_cast<Real>(0.01);
    const Real c = 1500;

    std::printf("       %10s %14s %16s\n", "r (m)", "RL (dB)", "fall from 100 m");
    const double ref = static_cast<double>(reverberation_level_db(
        sl, static_cast<Real>(20.0 * std::log10(100.0)), ss,
        ensonified_area_m2(100, phi, tau, c)));
    for (double r : {100.0, 316.23, 1000.0, 3162.3, 10000.0}) {
        const double tl = 20.0 * std::log10(r);
        const double rl = static_cast<double>(reverberation_level_db(
            sl, static_cast<Real>(tl), ss,
            ensonified_area_m2(static_cast<Real>(r), phi, tau, c)));
        const double fall = ref - rl;
        const double expect = 30.0 * std::log10(r / 100.0);
        std::printf("       %10.0f %14.3f %10.3f (%.3f)\n", r, rl, fall, expect);
        PT_CHECK_NEAR(fall, expect, pt::tol(1e-9, 1e-3));
    }
}

PT_TEST(volume_reverberation_falls_as_twenty_log_range) {
    // The volume grows as r^2, so it gives back 20 log10(r) against the
    // 40 log10(r) of two-way spreading. A decay measurement that shows 20 rather
    // than 30 says the water column is scattering, not the boundary.
    const Real sl = 200, sv = -70, omega = static_cast<Real>(0.05);
    const Real tau = static_cast<Real>(0.01), c = 1500;

    const double ref = static_cast<double>(reverberation_level_db(
        sl, static_cast<Real>(20.0 * std::log10(100.0)), sv,
        ensonified_volume_m3(100, omega, tau, c)));
    std::printf("       %10s %14s %16s\n", "r (m)", "RL (dB)", "fall from 100 m");
    for (double r : {100.0, 1000.0, 10000.0}) {
        const double tl = 20.0 * std::log10(r);
        const double rl = static_cast<double>(reverberation_level_db(
            sl, static_cast<Real>(tl), sv,
            ensonified_volume_m3(static_cast<Real>(r), omega, tau, c)));
        const double fall = ref - rl;
        std::printf("       %10.0f %14.3f %10.3f (%.3f)\n",
                    r, rl, fall, 20.0 * std::log10(r / 100.0));
        PT_CHECK_NEAR(fall, 20.0 * std::log10(r / 100.0), pt::tol(1e-9, 1e-3));
    }
}

// ---------------------------------------------------------------------------
// The result that decides sonar design
// ---------------------------------------------------------------------------

PT_TEST(source_level_cancels_from_the_echo_to_reverberation_ratio) {
    // EL - RL = TS - S_s - 10 log10(A). Source level and transmission loss are
    // absent because they raise the echo and the background together.
    //
    // The consequence is not a detail: against reverberation, a bigger
    // transmitter buys NOTHING. Only shrinking the ensonified area helps.
    const Real ts = 10, ss = -30, phi = static_cast<Real>(0.2);
    const Real tau = static_cast<Real>(0.01), c = 1500;
    const Real range = 2000;
    const Real area = ensonified_area_m2(range, phi, tau, c);

    const double ratio = static_cast<double>(echo_to_reverberation_ratio_db(ts, ss, area));
    std::printf("       TS %.0f dB, S_s %.0f dB, area %.0f m^2 -> E/R = %.2f dB\n",
                static_cast<double>(ts), static_cast<double>(ss),
                static_cast<double>(area), ratio);

    // Computed the long way, sweeping source level and transmission loss: the
    // ratio must not move at all.
    std::printf("       %10s %10s %12s %12s %10s\n", "SL", "TL", "EL", "RL", "EL-RL");
    for (double sl : {160.0, 200.0, 240.0}) {
        for (double tl : {40.0, 66.0, 90.0}) {
            const double el = sl - 2.0 * tl + static_cast<double>(ts);
            const double rl = static_cast<double>(reverberation_level_db(
                static_cast<Real>(sl), static_cast<Real>(tl), ss, area));
            std::printf("       %10.0f %10.0f %12.2f %12.2f %10.4f\n",
                        sl, tl, el, rl, el - rl);
            PT_CHECK_NEAR(el - rl, ratio, pt::tol(1e-9, 1e-3));
        }
    }

    // What DOES help: shorter pulse, narrower beam. Both shrink A.
    const double tenth_pulse = static_cast<double>(echo_to_reverberation_ratio_db(
        ts, ss, ensonified_area_m2(range, phi, tau / 10, c)));
    const double tenth_beam = static_cast<double>(echo_to_reverberation_ratio_db(
        ts, ss, ensonified_area_m2(range, phi / 10, tau, c)));
    std::printf("       pulse /10: %+.2f dB   beam /10: %+.2f dB   (both +10 dB)\n",
                tenth_pulse - ratio, tenth_beam - ratio);
    PT_CHECK_NEAR(tenth_pulse - ratio, 10.0, pt::tol(1e-9, 1e-3));
    PT_CHECK_NEAR(tenth_beam - ratio, 10.0, pt::tol(1e-9, 1e-3));

    // And that is the argument for pulse compression: a chirp's range
    // resolution is 1/B, not its length. A 20 ms pulse with 12 kHz of bandwidth
    // behaves as if it were 83 us long.
    const Real compressed_tau = static_cast<Real>(1.0 / 12000.0);
    const double compressed = static_cast<double>(echo_to_reverberation_ratio_db(
        ts, ss, ensonified_area_m2(range, phi, compressed_tau, c)));
    const double gain = compressed
                      - static_cast<double>(echo_to_reverberation_ratio_db(
                            ts, ss, ensonified_area_m2(range, phi,
                                                       static_cast<Real>(0.02), c)));
    std::printf("       20 ms pulse compressed to 1/B = 83 us: %+.1f dB of E/R\n", gain);
    PT_CHECK_NEAR(gain, 10.0 * std::log10(0.02 * 12000.0), pt::tol(1e-9, 1e-3));
    PT_CHECK(gain > 23.0);
}

PT_TEST(reverberation_limited_range_is_where_it_meets_the_noise) {
    // Inside this range the geometry is reverberation-limited and transmit
    // power does not help; outside it, it does.
    const Real sl = 200, ss = -30, phi = static_cast<Real>(0.2);
    const Real tau = static_cast<Real>(0.01), c = 1500;

    std::printf("       %12s %16s\n", "noise (dB)", "crossover (m)");
    double prev = 1e18;
    for (double nl : {40.0, 60.0, 80.0}) {
        const double r = static_cast<double>(reverberation_limited_range_m(
            sl, ss, phi, tau, c, static_cast<Real>(nl)));
        std::printf("       %12.0f %16.1f\n", nl, r);
        PT_CHECK(r > 0);
        PT_CHECK(r < prev);   // a louder ocean shortens the reverb-limited zone
        prev = r;

        // At the crossover the two levels must actually be equal.
        const double tl = 20.0 * std::log10(r);
        const double rl = static_cast<double>(reverberation_level_db(
            sl, static_cast<Real>(tl), ss,
            ensonified_area_m2(static_cast<Real>(r), phi, tau, c)));
        PT_CHECK_NEAR(rl, nl, 0.1);
    }

    // A quiet transmitter into a loud ocean is never reverberation-limited.
    PT_CHECK(reverberation_limited_range_m(100, ss, phi, tau, c, 120) == 0);
    PT_CHECK(reverberation_limited_range_m(sl, ss, phi, tau, c, 60, 100, 50) == 0);
}

// ---------------------------------------------------------------------------
// What it does to the detector
// ---------------------------------------------------------------------------

PT_TEST(reverberation_envelope_decays_then_meets_the_noise_floor) {
    ReverbProfile prof;
    prof.source_level_db = 200;
    prof.scattering_strength_db = -30;
    prof.noise_level_db = 60;
    prof.reference_level_db = 200;

    // Long enough to actually reach the floor. With these parameters the
    // crossover is at 5.3 km, i.e. 7.1 s of two-way travel, so a 2 s record
    // would only show the decay and never the transition.
    constexpr Real kFs = 2000;
    static std::array<Real, 24000> env{};   // 12 s
    const std::size_t n = reverberation_envelope(prof, kFs, static_cast<Real>(0.05), env);
    PT_CHECK(n == env.size());

    // Monotone decay until the ambient floor takes over, then flat.
    std::printf("       %10s %10s %14s\n", "t (s)", "r (m)", "envelope (dB)");
    for (std::size_t i : {std::size_t(0), std::size_t(2000), std::size_t(10000),
                          std::size_t(23999)}) {
        const double t = 0.05 + static_cast<double>(i) / 2000.0;
        std::printf("       %10.3f %10.0f %14.2f\n", t, 1500.0 * t / 2.0,
                    20.0 * std::log10(static_cast<double>(env[i])));
    }
    for (std::size_t i = 1; i < n; ++i) {
        PT_CHECK(env[i] <= env[i - 1] * static_cast<Real>(1.0001));
    }
    // The floor is the ambient level, not zero.
    const double floor_db = 20.0 * std::log10(static_cast<double>(env[n - 1]))
                          + static_cast<double>(prof.reference_level_db);
    std::printf("       settles at %.2f dB (ambient is %.0f)\n",
                floor_db, static_cast<double>(prof.noise_level_db));
    PT_CHECK(floor_db > static_cast<double>(prof.noise_level_db) - 0.5);
    PT_CHECK(floor_db < static_cast<double>(prof.noise_level_db) + 3.0);

    PT_CHECK(reverberation_envelope(prof, 0, 0, env) == 0);
    PT_CHECK(reverberation_envelope(prof, kFs, 0, std::span<Real>()) == 0);
}

PT_TEST(cfar_tracks_reverberation_where_a_fixed_threshold_cannot) {
    // The practical reason CA-CFAR exists. Against a background that falls by
    // tens of dB across one block, a fixed threshold is wrong everywhere: set
    // it for the near field and the far field is deaf, set it for the far field
    // and the near field is a wall of false alarms.
    constexpr Real        kFs  = 96000;
    constexpr std::size_t kFft = 8192;

    static PulseBank<4, kFft> bank(kFs);
    bank.clear();
    PulseSpec lfm;
    lfm.type = PulseType::LfmUp;
    lfm.f_start_hz = 8000;
    lfm.f_end_hz = 20000;
    lfm.duration_s = static_cast<Real>(0.01);
    bank.add(lfm);

    static AnalyzerScratch<kFft> scratch;
    static std::array<Real, kFft> block{};
    static std::array<Real, kFft> env{};
    std::array<PulseDescriptor, 32> pdw{};

    // A steeply decaying reverberation background across the block.
    ReverbProfile prof;
    prof.source_level_db = 200;
    prof.scattering_strength_db = -25;
    prof.pulse_length_s = static_cast<Real>(0.01);
    prof.noise_level_db = 40;
    prof.reference_level_db = 150;
    reverberation_envelope(prof, kFs, static_cast<Real>(0.02), env);

    const double first_db = 20.0 * std::log10(static_cast<double>(env[0]));
    const double last_db = 20.0 * std::log10(static_cast<double>(env[kFft - 1]));
    std::printf("       background falls %.1f dB across the block\n", first_db - last_db);
    PT_CHECK(first_db - last_db > 15.0);

    // Put a target well out in the decayed region, at a level a near-field
    // threshold would never reach down to.
    const std::size_t at = 6000;
    pt::Rng rng(31415);

    auto trial = [&](Real alpha, std::size_t guard, std::size_t train, bool with_target) {
        for (std::size_t i = 0; i < kFft; ++i) {
            block[i] = env[i] * static_cast<Real>(rng.normal());
        }
        if (with_target) {
            static std::array<Real, kFft> tmp{};
            const std::size_t m = render_real(lfm, kFs, tmp);
            const Real amp = env[at] * static_cast<Real>(6);
            for (std::size_t i = 0; i < m && at + i < kFft; ++i) block[at + i] += amp * tmp[i];
        }
        DetectorConfig cfg;
        cfg.cfar_guard = guard;
        cfg.cfar_train = train;
        cfg.threshold_alpha = alpha;
        cfg.dead_time_s = static_cast<Real>(0.004);
        return analyze_block(bank.view(), cfg, block, 0, scratch.view(), pdw);
    };

    // CA-CFAR: training cells short enough to track the decay.
    const Real alpha = cfar_alpha(256, static_cast<Real>(1e-5));
    std::size_t hits = 0, cfar_false = 0;
    for (int t = 0; t < 20; ++t) {
        const std::size_t n = trial(alpha, 64, 128, true);
        bool found = false;
        for (std::size_t i = 0; i < n; ++i) {
            const auto s = static_cast<std::size_t>(
                static_cast<double>(pdw[i].toa_s) * static_cast<double>(kFs) + 0.5);
            if (s + 2 >= at && s <= at + 2) found = true;
        }
        if (found) ++hits;
        cfar_false += trial(alpha, 64, 128, false);
    }
    std::printf("       CA-CFAR: %zu/20 detections, %zu false alarms in 20 empty blocks\n",
                hits, cfar_false);
    PT_CHECK(hits >= 18);
    PT_CHECK(cfar_false <= 4);

    // A fixed threshold set from the block's average power. The near field is
    // far above it and the far field far below, so it cannot be right anywhere.
    for (std::size_t i = 0; i < kFft; ++i) block[i] = env[i] * static_cast<Real>(rng.normal());
    double mean_power = 0;
    for (std::size_t i = 0; i < kFft; ++i) {
        mean_power += static_cast<double>(env[i]) * static_cast<double>(env[i]);
    }
    mean_power /= static_cast<double>(kFft);
    const double near_power = static_cast<double>(env[0]) * static_cast<double>(env[0]);
    const double far_power = static_cast<double>(env[kFft - 1]) * static_cast<double>(env[kFft - 1]);
    std::printf("       a single threshold at the block mean sits %.1f dB below the near\n"
                "       field and %.1f dB above the far field\n",
                10.0 * std::log10(near_power / mean_power),
                10.0 * std::log10(mean_power / far_power));
    // Both errors are large: that gap is exactly what CFAR removes.
    PT_CHECK(10.0 * std::log10(near_power / mean_power) > 3.0);
    PT_CHECK(10.0 * std::log10(mean_power / far_power) > 3.0);
}
