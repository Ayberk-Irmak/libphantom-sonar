// SPDX-License-Identifier: Apache-2.0
// Sound speed equations.
//
// The original strategy here was cross-validation: three equations derived from
// different data sets must agree inside their common validity box, so a mistyped
// coefficient breaks the agreement.
//
// THAT ARGUMENT IS WEAKER THAN IT LOOKS, and v0.11 found out how. The equations
// agree only to ~0.1 m/s, so the check can only see errors bigger than that --
// and the coefficient bug it missed was worth 0.016 m/s. Mutual agreement
// bounds how wrong you can be by how much your methods differ, which is not the
// same as being right.
//
// So the tests below now include the PUBLISHED UNESCO table: 220 values from
// Fofonoff & Millard (1983), tight to the table's own 0.05 m/s rounding.
#include "framework.hpp"

#include "phantom/sound_speed.hpp"
#include "data/unesco44_table.hpp"

using namespace phantom;
using namespace phantom::sound_speed;

PT_TEST(sound_speed_reference_point) {
    // T = 15 C, S = 35 PSU, surface. The textbook "about 1507 m/s" case.
    PT_CHECK_NEAR(medwin(15, 35, 0),       1506.8, 0.5);
    PT_CHECK_NEAR(mackenzie(15, 35, 0),    1506.7, 0.5);
    PT_CHECK_NEAR(chen_millero(15, 35, 0), 1506.7, 0.5);
}

PT_TEST(sound_speed_equations_agree) {
    // Common validity box of all three: 5..25 C, 34..36 PSU, 0..800 m.
    double worst_mw = 0.0;
    double worst_cm = 0.0;
    for (int ti = 5; ti <= 25; ++ti) {
        for (int si = 34; si <= 36; ++si) {
            for (int zi = 0; zi <= 800; zi += 50) {
                const auto t = static_cast<Real>(ti);
                const auto s = static_cast<Real>(si);
                const auto z = static_cast<Real>(zi);
                const double mk = static_cast<double>(mackenzie(t, s, z));
                const double mw = static_cast<double>(medwin(t, s, z));
                const double cm = static_cast<double>(unesco(t, s, z));
                const double dmw = std::fabs(mk - mw);
                const double dcm = std::fabs(mk - cm);
                if (dmw > worst_mw) worst_mw = dmw;
                if (dcm > worst_cm) worst_cm = dcm;
            }
        }
    }
    std::printf("       max |mackenzie - medwin|       = %.3f m/s\n", worst_mw);
    std::printf("       max |mackenzie - chen_millero| = %.3f m/s\n", worst_cm);
    PT_CHECK(worst_mw < 1.0);
    PT_CHECK(worst_cm < 1.0);
}

PT_TEST(sound_speed_monotonicity) {
    // Physically, c rises with temperature (well past 30 C), with salinity, and
    // with depth. Any sign error in a coefficient breaks one of these.
    for (int zi = 0; zi <= 4000; zi += 500) {
        const auto z = static_cast<Real>(zi);
        for (int ti = 2; ti < 29; ++ti) {
            PT_CHECK(mackenzie(static_cast<Real>(ti + 1), 35, z)
                   > mackenzie(static_cast<Real>(ti), 35, z));
        }
        for (int si = 30; si < 39; ++si) {
            PT_CHECK(mackenzie(10, static_cast<Real>(si + 1), z)
                   > mackenzie(10, static_cast<Real>(si), z));
        }
    }
    for (int zi = 0; zi < 8000; zi += 250) {
        PT_CHECK(mackenzie(4, 35, static_cast<Real>(zi + 250))
               > mackenzie(4, 35, static_cast<Real>(zi)));
    }
}

PT_TEST(sound_speed_deep_water) {
    // 5000 m, 4 C, 35 PSU. Deep-ocean measurements sit in the 1540-1560 band;
    // this is exactly where Medwin (valid to 1000 m) must NOT be trusted.
    const double deep = static_cast<double>(mackenzie(4, 35, 5000));
    std::printf("       mackenzie(4 C, 35 PSU, 5000 m) = %.2f m/s\n", deep);
    PT_CHECK(deep > 1530.0 && deep < 1575.0);
}

PT_TEST(depth_pressure_conversion) {
    // Roughly 1 bar per 10 m of seawater, plus compressibility.
    PT_CHECK_NEAR(depth_to_pressure_bar(0, 45),    0.0,   1e-6);
    PT_CHECK_NEAR(depth_to_pressure_bar(1000, 45), 101.0, 1.0);
    PT_CHECK_NEAR(depth_to_pressure_bar(5000, 45), 507.0, 5.0);
    // Latitude effect is real but small: gravity is ~0.5% stronger at the pole.
    PT_CHECK(depth_to_pressure_bar(1000, 90) > depth_to_pressure_bar(1000, 0));
    PT_CHECK(static_cast<double>(depth_to_pressure_bar(1000, 90))
           - static_cast<double>(depth_to_pressure_bar(1000, 0)) < 1.0);
}

PT_TEST(munk_profile_shape) {
    // The canonical benchmark profile: a minimum exactly on the axis.
    PT_CHECK_NEAR(munk(1300), 1500.0, 1e-9);
    for (int zi = 0; zi <= 5000; zi += 100) {
        if (zi == 1300) continue;
        PT_CHECK(static_cast<double>(munk(static_cast<Real>(zi))) > 1500.0);
    }
    std::printf("       munk(0 m) = %.3f  munk(5000 m) = %.3f m/s\n",
                static_cast<double>(munk(0)), static_cast<double>(munk(5000)));
}

PT_TEST(sound_speed_constexpr_paths) {
    // medwin/mackenzie must remain usable at compile time so that profiles can
    // be baked into flash on MCU targets.
    constexpr Real a = medwin(10, 35, 100);
    constexpr Real b = mackenzie(10, 35, 100);
    static_assert(a > 1480 && a < 1530, "medwin not constant-evaluable");
    static_assert(b > 1480 && b < 1530, "mackenzie not constant-evaluable");
    PT_CHECK_NEAR(a, b, 1.0);
}

// ---------------------------------------------------------------------------
// v0.11: against the PUBLISHED UNESCO table, not against ourselves
// ---------------------------------------------------------------------------

PT_TEST(chen_millero_1977_reproduces_the_official_unesco_check_value) {
    // UNESCO Technical Papers in Marine Science 44 (Fofonoff & Millard 1983),
    // p. 48 and the FORTRAN listing on p. 49:
    //   "CHECKVALUE: SVEL=1731.995 M/S, S=40 (PSS-78), T=40 DEG C, P=10000 DBAR"
    const double got = static_cast<double>(sound_speed::chen_millero_1977(
        static_cast<Real>(unesco44::kCheckTemperature),
        static_cast<Real>(unesco44::kCheckSalinity),
        static_cast<Real>(unesco44::kCheckPressure / 10.0)));   // dbar -> bar
    std::printf("       UNESCO 44 check value: published %.3f, computed %.4f\n",
                unesco44::kCheckSpeed, got);
    PT_CHECK_NEAR(got, unesco44::kCheckSpeed, pt::tol(1e-3, 0.05));
}

PT_TEST(chen_millero_1977_reproduces_all_220_published_table_values) {
    // The point of this test is that it compares against numbers PUBLISHED
    // ELSEWHERE. Every other sound-speed test in this file checks the equations
    // against each other, and they agree to ~0.1 m/s -- which is far too loose
    // to see a wrong coefficient in the last digit. This one is tight to the
    // table's own rounding.
    double worst = 0;
    double worst_s = 0, worst_t = 0, worst_p = 0;
    for (int i = 0; i < unesco44::kTableSize; ++i) {
        const unesco44::Entry& e = unesco44::kTable[i];
        const double got = static_cast<double>(sound_speed::chen_millero_1977(
            static_cast<Real>(e.temperature_c68),
            static_cast<Real>(e.salinity_psu),
            static_cast<Real>(e.pressure_dbar / 10.0)));
        const double d = std::fabs(got - e.speed_mps);
        if (d > worst) { worst = d; worst_s = e.salinity_psu;
                         worst_t = e.temperature_c68; worst_p = e.pressure_dbar; }
    }
    std::printf("       %d published values, worst |delta| = %.4f m/s"
                " (at S=%.0f T=%.0f p=%.0f dbar)\n",
                unesco44::kTableSize, worst, worst_s, worst_t, worst_p);
    std::printf("       the table is printed to 0.1 m/s, so 0.05 is its rounding half-width\n");
    // 0.05 in double. In float the polynomial itself loses more than the table
    // rounds, so the bound is the achievable precision rather than the table's.
    PT_CHECK(worst < pt::tol(0.05, 0.35));
}

PT_TEST(the_two_temperature_scales_are_a_real_difference_not_a_typo) {
    // v0.11 found a single chen_millero() carrying Wong & Zhu's coefficients
    // except for A02 and A03, which were the 1977 originals. The fix was to
    // implement BOTH equations properly. This test is what proves they are two
    // legitimate equations rather than one equation and one mistake: converting
    // the temperature scale makes them agree far better than they do raw.
    double worst_raw = 0, worst_converted = 0;
    for (int si = 0; si <= 40; si += 5) {
        for (int ti = 0; ti <= 30; ti += 2) {
            for (int pi = 0; pi <= 1000; pi += 100) {
                const Real s = static_cast<Real>(si);
                const Real p = static_cast<Real>(pi);
                const Real t90 = static_cast<Real>(ti);
                const double its90 = static_cast<double>(
                    sound_speed::chen_millero_its90(t90, s, p));
                // Wrong: same number fed to both, ignoring the scale.
                const double raw = static_cast<double>(
                    sound_speed::chen_millero_1977(t90, s, p));
                // Right: t68 = 1.00024 * t90.
                const double conv = static_cast<double>(
                    sound_speed::chen_millero_1977(sound_speed::t90_to_t68(t90), s, p));
                worst_raw = std::max(worst_raw, std::fabs(raw - its90));
                worst_converted = std::max(worst_converted, std::fabs(conv - its90));
            }
        }
    }
    std::printf("       1977 vs ITS-90, same T number      : %.4f m/s worst\n", worst_raw);
    std::printf("       1977 vs ITS-90, t68 = 1.00024*t90  : %.4f m/s worst\n", worst_converted);
    // Converting must help, and must land inside Wong & Zhu's own stated
    // "within 0.024 m/s" revision size for the UNESCO equation.
    PT_CHECK(worst_converted < worst_raw);
    PT_CHECK(worst_converted < pt::tol(0.03, 0.35));
}

PT_TEST(del_grosso_agrees_with_unesco_where_the_ocean_actually_exists) {
    // An INDEPENDENT equation: different laboratory, different data, different
    // functional form. Agreement here is evidence about the physics; agreement
    // between Medwin, Mackenzie and UNESCO is mostly evidence that they were
    // fitted to overlapping data.
    //
    // The useful result is that WHERE you compare them matters more than the
    // comparison. Their validity boxes are rectangles in (S, T, P), and a
    // rectangle contains combinations the ocean does not -- 26 C water under
    // 1000 bar exists in no sea, and neither equation was ever fitted to a
    // sample of it. So this reports the disagreement banded by depth, over
    // water that could actually be there, rather than one number over a box.

    // A crude but honest envelope: the deep ocean is 2-4 C everywhere, warm
    // water floats, the thermocline is a few hundred metres thick.
    auto plausible_max_temp_c = [](double pressure_bar) {
        return 4.0 + 26.0 * std::exp(-pressure_bar / 30.0);
    };

    struct Band { double max_bar; const char* label; double worst; };
    Band bands[] = {
        { 100,  "  0 - 1000 m  (most sonar work)", 0 },
        { 300,  "  0 - 3000 m  (deep channel)   ", 0 },
        { 600,  "  0 - 6000 m  (abyssal plain)  ", 0 },
        { 981,  "  0 - 9810 m  (Del Grosso max) ", 0 },
    };
    double worst_box = 0, box_s = 0, box_t = 0, box_p = 0;

    for (int si = 30; si <= 40; si += 2) {
        for (int ti = 0; ti <= 30; ti += 2) {
            for (int pi = 0; pi <= 1000; pi += 20) {
                const Real s = static_cast<Real>(si);
                const Real t = static_cast<Real>(ti);
                const Real p_bar = static_cast<Real>(pi);
                const double u = static_cast<double>(sound_speed::chen_millero_its90(t, s, p_bar));
                const double g = static_cast<double>(sound_speed::del_grosso(
                    t, s, sound_speed::bar_to_kgcm2(p_bar)));
                const double d = std::fabs(u - g);
                if (d > worst_box) { worst_box = d; box_s = si; box_t = ti; box_p = pi; }
                if (static_cast<double>(ti) > plausible_max_temp_c(static_cast<double>(pi))) continue;
                for (Band& band : bands) {
                    if (static_cast<double>(pi) <= band.max_bar) band.worst = std::max(band.worst, d);
                }
            }
        }
    }

    std::printf("       UNESCO vs Del Grosso, over water that could exist:\n");
    for (const Band& band : bands) {
        std::printf("       %s : %.3f m/s\n", band.label, band.worst);
    }
    std::printf("       full nominal box, no physical filter : %.3f m/s"
                " (S=%.0f T=%.0f p=%.0f bar)\n", worst_box, box_s, box_t, box_p);
    std::printf("       That corner is %.0f C water under %.0f bar. Neither equation was\n",
                box_t, box_p);
    std::printf("       fitted to such a sample and they extrapolate apart; it is not a bug\n");
    std::printf("       in either, and it is a reason to distrust both there.\n");

    // These bounds are REGRESSION bounds set from the measured values above
    // with a small margin, not predictions from theory -- there is no theory
    // that says how far two empirical fits should drift apart. What they are
    // for is catching a change in either implementation.
    //
    // The measured 0.41 m/s in the top kilometre is worth reading against Chen
    // & Millero's own quoted standard deviation of 0.19 m/s: two independent
    // equations differ by about twice the uncertainty each claims. That is the
    // real accuracy of "the speed of sound in seawater", and it is why
    // docs/validation.md refuses to quote ray paths to the metre.
    PT_CHECK(bands[0].worst < 0.5);
    PT_CHECK(bands[1].worst < 0.9);
    PT_CHECK(bands[2].worst < 1.0);
    PT_CHECK(bands[3].worst < 1.5);
    // And the unphysical corner is genuinely worse, by a lot. If this stops
    // being true, one of the two implementations has changed.
    PT_CHECK(worst_box > 3.0);
}

PT_TEST(the_library_default_is_the_modern_temperature_scale) {
    // unesco() and chen_millero() must route to the ITS-90 equation. This is a
    // one-line test guarding a decision that is otherwise invisible: data
    // collected since 1990 carries ITS-90 temperatures, and feeding those to
    // the 1977 equation is a systematic error, not a rounding one.
    const Real t = 12, s = 35, p = 300;
    PT_CHECK_NEAR(static_cast<double>(sound_speed::chen_millero(t, s, p)),
                  static_cast<double>(sound_speed::chen_millero_its90(t, s, p)),
                  pt::tol(1e-12, 1e-3));
    // And the scale conversions invert.
    PT_CHECK_NEAR(static_cast<double>(sound_speed::t68_to_t90(sound_speed::t90_to_t68(t))),
                  static_cast<double>(t), pt::tol(1e-12, 1e-4));
}
