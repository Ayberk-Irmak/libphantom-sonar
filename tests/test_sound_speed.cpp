// SPDX-License-Identifier: Apache-2.0
// Sound speed equations.
//
// There is no way to unit-test a 40-coefficient polynomial against itself, so
// the strategy is cross-validation: three equations derived independently from
// different data sets must agree inside their common validity box. A single
// mistyped coefficient shows up immediately as a broken agreement.
#include "framework.hpp"

#include "phantom/sound_speed.hpp"

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
