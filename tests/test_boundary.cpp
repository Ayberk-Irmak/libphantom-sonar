// SPDX-License-Identifier: Apache-2.0
// Surface and bottom reflection losses.
//
// The Rayleigh coefficient has three closed-form limits that between them pin
// the whole formula: the critical angle from the speed ratio, the impedance
// ratio at normal incidence, and unity below critical for a lossless bottom.
// Any one of them alone could be reproduced by a wrong expression; all three
// together could not.
#include "framework.hpp"

#include "phantom/boundary.hpp"
#include "phantom/eigenray.hpp"
#include "phantom/sound_speed.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr Real kCw = 1500;   // water sound speed used throughout

BottomProperties lossless(Real c2, Real m) {
    BottomProperties b;
    b.sound_speed_mps = c2;
    b.density_ratio = m;
    b.attenuation_db_per_wavelength = 0;
    return b;
}

}  // namespace

// ---------------------------------------------------------------------------
// Bottom
// ---------------------------------------------------------------------------

PT_TEST(bottom_critical_angle_is_the_speed_ratio) {
    // arccos(c1/c2). This single number decides which rays survive to long
    // range in shallow water, so it is worth pinning exactly.
    struct Case { double c2; double expect_deg; };
    const Case cases[] = {
        {1550.0, 14.59},   // fine sand
        {1650.0, 24.62},   // medium sand
        {1800.0, 33.56},   // coarse sand / gravel
    };
    std::printf("       %10s %14s %14s\n", "c2 (m/s)", "critical", "arccos(c1/c2)");
    for (const Case& c : cases) {
        const double got = static_cast<double>(
            rad2deg(bottom_critical_angle_rad(lossless(static_cast<Real>(c.c2), 2), kCw)));
        const double expect = std::acos(1500.0 / c.c2) * 180.0 / 3.14159265358979323846;
        std::printf("       %10.0f %13.3fd %13.3fd\n", c.c2, got, expect);
        PT_CHECK_NEAR(got, expect, pt::tol(1e-9, 1e-3));
        PT_CHECK_NEAR(got, c.expect_deg, 0.02);
    }

    // A bottom slower than the water has no critical angle: it leaks at every
    // grazing angle, which is why mud is acoustically so much worse than sand.
    PT_CHECK(bottom_critical_angle_rad(lossless(1450, 15), kCw) == 0);
    PT_CHECK(bottom_critical_angle_rad(lossless(1500, 15), kCw) == 0);
    PT_CHECK(bottom_critical_angle_rad(lossless(0, 2), kCw) == 0);
    PT_CHECK(bottom_critical_angle_rad(lossless(1650, 2), 0) == 0);
}

PT_TEST(bottom_reflection_is_total_below_the_critical_angle) {
    // With no attenuation, sub-critical incidence is total internal reflection:
    // |R| = 1 exactly. That is what traps energy in a shallow duct.
    const BottomProperties b = lossless(1650, 2);
    const double crit = static_cast<double>(bottom_critical_angle_rad(b, kCw));
    std::printf("       critical angle %.3f deg\n", crit * 180.0 / 3.14159265358979323846);

    for (double frac : {0.05, 0.3, 0.6, 0.9, 0.99}) {
        const auto th = static_cast<Real>(crit * frac);
        const double r = static_cast<double>(bottom_reflection_coefficient(b, th, kCw));
        PT_CHECK_NEAR(r, 1.0, pt::tol(1e-12, 1e-5));
    }

    // Just above it, energy starts leaving.
    const double r_above = static_cast<double>(
        bottom_reflection_coefficient(b, static_cast<Real>(crit * 1.05), kCw));
    std::printf("       |R| at 1.05x critical = %.4f (%.2f dB loss)\n",
                r_above, -20.0 * std::log10(r_above));
    PT_CHECK(r_above < 0.999);
}

PT_TEST(bottom_reflection_at_normal_incidence_is_the_impedance_ratio) {
    // At 90 degrees grazing, R reduces to (Z2 - Z1)/(Z2 + Z1) with Z = rho c.
    // A different limit of the same expression from the critical-angle one, so
    // the two together constrain it tightly.
    std::printf("       %8s %8s %14s %14s\n", "c2", "m", "|R| normal", "(Z2-Z1)/(Z2+Z1)");
    for (double c2 : {1450.0, 1650.0, 2000.0}) {
        for (double m : {1.2, 1.9, 2.5}) {
            const BottomProperties b = lossless(static_cast<Real>(c2), static_cast<Real>(m));
            const double got = static_cast<double>(
                bottom_reflection_coefficient(b, kHalfPi, kCw));
            const double z_ratio = m * c2 / 1500.0;
            const double expect = std::fabs((z_ratio - 1.0) / (z_ratio + 1.0));
            std::printf("       %8.0f %8.1f %14.6f %14.6f\n", c2, m, got, expect);
            PT_CHECK_NEAR(got, expect, pt::tol(1e-10, 1e-5));
        }
    }
}

PT_TEST(bottom_reflection_never_exceeds_unity) {
    // Energy conservation: a passive interface cannot return more than it
    // receives. This is also the test that pins the complex-square-root branch
    // choice, which is otherwise a coin flip that shows up as |R| > 1.
    std::printf("       sweeping c2, density, attenuation and angle...\n");
    double worst = 0.0;
    std::size_t checked = 0;
    for (double c2 : {1400.0, 1500.0, 1600.0, 1800.0, 2200.0}) {
        for (double m : {1.0, 1.5, 2.0, 2.5}) {
            for (double att : {0.0, 0.1, 0.5, 1.0, 3.0}) {
                BottomProperties b;
                b.sound_speed_mps = static_cast<Real>(c2);
                b.density_ratio = static_cast<Real>(m);
                b.attenuation_db_per_wavelength = static_cast<Real>(att);
                for (int i = 0; i <= 90; ++i) {
                    const auto th = deg2rad(static_cast<Real>(i));
                    const double r = static_cast<double>(
                        bottom_reflection_coefficient(b, th, kCw));
                    worst = std::max(worst, r);
                    PT_CHECK(r >= 0.0);
                    ++checked;
                }
            }
        }
    }
    std::printf("       %zu combinations, max |R| = %.12f\n", checked, worst);
    PT_CHECK(worst <= 1.0 + 1e-12);
}

PT_TEST(bottom_attenuation_leaks_below_the_critical_angle) {
    // The physically important consequence of a lossy sediment: sub-critical
    // rays are no longer perfectly trapped. Without this a shallow-water duct
    // would have unbounded range, which is not what the ocean does.
    const double crit = static_cast<double>(
        bottom_critical_angle_rad(lossless(1650, static_cast<Real>(1.9)), kCw));
    const auto th = static_cast<Real>(crit * 0.5);

    std::printf("       at half the critical angle (%.2f deg):\n",
                crit * 0.5 * 180.0 / 3.14159265358979323846);
    std::printf("       %14s %12s %14s\n", "att (dB/lambda)", "|R|", "loss/bounce");
    double prev_loss = -1.0;
    for (double att : {0.0, 0.1, 0.5, 1.0, 2.0}) {
        BottomProperties b;
        b.sound_speed_mps = 1650;
        b.density_ratio = static_cast<Real>(1.9);
        b.attenuation_db_per_wavelength = static_cast<Real>(att);
        const double r = static_cast<double>(bottom_reflection_coefficient(b, th, kCw));
        const double loss = static_cast<double>(bottom_loss_db(b, th, kCw));
        std::printf("       %14.1f %12.6f %12.4f dB\n", att, r, loss);
        PT_CHECK(loss >= prev_loss - 1e-9);   // monotone in attenuation
        prev_loss = loss;
    }
    PT_CHECK(prev_loss > 0.05);   // 2 dB/wavelength is not negligible

    // Over ten bounces even a fraction of a dB compounds into something that
    // decides whether a path survives.
    BottomProperties b;
    b.sound_speed_mps = 1650;
    b.density_ratio = static_cast<Real>(1.9);
    b.attenuation_db_per_wavelength = static_cast<Real>(0.8);
    const double per_bounce = static_cast<double>(bottom_loss_db(b, th, kCw));
    std::printf("       default sand at %.2f deg: %.3f dB/bounce, %.2f dB over 10\n",
                crit * 0.5 * 180.0 / 3.14159265358979323846, per_bounce, 10 * per_bounce);
    PT_CHECK(per_bounce > 0);
}

PT_TEST(bottom_loss_is_bounded_and_grows_past_critical) {
    const BottomProperties b;   // default sand
    const double crit = static_cast<double>(bottom_critical_angle_rad(b, kCw));
    std::printf("       default sediment, critical %.2f deg\n",
                crit * 180.0 / 3.14159265358979323846);
    std::printf("       %10s %12s\n", "grazing", "loss (dB)");
    for (int deg : {2, 5, 10, 20, 24, 40, 60, 90}) {
        const double loss = static_cast<double>(
            bottom_loss_db(b, deg2rad(static_cast<Real>(deg)), kCw));
        std::printf("       %9dd %12.3f\n", deg, loss);
        PT_CHECK(loss >= 0);
        PT_CHECK(loss <= 60.0);
    }
    // A steep ray loses far more per bounce than a grazing one -- which is why
    // steep paths die out and grazing paths carry shallow-water range.
    PT_CHECK(bottom_loss_db(b, deg2rad(static_cast<Real>(60)), kCw)
           > bottom_loss_db(b, deg2rad(static_cast<Real>(5)), kCw));
}

// ---------------------------------------------------------------------------
// Surface
// ---------------------------------------------------------------------------

PT_TEST(surface_is_perfect_when_flat) {
    SurfaceProperties s;
    s.rms_wave_height_m = 0;
    for (int deg : {1, 10, 45, 90}) {
        PT_CHECK_NEAR(surface_reflection_coefficient(s, deg2rad(static_cast<Real>(deg)),
                                                     10000, kCw), 1.0, 1e-12);
        PT_CHECK_NEAR(surface_loss_db(s, deg2rad(static_cast<Real>(deg)), 10000, kCw),
                      0.0, 1e-12);
    }
    // Degenerate inputs fall back to perfect reflection rather than a NaN.
    s.rms_wave_height_m = 1;
    PT_CHECK_NEAR(surface_reflection_coefficient(s, deg2rad(static_cast<Real>(10)), 0, kCw),
                  1.0, 1e-12);
    PT_CHECK_NEAR(surface_reflection_coefficient(s, deg2rad(static_cast<Real>(10)), 10000, 0),
                  1.0, 1e-12);
}

PT_TEST(surface_roughness_loss_follows_the_rayleigh_parameter) {
    // |R| = exp(-G^2/2) with G = 2 k sigma sin(theta). Checked against the
    // closed form across frequency, angle and wave height at once.
    SurfaceProperties s;
    s.rms_wave_height_m = static_cast<Real>(0.2);
    s.max_loss_db = 1000;   // uncapped, so the formula itself is under test

    std::printf("       sigma = 0.2 m\n");
    std::printf("       %10s %10s %10s %12s %12s\n",
                "f (Hz)", "graze", "G", "|R|", "exp(-G^2/2)");
    for (double f : {500.0, 2000.0, 8000.0}) {
        for (int deg : {2, 10, 30}) {
            const double th = static_cast<double>(deg) * 3.14159265358979323846 / 180.0;
            const double k = 2.0 * 3.14159265358979323846 * f / 1500.0;
            const double gamma = 2.0 * k * 0.2 * std::sin(th);
            const double expect = std::exp(-gamma * gamma / 2.0);
            const double got = static_cast<double>(surface_reflection_coefficient(
                s, static_cast<Real>(th), static_cast<Real>(f), kCw));
            if (f < 3000.0 || deg < 20) {
                std::printf("       %10.0f %9dd %10.3f %12.3e %12.3e\n",
                            f, deg, gamma, got, expect);
            }
            PT_CHECK_NEAR(got / std::max(expect, 1e-30), 1.0, pt::tol(1e-9, 1e-3));
        }
    }
}

PT_TEST(surface_loss_is_capped_because_the_diffuse_field_is_not_modelled) {
    // The coherent coefficient falls without limit, but the scattered energy
    // does not leave the ocean -- it goes into a diffuse field a ray model does
    // not carry. Reporting a 200 dB loss would be arithmetic, not physics.
    SurfaceProperties s;
    s.rms_wave_height_m = static_cast<Real>(0.5);
    s.max_loss_db = 30;

    const double raw = -20.0 * std::log10(static_cast<double>(
        surface_reflection_coefficient(s, deg2rad(static_cast<Real>(20)), 10000, kCw)));
    const double capped = static_cast<double>(
        surface_loss_db(s, deg2rad(static_cast<Real>(20)), 10000, kCw));
    std::printf("       0.5 m seas, 10 kHz, 20 deg: coherent loss %.0f dB, reported %.0f dB\n",
                raw, capped);
    PT_CHECK(raw > 100.0);
    PT_CHECK_NEAR(capped, 30.0, 1e-9);

    // Below the cap the value passes through untouched.
    s.rms_wave_height_m = static_cast<Real>(0.02);
    const double small = static_cast<double>(
        surface_loss_db(s, deg2rad(static_cast<Real>(5)), 2000, kCw));
    std::printf("       0.02 m ripple, 2 kHz, 5 deg: %.4f dB\n", small);
    PT_CHECK(small > 0);
    PT_CHECK(small < 30.0);
}

PT_TEST(wind_to_wave_height_is_pierson_moskowitz) {
    // H_1/3 = 0.0246 U^2, sigma = H_1/3 / 4.
    std::printf("       %10s %14s %14s\n", "wind (m/s)", "H_1/3 (m)", "sigma (m)");
    for (double u : {5.0, 10.0, 15.0, 20.0}) {
        const double sigma = static_cast<double>(
            wind_to_rms_wave_height_m(static_cast<Real>(u)));
        const double expect = 0.0246 * u * u / 4.0;
        std::printf("       %10.0f %14.3f %14.4f\n", u, 0.0246 * u * u, sigma);
        PT_CHECK_NEAR(sigma, expect, pt::tol(1e-12, 1e-6));
    }
    // Sea state 3 (about 10 m/s) gives roughly half a metre of RMS height,
    // which is the regime where surface loss starts to dominate at 10 kHz.
    PT_CHECK_NEAR(wind_to_rms_wave_height_m(10), 0.615, 0.01);
    PT_CHECK(wind_to_rms_wave_height_m(0) == 0);
    PT_CHECK(wind_to_rms_wave_height_m(-5) == 0);
}

// ---------------------------------------------------------------------------
// Integration with the eigenray path
// ---------------------------------------------------------------------------

PT_TEST(boundary_loss_uses_the_snell_invariant_angles) {
    // The grazing angle at each bounce is recovered from xi rather than tracked:
    // cos(theta) = xi * c at the boundary depth. Checked against a hand
    // computation for a known ray.
    Eigenray e;
    e.surface_bounces = 2;
    e.bottom_bounces = 3;
    // A ray whose grazing angle at 1500 m/s water is 10 degrees.
    const double th_water = 10.0 * 3.14159265358979323846 / 180.0;
    e.snell_invariant = static_cast<Real>(std::cos(th_water) / 1500.0);

    BoundaryModel model;
    model.surface.rms_wave_height_m = static_cast<Real>(0.1);
    model.bottom = BottomProperties{};

    const Real f = 5000;
    const double got = static_cast<double>(boundary_loss_db(e, model, 1500, 1500, f));

    const double expect =
        2.0 * static_cast<double>(surface_loss_db(model.surface, static_cast<Real>(th_water),
                                                  f, 1500))
      + 3.0 * static_cast<double>(bottom_loss_db(model.bottom, static_cast<Real>(th_water), 1500));
    std::printf("       2 surface + 3 bottom bounces at 10 deg: %.4f dB (expected %.4f)\n",
                got, expect);
    PT_CHECK_NEAR(got, expect, pt::tol(1e-9, 1e-4));

    // A ray that never touches a boundary pays nothing.
    Eigenray direct;
    direct.snell_invariant = e.snell_invariant;
    PT_CHECK(boundary_loss_db(direct, model, 1500, 1500, f) == 0);
    // An invalid invariant is not guessed at.
    Eigenray bad = e;
    bad.snell_invariant = 0;
    PT_CHECK(boundary_loss_db(bad, model, 1500, 1500, f) == 0);
}

PT_TEST(bounced_paths_are_now_quieter_than_the_direct_one) {
    // The point of the release. Before boundary losses existed, a path with
    // four bounces was reported at essentially the same level as the direct
    // path, which made shallow-water multipath look far too healthy.
    SoundSpeedProfile<64> p;
    fill_profile(p, 0, 200, 21, [](Real z) {
        return sound_speed::mackenzie(static_cast<Real>(15) - static_cast<Real>(0.01) * z, 35, z);
    });

    TraceConfig cfg;
    cfg.max_range_m = 20000;
    cfg.max_time_s = 60;
    cfg.bottom_depth_m = 200;

    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-25));
    search.angle_max_rad = deg2rad(static_cast<Real>(25));
    search.fan_count = 1201;

    static std::array<RayPoint, 8192> scratch;
    static std::array<Eigenray, 32> rays;
    const Real range = 3000;
    const std::size_t n = find_eigenrays(p.view(), 50, 120, range, cfg, search, scratch, rays);
    PT_CHECK(n >= 4);

    BoundaryModel model;
    model.surface.rms_wave_height_m = wind_to_rms_wave_height_m(8);   // ~8 m/s wind
    const Real f = 5000;
    const Real c_surface = speed_at(p.view(), 0);
    const Real c_bottom = speed_at(p.view(), 200);

    std::printf("       8 m/s wind (sigma %.3f m), 5 kHz, sand bottom\n",
                static_cast<double>(model.surface.rms_wave_height_m));
    std::printf("       %6s %6s %10s %12s %12s\n",
                "srf", "btm", "graze", "no bounds", "with bounds");

    double direct_total = -1;
    double worst_bounced = -1;
    std::size_t bounced = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (rays[i].near_caustic) continue;
        const double plain = static_cast<double>(
            transmission_loss_db(rays[i], range, f, speed_at(p.view(), 50)));
        const double total = static_cast<double>(total_transmission_loss_db(
            rays[i], range, f, speed_at(p.view(), 50), model, c_surface, c_bottom));
        const double graze = std::acos(std::min(1.0,
            static_cast<double>(rays[i].snell_invariant) * static_cast<double>(c_bottom)))
            * 180.0 / 3.14159265358979323846;
        std::printf("       %6u %6u %9.2fd %12.2f %12.2f\n",
                    rays[i].surface_bounces, rays[i].bottom_bounces, graze, plain, total);

        PT_CHECK(total >= plain - 1e-9);   // boundaries can only cost
        if (rays[i].surface_bounces + rays[i].bottom_bounces == 0) {
            direct_total = total;
            PT_CHECK_NEAR(total, plain, pt::tol(1e-9, 1e-4));   // and cost nothing here
        } else {
            ++bounced;
            worst_bounced = std::max(worst_bounced, total);
        }
    }
    PT_CHECK(bounced > 0);
    PT_CHECK(direct_total > 0);
    std::printf("       direct path %.2f dB, worst bounced %.2f dB\n",
                direct_total, worst_bounced);
    // Every bounced path now costs something the direct one does not.
    PT_CHECK(worst_bounced > direct_total);
}
