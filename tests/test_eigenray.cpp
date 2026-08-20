// SPDX-License-Identifier: Apache-2.0
// Eigenray search, transmission loss, and multipath echo synthesis.
//
// The load-bearing test here is the isovelocity one: in a uniform ocean with no
// boundaries, the eigenray is a straight line and the ray-tube spreading loss
// must reduce EXACTLY to 20 log10(R) with R the slant range. That is a closed
// form for a quantity computed through a Jacobian obtained by finite difference
// through the tracer, so it exercises the whole chain at once.
#include "framework.hpp"

#include "phantom/echo_synth.hpp"
#include "phantom/eigenray.hpp"
#include "phantom/sound_speed.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using namespace phantom;

namespace {

constexpr std::size_t kScratch = 8192;
std::array<RayPoint, kScratch> g_scratch;
std::array<Eigenray, 32>       g_rays;

TraceConfig deep_cfg(Real bottom_m) {
    TraceConfig cfg;
    cfg.max_range_m = 200000;
    cfg.max_time_s = 500;
    cfg.bottom_depth_m = bottom_m;
    return cfg;
}

}  // namespace

// ---------------------------------------------------------------------------
// Absorption
// ---------------------------------------------------------------------------

PT_TEST(thorp_absorption_matches_published_values) {
    // Thorp (1967). The published figures are the reason sonar designers pick
    // the band they do, so the decade behaviour matters more than any one point.
    // Tolerances reflect what Thorp actually is: a fit to measurements around
    // 4 C, quoted in the literature to about 10%. Values outside that band
    // would mean a mistyped coefficient; values inside it mean nothing more
    // than that the fit is a fit.
    struct Case { double f_hz; double expect_db_km; double tol; };
    const Case cases[] = {
        {100.0,    0.0042, 0.001},
        {1000.0,   0.069,  0.007},
        {10000.0,  1.19,   0.12},
        {50000.0,  17.5,   1.8},
        {100000.0, 34.1,   3.5},
    };
    std::printf("       %10s %14s %14s\n", "f (Hz)", "a (dB/km)", "published");
    for (const Case& c : cases) {
        const double got = static_cast<double>(
            thorp_absorption_db_per_km(static_cast<Real>(c.f_hz)));
        std::printf("       %10.0f %14.4f %14.4f\n", c.f_hz, got, c.expect_db_km);
        PT_CHECK_NEAR(got, c.expect_db_km, c.tol);
    }

    double prev = 0;
    for (int i = 2; i <= 60; ++i) {
        const double a = static_cast<double>(
            thorp_absorption_db_per_km(static_cast<Real>(i) * static_cast<Real>(500)));
        PT_CHECK(a > prev);
        prev = a;
    }
    PT_CHECK(thorp_absorption_db_per_km(0) == 0);
    PT_CHECK(thorp_absorption_db_per_km(-1) == 0);
}

// ---------------------------------------------------------------------------
// Eigenrays and spreading
// ---------------------------------------------------------------------------

PT_TEST(eigenray_isovelocity_is_the_straight_line) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(4000, 1500);

    const Real z_src = 1000, z_rcv = 1600, r = 5000;
    TraceConfig cfg = deep_cfg(4000);
    cfg.surface = BoundaryAction::Absorb;
    cfg.bottom = BoundaryAction::Absorb;

    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-40));
    search.angle_max_rad = deg2rad(static_cast<Real>(40));

    const std::size_t n = find_eigenrays(p.view(), z_src, z_rcv, r, cfg, search,
                                         g_scratch, g_rays);
    PT_CHECK(n == 1);
    if (n == 0) return;

    const double expect_angle = std::atan2(600.0, 5000.0);
    const double slant = std::sqrt(5000.0 * 5000.0 + 600.0 * 600.0);
    std::printf("       launch %.6f deg (geometric %.6f), path %.3f m (slant %.3f)\n",
                static_cast<double>(rad2deg(g_rays[0].launch_angle_rad)),
                expect_angle * 180.0 / 3.14159265358979323846,
                static_cast<double>(g_rays[0].path_length_m), slant);

    PT_CHECK_NEAR(g_rays[0].launch_angle_rad, expect_angle, pt::tol(1e-6, 1e-4));
    PT_CHECK_NEAR(g_rays[0].arrival_angle_rad, expect_angle, pt::tol(1e-6, 1e-4));
    PT_CHECK_REL(g_rays[0].path_length_m, slant, pt::tol(1e-7, 1e-4));
    PT_CHECK_REL(g_rays[0].travel_time_s, slant / 1500.0, pt::tol(1e-7, 1e-4));
    PT_CHECK(g_rays[0].surface_bounces == 0);
    PT_CHECK(g_rays[0].bottom_bounces == 0);

    // dz/dtheta_0 for a straight line at range r is r / cos^2(theta_0).
    const double expect_jac = 5000.0 / (std::cos(expect_angle) * std::cos(expect_angle));
    std::printf("       Jacobian %.3f m/rad (geometric %.3f)\n",
                static_cast<double>(g_rays[0].jacobian_m_per_rad), expect_jac);
    PT_CHECK_REL(g_rays[0].jacobian_m_per_rad, expect_jac, pt::tol(1e-4, 1e-2));
}

PT_TEST(spreading_loss_reduces_to_spherical_in_isovelocity_water) {
    // The closed form the whole module rests on: with no refraction the
    // ray-tube formula must give exactly 20 log10(slant range).
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(6000, 1500);

    TraceConfig cfg = deep_cfg(6000);
    cfg.surface = BoundaryAction::Absorb;
    cfg.bottom = BoundaryAction::Absorb;

    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-50));
    search.angle_max_rad = deg2rad(static_cast<Real>(50));

    std::printf("       %8s %8s %14s %14s %10s\n",
                "r (m)", "dz (m)", "TL ray-tube", "20log10(R)", "diff dB");
    double worst = 0.0;
    // Geometries chosen so every one is actually reachable: the ray must stay
    // inside the +/-50 degree fan AND the receiver must stay in the water. A
    // depth offset scaled to range fails both at the far end.
    struct Geom { double r; double dz; };
    const Geom cases[] = {
        {1000.0, 0.0},   {1000.0, 200.0},  {1000.0, 500.0},
        {5000.0, 0.0},   {5000.0, 1000.0}, {5000.0, 2500.0},
        {20000.0, 0.0},  {20000.0, 1000.0}, {20000.0, 3000.0},
    };
    for (const Geom& gm : cases) {
        {
            const double r = gm.r;
            const double dz = gm.dz;
            const Real z_src = 2500;
            const auto z_rcv = static_cast<Real>(2500.0 + dz);
            const std::size_t n = find_eigenrays(p.view(), z_src, z_rcv,
                                                 static_cast<Real>(r), cfg, search,
                                                 g_scratch, g_rays);
            if (n != 1) { PT_CHECK(n == 1); continue; }

            const double slant = std::sqrt(r * r + dz * dz);
            const double expect = 20.0 * std::log10(slant);
            const double got = static_cast<double>(g_rays[0].spreading_loss_db);
            std::printf("       %8.0f %8.0f %14.4f %14.4f %10.4f\n",
                        r, dz, got, expect, got - expect);
            worst = std::max(worst, std::fabs(got - expect));
        }
    }
    std::printf("       worst departure from spherical spreading: %.5f dB\n", worst);
    PT_CHECK(worst < pt::tol(0.01, 0.15));
}

PT_TEST(spreading_loss_closed_form_direct) {
    // The same formula fed analytic inputs rather than traced ones, so a bug in
    // the search cannot mask a bug in the loss, or the reverse.
    for (double deg : {0.0, 10.0, 30.0}) {
        const double th = deg * 3.14159265358979323846 / 180.0;
        const double r = 3000.0;
        const double jac = r / (std::cos(th) * std::cos(th));
        const double got = static_cast<double>(
            spreading_loss_db(static_cast<Real>(r), static_cast<Real>(th),
                              static_cast<Real>(th), static_cast<Real>(jac), 1500, 1500));
        const double expect = 20.0 * std::log10(r / std::cos(th));
        PT_CHECK_NEAR(got, expect, pt::tol(1e-9, 1e-4));
    }
    // Degenerate inputs return 0 rather than an infinity.
    PT_CHECK(spreading_loss_db(0, 0, 0, 1000, 1500, 1500) == 0);
    PT_CHECK(spreading_loss_db(1000, 0, 0, 0, 1500, 1500) == 0);           // caustic
    // cos(pi/2) is 6e-17, not 0, so this needs an epsilon inside the function
    // rather than a bare > 0 test -- otherwise it returns about +330 dB.
    PT_CHECK(spreading_loss_db(1000, 0, kHalfPi, 1000, 1500, 1500) == 0);  // vertical
    PT_CHECK(spreading_loss_db(1000, 0, 0, 1000, 0, 1500) == 0);
}

PT_TEST(eigenray_finds_the_multipath_structure) {
    // A shallow duct with reflecting boundaries: the direct path plus surface-
    // and bottom-reflected ones.
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

    const Real z_src = 50, z_rcv = 120, r = 3000;
    const std::size_t n = find_eigenrays(p.view(), z_src, z_rcv, r, cfg, search,
                                         g_scratch, g_rays);
    std::printf("       %zu eigenrays over 3 km in a 200 m duct\n", n);
    std::printf("       %10s %10s %10s %6s %6s %10s\n",
                "launch", "arrival", "t (ms)", "srf", "btm", "TL (dB)");
    for (std::size_t i = 0; i < n; ++i) {
        std::printf("       %9.3f%s %9.3f%s %10.3f %6u %6u %10.2f%s\n",
                    static_cast<double>(rad2deg(g_rays[i].launch_angle_rad)), "d",
                    static_cast<double>(rad2deg(g_rays[i].arrival_angle_rad)), "d",
                    static_cast<double>(g_rays[i].travel_time_s) * 1e3,
                    g_rays[i].surface_bounces, g_rays[i].bottom_bounces,
                    static_cast<double>(g_rays[i].spreading_loss_db),
                    g_rays[i].near_caustic ? "  [caustic]" : "");
    }
    PT_CHECK(n >= 4);

    const double straight = std::sqrt(3000.0 * 3000.0 + 70.0 * 70.0);
    for (std::size_t i = 0; i < n; ++i) {
        PT_CHECK(static_cast<double>(g_rays[i].path_length_m) >= straight - 1.0);
        PT_CHECK(static_cast<double>(g_rays[i].travel_time_s) > straight / 1520.0 - 1e-3);
        PT_CHECK(static_cast<double>(g_rays[i].travel_time_s) < 3.0);
        if (g_rays[i].surface_bounces + g_rays[i].bottom_bounces > 0) {
            PT_CHECK(static_cast<double>(g_rays[i].path_length_m) > straight);
        }
    }

    // Distinct, increasing launch angles: two "eigenrays" at the same angle
    // would be one root found twice.
    for (std::size_t i = 1; i < n; ++i) {
        PT_CHECK(g_rays[i].launch_angle_rad > g_rays[i - 1].launch_angle_rad);
        PT_CHECK(static_cast<double>(g_rays[i].launch_angle_rad)
               - static_cast<double>(g_rays[i - 1].launch_angle_rad) > 1e-6);
    }
}

PT_TEST(eigenray_search_rejects_the_unreachable) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1500);

    TraceConfig cfg = deep_cfg(1000);
    cfg.surface = BoundaryAction::Absorb;
    cfg.bottom = BoundaryAction::Absorb;
    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-2));
    search.angle_max_rad = deg2rad(static_cast<Real>(2));

    // +/-2 degrees over 1 km spans only +/-35 m; the receiver is 400 m below.
    PT_CHECK(find_eigenrays(p.view(), 500, 900, 1000, cfg, search, g_scratch, g_rays) == 0);

    PT_CHECK(find_eigenrays(p.view(), 500, 520, 0, cfg, search, g_scratch, g_rays) == 0);
    PT_CHECK(find_eigenrays(p.view(), 500, 520, 1000, cfg, search,
                            std::span<RayPoint>(), g_rays) == 0);
    PT_CHECK(find_eigenrays(p.view(), 500, 520, 1000, cfg, search, g_scratch,
                            std::span<Eigenray>()) == 0);
    EigenraySearch bad = search;
    bad.fan_count = 1;
    PT_CHECK(find_eigenrays(p.view(), 500, 520, 1000, cfg, bad, g_scratch, g_rays) == 0);
    bad = search;
    bad.angle_max_rad = bad.angle_min_rad;
    PT_CHECK(find_eigenrays(p.view(), 500, 520, 1000, cfg, bad, g_scratch, g_rays) == 0);
}

PT_TEST(transmission_loss_adds_absorption_over_the_path) {
    // TL = spreading + a(f) * path. At 100 kHz over 3 km the absorption term is
    // a hundred dB and dominates the geometry completely, which is the whole
    // reason high-frequency sonar is short-range.
    Eigenray e;
    e.launch_angle_rad = 0;
    e.arrival_angle_rad = 0;
    e.jacobian_m_per_rad = 3000;
    e.arrival_speed_mps = 1500;
    e.path_length_m = 3000;

    const double spread = 20.0 * std::log10(3000.0);
    std::printf("       spreading alone: %.2f dB over 3 km\n", spread);
    std::printf("       %10s %12s %12s\n", "f (Hz)", "absorption", "total TL");
    for (double f : {1000.0, 10000.0, 100000.0}) {
        const double tl = static_cast<double>(
            transmission_loss_db(e, 3000, static_cast<Real>(f), 1500));
        const double absorb = static_cast<double>(
            thorp_absorption_db_per_km(static_cast<Real>(f))) * 3.0;
        std::printf("       %10.0f %12.2f %12.2f\n", f, absorb, tl);
        PT_CHECK_NEAR(tl, spread + absorb, pt::tol(1e-6, 1e-3));
    }
    PT_CHECK(static_cast<double>(transmission_loss_db(e, 3000, 100000, 1500)) > 2 * spread);
    PT_CHECK_NEAR(transmission_loss_db(e, 3000, 1000, 1500), spread, 0.3);
}

PT_TEST(caustics_are_flagged_not_reported_as_levels) {
    // Where neighbouring rays cross, the tube collapses and ray theory predicts
    // infinite intensity. That is a failure of the method, not a property of
    // the ocean, so it is surfaced rather than returned as a number.
    PT_CHECK(spreading_loss_db(1000, 0, 0, 0, 1500, 1500) == 0);

    SoundSpeedProfile<128> p;
    fill_profile(p, 0, 5000, 101, [](Real z) { return sound_speed::munk(z); });
    TraceConfig cfg = deep_cfg(5000);
    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-14));
    search.angle_max_rad = deg2rad(static_cast<Real>(14));
    search.fan_count = 1401;
    search.caustic_jacobian_m_per_rad = 200;   // generous, to catch some

    const std::size_t n = find_eigenrays(p.view(), 1300, 1500, 40000, cfg, search,
                                         g_scratch, g_rays);
    std::size_t flagged = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (g_rays[i].near_caustic) {
            ++flagged;
            PT_CHECK(g_rays[i].spreading_loss_db == 0);
        } else {
            PT_CHECK(g_rays[i].spreading_loss_db > 0);
        }
    }
    std::printf("       %zu eigenrays at 40 km in the Munk channel, %zu near a caustic\n",
                n, flagged);
    PT_CHECK(n > 0);
}

// ---------------------------------------------------------------------------
// Multipath echo synthesis
// ---------------------------------------------------------------------------

PT_TEST(multipath_echoes_come_from_the_traced_paths) {
    // The coupling this release is about: the arrival structure of the reply is
    // computed from the ocean rather than chosen by hand.
    SoundSpeedProfile<64> p;
    fill_profile(p, 0, 200, 21, [](Real z) {
        return sound_speed::mackenzie(static_cast<Real>(15) - static_cast<Real>(0.01) * z, 35, z);
    });

    TraceConfig cfg;
    cfg.max_range_m = 20000;
    cfg.max_time_s = 60;
    cfg.bottom_depth_m = 200;

    EigenraySearch search;
    // A narrower fan than the structure test: a real projector has a limited
    // vertical beam, and steep paths travel far enough further to spread the
    // arrival structure over hundreds of milliseconds.
    search.angle_min_rad = deg2rad(static_cast<Real>(-10));
    search.angle_max_rad = deg2rad(static_cast<Real>(10));
    search.fan_count = 1201;

    const Real range = 3000;
    const std::size_t n = find_eigenrays(p.view(), 50, 120, range, cfg, search,
                                         g_scratch, g_rays);
    PT_CHECK(n >= 3);

    std::array<EchoSpec, 16> echoes{};
    const std::size_t m = echoes_from_eigenrays(
        std::span<const Eigenray>(g_rays.data(), n), range, 12000, 1510,
        static_cast<Real>(-3), echoes);
    std::printf("       %zu paths -> %zu echoes\n", n, m);
    PT_CHECK(m >= 3);

    std::printf("       %14s %14s\n", "delay (ms)", "level (dB)");
    bool zero_seen = false;
    double strongest = -1e9;
    for (std::size_t i = 0; i < m; ++i) {
        std::printf("       %14.4f %14.2f\n",
                    static_cast<double>(echoes[i].extra_delay_s) * 1e3,
                    static_cast<double>(echoes[i].target_strength_db));
        PT_CHECK(echoes[i].extra_delay_s >= 0);
        if (echoes[i].extra_delay_s == 0) zero_seen = true;
        // Multipath delay is a time, not a pretended range.
        PT_CHECK(echoes[i].range_offset_m == 0);
        strongest = std::max(strongest, static_cast<double>(echoes[i].target_strength_db));
    }
    PT_CHECK(zero_seen);
    PT_CHECK_NEAR(strongest, -3.0, pt::tol(1e-6, 1e-4));
    for (std::size_t i = 0; i < m; ++i) {
        PT_CHECK(static_cast<double>(echoes[i].target_strength_db) <= strongest + 1e-6);
    }

    double spread_ms = 0;
    for (std::size_t i = 0; i < m; ++i) {
        spread_ms = std::max(spread_ms, static_cast<double>(echoes[i].extra_delay_s) * 1e3);
    }
    // The number that matters operationally is the spread relative to the
    // pulse length: comparable to it and the paths smear one echo, much larger
    // and they resolve as separate targets.
    std::printf("       multipath spread %.2f ms two-way over a %.0f m path\n",
                spread_ms, static_cast<double>(range));
    std::printf("       (a 20 ms pulse would be smeared, not resolved, at this spread)\n");
    PT_CHECK(spread_ms > 0);
    PT_CHECK(spread_ms < 200.0);

    PT_CHECK(echoes_from_eigenrays(std::span<const Eigenray>(), range, 12000, 1510, 0, echoes) == 0);
    PT_CHECK(echoes_from_eigenrays(std::span<const Eigenray>(g_rays.data(), n), range,
                                   12000, 1510, 0, std::span<EchoSpec>()) == 0);
}

PT_TEST(multipath_echo_delays_survive_synthesis) {
    // Round trip: build echoes from paths, render them, and check the rendered
    // stream has energy at each expected delay. This is what ties the ray tracer
    // to the transmitted waveform.
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(200, 1500);

    TraceConfig cfg;
    cfg.max_range_m = 20000;
    cfg.max_time_s = 60;
    cfg.bottom_depth_m = 200;

    EigenraySearch search;
    search.angle_min_rad = deg2rad(static_cast<Real>(-20));
    search.angle_max_rad = deg2rad(static_cast<Real>(20));
    search.fan_count = 901;

    const Real range = 2000;
    const std::size_t n = find_eigenrays(p.view(), 50, 120, range, cfg, search,
                                         g_scratch, g_rays);
    PT_CHECK(n >= 2);

    std::array<EchoSpec, 16> echoes{};
    const std::size_t m = echoes_from_eigenrays(
        std::span<const Eigenray>(g_rays.data(), n), range, 12000, 1500, 0, echoes);
    PT_CHECK(m >= 2);

    PulseDescriptor pdw;
    pdw.type = PulseType::LfmUp;
    pdw.centre_freq_hz = 14000;
    pdw.bandwidth_hz = 12000;
    pdw.duration_s = static_cast<Real>(0.004);
    pdw.amplitude = 1;

    constexpr Real kFs = 96000;
    static std::array<Real, 65536> out{};
    const std::size_t need = swarm_length(pdw, std::span<const EchoSpec>(echoes.data(), m),
                                          kFs, 1500);
    PT_CHECK(need > 0);
    PT_CHECK(need <= out.size());
    const std::size_t written = synthesize_swarm(pdw, std::span<const EchoSpec>(echoes.data(), m),
                                                 kFs, 1500, out);
    PT_CHECK(written == need);

    for (std::size_t i = 0; i < m; ++i) {
        const auto at = static_cast<std::size_t>(
            static_cast<double>(echoes[i].extra_delay_s) * static_cast<double>(kFs));
        double peak = 0;
        for (std::size_t k = at; k < at + 64 && k < written; ++k) {
            peak = std::max(peak, std::fabs(static_cast<double>(out[k])));
        }
        PT_CHECK(peak > 0.05);
    }
    std::printf("       %zu paths rendered into %zu samples (%.2f ms)\n",
                m, written, static_cast<double>(written) / 96000.0 * 1e3);
}
