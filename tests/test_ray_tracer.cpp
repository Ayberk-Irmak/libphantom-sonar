// SPDX-License-Identifier: Apache-2.0
// Ray tracer verification.
//
// Every test here checks the tracer against a CLOSED-FORM result, not against a
// previously recorded output. Regression baselines only tell you the code did
// not change; these tell you it is right:
//
//   * isovelocity media  -> straight lines, exact geometry
//   * constant gradient  -> exact circular arcs, exact turning depth
//   * layer refinement   -> resampling a linear profile must change nothing
//   * Snell invariant    -> cos(theta)/c constant to machine precision
#include "framework.hpp"

#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

#include <array>

using namespace phantom;

namespace {

constexpr std::size_t kBuf = 4096;

// Largest deviation of cos(theta)/c from its launch value along a path.
double snell_drift(std::span<const RayPoint> path, double xi) {
    double worst = 0.0;
    for (const RayPoint& p : path) {
        const double v = std::cos(static_cast<double>(p.angle_rad))
                       / static_cast<double>(p.speed_mps);
        const double d = std::fabs(v - xi) / xi;
        if (d > worst) worst = d;
    }
    return worst;
}

}  // namespace

PT_TEST(ray_isovelocity_is_a_straight_line) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(2000, 1500);

    TraceConfig cfg;
    cfg.max_range_m = 5000;
    cfg.bottom_depth_m = 2000;

    std::array<RayPoint, kBuf> buf{};
    const Real angle = deg2rad(10);
    const TraceResult res = trace_ray(p.view(), 1000, angle, cfg, buf);

    PT_CHECK(res.status == TraceStatus::MaxRange);
    PT_CHECK(res.turning_points == 0);
    PT_CHECK(res.surface_bounces == 0 && res.bottom_bounces == 0);

    const RayPoint& end = buf[res.point_count - 1];
    const double tan10 = std::tan(10.0 * 3.14159265358979323846 / 180.0);
    const double cos10 = std::cos(10.0 * 3.14159265358979323846 / 180.0);
    PT_CHECK_NEAR(end.range_m, 5000.0, pt::tol(1e-9, 1e-2));
    PT_CHECK_NEAR(end.depth_m, 1000.0 + 5000.0 * tan10, pt::tol(1e-8, 1e-2));
    PT_CHECK_NEAR(end.time_s, (5000.0 / cos10) / 1500.0, pt::tol(1e-12, 1e-5));
    PT_CHECK_NEAR(end.angle_rad, static_cast<double>(angle), pt::tol(1e-12, 1e-6));
}

PT_TEST(ray_turning_depth_matches_snell) {
    // Linear profile: c(z) = 1500 + 0.02 z. A ray launched at 8 deg from 500 m
    // must turn exactly where c(z) = c(500) / cos(8 deg).
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(2000, 1540);

    TraceConfig cfg;
    cfg.max_range_m = 15000;
    cfg.bottom_depth_m = 2000;

    std::array<RayPoint, kBuf> buf{};
    const TraceResult res = trace_ray(p.view(), 500, deg2rad(8), cfg, buf);

    const double c0 = 1510.0;
    const double g  = 0.02;
    const double c_turn = c0 / std::cos(8.0 * 3.14159265358979323846 / 180.0);
    const double z_turn = 500.0 + (c_turn - c0) / g;

    PT_CHECK(res.turning_points == 1);
    bool seen = false;
    for (std::size_t i = 0; i < res.point_count; ++i) {
        if (std::fabs(static_cast<double>(buf[i].angle_rad)) < 1e-9) {
            PT_CHECK_NEAR(buf[i].depth_m, z_turn, pt::tol(1e-6, 5e-1));
            PT_CHECK_NEAR(buf[i].speed_mps, c_turn, pt::tol(1e-9, 1e-2));
            seen = true;
        }
    }
    PT_CHECK(seen);
    std::printf("       analytic turning depth = %.6f m\n", z_turn);
}

PT_TEST(ray_path_is_an_exact_circular_arc) {
    // The central claim of the design: in a constant-gradient layer the ray is
    // a circle of radius R = 1/(xi |g|) centred at the depth where the linear
    // sound speed extrapolates to zero. Verified point by point.
    SoundSpeedProfile<512> p;
    // 201 sample points of the SAME straight line, so the medium is unchanged
    // but the tracer is forced through 200 layer crossings.
    fill_profile(p, 0, 2000, 201, [](Real z) { return static_cast<Real>(1500) + static_cast<Real>(0.02) * z; });

    TraceConfig cfg;
    cfg.max_range_m = 15000;   // stops before the ray reaches the surface
    cfg.bottom_depth_m = 2000;

    std::array<RayPoint, kBuf> buf{};
    const Real angle = deg2rad(8);
    const TraceResult res = trace_ray(p.view(), 500, angle, cfg, buf);
    PT_CHECK(res.surface_bounces == 0 && res.bottom_bounces == 0);
    // ~75 crossings down to the vertex at 1242 m, then ~13 back up before the
    // range budget cuts the ray. Each one is an independent point on the circle.
    PT_CHECK(res.point_count > 60);

    const double g   = 0.02;
    const double c0  = 1510.0;
    const double xi  = static_cast<double>(res.snell_invariant);
    const double R   = 1.0 / (xi * std::fabs(g));
    const double z_c = 500.0 - c0 / g;                       // c extrapolates to 0
    const double r_c = 0.0 + R * std::sin(static_cast<double>(angle));  // g > 0

    double worst = 0.0;
    for (std::size_t i = 0; i < res.point_count; ++i) {
        const double dr = static_cast<double>(buf[i].range_m) - r_c;
        const double dz = static_cast<double>(buf[i].depth_m) - z_c;
        const double rel = std::fabs(std::sqrt(dr * dr + dz * dz) - R) / R;
        if (rel > worst) worst = rel;
    }
    std::printf("       arc radius %.1f m, worst radial error %.3g (relative)\n", R, worst);
    PT_CHECK(worst < pt::tol(1e-12, 1e-5));
    PT_CHECK(snell_drift(std::span<const RayPoint>(buf.data(), res.point_count), xi)
             < pt::tol(1e-13, 1e-5));
}

PT_TEST(ray_is_invariant_under_layer_refinement) {
    // Resampling a linear profile at 500x the resolution describes the exact
    // same ocean. If the arc bookkeeping is right, the answer cannot move.
    auto linear = [](Real z) { return static_cast<Real>(1500) + static_cast<Real>(0.018) * z; };

    SoundSpeedProfile<4> coarse;
    coarse.push(0, linear(0));
    coarse.push(3000, linear(3000));

    SoundSpeedProfile<2048> fine;
    fill_profile(fine, 0, 3000, 1001, linear);

    TraceConfig cfg;
    cfg.max_range_m = 40000;
    cfg.bottom_depth_m = 3000;

    std::array<RayPoint, kBuf> a{};
    std::array<RayPoint, kBuf> b{};
    const TraceResult ra = trace_ray(coarse.view(), 1000, deg2rad(6), cfg, a);
    const TraceResult rb = trace_ray(fine.view(), 1000, deg2rad(6), cfg, b);

    PT_CHECK(ra.status == TraceStatus::MaxRange);
    PT_CHECK(rb.status == TraceStatus::MaxRange);
    PT_CHECK(ra.turning_points == rb.turning_points);
    PT_CHECK(ra.surface_bounces == rb.surface_bounces);

    const RayPoint& ea = a[ra.point_count - 1];
    const RayPoint& eb = b[rb.point_count - 1];
    std::printf("       coarse: %zu pts, fine: %zu pts, depth delta = %.3g m\n",
                ra.point_count, rb.point_count,
                std::fabs(static_cast<double>(ea.depth_m - eb.depth_m)));
    PT_CHECK_NEAR(ea.range_m, eb.range_m, pt::tol(1e-9, 1e-2));
    PT_CHECK_NEAR(ea.depth_m, eb.depth_m, pt::tol(1e-7, 2.0));
    PT_CHECK_NEAR(ea.time_s, eb.time_s, pt::tol(1e-10, 1e-3));
}

PT_TEST(ray_boundary_reflections_are_counted) {
    // 45 deg in isovelocity water 1000 m deep from 500 m: pure sawtooth.
    // Bottom hits at r = 500, 2500, 4500, 6500, 8500; surface at 1500 ... 9500.
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1500);

    TraceConfig cfg;
    cfg.max_range_m = 10000;
    cfg.bottom_depth_m = 1000;

    std::array<RayPoint, kBuf> buf{};
    const TraceResult res = trace_ray(p.view(), 500, deg2rad(45), cfg, buf);

    PT_CHECK(res.bottom_bounces == 5);
    PT_CHECK(res.surface_bounces == 5);
    PT_CHECK(res.status == TraceStatus::MaxRange);

    const RayPoint& end = buf[res.point_count - 1];
    PT_CHECK_NEAR(end.range_m, 10000.0, pt::tol(1e-9, 1e-2));
    PT_CHECK_NEAR(end.depth_m, 500.0, pt::tol(1e-8, 1e-1));
    PT_CHECK_NEAR(end.time_s, (10000.0 * std::sqrt(2.0)) / 1500.0, pt::tol(1e-9, 1e-4));
}

PT_TEST(ray_absorbing_boundary_terminates) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1500);

    TraceConfig cfg;
    cfg.max_range_m = 10000;
    cfg.bottom_depth_m = 1000;
    cfg.bottom = BoundaryAction::Absorb;

    std::array<RayPoint, kBuf> buf{};
    const TraceResult res = trace_ray(p.view(), 500, deg2rad(45), cfg, buf);

    PT_CHECK(res.status == TraceStatus::Absorbed);
    PT_CHECK(res.bottom_bounces == 0);
    PT_CHECK_NEAR(buf[res.point_count - 1].range_m, 500.0, pt::tol(1e-9, 1e-2));
    PT_CHECK_NEAR(buf[res.point_count - 1].depth_m, 1000.0, pt::tol(1e-9, 1e-3));
}

PT_TEST(ray_range_budget_is_exact) {
    // Range truncation inverts the arc relation instead of overshooting and
    // clipping, so the final range lands on the budget to machine precision.
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });

    TraceConfig cfg;
    cfg.bottom_depth_m = 5000;

    std::array<RayPoint, kBuf> buf{};
    for (int a = -12; a <= 12; a += 3) {
        for (int rr = 7000; rr <= 43000; rr += 9000) {
            cfg.max_range_m = static_cast<Real>(rr);
            const TraceResult res = trace_ray(p.view(), 1300, deg2rad(static_cast<Real>(a)), cfg, buf);
            if (res.status != TraceStatus::MaxRange) continue;
            PT_CHECK_NEAR(res.final_range_m, static_cast<double>(rr), pt::tol(1e-8, 5.0));
            PT_CHECK_NEAR(buf[res.point_count - 1].range_m, static_cast<double>(rr), pt::tol(1e-8, 5.0));
        }
    }
}

PT_TEST(ray_snell_invariant_holds_in_munk_channel) {
    // The real ocean case: 501-layer Munk profile, a fan across the trapping
    // cone. cos(theta)/c must stay pinned for every ray, every point.
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });

    TraceConfig cfg;
    cfg.max_range_m = 100000;
    cfg.bottom_depth_m = 5000;

    std::array<RayPoint, kBuf> buf{};
    double worst = 0.0;
    std::size_t traced = 0;
    trace_fan(p.view(), 1300, deg2rad(-14), deg2rad(14), 29, cfg, buf,
              [&](std::size_t, Real, std::span<const RayPoint> path, const TraceResult& r) {
                  ++traced;
                  const double d = snell_drift(path, static_cast<double>(r.snell_invariant));
                  if (d > worst) worst = d;
              });
    std::printf("       %zu rays, worst Snell drift = %.3g (relative)\n", traced, worst);
    PT_CHECK(traced == 29);
    PT_CHECK(worst < pt::tol(1e-12, 1e-6));
}

PT_TEST(ray_trapped_rays_stay_inside_the_channel) {
    // A ray launched on the axis inside the trapping cone must never reach the
    // surface or the bottom; one launched outside it must.
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });

    TraceConfig cfg;
    cfg.max_range_m = 200000;
    cfg.bottom_depth_m = 5000;

    std::array<RayPoint, kBuf> buf{};
    const TraceResult trapped = trace_ray(p.view(), 1300, deg2rad(10), cfg, buf);
    PT_CHECK(trapped.surface_bounces == 0);
    PT_CHECK(trapped.bottom_bounces == 0);
    PT_CHECK(trapped.turning_points > 0);

    const TraceResult escaped = trace_ray(p.view(), 1300, deg2rad(20), cfg, buf);
    PT_CHECK(escaped.surface_bounces > 0);
}

PT_TEST(ray_rejects_degenerate_input) {
    SoundSpeedProfile<4> p;
    std::array<RayPoint, kBuf> buf{};
    TraceConfig cfg;

    // Empty profile.
    PT_CHECK(trace_ray(p.view(), 0, 0, cfg, buf).status == TraceStatus::Degenerate);

    p.push(0, 1500);
    p.push(1000, 1500);
    // Empty output span.
    PT_CHECK(trace_ray(p.view(), 0, 0, cfg, std::span<RayPoint>{}).status == TraceStatus::Degenerate);
    // Zero range budget.
    TraceConfig bad = cfg;
    bad.max_range_m = 0;
    PT_CHECK(trace_ray(p.view(), 0, 0, bad, buf).status == TraceStatus::Degenerate);
}

PT_TEST(ray_reports_buffer_exhaustion) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1500);

    TraceConfig cfg;
    cfg.max_range_m = 100000;
    cfg.bottom_depth_m = 1000;

    std::array<RayPoint, 8> tiny{};
    const TraceResult res = trace_ray(p.view(), 500, deg2rad(45), cfg, tiny);
    PT_CHECK(res.status == TraceStatus::BufferFull);
    PT_CHECK(res.point_count == tiny.size());
}

PT_TEST(ray_time_budget_is_never_exceeded) {
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });

    TraceConfig cfg;
    cfg.max_range_m = 500000;
    cfg.max_time_s = 20;
    cfg.bottom_depth_m = 5000;

    std::array<RayPoint, kBuf> buf{};
    const TraceResult res = trace_ray(p.view(), 1300, deg2rad(5), cfg, buf);
    PT_CHECK(res.status == TraceStatus::MaxTime);
    // The arc is cut exactly on the budget, not rounded back to the previous
    // layer boundary -- same guarantee the range budget gives.
    PT_CHECK_NEAR(res.final_time_s, 20.0, pt::tol(1e-12, 1e-4));
    PT_CHECK_NEAR(buf[res.point_count - 1].time_s, 20.0, pt::tol(1e-12, 1e-4));
    // 20 s of travel at ~1500 m/s cannot have gone much past 30 km.
    PT_CHECK(static_cast<double>(res.final_range_m) < 31000.0);
    PT_CHECK(static_cast<double>(res.final_range_m) > 28000.0);
}

PT_TEST(ray_time_budget_cuts_arcs_exactly) {
    // A coarse two-layer profile makes each arc kilometres long, so a naive
    // "stop before the step" policy would undershoot the budget badly.
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1520);

    TraceConfig cfg;
    cfg.max_range_m = 1000000;
    cfg.bottom_depth_m = 1000;

    std::array<RayPoint, kBuf> buf{};
    for (int a = -8; a <= 8; a += 2) {
        for (int ms = 500; ms <= 6000; ms += 1100) {
            cfg.max_time_s = static_cast<Real>(ms) / 1000;
            const TraceResult res = trace_ray(p.view(), 500, deg2rad(static_cast<Real>(a)), cfg, buf);
            PT_CHECK(res.status == TraceStatus::MaxTime);
            PT_CHECK_NEAR(res.final_time_s, static_cast<double>(cfg.max_time_s), pt::tol(1e-12, 1e-5));
            // Travel time and geometric path length must stay consistent: the
            // straight-line distance can never exceed c_max * t.
            const RayPoint& e = buf[res.point_count - 1];
            const double dist = std::sqrt(static_cast<double>(e.range_m * e.range_m)
                                        + static_cast<double>((e.depth_m - 500) * (e.depth_m - 500)));
            PT_CHECK(dist <= 1520.0 * static_cast<double>(cfg.max_time_s) + 1e-6);
        }
    }
}

PT_TEST(ray_vertical_and_horizontal_edge_cases) {
    SoundSpeedProfile<4> p;
    p.push(0, 1500);
    p.push(1000, 1520);

    TraceConfig cfg;
    cfg.max_range_m = 20000;
    cfg.max_time_s = 5;
    cfg.bottom_depth_m = 1000;

    std::array<RayPoint, kBuf> buf{};

    // Straight down: no horizontal progress, bounces between the boundaries.
    const TraceResult vert = trace_ray(p.view(), 500, deg2rad(90), cfg, buf);
    PT_CHECK_NEAR(vert.final_range_m, 0.0, pt::tol(1e-9, 1e-3));
    PT_CHECK(vert.bottom_bounces > 0);
    PT_CHECK(static_cast<double>(vert.final_time_s) <= 5.0);

    // Straight up from a source already at the surface: reflects on launch.
    const TraceResult up = trace_ray(p.view(), 0, deg2rad(-90), cfg, buf);
    PT_CHECK(up.surface_bounces >= 1);

    // Exactly horizontal launch inside a positive gradient: the ray must
    // refract upward, toward the lower sound speed. This is the case that makes
    // naive tracers stall, because sin(theta) is exactly zero. The first arc is
    // ~8.7 km long, so give it a time budget that lets it complete.
    cfg.max_time_s = 60;
    const TraceResult horiz = trace_ray(p.view(), 500, 0, cfg, buf);
    PT_CHECK(horiz.point_count >= 2);
    PT_CHECK(horiz.status != TraceStatus::Degenerate);
    PT_CHECK(static_cast<double>(buf[horiz.point_count - 1].depth_m) < 500.0);
    PT_CHECK(horiz.surface_bounces >= 1);
}
