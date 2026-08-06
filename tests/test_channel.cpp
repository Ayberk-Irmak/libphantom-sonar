// SPDX-License-Identifier: Apache-2.0
#include "framework.hpp"

#include "phantom/channel.hpp"
#include "phantom/coverage.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

#include <array>

using namespace phantom;

PT_TEST(sofar_axis_found_in_munk_profile) {
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });
    const ChannelInfo info = analyze_sofar(p.view());

    PT_CHECK(info.found);
    PT_CHECK_NEAR(info.axis_depth_m, 1300.0, 10.0);
    PT_CHECK_NEAR(info.axis_speed_mps, 1500.0, 1e-6);

    // The Munk profile is capped by the surface, not the sea floor.
    PT_CHECK(info.limited_by_surface);
    PT_CHECK_NEAR(info.limiting_speed_mps, 1548.52, 0.05);

    // Trapping cone half-angle: acos(1500 / 1548.52) = 14.37 deg.
    const double deg = static_cast<double>(rad2deg(info.max_trapped_angle_rad));
    std::printf("       axis %.1f m, trapping cone +/- %.3f deg, conjugates %.0f..%.0f m\n",
                static_cast<double>(info.axis_depth_m), deg,
                static_cast<double>(info.upper_conjugate_m),
                static_cast<double>(info.lower_conjugate_m));
    PT_CHECK_NEAR(deg, 14.37, 0.15);

    // Upper conjugate is the surface itself; lower sits near 4800 m.
    PT_CHECK_NEAR(info.upper_conjugate_m, 0.0, 5.0);
    PT_CHECK_NEAR(info.lower_conjugate_m, 4800.0, 40.0);
}

PT_TEST(sofar_trapping_angle_predicts_ray_behaviour) {
    // The whole point of the analysis: the reported cone must actually match
    // what the ray tracer does. Just inside -> trapped. Just outside -> escapes.
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });
    const ChannelInfo info = analyze_sofar(p.view());
    const Real cone = info.max_trapped_angle_rad;

    TraceConfig cfg;
    cfg.max_range_m = 250000;
    cfg.bottom_depth_m = 5000;

    std::array<RayPoint, 8192> buf{};
    const TraceResult inside  = trace_ray(p.view(), info.axis_depth_m,
                                          cone - deg2rad(static_cast<Real>(0.5)), cfg, buf);
    const TraceResult outside = trace_ray(p.view(), info.axis_depth_m,
                                          cone + deg2rad(static_cast<Real>(0.5)), cfg, buf);

    PT_CHECK(inside.surface_bounces == 0);
    PT_CHECK(inside.bottom_bounces == 0);
    PT_CHECK(outside.surface_bounces > 0);
}

PT_TEST(sofar_rejects_monotonic_profiles) {
    // A profile whose minimum is at an endpoint has no duct: there is no upper
    // wall to refract energy back down.
    SoundSpeedProfile<8> rising;
    rising.push(0, 1500);
    rising.push(1000, 1520);
    rising.push(2000, 1540);
    PT_CHECK(!analyze_sofar(rising.view()).found);

    SoundSpeedProfile<8> falling;
    falling.push(0, 1540);
    falling.push(1000, 1520);
    falling.push(2000, 1500);
    PT_CHECK(!analyze_sofar(falling.view()).found);
}

PT_TEST(surface_duct_detection) {
    // Warm mixed layer over a thermocline: speed rises to 80 m, then falls.
    SoundSpeedProfile<16> p;
    p.push(0, 1500);
    p.push(80, 1508);
    p.push(300, 1490);
    p.push(1500, 1485);
    p.push(3000, 1520);
    PT_CHECK_NEAR(surface_duct_depth(p.view()), 80.0, 1e-9);

    SoundSpeedProfile<8> nodust;
    nodust.push(0, 1520);
    nodust.push(500, 1490);
    nodust.push(3000, 1520);
    PT_CHECK(surface_duct_depth(nodust.view()) < 0);
}

PT_TEST(critical_angle_matches_snell) {
    SoundSpeedProfile<8> p;
    p.push(0, 1500);
    p.push(1000, 1520);
    const ProfileView v = p.view();
    // A ray from 0 m turning at 1000 m: cos(theta) = 1500/1520.
    PT_CHECK_NEAR(critical_angle(v, 0, 1000), std::acos(1500.0 / 1520.0), pt::tol(1e-12, 1e-6));
    // No turning possible when the sound speed does not increase.
    PT_CHECK_NEAR(critical_angle(v, 1000, 0), 0.0, 1e-12);
}

PT_TEST(coverage_grid_finds_shadow_zones) {
    // Downward-refracting water (summer thermocline) from a shallow source
    // leaves a large unensonified wedge under the source. That wedge is the
    // classic shadow zone, and it must show up as empty cells.
    SoundSpeedProfile<512> p;
    fill_profile(p, 0, 1000, 201, [](Real z) {
        return sound_speed::mackenzie(static_cast<Real>(20) - static_cast<Real>(0.018) * z,
                                      35, z);
    });

    TraceConfig cfg;
    cfg.max_range_m = 12000;
    cfg.bottom_depth_m = 1000;
    cfg.bottom = BoundaryAction::Absorb;
    cfg.surface = BoundaryAction::Absorb;

    CoverageGrid<240, 100> grid(12000, 0, 1000);
    CoverageView view = grid.view();
    coverage_clear(view);
    PT_CHECK(view.valid());

    std::array<RayPoint, 4096> buf{};
    trace_fan(p.view(), 50, deg2rad(-12), deg2rad(12), 241, cfg, buf,
              [&](std::size_t, Real, std::span<const RayPoint> path, const TraceResult&) {
                  coverage_mark(view, path);
              });

    const double shadow = static_cast<double>(coverage_shadow_fraction(view));
    std::printf("       shadow fraction = %.1f%% of the 240x100 grid\n", shadow * 100.0);
    PT_CHECK(shadow > 0.05);   // there IS a shadow zone
    PT_CHECK(shadow < 0.95);   // ... and the fan did ensonify something

    // Right next to the source must be lit; the far bottom corner must not be.
    PT_CHECK(coverage_at(view, 100, 50) > 0);
    // Queries outside the grid are shadow by definition, never a crash.
    PT_CHECK(coverage_at(view, -1, 50) == 0);
    PT_CHECK(coverage_at(view, 1e9, 50) == 0);
    PT_CHECK(coverage_at(view, 100, 1e9) == 0);
}

PT_TEST(coverage_clear_resets_every_cell) {
    CoverageGrid<16, 8> grid(1000, 0, 100);
    CoverageView view = grid.view();
    std::array<RayPoint, 2> path{RayPoint{0, 50, 0, 0, 1500},
                                 RayPoint{1000, 50, 0, static_cast<Real>(0.6), 1500}};
    coverage_mark(view, path);
    PT_CHECK(static_cast<double>(coverage_shadow_fraction(view)) < 1.0);
    coverage_clear(view);
    PT_CHECK_NEAR(coverage_shadow_fraction(view), 1.0, 1e-12);
}
