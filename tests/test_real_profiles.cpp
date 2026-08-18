// SPDX-License-Identifier: Apache-2.0
// Real ocean profiles.
//
// Every profile the library had been tested against until v0.11 was analytic:
// Munk, a linear gradient, an isovelocity slab. Those are smooth, and they are
// smooth in a way real water is not. This file runs the same machinery over six
// World Ocean Atlas 2023 profiles and asserts what the data actually shows --
// including where the library's assumptions are least defensible.
#include "framework.hpp"

#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

#include "data/woa_profiles.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace phantom;

namespace {

using RealProfile = SoundSpeedProfile<128>;

// Builds a sound-speed profile from a site's T/S levels using the library
// default (UNESCO / ITS-90), with the site's own latitude in the depth-to-
// pressure conversion.
bool build(const woa23::Site& site, RealProfile& out) {
    out.clear();
    for (std::size_t i = 0; i < site.level_count; ++i) {
        const woa23::Level& lv = site.levels[i];
        const Real z = static_cast<Real>(lv.depth_m);
        const Real c = sound_speed::unesco(static_cast<Real>(lv.temperature_c),
                                           static_cast<Real>(lv.salinity_psu),
                                           z, static_cast<Real>(site.latitude_deg));
        if (!out.push(z, c)) return false;
    }
    return true;
}

struct Axis { Real depth_m; Real speed_mps; std::size_t index; };

Axis find_axis(const RealProfile& p) {
    const ProfileView v = p.view();
    Axis a{v.depth_m[0], v.speed_mps[0], 0};
    for (std::size_t i = 1; i < v.point_count(); ++i) {
        if (v.speed_mps[i] < a.speed_mps) a = Axis{v.depth_m[i], v.speed_mps[i], i};
    }
    return a;
}

Real speed_i(const RealProfile& p, std::size_t i) { return p.view().speed_mps[i]; }
Real depth_i(const RealProfile& p, std::size_t i) { return p.view().depth_m[i]; }

}  // namespace

PT_TEST(every_real_profile_loads_and_is_physical) {
    PT_CHECK(woa23::kSiteCount == 6);
    for (std::size_t k = 0; k < woa23::kSiteCount; ++k) {
        const woa23::Site& site = woa23::kSites[k];
        RealProfile p;
        PT_CHECK(build(site, p));
        PT_CHECK(p.size() == site.level_count);

        // Depths strictly increasing, speeds inside the range any ocean can
        // produce. 1400-1600 m/s covers the Black Sea's near-fresh surface at
        // one end and the Levantine's hot salty water at the other.
        for (std::size_t i = 0; i < p.size(); ++i) {
            if (i > 0) PT_CHECK(depth_i(p, i) > depth_i(p, i - 1));
            PT_CHECK(speed_i(p, i) > 1400);
            PT_CHECK(speed_i(p, i) < 1600);
        }
    }
}

PT_TEST(real_profiles_have_the_structure_their_descriptions_claim) {
    // The descriptions shipped with the data are checked against the data.
    // An earlier draft of tools/fetch_woa_profiles.py called the Norwegian Sea
    // upward-refracting and put the North Atlantic axis at 1100 m; both were
    // wrong, and nothing would have caught it without a test like this.
    std::printf("       %-11s %7s %9s %9s %9s %9s\n",
                "site", "z_max", "c_surface", "c_axis", "z_axis", "excess");
    for (std::size_t k = 0; k < woa23::kSiteCount; ++k) {
        const woa23::Site& site = woa23::kSites[k];
        RealProfile p;
        PT_CHECK(build(site, p));
        const Axis a = find_axis(p);
        const double excess = static_cast<double>(speed_i(p, 0) - a.speed_mps);
        std::printf("       %-11s %7.0f %9.2f %9.2f %9.0f %8.2f m/s\n",
                    site.name, static_cast<double>(depth_i(p, p.size() - 1)),
                    static_cast<double>(speed_i(p, 0)),
                    static_cast<double>(a.speed_mps),
                    static_cast<double>(a.depth_m), excess);
    }

    auto site_by_name = [](const char* want) -> const woa23::Site& {
        for (std::size_t k = 0; k < woa23::kSiteCount; ++k) {
            const char* n = woa23::kSites[k].name;
            std::size_t i = 0;
            while (n[i] != '\0' && want[i] != '\0' && n[i] == want[i]) ++i;
            if (n[i] == '\0' && want[i] == '\0') return woa23::kSites[k];
        }
        // Never reached unless a site is renamed; the caller's PT_CHECKs on the
        // returned profile will fail loudly rather than silently pass.
        return woa23::kSites[0];
    };

    // North Atlantic: the textbook deep channel, axis near 950 m.
    {
        RealProfile p;
        PT_CHECK(build(site_by_name("n-atlantic"), p));
        const Axis a = find_axis(p);
        PT_CHECK_NEAR(static_cast<double>(a.depth_m), 950.0, 150.0);
        // A real channel: the axis is a genuine interior minimum, well below
        // both the surface and the bottom.
        PT_CHECK(a.speed_mps < speed_i(p, 0) - 10);
        PT_CHECK(a.speed_mps < speed_i(p, p.size() - 1) - 10);
    }

    // Black Sea: axis very shallow, because salinity not temperature sets it.
    {
        RealProfile p;
        PT_CHECK(build(site_by_name("black-sea"), p));
        const Axis a = find_axis(p);
        PT_CHECK(a.depth_m < 150);
    }

    // Norwegian Sea: near-surface water refracts DOWNWARD -- sound speed falls
    // with depth over the top hundred metres, so a surface source has a shadow
    // zone beneath it rather than a duct.
    {
        RealProfile p;
        PT_CHECK(build(site_by_name("norwegian"), p));
        PT_CHECK(speed_i(p, 4) < speed_i(p, 0));
        const Axis a = find_axis(p);
        PT_CHECK(a.depth_m > 500);
    }
}

PT_TEST(the_black_sea_is_where_ignoring_salinity_stops_being_safe) {
    // Most sound-speed work treats salinity as a small correction around 35 PSU,
    // and in most of the ocean that is true. The Black Sea is where it is not:
    // the surface is near-fresh. This measures the error a caller would make by
    // assuming 35 PSU, at every site, so the claim is quantified rather than
    // asserted.
    std::printf("       error from assuming S = 35 PSU instead of the measured value:\n");
    std::printf("       %-11s %12s %12s\n", "site", "surface", "worst");
    double black_sea_worst = 0, others_worst = 0;
    for (std::size_t k = 0; k < woa23::kSiteCount; ++k) {
        const woa23::Site& site = woa23::kSites[k];
        double worst = 0, surface = 0;
        for (std::size_t i = 0; i < site.level_count; ++i) {
            const woa23::Level& lv = site.levels[i];
            const Real z = static_cast<Real>(lv.depth_m);
            const Real lat = static_cast<Real>(site.latitude_deg);
            const double truth = static_cast<double>(sound_speed::unesco(
                static_cast<Real>(lv.temperature_c), static_cast<Real>(lv.salinity_psu), z, lat));
            const double assumed = static_cast<double>(sound_speed::unesco(
                static_cast<Real>(lv.temperature_c), static_cast<Real>(35), z, lat));
            const double d = std::fabs(truth - assumed);
            if (i == 0) surface = d;
            worst = std::max(worst, d);
        }
        std::printf("       %-11s %10.2f m/s %10.2f m/s\n", site.name, surface, worst);
        const char* n = site.name;
        const bool is_black = n[0] == 'b' && n[1] == 'l';
        if (is_black) black_sea_worst = worst;
        else others_worst = std::max(others_worst, worst);
    }
    std::printf("       Black Sea %.1f m/s worst, vs %.1f m/s at every other site.\n",
                black_sea_worst, others_worst);
    std::printf("       That is %.1f%% of the sound speed. For scale, the whole excess\n",
                100.0 * black_sea_worst / 1500.0);
    std::printf("       depth of the Black Sea channel measured above is 24 m/s, so the\n");
    std::printf("       error is the same size as the feature -- assume 35 PSU there and\n");
    std::printf("       the channel you trace is not the channel that exists.\n");
    PT_CHECK(black_sea_worst > 15.0);
    PT_CHECK(others_worst < 8.0);
    PT_CHECK(black_sea_worst > 2.5 * others_worst);
}

PT_TEST(rays_trace_through_real_water) {
    // The tracer was built and verified on Munk, which is analytic and smooth.
    // Real profiles are piecewise linear between standard levels, with a kink at
    // every one of them. This checks the tracer survives that.
    //
    // TWO checks, because the obvious one is weaker than it looks.
    //
    // Snell's invariant xi = cos(theta)/c is the natural thing to measure, and
    // it comes out at EXACTLY zero drift here -- which should be suspicious
    // rather than reassuring. The tracer carries xi and reconstructs theta from
    // it as theta = +/- acos(xi*c), so measuring cos(theta)/c largely measures
    // whether acos and cos round-trip. It does catch a ray that jumps layers
    // wrongly, but it cannot catch a ray that turns in the wrong PLACE.
    //
    // The turning depth can. At a turning point theta = 0, so the local sound
    // speed must equal c_source / cos(theta_launch) -- a statement about where
    // the profile lookup put the ray, independent of how theta is stored.
    for (std::size_t k = 0; k < woa23::kSiteCount; ++k) {
        const woa23::Site& site = woa23::kSites[k];
        RealProfile p;
        PT_CHECK(build(site, p));
        const ProfileView v = p.view();
        const Axis a = find_axis(p);
        const Real c_src = speed_at(v, a.depth_m);

        TraceConfig cfg;
        cfg.max_time_s = 20;
        cfg.max_range_m = 50000;
        cfg.surface = BoundaryAction::Absorb;
        cfg.bottom  = BoundaryAction::Absorb;

        double worst_drift = 0;
        double worst_turn_err = 0;
        std::size_t traced = 0, turns_checked = 0;
        for (int deg = -10; deg <= 10; deg += 2) {
            if (deg == 0) continue;              // no turning point to check
            static std::array<RayPoint, 8192> pts{};
            const Real launch = deg2rad(static_cast<Real>(deg));
            const TraceResult r = trace_ray(v, a.depth_m, launch, cfg, pts);
            if (r.point_count < 3) continue;
            ++traced;

            const Real xi0 = std::cos(pts[0].angle_rad) / pts[0].speed_mps;
            for (std::size_t i = 1; i < r.point_count; ++i) {
                const Real xi = std::cos(pts[i].angle_rad) / pts[i].speed_mps;
                worst_drift = std::max(worst_drift,
                    std::fabs(static_cast<double>((xi - xi0) / xi0)));
            }

            // A turning point is a point the tracer emits with angle exactly
            // zero -- it cuts the arc there rather than stepping past it.
            //
            // Detecting it as "the angle changed sign across i" does NOT work,
            // and getting that wrong is what the first version of this test did:
            // because the turn itself sits at zero, and zero is not > 0, the
            // point BEFORE the turn also looks like a sign change. It then
            // measured the sound speed one layer short of the turn and reported
            // errors up to 4.3 m/s that were entirely its own.
            const Real predicted = c_src / std::cos(launch);
            for (std::size_t i = 1; i + 1 < r.point_count; ++i) {
                if (std::fabs(static_cast<double>(pts[i].angle_rad)) > 1e-7) continue;
                // Skip boundary reflections: those reverse without turning.
                if (pts[i].depth_m <= v.min_depth() + 1) continue;
                if (pts[i].depth_m >= v.max_depth() - 1) continue;
                ++turns_checked;
                worst_turn_err = std::max(worst_turn_err,
                    std::fabs(static_cast<double>(speed_at(v, pts[i].depth_m) - predicted)));
            }
        }
        std::printf("       %-11s %2zu rays, Snell drift %.1e (near-tautological),"
                    " %2zu turns, worst c_turn error %.3f m/s\n",
                    site.name, traced, worst_drift, turns_checked, worst_turn_err);
        PT_CHECK(traced >= 8);
        PT_CHECK(worst_drift < pt::tol(1e-12, 1e-5));
        PT_CHECK(turns_checked >= 1);
        // Tight, because the tracer places the turn EXACTLY: it cuts the arc
        // at theta = 0 rather than stepping past and interpolating. There is no
        // profile-resolution term to allow for.
        PT_CHECK(worst_turn_err < pt::tol(1e-9, 0.05));
    }
}
