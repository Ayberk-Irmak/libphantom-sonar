// SPDX-License-Identifier: Apache-2.0
//
// The canonical benchmark: a fan of rays through the Munk deep sound channel.
// Writes CSV that tools/plot_rays.py turns into the figure in the README.
//
//   ./munk_simulation [output_dir]
#include "phantom/channel.hpp"
#include "phantom/coverage.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string>

using namespace phantom;

namespace {

constexpr std::size_t kProfilePoints = 501;
constexpr Real        kBottomDepth   = 5000;
constexpr Real        kMaxRange      = 100000;   // 100 km
// Two fans, on purpose. A sparse one is what you can actually read in a figure;
// a dense one is what you need before calling any cell a "shadow zone". Quoting
// a shadow fraction computed from a sparse fan is measuring your ray spacing,
// not the ocean.
constexpr std::size_t kPlotRays      = 61;
constexpr std::size_t kCoverageRays  = 721;
constexpr std::size_t kMaxPoints     = 8192;

constexpr std::size_t kGridR = 500;
constexpr std::size_t kGridZ = 250;

// Statically allocated. The library never allocates; neither does this example.
SoundSpeedProfile<kProfilePoints + 1> g_profile;
std::array<RayPoint, kMaxPoints>      g_scratch;
CoverageGrid<kGridR, kGridZ>          g_grid(kMaxRange, 0, kBottomDepth);

bool write_profile_csv(const std::string& path, const ProfileView& svp) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr) return false;
    std::fprintf(f, "depth_m,munk_mps,mackenzie_mps,medwin_mps\n");
    for (std::size_t i = 0; i < svp.point_count(); ++i) {
        const Real z = svp.depth_m[i];
        // A representative deep-ocean T/S column, for equation comparison only.
        const Real t = static_cast<Real>(4) + static_cast<Real>(16) * std::exp(-z / static_cast<Real>(300));
        std::fprintf(f, "%.3f,%.6f,%.6f,%.6f\n",
                     static_cast<double>(z),
                     static_cast<double>(svp.speed_mps[i]),
                     static_cast<double>(sound_speed::mackenzie(t, 35, z)),
                     static_cast<double>(sound_speed::medwin(t, 35, z)));
    }
    std::fclose(f);
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string out_dir = (argc > 1) ? std::string(argv[1]) : std::string(".");
    // Optional second argument overrides the dense fan size, so the shadow
    // fraction can be checked for convergence (see docs/validation.md).
    std::size_t coverage_rays = kCoverageRays;
    if (argc > 2) {
        const long n = std::strtol(argv[2], nullptr, 10);
        if (n >= 2 && n <= 100000) coverage_rays = static_cast<std::size_t>(n);
    }

    // ---- 1. Ocean model -----------------------------------------------------
    if (!fill_profile(g_profile, 0, kBottomDepth, kProfilePoints,
                      [](Real z) { return sound_speed::munk(z); })) {
        std::fprintf(stderr, "failed to build profile\n");
        return 1;
    }
    const ProfileView svp = g_profile.view();

    std::printf("libphantom-sonar %s  --  Munk deep sound channel simulation\n\n",
                PHANTOM_VERSION_STRING);
    std::printf("Ocean model\n");
    std::printf("  profile points     : %zu (0 - %.0f m)\n",
                svp.point_count(), static_cast<double>(kBottomDepth));
    std::printf("  surface speed      : %.3f m/s\n", static_cast<double>(svp.speed_mps.front()));
    std::printf("  bottom speed       : %.3f m/s\n", static_cast<double>(svp.speed_mps.back()));

    // ---- 2. Channel analysis ------------------------------------------------
    const ChannelInfo ch = analyze_sofar(svp);
    std::printf("\nDeep sound channel\n");
    if (ch.found) {
        std::printf("  axis depth         : %.1f m\n", static_cast<double>(ch.axis_depth_m));
        std::printf("  axis speed         : %.3f m/s\n", static_cast<double>(ch.axis_speed_mps));
        std::printf("  limited by         : %s (%.3f m/s)\n",
                    ch.limited_by_surface ? "surface" : "sea floor",
                    static_cast<double>(ch.limiting_speed_mps));
        std::printf("  trapping cone      : +/- %.3f deg\n",
                    static_cast<double>(rad2deg(ch.max_trapped_angle_rad)));
        std::printf("  conjugate depths   : %.0f m .. %.0f m\n",
                    static_cast<double>(ch.upper_conjugate_m),
                    static_cast<double>(ch.lower_conjugate_m));
    } else {
        std::printf("  none detected\n");
    }

    // ---- 3. Ray fan ---------------------------------------------------------
    TraceConfig cfg;
    cfg.max_range_m = kMaxRange;
    cfg.max_time_s = 200;
    cfg.bottom_depth_m = kBottomDepth;

    CoverageView grid = g_grid.view();
    coverage_clear(grid);

    const std::string ray_path = out_dir + "/rays.csv";
    std::FILE* rf = std::fopen(ray_path.c_str(), "w");
    if (rf == nullptr) {
        std::fprintf(stderr, "cannot write %s\n", ray_path.c_str());
        return 1;
    }
    std::fprintf(rf, "ray_id,launch_deg,range_m,depth_m,angle_deg,time_s,speed_mps\n");

    std::size_t total_points = 0;
    std::size_t total_steps  = 0;

    const Real source_depth = ch.found ? ch.axis_depth_m : static_cast<Real>(1000);

    trace_fan(svp, source_depth, deg2rad(-15), deg2rad(15), kPlotRays, cfg, g_scratch,
              [&](std::size_t id, Real launch, std::span<const RayPoint> path,
                  const TraceResult&) {
                  for (const RayPoint& p : path) {
                      std::fprintf(rf, "%zu,%.4f,%.4f,%.4f,%.6f,%.9f,%.4f\n",
                                   id, static_cast<double>(rad2deg(launch)),
                                   static_cast<double>(p.range_m),
                                   static_cast<double>(p.depth_m),
                                   static_cast<double>(rad2deg(p.angle_rad)),
                                   static_cast<double>(p.time_s),
                                   static_cast<double>(p.speed_mps));
                  }
              });
    std::fclose(rf);

    // Dense fan, timed with no I/O in the loop: this is the real hot path.
    const auto t0 = std::chrono::steady_clock::now();
    trace_fan(svp, source_depth, deg2rad(-15), deg2rad(15), coverage_rays, cfg, g_scratch,
              [&](std::size_t, Real, std::span<const RayPoint> path, const TraceResult& res) {
                  total_points += path.size();
                  total_steps  += res.steps;
                  coverage_mark(grid, path);
              });
    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    std::printf("\nRay fan\n");
    std::printf("  plotted rays       : %zu over +/- 15 deg, to %.0f km\n",
                kPlotRays, static_cast<double>(kMaxRange) / 1000.0);
    std::printf("  coverage rays      : %zu (%.4f deg spacing)\n",
                coverage_rays, 30.0 / static_cast<double>(coverage_rays - 1));
    std::printf("  arc steps          : %zu\n", total_steps);
    std::printf("  polyline points    : %zu\n", total_points);
    std::printf("  trace + grid stamp : %.2f ms total, %.1f ns per arc step\n",
                static_cast<double>(elapsed_ns) / 1e6,
                total_steps ? static_cast<double>(elapsed_ns) / static_cast<double>(total_steps) : 0.0);

    // ---- 4. Coverage / shadow zones ----------------------------------------
    const Real shadow = coverage_shadow_fraction(grid);
    std::printf("\nEnsonification\n");
    std::printf("  grid               : %zu x %zu cells (%.0f m x %.0f m each)\n",
                kGridR, kGridZ,
                static_cast<double>(kMaxRange) / static_cast<double>(kGridR),
                static_cast<double>(kBottomDepth) / static_cast<double>(kGridZ));
    std::printf("  shadow fraction    : %.1f %%  (from the dense fan; this number\n"
                "                       is only meaningful once it stops moving as\n"
                "                       you add rays -- see docs/validation.md)\n",
                static_cast<double>(shadow) * 100.0);

    const std::string cov_path = out_dir + "/coverage.csv";
    if (std::FILE* cf = std::fopen(cov_path.c_str(), "w")) {
        std::fprintf(cf, "range_m,depth_m,hits\n");
        for (std::size_t iz = 0; iz < kGridZ; ++iz) {
            for (std::size_t ir = 0; ir < kGridR; ++ir) {
                const double r = (static_cast<double>(ir) + 0.5)
                               * static_cast<double>(kMaxRange) / static_cast<double>(kGridR);
                const double z = (static_cast<double>(iz) + 0.5)
                               * static_cast<double>(kBottomDepth) / static_cast<double>(kGridZ);
                std::fprintf(cf, "%.1f,%.1f,%u\n", r, z,
                             static_cast<unsigned>(grid.cells[iz * kGridR + ir]));
            }
        }
        std::fclose(cf);
    }

    write_profile_csv(out_dir + "/profile.csv", svp);

    std::printf("\nWrote %s/rays.csv, %s/coverage.csv, %s/profile.csv\n",
                out_dir.c_str(), out_dir.c_str(), out_dir.c_str());
    std::printf("Plot with:  python3 tools/plot_rays.py %s\n", out_dir.c_str());
    return 0;
}
