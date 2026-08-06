// SPDX-License-Identifier: Apache-2.0
//
// Hot-path timing. The published budget is stated PER ARC STEP, not per "cycle"
// -- a full 100 km trace crosses thousands of layers and obviously cannot fit
// in a sub-microsecond budget. Per-step is the number that actually constrains
// a real-time control loop.
#include "phantom/channel.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <span>
#include <vector>

using namespace phantom;

namespace {

using Clock = std::chrono::steady_clock;

SoundSpeedProfile<1024>         g_munk;
std::array<RayPoint, 16384>     g_buf;

double median(std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    return v.empty() ? 0.0 : v[v.size() / 2];
}

void report(const char* name, double ns_per_op, const char* unit) {
    std::printf("  %-34s %9.1f ns / %s\n", name, ns_per_op, unit);
}

}  // namespace

int main() {
    std::printf("libphantom-sonar %s benchmark  |  Real = %s\n",
                PHANTOM_VERSION_STRING, sizeof(Real) == 4 ? "float" : "double");
    std::printf("--------------------------------------------------------------\n");

    fill_profile(g_munk, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });
    const ProfileView svp = g_munk.view();

    // ---- sound speed equations ---------------------------------------------
    {
        constexpr int kN = 2000000;
        volatile Real sink = 0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kN; ++i) {
            const auto z = static_cast<Real>(i % 5000);
            sink = sound_speed::mackenzie(10, 35, z);
        }
        const auto t1 = Clock::now();
        (void)sink;
        report("mackenzie()",
               static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / kN,
               "call");
    }
    {
        constexpr int kN = 2000000;
        volatile Real sink = 0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kN; ++i) {
            const auto z = static_cast<Real>(i % 5000);
            sink = sound_speed::chen_millero(10, 35, z / 10);
        }
        const auto t1 = Clock::now();
        (void)sink;
        report("chen_millero()  [UNESCO]",
               static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / kN,
               "call");
    }

    // ---- profile lookup -----------------------------------------------------
    {
        constexpr int kN = 5000000;
        volatile Real sink = 0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kN; ++i) {
            sink = speed_at(svp, static_cast<Real>(i % 5000));
        }
        const auto t1 = Clock::now();
        (void)sink;
        report("speed_at()  [501-layer binary search]",
               static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / kN,
               "call");
    }

    // ---- the ray tracing hot path ------------------------------------------
    {
        TraceConfig cfg;
        cfg.max_range_m = 100000;
        cfg.bottom_depth_m = 5000;

        constexpr int kReps = 200;
        std::vector<double> per_step;
        per_step.reserve(kReps);
        std::size_t last_steps = 0;
        double total_ns = 0;
        std::size_t total_steps = 0;

        for (int rep = 0; rep < kReps; ++rep) {
            const Real angle = deg2rad(static_cast<Real>(-14 + (rep % 29)));
            const auto t0 = Clock::now();
            const TraceResult res = trace_ray(svp, 1300, angle, cfg, g_buf);
            const auto t1 = Clock::now();
            const auto ns = static_cast<double>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
            if (res.steps > 0) {
                per_step.push_back(ns / static_cast<double>(res.steps));
                total_ns += ns;
                total_steps += res.steps;
                last_steps = res.steps;
            }
        }

        std::printf("\n  100 km Munk trace: %zu arc steps per ray\n", last_steps);
        report("trace_ray()  median", median(per_step), "arc step");
        report("trace_ray()  mean",
               total_steps ? total_ns / static_cast<double>(total_steps) : 0.0, "arc step");
        std::printf("  %-34s %9.1f us / ray\n", "trace_ray()  full 100 km ray",
                    total_ns / static_cast<double>(per_step.size()) / 1000.0);
    }

    // ---- channel analysis ---------------------------------------------------
    {
        constexpr int kN = 200000;
        volatile Real sink = 0;
        const auto t0 = Clock::now();
        for (int i = 0; i < kN; ++i) {
            sink = analyze_sofar(svp).max_trapped_angle_rad;
        }
        const auto t1 = Clock::now();
        (void)sink;
        report("analyze_sofar()  [501 points]",
               static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / kN,
               "call");
    }

    std::printf("--------------------------------------------------------------\n");
    std::printf("Reported on this machine only. Re-run on your target before\n"
                "quoting any of these numbers as a real-time guarantee.\n");
    return 0;
}
