// SPDX-License-Identifier: Apache-2.0
#include "phantom/coverage.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);

struct Cell {
    std::size_t ir;
    std::size_t iz;
    bool inside;
};

Cell locate(const CoverageView& g, Real range_m, Real depth_m) noexcept {
    Cell out{0, 0, false};
    if (range_m < kZero || range_m > g.range_max_m) return out;
    if (depth_m < g.depth_min_m || depth_m > g.depth_max_m) return out;

    const Real fr = range_m / g.range_max_m;
    const Real fz = (depth_m - g.depth_min_m) / (g.depth_max_m - g.depth_min_m);

    auto ir = static_cast<std::size_t>(fr * static_cast<Real>(g.nr));
    auto iz = static_cast<std::size_t>(fz * static_cast<Real>(g.nz));
    if (ir >= g.nr) ir = g.nr - 1;
    if (iz >= g.nz) iz = g.nz - 1;
    return Cell{ir, iz, true};
}

inline void stamp(const CoverageView& g, const Cell& cell) noexcept {
    if (!cell.inside) return;
    std::uint16_t& v = g.cells[cell.iz * g.nr + cell.ir];
    if (v != 0xFFFFu) ++v;  // saturate rather than wrap
}

}  // namespace

void coverage_clear(const CoverageView& grid) noexcept {
    if (!grid.valid()) return;
    for (std::uint16_t& c : grid.cells) c = 0;
}

void coverage_mark(const CoverageView& grid, std::span<const RayPoint> path) noexcept {
    if (!grid.valid() || path.empty()) return;

    const Real dr_cell = grid.range_max_m / static_cast<Real>(grid.nr);
    const Real dz_cell = (grid.depth_max_m - grid.depth_min_m) / static_cast<Real>(grid.nz);

    stamp(grid, locate(grid, path[0].range_m, path[0].depth_m));

    for (std::size_t i = 1; i < path.size(); ++i) {
        const RayPoint& a = path[i - 1];
        const RayPoint& b = path[i];
        const Real dr = b.range_m - a.range_m;
        const Real dz = b.depth_m - a.depth_m;

        // Half-cell sampling: guarantees no cell along the segment is skipped.
        const Real steps_r = std::fabs(dr) / (dr_cell * static_cast<Real>(0.5));
        const Real steps_z = std::fabs(dz) / (dz_cell * static_cast<Real>(0.5));
        const Real steps_f = (steps_r > steps_z) ? steps_r : steps_z;

        auto steps = static_cast<std::size_t>(steps_f) + 1;
        // The polyline comes from the ray tracer, so segments span at most one
        // layer; this cap only guards against a pathological caller-built path.
        constexpr std::size_t kMaxSubsteps = 4096;
        if (steps > kMaxSubsteps) steps = kMaxSubsteps;

        for (std::size_t s = 1; s <= steps; ++s) {
            const Real f = static_cast<Real>(s) / static_cast<Real>(steps);
            stamp(grid, locate(grid, a.range_m + dr * f, a.depth_m + dz * f));
        }
    }
}

std::uint16_t coverage_at(const CoverageView& grid, Real range_m, Real depth_m) noexcept {
    if (!grid.valid()) return 0;
    const Cell cell = locate(grid, range_m, depth_m);
    if (!cell.inside) return 0;
    return grid.cells[cell.iz * grid.nr + cell.ir];
}

Real coverage_shadow_fraction(const CoverageView& grid) noexcept {
    if (!grid.valid()) return kZero;
    std::size_t empty = 0;
    for (const std::uint16_t v : grid.cells) {
        if (v == 0) ++empty;
    }
    return static_cast<Real>(empty) / static_cast<Real>(grid.cells.size());
}

}  // namespace phantom
