// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — ensonification grid and shadow zone extraction.
//
// A shadow zone is not a property of a single ray; it is the absence of every
// ray. So it is computed the honest way: fan out a dense set of rays, stamp the
// range-depth cells each one crosses, then report the cells nothing reached.
// The grid is caller-owned and fixed size -- no allocation in this library.
#ifndef PHANTOM_COVERAGE_HPP
#define PHANTOM_COVERAGE_HPP

#include "phantom/types.hpp"

#include <array>
#include <span>

namespace phantom {

// Non-owning view of a range-depth hit-count grid, row-major: index = iz*nr+ir.
struct CoverageView {
    std::span<std::uint16_t> cells;
    std::size_t nr = 0;
    std::size_t nz = 0;
    Real range_max_m = 0;
    Real depth_min_m = 0;
    Real depth_max_m = 0;

    [[nodiscard]] constexpr bool valid() const noexcept {
        return nr > 0 && nz > 0 && cells.size() == nr * nz
            && range_max_m > 0 && depth_max_m > depth_min_m;
    }
};

void coverage_clear(const CoverageView& grid) noexcept;

// Stamps every cell the ray polyline passes through. Segments are subdivided at
// sub-cell resolution, so an arc that crosses a cell corner is not missed.
void coverage_mark(const CoverageView& grid, std::span<const RayPoint> path) noexcept;

// Hit count at a physical location; 0 means shadow. Out-of-grid returns 0.
[[nodiscard]] std::uint16_t coverage_at(const CoverageView& grid, Real range_m, Real depth_m) noexcept;

// Fraction of cells with zero hits, in [0, 1].
[[nodiscard]] Real coverage_shadow_fraction(const CoverageView& grid) noexcept;

// Owning fixed-size grid. NR * NZ * 2 bytes of storage, statically allocated.
template <std::size_t NR, std::size_t NZ>
class CoverageGrid {
    static_assert(NR > 0 && NZ > 0, "grid dimensions must be non-zero");

  public:
    CoverageGrid(Real range_max_m, Real depth_min_m, Real depth_max_m) noexcept
        : range_max_(range_max_m), depth_min_(depth_min_m), depth_max_(depth_max_m) {}

    [[nodiscard]] CoverageView view() noexcept {
        return CoverageView{std::span<std::uint16_t>(cells_.data(), cells_.size()),
                            NR, NZ, range_max_, depth_min_, depth_max_};
    }

    [[nodiscard]] static constexpr std::size_t cell_count() noexcept { return NR * NZ; }

  private:
    std::array<std::uint16_t, NR * NZ> cells_{};
    Real range_max_;
    Real depth_min_;
    Real depth_max_;
};

}  // namespace phantom

#endif  // PHANTOM_COVERAGE_HPP
