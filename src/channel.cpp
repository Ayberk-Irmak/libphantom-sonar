// SPDX-License-Identifier: Apache-2.0
#include "phantom/channel.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// Depth at which the profile first reaches `target` speed, scanning from
// `from_index` toward `to_index` (either direction). Returns `fallback` if the
// target is never crossed.
Real crossing_depth(const ProfileView& svp, std::size_t from_index, std::size_t to_index,
                    Real target, Real fallback) noexcept {
    const bool upward = to_index < from_index;
    std::size_t i = from_index;
    while (i != to_index) {
        const std::size_t j = upward ? i - 1 : i + 1;
        const Real c0 = svp.speed_mps[i];
        const Real c1 = svp.speed_mps[j];
        if ((c0 - target) * (c1 - target) <= kZero && c1 != c0) {
            const Real f = (target - c0) / (c1 - c0);
            return svp.depth_m[i] + f * (svp.depth_m[j] - svp.depth_m[i]);
        }
        i = j;
    }
    return fallback;
}

}  // namespace

ChannelInfo analyze_sofar(const ProfileView& svp) noexcept {
    ChannelInfo info{};
    if (!svp.valid()) return info;

    const std::size_t n = svp.point_count();

    std::size_t axis = 0;
    for (std::size_t i = 1; i < n; ++i) {
        if (svp.speed_mps[i] < svp.speed_mps[axis]) axis = i;
    }

    info.axis_depth_m      = svp.depth_m[axis];
    info.axis_speed_mps    = svp.speed_mps[axis];
    info.surface_speed_mps = svp.speed_mps[0];
    info.bottom_speed_mps  = svp.speed_mps[n - 1];

    // A minimum at either endpoint is not a duct: there is no upper (or lower)
    // wall to refract energy back, so nothing is trapped.
    info.found = (axis > 0) && (axis + 1 < n);
    if (!info.found) return info;

    info.limited_by_surface = info.surface_speed_mps <= info.bottom_speed_mps;
    info.limiting_speed_mps = info.limited_by_surface ? info.surface_speed_mps
                                                      : info.bottom_speed_mps;

    const Real ratio = info.axis_speed_mps / info.limiting_speed_mps;
    info.max_trapped_angle_rad = (ratio < kOne) ? std::acos(ratio) : kZero;

    info.upper_conjugate_m = crossing_depth(svp, axis, 0, info.limiting_speed_mps, svp.min_depth());
    info.lower_conjugate_m = crossing_depth(svp, axis, n - 1, info.limiting_speed_mps, svp.max_depth());
    return info;
}

Real critical_angle(const ProfileView& svp, Real depth_m, Real turn_depth_m) noexcept {
    if (!svp.valid()) return kZero;
    const Real c0 = speed_at(svp, depth_m);
    const Real c1 = speed_at(svp, turn_depth_m);
    if (!(c0 > kZero) || !(c1 > c0)) return kZero;
    return std::acos(c0 / c1);
}

Real surface_duct_depth(const ProfileView& svp) noexcept {
    if (!svp.valid()) return static_cast<Real>(-1);
    const std::size_t n = svp.point_count();
    // Walk down while sound speed rises; the first turn-over is the duct base.
    std::size_t i = 1;
    while (i < n && svp.speed_mps[i] > svp.speed_mps[i - 1]) ++i;
    if (i == 1 || i >= n) return static_cast<Real>(-1);  // no duct, or isothermal to the bottom
    return svp.depth_m[i - 1];
}

}  // namespace phantom
