// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — deep sound channel (SOFAR) and duct analysis.
#ifndef PHANTOM_CHANNEL_HPP
#define PHANTOM_CHANNEL_HPP

#include "phantom/profile.hpp"
#include "phantom/types.hpp"

namespace phantom {

// Result of a deep sound channel analysis.
//
// A ray launched on the channel axis stays trapped while its turning speed
// c_axis / cos(theta) stays below the lower of the two speeds that cap the
// duct. Hence:
//
//   max_trapped_angle = acos( c_axis / c_limit )
//
// The conjugate depths are where c(z) == c_limit above and below the axis;
// energy launched inside the trapping cone never leaves that depth band, which
// is precisely the band an AUV wants to sit in (or avoid).
struct ChannelInfo {
    bool found = false;                       // a genuine interior minimum exists
    Real axis_depth_m = 0;
    Real axis_speed_mps = 0;
    Real surface_speed_mps = 0;
    Real bottom_speed_mps = 0;
    Real limiting_speed_mps = 0;              // min(surface, bottom) speed
    Real upper_conjugate_m = 0;               // depth above axis where c == limit
    Real lower_conjugate_m = 0;               // depth below axis where c == limit
    Real max_trapped_angle_rad = 0;           // half-aperture of the trapping cone
    bool limited_by_surface = false;          // which boundary caps the duct
};

// Locates the deep sound channel. O(n) over the profile, no allocation.
[[nodiscard]] ChannelInfo analyze_sofar(const ProfileView& svp) noexcept;

// Critical grazing angle at `depth_m` for a ray that must turn before reaching
// `turn_depth_m`. Returns 0 when no such ray exists (sound speed does not
// increase enough to refract it back).
[[nodiscard]] Real critical_angle(const ProfileView& svp, Real depth_m, Real turn_depth_m) noexcept;

// Depth of the shallowest sound speed maximum below the surface, i.e. the base
// of a surface duct. Returns a negative value when there is no surface duct.
[[nodiscard]] Real surface_duct_depth(const ProfileView& svp) noexcept;

}  // namespace phantom

#endif  // PHANTOM_CHANNEL_HPP
