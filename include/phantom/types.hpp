// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — core value types. POD only, no allocation, no virtuals.
#ifndef PHANTOM_TYPES_HPP
#define PHANTOM_TYPES_HPP

#include "phantom/config.hpp"

namespace phantom {

inline constexpr Real kPi      = static_cast<Real>(3.14159265358979323846);
inline constexpr Real kHalfPi  = kPi / static_cast<Real>(2);
inline constexpr Real kDeg2Rad = kPi / static_cast<Real>(180);
inline constexpr Real kRad2Deg = static_cast<Real>(180) / kPi;

constexpr Real deg2rad(Real degrees) noexcept { return degrees * kDeg2Rad; }

// Wraps an angle into [-pi, pi]. Needed wherever a bearing difference is taken:
// without it a track sitting near +/-180 degrees produces a 2*pi innovation and
// the filter diverges on the first update.
constexpr Real wrap_pi(Real a) noexcept {
    while (a > kPi) a -= static_cast<Real>(2) * kPi;
    while (a < -kPi) a += static_cast<Real>(2) * kPi;
    return a;
}
constexpr Real rad2deg(Real radians) noexcept { return radians * kRad2Deg; }

constexpr Real clamp(Real v, Real lo, Real hi) noexcept {
    return v < lo ? lo : (v > hi ? hi : v);
}

// One sample along a traced ray. Depth is positive downward from the surface;
// angle is measured from the horizontal, positive downgoing.
struct RayPoint {
    Real range_m   = 0;  // horizontal range from the source
    Real depth_m   = 0;  // positive down
    Real angle_rad = 0;  // grazing angle, positive downgoing
    Real time_s    = 0;  // one-way travel time from the source
    Real speed_mps = 0;  // local sound speed at this point
};

enum class TraceStatus : std::uint8_t {
    MaxRange = 0,   // stopped on the horizontal range limit (normal exit)
    MaxTime,        // stopped on the travel-time limit
    BufferFull,     // caller's RayPoint span ran out
    Absorbed,       // lost at a boundary with BoundaryAction::Absorb
    MaxBounces,     // exceeded the configured boundary-interaction budget
    Degenerate,     // ill-posed input (bad profile, zero sound speed, ...)
};

enum class BoundaryAction : std::uint8_t { Reflect = 0, Absorb = 1 };

struct TraceConfig {
    Real max_range_m     = static_cast<Real>(50000);
    Real max_time_s      = static_cast<Real>(120);
    Real surface_depth_m = static_cast<Real>(0);
    // 0 means "use the deepest point of the sound speed profile".
    Real bottom_depth_m  = static_cast<Real>(0);
    BoundaryAction surface = BoundaryAction::Reflect;
    BoundaryAction bottom  = BoundaryAction::Reflect;
    std::uint32_t max_bounces = 200;
};

struct TraceResult {
    std::size_t   point_count     = 0;
    std::uint32_t surface_bounces = 0;
    std::uint32_t bottom_bounces  = 0;
    std::uint32_t turning_points  = 0;
    std::uint32_t steps           = 0;  // arc segments integrated
    Real          final_range_m   = 0;
    Real          final_time_s    = 0;
    // Arc length travelled, exact: each constant-gradient segment contributes
    // R*|dtheta| with R = 1/(xi|g|), and each isovelocity segment |dz/sin|.
    // Needed for absorption, which is quoted per unit path length rather than
    // per unit range.
    Real          path_length_m   = 0;
    Real          snell_invariant = 0;  // xi = cos(theta)/c, constant by construction
    TraceStatus   status          = TraceStatus::Degenerate;
};

}  // namespace phantom

#endif  // PHANTOM_TYPES_HPP
