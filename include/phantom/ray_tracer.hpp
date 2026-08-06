// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — analytic constant-gradient ray tracer.
//
// Within a layer of constant sound speed gradient g = dc/dz, an acoustic ray is
// EXACTLY a circular arc. There is no integration error to trade against step
// size: one arc per layer crossing, closed form, and the result is correct to
// machine precision. See docs/math_spec.md for the derivation.
//
//   Snell invariant    xi = cos(theta) / c(z)          (constant along a ray)
//   arc radius         R  = 1 / (xi * |g|)
//   range increment    dr = (sin(theta1) - sin(theta2)) / (xi * g)
//   time increment     dt = ln[ c2 (1 + sin theta1) / (c1 (1 + sin theta2)) ] / g
//   turning depth      where c(z) = 1 / xi
#ifndef PHANTOM_RAY_TRACER_HPP
#define PHANTOM_RAY_TRACER_HPP

#include "phantom/profile.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom {

// Traces a single ray from (range 0, `source_depth_m`) launched at
// `launch_angle_rad` from the horizontal (positive = downgoing).
//
// Writes the polyline into `out` and returns how far it got. Never allocates,
// never throws, no virtual dispatch. The first written point is the source.
TraceResult trace_ray(const ProfileView& svp,
                      Real source_depth_m,
                      Real launch_angle_rad,
                      const TraceConfig& cfg,
                      std::span<RayPoint> out) noexcept;

// Traces `count` rays fanned evenly over [angle_begin_rad, angle_end_rad],
// reusing one scratch buffer. `sink(index, angle_rad, path, result)` is invoked
// per ray with a span valid only for the duration of the call.
template <typename Sink>
void trace_fan(const ProfileView& svp,
               Real source_depth_m,
               Real angle_begin_rad,
               Real angle_end_rad,
               std::size_t count,
               const TraceConfig& cfg,
               std::span<RayPoint> scratch,
               Sink&& sink) {
    if (count == 0) return;
    const Real step = (count == 1)
                        ? static_cast<Real>(0)
                        : (angle_end_rad - angle_begin_rad) / static_cast<Real>(count - 1);
    for (std::size_t i = 0; i < count; ++i) {
        const Real angle = angle_begin_rad + step * static_cast<Real>(i);
        const TraceResult r = trace_ray(svp, source_depth_m, angle, cfg, scratch);
        sink(i, angle, scratch.subspan(0, r.point_count), r);
    }
}

}  // namespace phantom

#endif  // PHANTOM_RAY_TRACER_HPP
