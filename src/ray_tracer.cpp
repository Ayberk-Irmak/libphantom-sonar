// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — analytic constant-gradient ray tracing kernel.
//
// Non-templated on purpose: this translation unit is compiled once and every
// caller shares it, which keeps flash usage flat on MCU targets and keeps the
// hot loop in one place for profiling.
#include "phantom/ray_tracer.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwoLocal = static_cast<Real>(2);

// |dc/dz| below this is treated as an isovelocity layer. At 1e-12 s^-1 a ray
// would need >1e9 m of range to deviate one metre, so the branch is safe.
constexpr Real kGradEps  = static_cast<Real>(1e-12);
constexpr Real kXiEps    = static_cast<Real>(1e-15);  // xi == 0 <=> vertical ray
constexpr Real kSinEps   = static_cast<Real>(1e-12);  // |sin(theta)| == 0 <=> at a vertex
constexpr Real kDepthEps = static_cast<Real>(1e-9);   // boundary snapping, metres

constexpr std::uint32_t kMaxStalls = 4;  // consecutive zero-length steps tolerated

inline Real safe_sqrt(Real v) noexcept { return v > kZero ? std::sqrt(v) : kZero; }

// Travel time across one constant-gradient arc.
//
//   dt = (1/g) * ln[ c2 (1 + sin1) / (c1 (1 + sin2)) ]        (form A)
//      = (1/g) * ln[ c1 (1 - sin2) / (c2 (1 - sin1)) ]        (form B)
//
// The two are algebraically identical via (1+sin)/cos == cos/(1-sin). Form A
// loses precision as sin -> -1 (steep upgoing), form B as sin -> +1. Within a
// single layer both sines share a sign, so their sum picks the stable branch.
inline Real arc_time(Real c1, Real sin1, Real c2, Real sin2, Real g) noexcept {
    if (sin1 + sin2 >= kZero) {
        return std::log((c2 * (kOne + sin1)) / (c1 * (kOne + sin2))) / g;
    }
    return std::log((c1 * (kOne - sin2)) / (c2 * (kOne - sin1))) / g;
}

// Inverse of arc_time: the grazing angle reached after travelling for `dt`.
//
// With X(theta) = (1 + sin theta) / cos theta = tan(pi/4 + theta/2), the time
// relation above is simply X2 = X1 * exp(-g dt), so the angle inverts in closed
// form. This is what lets a travel-time budget cut an arc exactly, instead of
// stopping short at the previous layer boundary.
inline Real angle_after_time(Real sin1, Real g, Real dt) noexcept {
    const Real theta1 = std::asin(clamp(sin1, -kOne, kOne));
    const Real x1     = std::tan(kPi / static_cast<Real>(4) + theta1 / static_cast<Real>(2));
    const Real x2     = x1 * std::exp(-g * dt);
    return static_cast<Real>(2) * (std::atan(x2) - kPi / static_cast<Real>(4));
}

}  // namespace

TraceResult trace_ray(const ProfileView& svp,
                      Real source_depth_m,
                      Real launch_angle_rad,
                      const TraceConfig& cfg,
                      std::span<RayPoint> out) noexcept {
    TraceResult res{};
    res.status = TraceStatus::Degenerate;

    if (!svp.valid() || out.empty()) return res;

    // Clamp the working water column to the intersection of the configured
    // boundaries and the profile support, so speed_at() never extrapolates.
    Real z_top = cfg.surface_depth_m;
    Real z_bot = (cfg.bottom_depth_m > kZero) ? cfg.bottom_depth_m : svp.max_depth();
    if (z_top < svp.min_depth()) z_top = svp.min_depth();
    if (z_bot > svp.max_depth()) z_bot = svp.max_depth();
    if (!(z_bot > z_top)) return res;
    if (!(cfg.max_range_m > kZero) || !(cfg.max_time_s > kZero)) return res;

    Real z = clamp(source_depth_m, z_top, z_bot);
    Real c = speed_at(svp, z);
    if (!(c > kZero)) return res;

    const Real theta = clamp(launch_angle_rad, -kHalfPi, kHalfPi);
    Real sin_th = std::sin(theta);
    const Real xi = std::cos(theta) / c;  // Snell invariant, constant hereafter
    res.snell_invariant = xi;

    // A ray launched into a boundary reflects before it takes its first step.
    if (z <= z_top + kDepthEps && sin_th < kZero) {
        if (cfg.surface == BoundaryAction::Absorb) { res.status = TraceStatus::Absorbed; return res; }
        sin_th = -sin_th;
        ++res.surface_bounces;
    } else if (z >= z_bot - kDepthEps && sin_th > kZero) {
        if (cfg.bottom == BoundaryAction::Absorb) { res.status = TraceStatus::Absorbed; return res; }
        sin_th = -sin_th;
        ++res.bottom_bounces;
    }

    Real r = kZero;
    Real t = kZero;
    std::uint32_t bounces = 0;
    std::uint32_t stalls  = 0;

    auto emit = [&](Real rr, Real zz, Real ss, Real tt, Real cc) noexcept -> bool {
        if (res.point_count >= out.size()) return false;
        RayPoint& p = out[res.point_count++];
        p.range_m   = rr;
        p.depth_m   = zz;
        p.angle_rad = std::atan2(ss, xi * cc);
        p.time_s    = tt;
        p.speed_mps = cc;
        return true;
    };

    // A ray pinned at a sound speed minimum, or horizontal inside an
    // isovelocity layer, propagates along the horizontal. Both are physical;
    // both would otherwise produce zero-length steps forever.
    auto finish_horizontal = [&]() noexcept {
        const Real r_left = cfg.max_range_m - r;
        const Real dt     = r_left / c;
        if (t + dt > cfg.max_time_s) {
            const Real dt_cap = cfg.max_time_s - t;
            res.path_length_m += dt_cap * c;
            r += dt_cap * c;
            t  = cfg.max_time_s;
            res.status = TraceStatus::MaxTime;
        } else {
            res.path_length_m += r_left;
            r = cfg.max_range_m;
            t += dt;
            res.status = TraceStatus::MaxRange;
        }
        if (!emit(r, z, kZero, t, c)) res.status = TraceStatus::BufferFull;
    };

    // Cuts the current arc exactly on the travel-time budget. Mirrors the range
    // clip below so that both budgets are honoured to machine precision rather
    // than rounded back to the previous layer boundary.
    auto clip_to_time = [&](Real g_layer, bool iso_layer) noexcept {
        const Real dt_left = cfg.max_time_s - t;
        const Real cos_th  = xi * c;
        Real z_end;
        Real c_end;
        Real sin_end;
        Real dr_end;
        if (iso_layer) {
            const Real path = c * dt_left;
            dr_end  = path * cos_th;
            z_end   = z + path * sin_th;
            c_end   = c;
            sin_end = sin_th;
        } else if (cos_th <= static_cast<Real>(1e-9)) {
            // Vertical: the angle never changes, only the depth.
            const Real sgn = (sin_th >= kZero) ? kOne : -kOne;
            c_end   = c * std::exp(sgn * g_layer * dt_left);
            z_end   = z + (c_end - c) / g_layer;
            sin_end = sin_th;
            dr_end  = kZero;
        } else {
            const Real theta2 = angle_after_time(sin_th, g_layer, dt_left);
            sin_end = std::sin(theta2);
            c_end   = std::cos(theta2) / xi;
            z_end   = z + (c_end - c) / g_layer;
            dr_end  = (sin_th - sin_end) / (xi * g_layer);
        }
        r += (dr_end > kZero) ? dr_end : kZero;
        res.path_length_m += dt_left * ((c + c_end) / kTwoLocal);
        t  = cfg.max_time_s;
        z  = clamp(z_end, z_top, z_bot);
        c  = c_end;
        sin_th = sin_end;
        ++res.steps;
        if (!emit(r, z, sin_th, t, c)) res.status = TraceStatus::BufferFull;
        else                           res.status = TraceStatus::MaxTime;
    };

    if (!emit(r, z, sin_th, t, c)) {
        res.status = TraceStatus::BufferFull;
        return res;
    }
    res.status = TraceStatus::MaxRange;

    for (;;) {
        if (r >= cfg.max_range_m) { res.status = TraceStatus::MaxRange; break; }

        // --- 1. Travel direction and the layer we are about to cross --------
        bool down;
        if (sin_th > kSinEps) {
            down = true;
        } else if (sin_th < -kSinEps) {
            down = false;
        } else {
            // At a vertex the ray curves toward the lower sound speed.
            const std::size_t below = find_layer(svp, z, true);
            down = (svp.gradient[below] < kZero);
        }

        const std::size_t seg = find_layer(svp, z, down);
        const Real g   = svp.gradient[seg];
        const bool iso = std::fabs(g) <= kGradEps;

        Real z_target = down ? svp.depth_m[seg + 1] : svp.depth_m[seg];
        if (down && z_target > z_bot)  z_target = z_bot;
        if (!down && z_target < z_top) z_target = z_top;

        // --- 2. Endpoint of the arc, or the turning point if it comes first --
        Real c_target   = c + g * (z_target - z);
        Real cos_target = xi * c_target;
        bool turning    = false;

        if (cos_target >= kOne) {
            if (iso) { finish_horizontal(); break; }  // exactly horizontal, isovelocity
            const Real c_turn = kOne / xi;            // xi > 0 whenever cos_target >= 1
            z_target   = z + (c_turn - c) / g;
            c_target   = c_turn;
            cos_target = kOne;
            turning    = true;
            if (std::fabs(z_target - z) < kDepthEps) {
                // Vertex sitting on a sound speed minimum: an axial ray.
                finish_horizontal();
                break;
            }
        }

        Real sin_target;
        if (turning) {
            sin_target = kZero;
        } else {
            sin_target = safe_sqrt(kOne - cos_target * cos_target);
            if (!down) sin_target = -sin_target;
        }

        // --- 3. Closed-form range and time increments -----------------------
        Real dr;
        Real dt;
        Real ds;   // arc length of this segment
        if (!iso) {
            dr = (xi > kXiEps) ? (sin_th - sin_target) / (xi * g) : kZero;
            dt = arc_time(c, sin_th, c_target, sin_target, g);
            // On a circular arc of radius R = 1/(xi|g|) the length is R|dtheta|.
            // For a vertical ray xi -> 0 and R diverges, but the path is then
            // simply the depth change.
            if (xi > kXiEps) {
                const Real th1 = std::asin(clamp(sin_th, -kOne, kOne));
                const Real th2 = std::asin(clamp(sin_target, -kOne, kOne));
                ds = std::fabs(th1 - th2) / (xi * std::fabs(g));
            } else {
                ds = std::fabs(z_target - z);
            }
        } else {
            const Real dz     = z_target - z;
            const Real abs_sn = std::fabs(sin_th);
            dr = (abs_sn > kSinEps) ? dz * (xi * c) / sin_th : kZero;
            dt = (abs_sn > kSinEps) ? std::fabs(dz) / (c * abs_sn) : kZero;
            ds = (abs_sn > kSinEps) ? std::fabs(dz) / abs_sn : kZero;
        }
        if (dr < kZero) dr = kZero;  // guards against catastrophic cancellation

        // --- 4. Clip against the range budget, exactly ----------------------
        const Real r_left = cfg.max_range_m - r;
        if (dr > r_left) {
            Real z_end;
            Real c_end;
            Real sin_end;
            if (!iso && xi > kXiEps) {
                sin_end = clamp(sin_th - xi * g * r_left, -kOne, kOne);
                c_end   = safe_sqrt(kOne - sin_end * sin_end) / xi;
                z_end   = z + (c_end - c) / g;
            } else {
                const Real cos_th = xi * c;
                z_end   = z + r_left * (sin_th / cos_th);
                c_end   = c;
                sin_end = sin_th;
            }
            const Real dt_end = iso ? (r_left / (c * xi * c))
                                    : arc_time(c, sin_th, c_end, sin_end, g);
            // Both budgets can bind within the same arc; the earlier one wins.
            if (t + dt_end > cfg.max_time_s) { clip_to_time(g, iso); break; }
            r = cfg.max_range_m;
            t += dt_end;
            if (!iso && xi > kXiEps) {
                const Real th1 = std::asin(clamp(sin_th, -kOne, kOne));
                const Real th2 = std::asin(clamp(sin_end, -kOne, kOne));
                res.path_length_m += std::fabs(th1 - th2) / (xi * std::fabs(g));
            } else {
                const Real cos_th = xi * c;
                res.path_length_m += (cos_th > kZero) ? r_left / cos_th
                                                      : std::fabs(z_end - z);
            }
            z = clamp(z_end, z_top, z_bot);
            c = c_end;
            sin_th = sin_end;
            ++res.steps;
            if (!emit(r, z, sin_th, t, c)) res.status = TraceStatus::BufferFull;
            else                           res.status = TraceStatus::MaxRange;
            break;
        }

        // --- 5. Clip against the time budget, exactly -----------------------
        if (t + dt > cfg.max_time_s) { clip_to_time(g, iso); break; }

        // --- 6. Commit ------------------------------------------------------
        const bool stalled = (dr < kDepthEps) && (std::fabs(z_target - z) < kDepthEps);
        r += dr;
        t += dt;
        res.path_length_m += (ds > kZero) ? ds : kZero;
        z  = z_target;
        c  = c_target;
        sin_th = sin_target;
        ++res.steps;
        if (turning) ++res.turning_points;

        stalls = stalled ? (stalls + 1) : 0;
        if (stalls >= kMaxStalls) { res.status = TraceStatus::Degenerate; break; }

        // --- 7. Boundary interaction ----------------------------------------
        if (!turning) {
            if (z <= z_top + kDepthEps) {
                z = z_top;
                if (cfg.surface == BoundaryAction::Absorb) {
                    res.status = TraceStatus::Absorbed;
                    emit(r, z, sin_th, t, c);
                    break;
                }
                sin_th = -sin_th;
                ++res.surface_bounces;
                ++bounces;
            } else if (z >= z_bot - kDepthEps) {
                z = z_bot;
                if (cfg.bottom == BoundaryAction::Absorb) {
                    res.status = TraceStatus::Absorbed;
                    emit(r, z, sin_th, t, c);
                    break;
                }
                sin_th = -sin_th;
                ++res.bottom_bounces;
                ++bounces;
            }
        }

        if (!emit(r, z, sin_th, t, c)) { res.status = TraceStatus::BufferFull; break; }
        if (bounces > cfg.max_bounces) { res.status = TraceStatus::MaxBounces; break; }
    }

    res.final_range_m = r;
    res.final_time_s  = t;
    return res;
}

}  // namespace phantom
