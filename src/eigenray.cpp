// SPDX-License-Identifier: Apache-2.0
#include "phantom/eigenray.hpp"

#include <cmath>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// One probe: trace to the receiver range and report where the ray got to.
struct Probe {
    bool valid = false;
    Real depth_m = 0;
    Real angle_rad = 0;
    Real time_s = 0;
    Real path_m = 0;
    Real speed_mps = 0;
    std::uint32_t surface_bounces = 0;
    std::uint32_t bottom_bounces = 0;
};

Probe probe_at(const ProfileView& svp, Real source_depth_m, Real angle_rad,
               const TraceConfig& cfg, std::span<RayPoint> scratch) noexcept {
    Probe p;
    const TraceResult r = trace_ray(svp, source_depth_m, angle_rad, cfg, scratch);
    // Only a ray that actually reached the range budget has a defined depth
    // there. One that ran out of time, buffer or water has not arrived, and
    // treating its last point as an arrival would invent a path.
    if (r.status != TraceStatus::MaxRange || r.point_count == 0) return p;
    const RayPoint& end = scratch[r.point_count - 1];
    p.valid = true;
    p.depth_m = end.depth_m;
    p.angle_rad = end.angle_rad;
    p.time_s = end.time_s;
    p.path_m = r.path_length_m;
    p.speed_mps = end.speed_mps;
    p.surface_bounces = r.surface_bounces;
    p.bottom_bounces = r.bottom_bounces;
    return p;
}

}  // namespace

Real thorp_absorption_db_per_km(Real frequency_hz) noexcept {
    if (!(frequency_hz > kZero)) return kZero;
    const Real f = frequency_hz / static_cast<Real>(1000);   // kHz
    const Real f2 = f * f;
    return static_cast<Real>(0.11) * f2 / (kOne + f2)
         + static_cast<Real>(44) * f2 / (static_cast<Real>(4100) + f2)
         + static_cast<Real>(2.75e-4) * f2
         + static_cast<Real>(0.003);
}

Real spreading_loss_db(Real range_m,
                       Real launch_angle_rad,
                       Real arrival_angle_rad,
                       Real jacobian_m_per_rad,
                       Real source_speed_mps,
                       Real arrival_speed_mps) noexcept {
    if (!(range_m > kZero) || !(source_speed_mps > kZero) || !(arrival_speed_mps > kZero)) {
        return kZero;
    }
    const Real j = std::fabs(jacobian_m_per_rad);
    const Real cos_rcv = std::fabs(std::cos(arrival_angle_rad));
    // An epsilon rather than > 0: cos(pi/2) evaluates to 6e-17 in double and
    // -4e-8 in float, never zero, and dividing by it would return a level of
    // some hundreds of dB instead of reporting that a vertically arriving ray
    // has no horizontal tube. 1e-6 rad of grazing angle is far below anything
    // physically meaningful and clears both precisions.
    constexpr Real kCosEps = static_cast<Real>(1e-6);
    if (!(j > kZero) || !(cos_rcv > kCosEps)) return kZero;   // caustic or vertical

    const Real intensity = (arrival_speed_mps * std::fabs(std::cos(launch_angle_rad)))
                         / (source_speed_mps * range_m * cos_rcv * j);
    if (!(intensity > kZero)) return kZero;
    return -static_cast<Real>(10) * std::log10(intensity);
}

Real transmission_loss_db(const Eigenray& ray, Real range_m,
                          Real frequency_hz, Real source_speed_mps) noexcept {
    const Real spread = spreading_loss_db(range_m, ray.launch_angle_rad,
                                          ray.arrival_angle_rad, ray.jacobian_m_per_rad,
                                          source_speed_mps, ray.arrival_speed_mps);
    const Real absorption = thorp_absorption_db_per_km(frequency_hz)
                          * (ray.path_length_m / static_cast<Real>(1000));
    return spread + absorption;
}

std::size_t find_eigenrays(const ProfileView& svp,
                           Real source_depth_m,
                           Real receiver_depth_m,
                           Real range_m,
                           const TraceConfig& cfg,
                           const EigenraySearch& search,
                           std::span<RayPoint> scratch,
                           std::span<Eigenray> out) noexcept {
    if (!svp.valid() || scratch.empty() || out.empty()) return 0;
    if (!(range_m > kZero) || search.fan_count < 2) return 0;
    if (!(search.angle_max_rad > search.angle_min_rad)) return 0;

    // Trace exactly to the receiver. The tracer honours a range budget to
    // machine precision, so the final point IS the arrival -- no interpolation
    // and none of the chord-versus-arc error that resampling would introduce.
    TraceConfig probe_cfg = cfg;
    probe_cfg.max_range_m = range_m;

    const Real d_angle = (search.angle_max_rad - search.angle_min_rad)
                       / static_cast<Real>(search.fan_count - 1);
    const Real c_src = speed_at(svp, source_depth_m);

    std::size_t written = 0;
    Real prev_angle = search.angle_min_rad;
    Probe prev = probe_at(svp, source_depth_m, prev_angle, probe_cfg, scratch);

    for (std::size_t i = 1; i < search.fan_count && written < out.size(); ++i) {
        const Real angle = search.angle_min_rad + d_angle * static_cast<Real>(i);
        const Probe cur = probe_at(svp, source_depth_m, angle, probe_cfg, scratch);

        if (prev.valid && cur.valid) {
            const Real f_lo = prev.depth_m - receiver_depth_m;
            const Real f_hi = cur.depth_m - receiver_depth_m;

            // A sign change brackets an eigenray. Equal signs with one exactly
            // zero is handled by the <= 0 test; two consecutive exact zeros
            // would be a degenerate fan and are not chased.
            if ((f_lo <= kZero && f_hi >= kZero) || (f_lo >= kZero && f_hi <= kZero)) {
                Real lo = prev_angle;
                Real hi = angle;
                Real f_lo_v = f_lo;
                Probe mid;
                Real mid_angle = (lo + hi) / static_cast<Real>(2);

                for (std::size_t k = 0; k < search.refine_steps; ++k) {
                    mid_angle = (lo + hi) / static_cast<Real>(2);
                    mid = probe_at(svp, source_depth_m, mid_angle, probe_cfg, scratch);
                    if (!mid.valid) break;
                    const Real f_mid = mid.depth_m - receiver_depth_m;
                    if ((f_lo_v <= kZero && f_mid <= kZero)
                        || (f_lo_v >= kZero && f_mid >= kZero)) {
                        lo = mid_angle;
                        f_lo_v = f_mid;
                    } else {
                        hi = mid_angle;
                    }
                }

                if (mid.valid
                    && std::fabs(mid.depth_m - receiver_depth_m) <= search.depth_tolerance_m) {
                    // The Jacobian dz/dtheta_0, by central difference over a
                    // step small enough to be local but large enough not to be
                    // dominated by the tracer's own rounding.
                    const Real h = d_angle / static_cast<Real>(8);
                    const Probe a = probe_at(svp, source_depth_m, mid_angle - h, probe_cfg, scratch);
                    const Probe b = probe_at(svp, source_depth_m, mid_angle + h, probe_cfg, scratch);

                    Eigenray& e = out[written];
                    e.launch_angle_rad = mid_angle;
                    e.arrival_angle_rad = mid.angle_rad;
                    e.travel_time_s = mid.time_s;
                    e.path_length_m = mid.path_m;
                    e.arrival_speed_mps = mid.speed_mps;
                    e.surface_bounces = mid.surface_bounces;
                    e.bottom_bounces = mid.bottom_bounces;

                    if (a.valid && b.valid) {
                        e.jacobian_m_per_rad =
                            (b.depth_m - a.depth_m) / (static_cast<Real>(2) * h);
                    } else {
                        e.jacobian_m_per_rad = kZero;
                    }
                    e.near_caustic = std::fabs(e.jacobian_m_per_rad)
                                   < search.caustic_jacobian_m_per_rad;
                    e.spreading_loss_db =
                        e.near_caustic
                            ? kZero
                            : spreading_loss_db(range_m, e.launch_angle_rad, e.arrival_angle_rad,
                                                e.jacobian_m_per_rad, c_src, e.arrival_speed_mps);
                    ++written;
                }
            }
        }

        prev = cur;
        prev_angle = angle;
    }

    return written;
}

}  // namespace phantom
