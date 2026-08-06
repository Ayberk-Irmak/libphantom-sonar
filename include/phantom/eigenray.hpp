// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — eigenrays and transmission loss.
//
// This is where the two halves of the library meet. Up to now the ocean model
// and the pulse chain shared only a sound speed; here the ray tracer supplies
// the delays and levels that the echo synthesiser transmits, so a multipath
// arrival structure comes out of the physics rather than out of a parameter.
//
// An eigenray is a launch angle whose ray passes through a given receiver. The
// search is exact rather than interpolated: the tracer already stops on a range
// budget to machine precision, so setting that budget to the receiver range and
// reading the final depth gives z(theta_0) with no chord-versus-arc error. What
// remains is a root find on z(theta_0) - z_receiver.
#ifndef PHANTOM_EIGENRAY_HPP
#define PHANTOM_EIGENRAY_HPP

#include "phantom/boundary.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/types.hpp"

#include <span>

namespace phantom {

struct Eigenray {
    Real launch_angle_rad = 0;
    Real arrival_angle_rad = 0;
    Real travel_time_s = 0;
    Real path_length_m = 0;
    Real arrival_speed_mps = 0;

    // |dz/dtheta_0| at the receiver, in metres per radian. This is the ray-tube
    // Jacobian: the whole of geometric spreading is in it.
    Real jacobian_m_per_rad = 0;

    // Spreading loss alone, dB re 1 m. Absorption is not included -- it depends
    // on frequency, which an eigenray does not have.
    Real spreading_loss_db = 0;

    std::uint32_t surface_bounces = 0;
    std::uint32_t bottom_bounces = 0;

    // xi = cos(theta)/c, constant along the ray. In a range-independent ocean
    // this recovers the grazing angle at ANY depth, so the angle at each
    // boundary bounce is exact without having been tracked:
    //     cos(theta_boundary) = xi * c(z_boundary)
    Real snell_invariant = 0;

    // A caustic is where neighbouring rays cross and the tube collapses. Ray
    // theory predicts infinite intensity there, which is a known failure of the
    // method rather than a property of the ocean. When this is set the level is
    // not trustworthy and the caller should say so rather than plot it.
    bool near_caustic = false;
};

// Absorption in seawater, Thorp (1967), dB per kilometre:
//
//   a(f) = 0.11 f^2/(1 + f^2) + 44 f^2/(4100 + f^2) + 2.75e-4 f^2 + 0.003
//
// with f in kHz. The four terms are boric acid relaxation, magnesium sulphate
// relaxation, pure-water viscosity, and a low-frequency floor. Valid roughly
// 100 Hz to 1 MHz; it is a fit to measurements at about 4 degrees C and does
// not carry temperature or depth dependence, so treat it as within ~10%.
//
// The scale of it is what matters for design: 0.07 dB/km at 1 kHz, 1.2 dB/km at
// 10 kHz, 34 dB/km at 100 kHz. Absorption is why long-range sonar is
// low-frequency and why a tank experiment can use 200 kHz without noticing.
[[nodiscard]] Real thorp_absorption_db_per_km(Real frequency_hz) noexcept;

// Spreading loss for one ray tube, dB re 1 m:
//
//   TL = -10 log10[ c_rcv cos(theta_0) / (c_src * r * cos(theta_rcv) * |dz/dtheta_0|) ]
//
// Derived by conserving power in the tube: the source radiates
// P cos(theta_0) dtheta_0 / 2 into the fan element, and at range r that power
// crosses an annulus of circumference 2 pi r and cross-section dz cos(theta).
// In an isovelocity ocean this reduces exactly to 20 log10(R) with R the slant
// range, which is what the test suite checks it against.
[[nodiscard]] Real spreading_loss_db(Real range_m,
                                     Real launch_angle_rad,
                                     Real arrival_angle_rad,
                                     Real jacobian_m_per_rad,
                                     Real source_speed_mps,
                                     Real arrival_speed_mps) noexcept;

// Total one-way transmission loss for a path at a given frequency:
// spreading plus Thorp absorption over the arc length. Boundary losses are NOT
// included -- see boundary_loss_db and total_transmission_loss_db.
[[nodiscard]] Real transmission_loss_db(const Eigenray& ray, Real range_m,
                                        Real frequency_hz, Real source_speed_mps) noexcept;

// Loss from every boundary interaction along the path, dB.
//
// The grazing angle at each bounce comes from the Snell invariant rather than
// from having recorded it, which is exact in a range-independent ocean:
// cos(theta) = xi * c at the boundary depth. Each surface bounce costs
// surface_loss_db and each bottom bounce bottom_loss_db, at that angle.
//
// Sub-critical bottom bounces cost almost nothing and super-critical ones cost
// several dB apiece, so a path with four bounces can be anywhere from
// negligible to gone depending only on its launch angle.
[[nodiscard]] Real boundary_loss_db(const Eigenray& ray,
                                    const BoundaryModel& model,
                                    Real surface_speed_mps,
                                    Real bottom_speed_mps,
                                    Real frequency_hz) noexcept;

// Spreading + absorption + boundaries: what a path actually costs.
[[nodiscard]] Real total_transmission_loss_db(const Eigenray& ray, Real range_m,
                                              Real frequency_hz,
                                              Real source_speed_mps,
                                              const BoundaryModel& model,
                                              Real surface_speed_mps,
                                              Real bottom_speed_mps) noexcept;

struct EigenraySearch {
    Real angle_min_rad = deg2rad(static_cast<Real>(-45));
    Real angle_max_rad = deg2rad(static_cast<Real>(45));
    // Fan resolution. Two eigenrays closer together in launch angle than one
    // fan step are seen as one; a fan too coarse silently loses paths, so this
    // is the knob that decides completeness.
    std::size_t fan_count = 721;
    // Bisection iterations used to polish each bracketed root.
    std::size_t refine_steps = 40;

    // Depth tolerance at the receiver, metres. A root that cannot be polished
    // to inside this is DISCARDED -- so a tolerance tighter than the build can
    // achieve does not degrade the answer, it silently removes paths.
    //
    // Measured on a 200 m duct at 3 km, 14 eigenrays present:
    //
    //   tolerance   double   float
    //     0.001 m      14       0
    //     0.010 m      14       3
    //     0.050 m      14      11
    //     0.100 m      14      14
    //
    // Double precision resolves the receiver to well under a millimetre; a
    // single-precision build cannot do better than about 0.1 m at this range,
    // because the tracer's own depth error there is larger than that. The
    // default follows the build for exactly this reason.
#if defined(PHANTOM_REAL_FLOAT)
    Real depth_tolerance_m = static_cast<Real>(0.1);
#else
    Real depth_tolerance_m = static_cast<Real>(0.01);
#endif
    // |dz/dtheta_0| below this counts as a caustic.
    Real caustic_jacobian_m_per_rad = static_cast<Real>(1);
};

// Finds eigenrays from (0, source_depth_m) to (range_m, receiver_depth_m).
//
// `scratch` is the tracer's per-ray point buffer; nothing is allocated. Results
// are written in increasing launch angle. Returns how many were written, which
// may be fewer than found if `out` fills.
std::size_t find_eigenrays(const ProfileView& svp,
                           Real source_depth_m,
                           Real receiver_depth_m,
                           Real range_m,
                           const TraceConfig& cfg,
                           const EigenraySearch& search,
                           std::span<RayPoint> scratch,
                           std::span<Eigenray> out) noexcept;

}  // namespace phantom

#endif  // PHANTOM_EIGENRAY_HPP
