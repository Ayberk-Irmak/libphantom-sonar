// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — sound speed in seawater.
//
// Three independent equations are provided on purpose. They disagree by
// O(0.1 m/s) inside their common validity box, and that disagreement is used as
// a self-check in the unit tests: a single mistyped coefficient breaks the
// agreement immediately.
//
//   medwin()       Medwin (1975), 6 terms. Fast, constexpr, valid to ~1000 m.
//   mackenzie()    Mackenzie (1981), 9 terms. Valid to 8000 m. DEFAULT.
//   chen_millero() Chen & Millero (1977) — the actual UNESCO algorithm.
//                  Takes PRESSURE, not depth. Reference / validation path.
//
// A note on a common error: the six-term polynomial that circulates as "the
// UNESCO equation" is Medwin's simplification, not UNESCO. UNESCO is
// Chen-Millero. Medwin's depth term is linear (0.016 z) and its stated validity
// stops at 1000 m, which is at or above the deep sound channel axis at most
// latitudes -- so Medwin is the wrong tool for SOFAR work. Use mackenzie().
#ifndef PHANTOM_SOUND_SPEED_HPP
#define PHANTOM_SOUND_SPEED_HPP

#include "phantom/types.hpp"

#include <cmath>

namespace phantom::sound_speed {

// ---------------------------------------------------------------------------
// Medwin (1975). T in degrees C, S in PSU, z in metres.
// Validity: 0 <= T <= 35, 0 <= S <= 45, 0 <= z <= 1000.
// ---------------------------------------------------------------------------
constexpr Real medwin(Real temperature_c, Real salinity_psu, Real depth_m) noexcept {
    const Real t = temperature_c;
    return static_cast<Real>(1449.2)
         + static_cast<Real>(4.6) * t
         - static_cast<Real>(0.055) * t * t
         + static_cast<Real>(0.00029) * t * t * t
         + (static_cast<Real>(1.34) - static_cast<Real>(0.010) * t)
             * (salinity_psu - static_cast<Real>(35))
         + static_cast<Real>(0.016) * depth_m;
}

// ---------------------------------------------------------------------------
// Mackenzie (1981), nine-term. T in degrees C, S in PSU, D in metres.
// Validity: 2 <= T <= 30, 25 <= S <= 40, 0 <= D <= 8000.
// This is the library default: it carries a quadratic depth term and a
// T*D^3 cross term, so it stays honest through the deep sound channel.
// ---------------------------------------------------------------------------
constexpr Real mackenzie(Real temperature_c, Real salinity_psu, Real depth_m) noexcept {
    const Real t  = temperature_c;
    const Real ds = salinity_psu - static_cast<Real>(35);
    const Real d  = depth_m;
    return static_cast<Real>(1448.96)
         + static_cast<Real>(4.591) * t
         - static_cast<Real>(5.304e-2) * t * t
         + static_cast<Real>(2.374e-4) * t * t * t
         + static_cast<Real>(1.340) * ds
         + static_cast<Real>(1.630e-2) * d
         + static_cast<Real>(1.675e-7) * d * d
         - static_cast<Real>(1.025e-2) * t * ds
         - static_cast<Real>(7.139e-13) * t * d * d * d;
}

// ---------------------------------------------------------------------------
// Chen & Millero (1977) — UNESCO algorithm.
// T in degrees C, S in PSU, P in BAR (gauge, i.e. 0 at the surface).
// Validity: 0 <= T <= 40, 0 <= S <= 40, 0 <= P <= 1000 bar.
//   c = Cw(T,P) + A(T,P)*S + B(T,P)*S^1.5 + D(T,P)*S^2
//
// Not constexpr: the S^1.5 term needs sqrt, which C++20 does not make constant
// evaluable (C++26 does). medwin() and mackenzie() are polynomial and stay
// constexpr, so compile-time profile generation still works with those.
// ---------------------------------------------------------------------------
inline Real chen_millero(Real temperature_c, Real salinity_psu, Real pressure_bar) noexcept {
    const Real t = temperature_c;
    const Real p = pressure_bar;
    const Real t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
    const Real p2 = p * p, p3 = p2 * p;

    // Pure-water term.
    const Real cw0 = static_cast<Real>(1402.388)
                   + static_cast<Real>(5.03830) * t
                   - static_cast<Real>(5.81090e-2) * t2
                   + static_cast<Real>(3.3432e-4) * t3
                   - static_cast<Real>(1.47797e-6) * t4
                   + static_cast<Real>(3.1419e-9) * t5;
    const Real cw1 = static_cast<Real>(0.153563)
                   + static_cast<Real>(6.8999e-4) * t
                   - static_cast<Real>(8.1829e-6) * t2
                   + static_cast<Real>(1.3632e-7) * t3
                   - static_cast<Real>(6.1260e-10) * t4;
    const Real cw2 = static_cast<Real>(3.1260e-5)
                   - static_cast<Real>(1.7111e-6) * t
                   + static_cast<Real>(2.5986e-8) * t2
                   - static_cast<Real>(2.5353e-10) * t3
                   + static_cast<Real>(1.0415e-12) * t4;
    const Real cw3 = -static_cast<Real>(9.7729e-9)
                   + static_cast<Real>(3.8513e-10) * t
                   - static_cast<Real>(2.3654e-12) * t2;
    const Real cw = cw0 + cw1 * p + cw2 * p2 + cw3 * p3;

    // Linear salinity term.
    const Real a0 = static_cast<Real>(1.389)
                  - static_cast<Real>(1.262e-2) * t
                  + static_cast<Real>(7.164e-5) * t2
                  + static_cast<Real>(2.006e-6) * t3
                  - static_cast<Real>(3.21e-8) * t4;
    const Real a1 = static_cast<Real>(9.4742e-5)
                  - static_cast<Real>(1.2583e-5) * t
                  - static_cast<Real>(6.4928e-8) * t2
                  + static_cast<Real>(1.0515e-8) * t3
                  - static_cast<Real>(2.0142e-10) * t4;
    const Real a2 = -static_cast<Real>(3.9064e-7)
                  + static_cast<Real>(9.1061e-9) * t
                  - static_cast<Real>(1.6009e-10) * t2
                  + static_cast<Real>(7.994e-12) * t3;
    const Real a3 = static_cast<Real>(1.100e-10)
                  + static_cast<Real>(6.651e-12) * t
                  - static_cast<Real>(3.391e-13) * t2;
    const Real a = a0 + a1 * p + a2 * p2 + a3 * p3;

    // S^1.5 term.
    const Real b0 = -static_cast<Real>(1.922e-2) - static_cast<Real>(4.42e-5) * t;
    const Real b1 = static_cast<Real>(7.3637e-5) + static_cast<Real>(1.7950e-7) * t;
    const Real b  = b0 + b1 * p;

    // S^2 term.
    const Real d = static_cast<Real>(1.727e-3) - static_cast<Real>(7.9836e-6) * p;

    const Real s   = salinity_psu;
    const Real s15 = (s > 0) ? s * std::sqrt(s) : static_cast<Real>(0);
    return cw + a * s + b * s15 + d * s * s;
}

// ---------------------------------------------------------------------------
// Depth <-> pressure, Leroy & Parthiot (1998) international formula.
// Returns gauge pressure in bar. Latitude matters at the 0.5 bar level over
// 1000 m, which is ~0.008 m/s of sound speed -- small, but free to get right.
// ---------------------------------------------------------------------------
inline Real depth_to_pressure_bar(Real depth_m, Real latitude_deg = static_cast<Real>(45)) noexcept {
    const Real z   = depth_m;
    const Real phi = latitude_deg * kDeg2Rad;
    const Real sp  = std::sin(phi);
    // h45: pressure in MPa at 45 degrees latitude.
    const Real h45 = static_cast<Real>(1.00818e-2) * z
                   + static_cast<Real>(2.465e-8)  * z * z
                   - static_cast<Real>(1.25e-13)  * z * z * z
                   + static_cast<Real>(2.8e-19)   * z * z * z * z;
    const Real gphi = static_cast<Real>(9.7803) * (static_cast<Real>(1) + static_cast<Real>(5.3e-3) * sp * sp);
    const Real k = (gphi - static_cast<Real>(2e-5) * z)
                 / (static_cast<Real>(9.80612) - static_cast<Real>(2e-5) * z);
    return static_cast<Real>(10) * h45 * k;  // MPa -> bar
}

// Convenience wrapper: UNESCO sound speed from depth instead of pressure.
inline Real unesco(Real temperature_c, Real salinity_psu, Real depth_m,
                   Real latitude_deg = static_cast<Real>(45)) noexcept {
    return chen_millero(temperature_c, salinity_psu, depth_to_pressure_bar(depth_m, latitude_deg));
}

// ---------------------------------------------------------------------------
// Munk (1974) canonical deep sound channel profile. The standard benchmark
// case for ray tracers -- Bellhop ships the same profile as MunkB_ray.
//   c(z) = c1 * [1 + eps * (eta + exp(-eta) - 1)],  eta = 2 (z - z1) / B
// ---------------------------------------------------------------------------
inline Real munk(Real depth_m,
                 Real axis_depth_m  = static_cast<Real>(1300),
                 Real axis_speed_mps = static_cast<Real>(1500),
                 Real epsilon        = static_cast<Real>(7.37e-3),
                 Real scale_m        = static_cast<Real>(1300)) noexcept {
    const Real eta = static_cast<Real>(2) * (depth_m - axis_depth_m) / scale_m;
    return axis_speed_mps
         * (static_cast<Real>(1) + epsilon * (eta + std::exp(-eta) - static_cast<Real>(1)));
}

}  // namespace phantom::sound_speed

#endif  // PHANTOM_SOUND_SPEED_HPP
