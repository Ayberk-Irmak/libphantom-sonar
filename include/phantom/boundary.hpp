// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — surface and bottom reflection losses.
//
// Until now every boundary in this library reflected perfectly. That made
// bounced paths louder than they are, and in shallow water most paths bounce,
// so it was the largest remaining lie in the transmission loss.
//
// Two mechanisms, and they fail in opposite directions:
//
//   BOTTOM   a fluid-fluid interface. Below a critical grazing angle the wave
//            is totally reflected and the water column traps it; above it,
//            energy leaks into the sediment. The critical angle is what decides
//            which rays survive to long range in shallow water, and it depends
//            only on the speed ratio.
//
//   SURFACE  a pressure-release boundary, perfect when flat. Roughness scatters
//            the coherent (specular) component away exponentially in
//            frequency, angle and wave height -- so at any real sea state a
//            high-frequency grazing path loses its specular return entirely.
//
// Both are plane-wave, flat-interface results applied per bounce. That is the
// standard ray-theory treatment and it is an approximation; the limits are
// stated at each function.
#ifndef PHANTOM_BOUNDARY_HPP
#define PHANTOM_BOUNDARY_HPP

#include "phantom/types.hpp"

namespace phantom {

// A fluid sediment. Defaults are a medium sand: faster and denser than water,
// with the modest attenuation such sediments show.
struct BottomProperties {
    Real sound_speed_mps = 1650;               // compressional speed in sediment
    Real density_ratio = static_cast<Real>(1.9);  // rho_sediment / rho_water
    // Attenuation inside the sediment, dB per wavelength. This is the quantity
    // sediment tables publish, and it is what makes the reflection lossy even
    // below the critical angle -- with it set to zero, sub-critical rays are
    // trapped forever and shallow-water range becomes unbounded.
    Real attenuation_db_per_wavelength = static_cast<Real>(0.8);
};

struct SurfaceProperties {
    // RMS wave height. Use wind_to_rms_wave_height_m() if you have wind instead.
    Real rms_wave_height_m = 0;
    // The coherent coefficient falls off without limit, but the energy does not
    // vanish -- it goes into the diffuse field, which a ray model does not
    // carry. Capping the reported loss acknowledges that rather than returning
    // a path 200 dB down that would not be missing in the water.
    Real max_loss_db = 30;
};

// Critical grazing angle, arccos(c_water / c_bottom), radians.
//
// Zero when the bottom is slower than the water: a slow bottom has no critical
// angle and leaks at every angle, which is why a mud bottom is acoustically far
// worse than sand.
[[nodiscard]] Real bottom_critical_angle_rad(const BottomProperties& bottom,
                                             Real water_speed_mps) noexcept;

// Rayleigh plane-wave reflection coefficient magnitude at the water-sediment
// interface, for a grazing angle measured from the interface.
//
//   R = (m nu sin(th) - sin(th_2)) / (m nu sin(th) + sin(th_2))
//   sin(th_2) = sqrt(1 - nu^2 cos^2(th)),   nu = c2/c1,  m = rho2/rho1
//
// With attenuation, c2 is complex and so is R, so |R| < 1 at every angle.
// At normal incidence this reduces to the impedance ratio form
// (Z2 - Z1)/(Z2 + Z1); below the critical angle and with no attenuation it is
// exactly 1.
//
// Plane-wave and flat-interface: it ignores the finite beam width that produces
// beam displacement near the critical angle, and any sediment layering.
[[nodiscard]] Real bottom_reflection_coefficient(const BottomProperties& bottom,
                                                 Real grazing_angle_rad,
                                                 Real water_speed_mps) noexcept;

// The same, expressed as a positive loss in dB. Capped so a total-transmission
// geometry does not return an infinity.
[[nodiscard]] Real bottom_loss_db(const BottomProperties& bottom,
                                  Real grazing_angle_rad,
                                  Real water_speed_mps) noexcept;

// RMS wave height from wind speed, via the Pierson-Moskowitz fully developed
// sea: H_1/3 = 0.0246 U^2 with U in m/s, and sigma = H_1/3 / 4.
//
// "Fully developed" is doing real work in that sentence -- it assumes the wind
// has blown long enough over enough fetch. A rising wind gives a lower sea than
// this predicts.
[[nodiscard]] Real wind_to_rms_wave_height_m(Real wind_speed_mps) noexcept;

// Coherent (specular) reflection coefficient magnitude at a rough
// pressure-release surface:
//
//   |R| = exp(-G^2 / 2),   G = 2 k sigma sin(theta),  k = 2 pi f / c
//
// G is the Rayleigh roughness parameter. This is the SPECULAR component only:
// the scattered energy is not destroyed, it goes into the diffuse field, which
// a ray model does not represent. So this over-predicts loss for a real
// receiver, which is why SurfaceProperties carries a cap.
[[nodiscard]] Real surface_reflection_coefficient(const SurfaceProperties& surface,
                                                  Real grazing_angle_rad,
                                                  Real frequency_hz,
                                                  Real water_speed_mps) noexcept;

// The same as a positive loss in dB, capped at `surface.max_loss_db`.
[[nodiscard]] Real surface_loss_db(const SurfaceProperties& surface,
                                   Real grazing_angle_rad,
                                   Real frequency_hz,
                                   Real water_speed_mps) noexcept;

struct BoundaryModel {
    BottomProperties bottom;
    SurfaceProperties surface;
};

}  // namespace phantom

#endif  // PHANTOM_BOUNDARY_HPP
