// SPDX-License-Identifier: Apache-2.0
#include "phantom/boundary.hpp"

#include <cmath>
#include <complex>

namespace phantom {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);
constexpr Real kTwo  = static_cast<Real>(2);

// A reflection loss beyond this is not physically distinguishable from "the
// path is gone", and returning an infinity would poison every sum it enters.
constexpr Real kMaxLossDb = static_cast<Real>(60);

// Attenuation of alpha dB per wavelength corresponds to a fractional imaginary
// part of the wavenumber of alpha / (40 pi log10 e) = alpha / 54.58:
//   alpha_dB/lambda = 20 log10(e) * alpha_np * lambda = 8.686 * alpha_np * lambda
//   k = (w/c)(1 - i alpha_np c / w) = (w/c)(1 - i alpha / (8.686 * 2 pi))
constexpr Real kDbPerWavelengthToDelta = static_cast<Real>(54.5751);

}  // namespace

Real bottom_critical_angle_rad(const BottomProperties& bottom,
                               Real water_speed_mps) noexcept {
    if (!(water_speed_mps > kZero) || !(bottom.sound_speed_mps > kZero)) return kZero;
    const Real ratio = water_speed_mps / bottom.sound_speed_mps;
    // A bottom slower than the water has no critical angle at all: there is no
    // total internal reflection and it leaks at every grazing angle.
    if (!(ratio < kOne)) return kZero;
    return std::acos(ratio);
}

Real bottom_reflection_coefficient(const BottomProperties& bottom,
                                   Real grazing_angle_rad,
                                   Real water_speed_mps) noexcept {
    if (!(water_speed_mps > kZero) || !(bottom.sound_speed_mps > kZero)) return kOne;
    if (!(bottom.density_ratio > kZero)) return kOne;

    const Real th = std::fabs(clamp(grazing_angle_rad, -kHalfPi, kHalfPi));
    const Real sin_th = std::sin(th);
    const Real cos_th = std::cos(th);

    // Complex sediment speed carries the attenuation. c2' = c2 / (1 - i delta)
    // so that the wave decays as it propagates into the sediment.
    const Real delta = (bottom.attenuation_db_per_wavelength > kZero)
                     ? bottom.attenuation_db_per_wavelength / kDbPerWavelengthToDelta
                     : kZero;
    using C = std::complex<Real>;
    const C c2 = C(bottom.sound_speed_mps, kZero) / C(kOne, -delta);
    const C nu = c2 / C(water_speed_mps, kZero);

    // sin(theta_2) from Snell in the sediment. Beyond the critical angle the
    // argument goes negative and the root is imaginary: an evanescent wave.
    C sin_th2 = std::sqrt(C(kOne, kZero) - nu * nu * C(cos_th * cos_th, kZero));

    // Branch choice. The transmitted wave must decay into the sediment, not
    // grow, and the principal square root does not always give that. Rather
    // than reason about sign conventions, take the branch that conserves
    // energy: the other one yields |R| > 1, which no passive interface can do.
    const C num = C(bottom.density_ratio, kZero) * nu * C(sin_th, kZero) - sin_th2;
    const C den = C(bottom.density_ratio, kZero) * nu * C(sin_th, kZero) + sin_th2;
    if (std::abs(den) <= kZero) return kOne;
    Real r = std::abs(num / den);
    if (r > kOne) {
        sin_th2 = -sin_th2;
        const C num2 = C(bottom.density_ratio, kZero) * nu * C(sin_th, kZero) - sin_th2;
        const C den2 = C(bottom.density_ratio, kZero) * nu * C(sin_th, kZero) + sin_th2;
        if (std::abs(den2) > kZero) {
            const Real r2 = std::abs(num2 / den2);
            if (r2 <= kOne) r = r2;
        }
    }
    return (r > kOne) ? kOne : r;
}

Real bottom_loss_db(const BottomProperties& bottom,
                    Real grazing_angle_rad,
                    Real water_speed_mps) noexcept {
    const Real r = bottom_reflection_coefficient(bottom, grazing_angle_rad, water_speed_mps);
    if (!(r > kZero)) return kMaxLossDb;
    const Real loss = -static_cast<Real>(20) * std::log10(r);
    return (loss > kMaxLossDb) ? kMaxLossDb : ((loss > kZero) ? loss : kZero);
}

Real wind_to_rms_wave_height_m(Real wind_speed_mps) noexcept {
    if (!(wind_speed_mps > kZero)) return kZero;
    // Pierson-Moskowitz fully developed sea, then RMS = H_1/3 / 4.
    const Real h_significant = static_cast<Real>(0.0246) * wind_speed_mps * wind_speed_mps;
    return h_significant / static_cast<Real>(4);
}

Real surface_reflection_coefficient(const SurfaceProperties& surface,
                                    Real grazing_angle_rad,
                                    Real frequency_hz,
                                    Real water_speed_mps) noexcept {
    if (!(surface.rms_wave_height_m > kZero)) return kOne;   // flat: perfect
    if (!(frequency_hz > kZero) || !(water_speed_mps > kZero)) return kOne;

    const Real th = std::fabs(clamp(grazing_angle_rad, -kHalfPi, kHalfPi));
    const Real k = kTwo * kPi * frequency_hz / water_speed_mps;
    const Real gamma = kTwo * k * surface.rms_wave_height_m * std::sin(th);
    // exp(-G^2/2) underflows long before it stops mattering; clamp the exponent
    // so the result is a small number rather than a denormal or a zero that the
    // caller then takes a log of.
    const Real arg = gamma * gamma / kTwo;
    if (arg > static_cast<Real>(80)) return static_cast<Real>(1e-35);
    return std::exp(-arg);
}

Real surface_loss_db(const SurfaceProperties& surface,
                     Real grazing_angle_rad,
                     Real frequency_hz,
                     Real water_speed_mps) noexcept {
    const Real r = surface_reflection_coefficient(surface, grazing_angle_rad,
                                                  frequency_hz, water_speed_mps);
    const Real cap = (surface.max_loss_db > kZero) ? surface.max_loss_db : kMaxLossDb;
    if (!(r > kZero)) return cap;
    const Real loss = -static_cast<Real>(20) * std::log10(r);
    return (loss > cap) ? cap : ((loss > kZero) ? loss : kZero);
}

}  // namespace phantom
