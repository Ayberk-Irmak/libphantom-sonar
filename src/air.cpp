// SPDX-License-Identifier: Apache-2.0
#include "phantom/air.hpp"

#include <cmath>

namespace phantom::air {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// ISO 9613-1:1993 reference conditions.
constexpr Real kPr  = static_cast<Real>(101.325);   // kPa
constexpr Real kT0  = static_cast<Real>(293.15);    // K, 20 C
constexpr Real kT01 = static_cast<Real>(273.16);    // K, triple point

Real kelvin(Real temperature_c) noexcept {
    return temperature_c + static_cast<Real>(273.15);
}

}  // namespace

Real water_vapour_molar_percent(Real temperature_c, Real relative_humidity_percent,
                                Real pressure_kpa) noexcept {
    if (!(pressure_kpa > kZero)) return kZero;
    const Real t = kelvin(temperature_c);
    // ISO 9613-1 annex B: saturation vapour pressure relative to the reference.
    const Real psat_over_pr = std::pow(static_cast<Real>(10),
        -static_cast<Real>(6.8346) * std::pow(kT01 / t, static_cast<Real>(1.261))
        + static_cast<Real>(4.6151));
    return relative_humidity_percent * psat_over_pr * kPr / pressure_kpa;
}

Real oxygen_relaxation_hz(Real temperature_c, Real relative_humidity_percent,
                          Real pressure_kpa) noexcept {
    const Real h = water_vapour_molar_percent(temperature_c, relative_humidity_percent, pressure_kpa);
    const Real pa_pr = pressure_kpa / kPr;
    // ISO 9613-1 equation (3).
    return pa_pr * (static_cast<Real>(24)
        + static_cast<Real>(4.04e4) * h * (static_cast<Real>(0.02) + h)
          / (static_cast<Real>(0.391) + h));
}

Real nitrogen_relaxation_hz(Real temperature_c, Real relative_humidity_percent,
                            Real pressure_kpa) noexcept {
    const Real h = water_vapour_molar_percent(temperature_c, relative_humidity_percent, pressure_kpa);
    const Real pa_pr = pressure_kpa / kPr;
    const Real tr = kelvin(temperature_c) / kT0;
    // ISO 9613-1 equation (4).
    return pa_pr * std::pow(tr, static_cast<Real>(-0.5))
         * (static_cast<Real>(9) + static_cast<Real>(280) * h
            * std::exp(-static_cast<Real>(4.170)
                       * (std::pow(tr, -kOne / static_cast<Real>(3)) - kOne)));
}

Real absorption_db_per_km(Real frequency_hz, Real temperature_c,
                          Real relative_humidity_percent, Real pressure_kpa) noexcept {
    if (!(frequency_hz > kZero) || !(pressure_kpa > kZero)) return kZero;
    const Real t = kelvin(temperature_c);
    const Real tr = t / kT0;
    const Real pa_pr = pressure_kpa / kPr;
    const Real frO = oxygen_relaxation_hz(temperature_c, relative_humidity_percent, pressure_kpa);
    const Real frN = nitrogen_relaxation_hz(temperature_c, relative_humidity_percent, pressure_kpa);
    const Real f2 = frequency_hz * frequency_hz;

    // ISO 9613-1 equation (5). The first term is classical absorption --
    // viscosity and heat conduction -- and is negligible below about 100 kHz;
    // the other two are the molecular relaxation that actually sets the loss.
    const Real classical = static_cast<Real>(1.84e-11) * (kOne / pa_pr)
                         * std::pow(tr, static_cast<Real>(0.5));
    const Real relax = std::pow(tr, static_cast<Real>(-2.5))
        * (static_cast<Real>(0.01275) * std::exp(-static_cast<Real>(2239.1) / t)
             / (frO + f2 / frO)
         + static_cast<Real>(0.1068) * std::exp(-static_cast<Real>(3352.0) / t)
             / (frN + f2 / frN));

    // The standard gives dB/m; the table, and everyone quoting it, uses dB/km.
    return static_cast<Real>(8.686) * f2 * (classical + relax) * static_cast<Real>(1000);
}

Real midband_frequency_hz(int k) noexcept {
    return static_cast<Real>(1000) * std::pow(static_cast<Real>(10),
                                              static_cast<Real>(k) / static_cast<Real>(10));
}

Real sound_speed(Real temperature_c, Real relative_humidity_percent,
                 Real pressure_kpa) noexcept {
    // Dry air: c = c0 * sqrt(T / T_ice), with c0 = 331.45 m/s at 0 C.
    const Real t = kelvin(temperature_c);
    const Real dry = static_cast<Real>(331.45)
                   * std::sqrt(t / static_cast<Real>(273.15));
    // Humidity raises the speed: water vapour (18 g/mol) is lighter than the
    // air it displaces (29 g/mol), so the mixture's mean molar mass falls.
    // First order in the vapour mole fraction, which is enough at sea level.
    const Real h = water_vapour_molar_percent(temperature_c, relative_humidity_percent,
                                              pressure_kpa) / static_cast<Real>(100);
    return dry * (kOne + static_cast<Real>(0.16) * h);
}

Real density_kgm3(Real temperature_c, Real relative_humidity_percent,
                  Real pressure_kpa) noexcept {
    const Real t = kelvin(temperature_c);
    const Real h = water_vapour_molar_percent(temperature_c, relative_humidity_percent,
                                              pressure_kpa) / static_cast<Real>(100);
    // Mean molar mass of the mixture, kg/mol, then the ideal gas law.
    const Real m = (kOne - h) * static_cast<Real>(0.0289645)
                 + h * static_cast<Real>(0.018015);
    return pressure_kpa * static_cast<Real>(1000) * m / (static_cast<Real>(8.31446) * t);
}

Real impedance_rayl(Real temperature_c, Real relative_humidity_percent,
                    Real pressure_kpa) noexcept {
    return density_kgm3(temperature_c, relative_humidity_percent, pressure_kpa)
         * sound_speed(temperature_c, relative_humidity_percent, pressure_kpa);
}

Real doppler_scale(Real closing_speed_mps, Real temperature_c) noexcept {
    const Real c = sound_speed(temperature_c);
    return (c > kZero) ? closing_speed_mps / c : kZero;
}

}  // namespace phantom::air
