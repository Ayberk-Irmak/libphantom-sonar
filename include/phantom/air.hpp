// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — acoustics in air.
//
// WHY AIR IS IN AN UNDERWATER LIBRARY.
//
// Because it is the only medium most people can test in. A hydrophone, a
// projector and an amplifier cost several hundred euros and need water to put
// them in; a speaker and a microphone are already on the desk. The signal
// processing does not care which medium it is in -- a matched filter, a Doppler
// bank, a CFAR detector and a tracker are the same code at 340 m/s as at
// 1500 m/s -- so an air bench exercises the whole chain end to end before any
// wet hardware is bought.
//
// WHAT ACTUALLY CHANGES. Three things, and they are the reason this header
// exists rather than a note saying "set c = 343":
//
//   SOUND SPEED   ~343 m/s instead of ~1500, so every range is 4.4x shorter for
//                 the same travel time, and every Doppler scaling is 4.4x
//                 LARGER for the same closing speed. Air is the harsher Doppler
//                 environment, which makes it a good place to test the
//                 time-scaling machinery of the comm module.
//
//   ABSORPTION    dominated by molecular relaxation of oxygen and nitrogen, and
//                 therefore strongly dependent on HUMIDITY -- not a small
//                 correction. At 10 C and 4 kHz, going from 10% to 20% relative
//                 humidity nearly doubles the loss. Nothing in seawater behaves
//                 like that.
//
//   IMPEDANCE     ~415 rayl against seawater's ~1.5e6, so a target's echo
//                 strength and every boundary reflection coefficient differ by
//                 orders of magnitude. This header does not model targets; it
//                 gives the impedance so a caller can.
#ifndef PHANTOM_AIR_HPP
#define PHANTOM_AIR_HPP

#include "phantom/types.hpp"

namespace phantom::air {

// Standard sea-level pressure, kPa.
inline constexpr Real kStandardPressureKpa = static_cast<Real>(101.325);

// Speed of sound in air, m/s.
//
// The dry-air part is the textbook square-root law about the 0 C value of
// 331.45 m/s; the humidity term is a first-order correction, since water vapour
// is lighter than the air it displaces and so raises the speed slightly.
//
// This is NOT Cramer (1993), which carries CO2 concentration and a pressure
// dependence and is accurate to ~300 ppm. Over 0-40 C at sea level the
// difference is a few tenths of a m/s -- about 0.1% -- which is smaller than
// the temperature uncertainty of any bench measurement. The tests state the
// achieved agreement rather than claiming the standard.
[[nodiscard]] Real sound_speed(Real temperature_c,
                               Real relative_humidity_percent = static_cast<Real>(0),
                               Real pressure_kpa = kStandardPressureKpa) noexcept;

// Molar concentration of water vapour, in percent, from relative humidity.
// This is the `h` that drives both relaxation frequencies below, and the
// quantity the ISO standard is actually written in terms of.
[[nodiscard]] Real water_vapour_molar_percent(Real temperature_c, Real relative_humidity_percent,
                                              Real pressure_kpa = kStandardPressureKpa) noexcept;

// Oxygen and nitrogen relaxation frequencies, Hz. ISO 9613-1:1993 equations
// (3) and (4). Exposed because they are the whole shape of the absorption
// curve: below f_rO absorption rises as f^2, between the two it flattens, and
// above f_rN it is f^2 again with a different constant.
[[nodiscard]] Real oxygen_relaxation_hz(Real temperature_c, Real relative_humidity_percent,
                                        Real pressure_kpa = kStandardPressureKpa) noexcept;
[[nodiscard]] Real nitrogen_relaxation_hz(Real temperature_c, Real relative_humidity_percent,
                                          Real pressure_kpa = kStandardPressureKpa) noexcept;

// Pure-tone atmospheric absorption, dB/km. ISO 9613-1:1993 equation (5).
//
// Verified against Table 1 of the standard itself -- see docs/validation.md.
// The standard's stated accuracy is +/- 10% over 0.05% to 5% molar water
// vapour, -20 to +50 C, which is most weather.
//
// A NOTE ON FREQUENCY. ISO 9613-1 note 5 records that Table 1 was computed at
// the EXACT midband frequencies 1000 * 10^(k/10), not the preferred nominal
// values printed in its own row headings. Comparing against the table with
// f = 50 rather than f = 50.119 introduces a visible error, and it is not the
// implementation's.
[[nodiscard]] Real absorption_db_per_km(Real frequency_hz, Real temperature_c,
                                        Real relative_humidity_percent,
                                        Real pressure_kpa = kStandardPressureKpa) noexcept;

// Exact one-third-octave midband frequency, 1000 * 10^(k/10). k = -13 gives the
// 50 Hz band, k = 0 gives 1 kHz, k = +10 gives 10 kHz.
[[nodiscard]] Real midband_frequency_hz(int k) noexcept;

// Characteristic acoustic impedance, rayl (Pa*s/m). About 415 at 20 C, against
// about 1.5e6 for seawater -- a factor of 3600, which is why an air bench says
// nothing quantitative about target strength even though it exercises every
// line of the signal processing.
[[nodiscard]] Real impedance_rayl(Real temperature_c,
                                  Real relative_humidity_percent = static_cast<Real>(0),
                                  Real pressure_kpa = kStandardPressureKpa) noexcept;

// Air density, kg/m^3.
[[nodiscard]] Real density_kgm3(Real temperature_c, Real relative_humidity_percent,
                                Real pressure_kpa = kStandardPressureKpa) noexcept;

// The Doppler time-scaling factor v/c for a closing speed, in air.
//
// Worth its own function because the number is startling next to water: 1 m/s
// in air is 2.9e-3, where in water it is 6.7e-4. Every chip-slip and
// Doppler-tolerance limit in the comm module tightens by 4.4x.
[[nodiscard]] Real doppler_scale(Real closing_speed_mps, Real temperature_c) noexcept;

}  // namespace phantom::air

#endif  // PHANTOM_AIR_HPP
