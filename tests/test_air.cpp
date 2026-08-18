// SPDX-License-Identifier: Apache-2.0
// Acoustics in air.
//
// Air is here because it is the medium anyone can test in: a speaker and a
// microphone are already on the desk, where a hydrophone and a projector are
// several hundred euros and need water. The DSP does not care which medium it
// is in, so an air bench exercises the whole chain before any wet hardware is
// bought.
//
// The absorption model is verified against ISO 9613-1:1993 Table 1 -- the
// published standard, not a second implementation of the same equations.
#include "framework.hpp"

#include "phantom/air.hpp"
#include "phantom/comm.hpp"
#include "phantom/sound_speed.hpp"

#include "data/iso9613_table.hpp"

#include <cmath>
#include <cstdio>

using namespace phantom;

PT_TEST(atmospheric_absorption_matches_the_published_iso_table) {
    // 105 values from ISO 9613-1:1993 Table 1 (g), 10 C, one standard
    // atmosphere. The table is printed to three significant figures, so that is
    // all a test may demand of it -- and all it needs to demand.
    double worst = 0;
    double worst_f = 0, worst_rh = 0, worst_pub = 0, worst_got = 0;
    for (int i = 0; i < iso9613::kTableSize; ++i) {
        const iso9613::Entry& e = iso9613::kTable[i];
        // Note 5 of the standard: the coefficients were computed at the EXACT
        // midband frequencies, not the nominal ones in the row headings.
        const Real f = air::midband_frequency_hz(e.band_index);
        const double got = static_cast<double>(air::absorption_db_per_km(
            f, static_cast<Real>(e.temperature_c),
            static_cast<Real>(e.relative_humidity_pct)));
        const double rel = std::fabs(got - e.attenuation_db_per_km) / e.attenuation_db_per_km;
        if (rel > worst) {
            worst = rel;
            worst_f = e.nominal_frequency_hz;
            worst_rh = e.relative_humidity_pct;
            worst_pub = e.attenuation_db_per_km;
            worst_got = got;
        }
    }
    std::printf("       %d published ISO 9613-1 values, worst relative error %.3f%%\n",
                iso9613::kTableSize, 100.0 * worst);
    std::printf("       worst at %.0f Hz, %.0f%% RH: published %.4g, computed %.4g dB/km\n",
                worst_f, worst_rh, worst_pub, worst_got);
    PT_CHECK(worst < pt::tol(0.006, 0.02));

    // Using the NOMINAL frequency instead of the midband one is the error the
    // standard's note 5 exists to prevent. It is small but systematic, and this
    // measures it rather than leaving the note as a claim.
    double worst_nominal = 0;
    for (int i = 0; i < iso9613::kTableSize; ++i) {
        const iso9613::Entry& e = iso9613::kTable[i];
        const double got = static_cast<double>(air::absorption_db_per_km(
            static_cast<Real>(e.nominal_frequency_hz),
            static_cast<Real>(e.temperature_c),
            static_cast<Real>(e.relative_humidity_pct)));
        worst_nominal = std::max(worst_nominal,
            std::fabs(got - e.attenuation_db_per_km) / e.attenuation_db_per_km);
    }
    std::printf("       same table with NOMINAL frequencies: worst error %.3f%% (%.1fx worse)\n",
                100.0 * worst_nominal, worst_nominal / worst);
    PT_CHECK(worst_nominal > worst * 2);
}

PT_TEST(absorption_in_air_is_driven_by_humidity_in_a_way_water_never_is) {
    // The structural difference between the two media, and the reason air needs
    // its own model rather than a different constant. Seawater absorption
    // varies with temperature and salinity by tens of percent; air absorption
    // varies with HUMIDITY by factors.
    std::printf("       absorption at 10 C, dB/km:\n");
    std::printf("       %8s %10s %10s %10s %10s\n", "freq", "10%RH", "20%RH", "50%RH", "100%RH");
    double ratio_4k = 0;
    for (int k : {-3, 0, 3, 6, 9}) {
        const Real f = air::midband_frequency_hz(k);
        const double a10 = static_cast<double>(air::absorption_db_per_km(f, 10, 10));
        const double a20 = static_cast<double>(air::absorption_db_per_km(f, 10, 20));
        const double a50 = static_cast<double>(air::absorption_db_per_km(f, 10, 50));
        const double a100 = static_cast<double>(air::absorption_db_per_km(f, 10, 100));
        std::printf("       %8.0f %10.2f %10.2f %10.2f %10.2f\n",
                    static_cast<double>(f), a10, a20, a50, a100);
        if (k == 6) ratio_4k = a20 / a10;
    }
    std::printf("       At 4 kHz, 10%% -> 20%% humidity multiplies the loss by %.2f.\n", ratio_4k);
    std::printf("       Nothing in seawater behaves like that, and it is why a bench\n");
    std::printf("       measurement in air must record the humidity or mean nothing.\n");
    PT_CHECK(ratio_4k > 1.5);

    // Non-monotonic in humidity: at low frequency more moisture REDUCES the
    // loss, at high frequency it increases it. A model that only knew "damp air
    // absorbs more" would be wrong half the time.
    const Real f500 = air::midband_frequency_hz(-3);
    PT_CHECK(air::absorption_db_per_km(f500, 10, 50) < air::absorption_db_per_km(f500, 10, 10));
    const Real f8k = air::midband_frequency_hz(9);
    PT_CHECK(air::absorption_db_per_km(f8k, 10, 50) > air::absorption_db_per_km(f8k, 10, 10));
}

PT_TEST(sound_speed_in_air_hits_its_textbook_anchors) {
    // Not a primary standard: this is the square-root law about 331.45 m/s at
    // 0 C, and the test says how close it gets rather than claiming Cramer.
    const double c0 = static_cast<double>(air::sound_speed(0, 0));
    const double c20 = static_cast<double>(air::sound_speed(20, 0));
    std::printf("       dry air:  0 C -> %.2f m/s (textbook 331.45)\n", c0);
    std::printf("                20 C -> %.2f m/s (textbook 343.2)\n", c20);
    PT_CHECK_NEAR(c0, 331.45, 0.05);
    PT_CHECK_NEAR(c20, 343.2, 0.5);

    // Humidity raises it slightly -- water vapour is lighter than the air it
    // displaces. A small effect, and in the right direction.
    const double c20_humid = static_cast<double>(air::sound_speed(20, 100));
    std::printf("       20 C at 100%% RH -> %.2f m/s (+%.2f m/s over dry)\n",
                c20_humid, c20_humid - c20);
    PT_CHECK(c20_humid > c20);
    PT_CHECK(c20_humid - c20 < 2.0);

    // Impedance, the number that says an air bench cannot measure target
    // strength however well it exercises the processing.
    const double z_air = static_cast<double>(air::impedance_rayl(20, 50));
    std::printf("       impedance: air %.0f rayl, seawater ~1.5e6 -> a factor of %.0f\n",
                z_air, 1.5e6 / z_air);
    PT_CHECK(z_air > 380 && z_air < 430);
}

PT_TEST(air_is_the_harsher_doppler_environment_by_a_factor_of_four) {
    // The reason an air bench is a GOOD test of the comm module rather than a
    // weak one: with c 4.4x smaller, every Doppler time-scaling is 4.4x larger
    // for the same closing speed, so a code length that is comfortable in water
    // falls apart in air.
    const double s_air = static_cast<double>(air::doppler_scale(1, 20));
    const double s_water = 1.0 / 1500.0;
    std::printf("       1 m/s closing: v/c = %.2e in air, %.2e in water (%.2fx)\n",
                s_air, s_water, s_air / s_water);
    PT_CHECK_NEAR(s_air / s_water, 4.37, 0.1);

    // What that costs a spread-spectrum link, using the comm module's own
    // chip-slip model with the air sound speed substituted.
    comm::DsssConfig cfg;
    std::printf("       chip slip at 1 m/s:\n");
    std::printf("       %8s %12s %12s\n", "chips", "water", "air");
    for (std::size_t chips : {std::size_t(127), std::size_t(511), std::size_t(2047)}) {
        cfg.chips_per_bit = chips;
        std::printf("       %8zu %12.3f %12.3f\n", chips,
                    static_cast<double>(comm::chip_slip(cfg, 1, static_cast<Real>(1500))),
                    static_cast<double>(comm::chip_slip(cfg, 1, air::sound_speed(20))));
    }
    cfg.chips_per_bit = 511;
    const Real slip_air = comm::chip_slip(cfg, 1, air::sound_speed(20));
    std::printf("       A 511-chip code slips %.2f chips per bit at walking pace in air.\n",
                static_cast<double>(slip_air));
    std::printf("       If the Doppler machinery survives that, water is the easy case.\n");
    PT_CHECK(static_cast<double>(slip_air) > 1.0);
}

PT_TEST(relaxation_frequencies_bracket_the_absorption_curve_shape) {
    // The two relaxation frequencies are the whole shape of the curve, and
    // exposing them is what lets a caller understand a measurement rather than
    // just look it up.
    std::printf("       %6s %10s %14s %14s\n", "T (C)", "RH %", "f_rO (Hz)", "f_rN (Hz)");
    for (int rh : {10, 50, 100}) {
        const double fo = static_cast<double>(air::oxygen_relaxation_hz(10, static_cast<Real>(rh)));
        const double fn = static_cast<double>(air::nitrogen_relaxation_hz(10, static_cast<Real>(rh)));
        std::printf("       %6d %10d %14.1f %14.1f\n", 10, rh, fo, fn);
        // Oxygen relaxes far above nitrogen at every humidity, which is what
        // makes the mid-band plateau in the table.
        PT_CHECK(fo > fn);
    }
    // Both rise with humidity -- more water vapour makes relaxation faster,
    // which is why damp air is quieter at low frequency and louder at high.
    PT_CHECK(air::oxygen_relaxation_hz(10, 100) > air::oxygen_relaxation_hz(10, 10));
    PT_CHECK(air::nitrogen_relaxation_hz(10, 100) > air::nitrogen_relaxation_hz(10, 10));

    // Degenerate inputs are refused rather than producing a NaN downstream.
    PT_CHECK(air::absorption_db_per_km(0, 20, 50) == 0);
    PT_CHECK(air::absorption_db_per_km(1000, 20, 50, 0) == 0);
}
