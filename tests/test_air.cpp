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
#include "phantom/bench.hpp"
#include "phantom/fft.hpp"
#include "phantom/comm.hpp"
#include "phantom/sound_speed.hpp"

#include "data/iso9613_table.hpp"

#include <algorithm>
#include <array>
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

// ---------------------------------------------------------------------------
// v0.16: turning a recording into a measurement
// ---------------------------------------------------------------------------

namespace {

constexpr std::size_t kBenchFft = 4096;
std::array<Complex, kBenchFft / 2> g_tw{};
std::array<Complex, 2 * kBenchFft> g_bench_work{};
std::array<Real, kBenchFft> g_active{};
std::array<Real, kBenchFft> g_silent{};
std::array<Real, kBenchFft> g_probe{};
std::array<Real, kBenchFft> g_h{};

FftView bench_fft() {
    fft_init(g_tw);
    return FftView{std::span<const Complex>(g_tw.data(), g_tw.size()), kBenchFft};
}

// A 20 ms LFM chirp, the bench probe.
void make_probe(std::span<Real> out, std::size_t n, Real fs) {
    const Real T = static_cast<Real>(n) / fs;
    for (std::size_t i = 0; i < n; ++i) {
        const Real t = static_cast<Real>(i) / fs;
        const Real ph = static_cast<Real>(2) * kPi
                      * (static_cast<Real>(2000) * t
                         + static_cast<Real>(0.5) * static_cast<Real>(6000) / T * t * t);
        out[i] = std::sin(ph);
    }
}

}  // namespace

PT_TEST(a_dead_channel_is_recognised_before_it_produces_a_number) {
    // The test that exists because v0.15's bench could not run.
    //
    // The dead channel did not look dead. It gave 13000 distinct sample values,
    // a healthy RMS, and a matched filter peak 44 dB above background -- all of
    // which a working microphone also gives, and all of which were electrical
    // noise and a startup transient. Only a CONTROLLED comparison settled it.
    const Real fs = 48000;
    const FftView fft = bench_fft();
    PT_CHECK(fft.valid());
    bench::QualifyConfig cfg;
    cfg.sample_rate_hz = fs;

    pt::Rng rng(1616);
    const std::size_t n_probe = 960;   // 20 ms
    make_probe(g_probe, n_probe, fs);

    auto fill = [&](std::span<Real> dst, double probe_gain, double noise) {
        for (std::size_t i = 0; i < dst.size(); ++i) {
            dst[i] = static_cast<Real>(noise * rng.normal());
        }
        for (std::size_t i = 0; i < n_probe && i + 500 < dst.size(); ++i) {
            dst[i + 500] += static_cast<Real>(probe_gain) * g_probe[i];
        }
    };

    // 1. A working channel: the probe arrives.
    fill(g_active, 0.30, 0.02);
    fill(g_silent, 0.0, 0.02);
    bench::QualifyResult r;
    PT_CHECK(bench::qualify_channel(fft, g_active, g_silent, cfg, g_bench_work, r));
    std::printf("       working channel : active %.1f dB, silent %.1f dB, excess %+.2f dB -> %s\n",
                static_cast<double>(r.active_db), static_cast<double>(r.silent_db),
                static_cast<double>(r.excess_db), r.usable ? "USABLE" : "unusable");
    PT_CHECK(r.usable);
    PT_CHECK(r.excess_db > 6);

    // 2. A dead channel: noise only, whatever the playback did. This is the
    //    real machine's behaviour, where the measured excess was -4.60 dB.
    fill(g_active, 0.0, 0.02);
    fill(g_silent, 0.0, 0.02);
    PT_CHECK(bench::qualify_channel(fft, g_active, g_silent, cfg, g_bench_work, r));
    std::printf("       dead channel    : excess %+.2f dB -> %s\n",
                static_cast<double>(r.excess_db), r.usable ? "USABLE" : "unusable");
    PT_CHECK(!r.usable);
    PT_CHECK(std::fabs(static_cast<double>(r.excess_db)) < 3.0);

    // 3. A clipping channel: the probe IS there, loudly, but every level it
    //    reports is fiction. Excess alone would pass this.
    fill(g_active, 3.0, 0.02);
    for (Real& v : g_active) v = std::clamp(v, static_cast<Real>(-1), static_cast<Real>(1));
    fill(g_silent, 0.0, 0.02);
    PT_CHECK(bench::qualify_channel(fft, g_active, g_silent, cfg, g_bench_work, r));
    std::printf("       clipping channel: excess %+.2f dB, %.2f%% clipped -> %s\n",
                static_cast<double>(r.excess_db),
                100.0 * static_cast<double>(r.clipped_fraction),
                r.usable ? "USABLE" : "unusable");
    PT_CHECK(r.excess_db > 6);          // it would pass on excess alone...
    PT_CHECK(!r.usable);                // ...and is rejected anyway
}

PT_TEST(two_distances_cancel_the_sound_card_latency) {
    // A single delay measurement over a desk is 97% sound-card buffering. The
    // difference of two is not.
    const double c_true = 343.37;
    const double latency = 0.032;       // 32 ms, a typical figure
    const double d1 = 0.40, d2 = 1.60;
    const double t1 = d1 / c_true + latency;
    const double t2 = d2 / c_true + latency;

    std::printf("       truth: c = %.2f m/s, fixed latency = %.1f ms\n", c_true, 1000 * latency);
    std::printf("       measured delays: %.2f ms at %.2f m, %.2f ms at %.2f m\n",
                1000 * t1, d1, 1000 * t2, d2);
    std::printf("       naive single-shot speed d1/t1 = %.1f m/s (out by %.0fx)\n",
                d1 / t1, c_true / (d1 / t1));

    const bench::TwoDistanceResult r = bench::two_distance_solve(
        static_cast<Real>(d1), static_cast<Real>(t1),
        static_cast<Real>(d2), static_cast<Real>(t2), static_cast<Real>(1e-5));
    PT_CHECK(r.valid);
    std::printf("       two-distance solve: c = %.2f m/s, latency = %.2f ms\n",
                static_cast<double>(r.sound_speed_mps), 1000 * static_cast<double>(r.fixed_latency_s));
    PT_CHECK_NEAR(static_cast<double>(r.sound_speed_mps), c_true, pt::tol(1e-6, 0.5));
    PT_CHECK_NEAR(static_cast<double>(r.fixed_latency_s), latency, pt::tol(1e-9, 1e-4));

    // Two distances too close together cannot resolve anything, and the
    // function says so rather than dividing by noise.
    const bench::TwoDistanceResult bad = bench::two_distance_solve(
        static_cast<Real>(1.0), static_cast<Real>(t1),
        static_cast<Real>(1.001), static_cast<Real>(t1 + 3e-6),
        static_cast<Real>(1e-5));
    PT_CHECK(!bad.valid);
}

PT_TEST(the_impulse_response_recovers_a_known_multipath) {
    // Deconvolution against a channel whose arrivals are known exactly, so the
    // estimator is checked rather than admired.
    const Real fs = 48000;
    const FftView fft = bench_fft();
    const std::size_t n_probe = 960;
    make_probe(g_probe, n_probe, fs);

    struct Path { std::size_t delay; double gain; };
    const Path paths[] = {{200, 1.0}, {320, 0.5}, {512, 0.25}};

    pt::Rng rng(99);
    for (Real& v : g_active) v = static_cast<Real>(0.001 * rng.normal());
    for (const Path& p : paths) {
        for (std::size_t i = 0; i < n_probe && p.delay + i < kBenchFft; ++i) {
            g_active[p.delay + i] += static_cast<Real>(p.gain) * g_probe[i];
        }
    }

    const std::size_t n = bench::estimate_impulse_response(
        fft, std::span<const Real>(g_probe.data(), n_probe), g_active,
        static_cast<Real>(1e-3), g_bench_work, g_h);
    PT_CHECK(n == kBenchFft);

    std::array<bench::Arrival, 8> arrivals{};
    // Separation of 1 ms: several times the probe's own 1/6000 s resolution,
    // and well under the 2.5 ms between the paths above.
    const std::size_t found = bench::find_arrivals(
        std::span<const Real>(g_h.data(), n), fs,
        static_cast<Real>(0.001), static_cast<Real>(-20), arrivals);

    std::printf("       three paths at %zu/%zu/%zu samples, gains 1.0/0.5/0.25:\n",
                paths[0].delay, paths[1].delay, paths[2].delay);
    std::printf("       %10s %12s %12s %12s\n", "found", "sample", "truth", "rel dB");
    PT_CHECK(found >= 3);
    for (std::size_t i = 0; i < found && i < 3; ++i) {
        const double s = static_cast<double>(arrivals[i].delay_s) * static_cast<double>(fs);
        std::printf("       %10zu %12.1f %12zu %12.2f\n", i, s, paths[i].delay,
                    static_cast<double>(arrivals[i].relative_db));
        PT_CHECK_NEAR(s, static_cast<double>(paths[i].delay), 2.0);
    }
    // Amplitudes in the right ratio: -6 dB and -12 dB.
    PT_CHECK_NEAR(static_cast<double>(arrivals[1].relative_db), -6.02, 1.5);
    PT_CHECK_NEAR(static_cast<double>(arrivals[2].relative_db), -12.04, 1.5);

    std::printf("       Regularisation matters: without it the deconvolution divides\n");
    std::printf("       by the probe's near-zero out-of-band spectrum and amplifies\n");
    std::printf("       noise without bound -- while still looking like a response.\n");
}
