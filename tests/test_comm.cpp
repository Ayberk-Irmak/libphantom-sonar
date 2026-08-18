// SPDX-License-Identifier: Apache-2.0
// Spread-spectrum acoustic communication.
//
// The claims here are all quantitative and all checkable against closed forms:
// a maximal-length code has period 2^n-1 and off-peak autocorrelation -1/N,
// despreading buys 10 log10(N) dB, and an HFM preamble survives a Doppler that
// destroys an LFM one. Nothing below is verified against a recorded baseline.
#include "framework.hpp"

#include "phantom/comm.hpp"
#include "phantom/matched_filter.hpp"
#include "phantom/waveform.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace phantom;
using namespace phantom::comm;

namespace {
std::array<Real, 40000> g_chips{};
std::array<Real, 650000> g_signal{};
std::array<Real, 4096> g_soft{};
std::array<std::uint8_t, 4096> g_bits{};
std::array<std::uint8_t, 4096> g_rx{};
}  // namespace

PT_TEST(m_sequences_really_are_maximal_length) {
    // The tap table is a table, and a table can be mistyped. A non-primitive
    // polynomial gives a SHORT cycle, and a short cycle would quietly destroy
    // the processing gain while everything still appeared to run -- so the
    // period is measured rather than trusted.
    std::printf("       %7s %9s %11s %14s\n", "degree", "length", "balance", "worst off-peak");
    for (std::uint32_t d = kMinLfsrDegree; d <= kMaxLfsrDegree; ++d) {
        const std::size_t n = msequence_length(d);
        PT_CHECK(n == (std::size_t{1} << d) - 1);
        PT_CHECK(generate_msequence(d, 1u, g_chips) == n);

        // Balance: exactly one more +1 than -1 over a full period. This is a
        // property only a MAXIMAL sequence has, and it is O(n).
        long sum = 0;
        for (std::size_t i = 0; i < n; ++i) sum += (g_chips[i] > 0) ? 1 : -1;
        PT_CHECK(sum == 1 || sum == -1);

        // Off-peak periodic autocorrelation is exactly -1 at every shift --
        // the property the processing gain rests on, and definitive proof of
        // maximal length. Exhaustive up to 2047 chips; beyond that it is
        // O(n^2) and a spread of shifts is enough to catch a broken table.
        const std::size_t stride = (n <= 2047) ? 1 : (n / 512);
        double worst = 0;
        for (std::size_t shift = stride; shift < n; shift += stride) {
            double acc = 0;
            for (std::size_t i = 0; i < n; ++i) {
                acc += static_cast<double>(g_chips[i]) * static_cast<double>(g_chips[(i + shift) % n]);
            }
            worst = std::max(worst, std::fabs(acc + 1.0));
        }
        std::printf("       %7u %9zu %11ld %14.1e\n", d, n, sum, worst);
        PT_CHECK(worst < 1e-9);
    }
    // A zero seed is the LFSR's fixed point and must be refused, not run.
    PT_CHECK(generate_msequence(10, 0u, g_chips) == 0);
    PT_CHECK(generate_msequence(4, 1u, g_chips) == 0);
    PT_CHECK(generate_msequence(16, 1u, g_chips) == 0);
}

namespace {

// Output SNR of the correlator: its mean over its scatter, rectified by the
// known transmitted symbol. A direct measurement, far sharper than estimating a
// bit error rate from a few hundred bits.
double measure_output_snr_db(const DsssConfig& cfg, std::size_t chips,
                             double amplitude, double noise_sigma,
                             double tone_hz, std::uint64_t seed) {
    PT_CHECK(generate_msequence(9, 0x1FFu, g_chips) > 0);
    constexpr std::size_t kNBits = 100;
    pt::Rng rng(seed);
    for (std::size_t i = 0; i < kNBits; ++i) {
        g_bits[i] = static_cast<std::uint8_t>(rng.uniform01() > 0.5 ? 1 : 0);
    }
    const std::size_t n = dsss_modulate(cfg, std::span<const std::uint8_t>(g_bits.data(), kNBits),
                                        std::span<const Real>(g_chips.data(), chips), g_signal);
    PT_CHECK(n > 0);
    for (std::size_t i = 0; i < n; ++i) {
        double interference = noise_sigma * rng.normal();
        if (tone_hz > 0) {
            interference += noise_sigma * std::cos(
                2.0 * 3.14159265358979323846 * tone_hz * static_cast<double>(i)
                / static_cast<double>(cfg.sample_rate_hz) + 0.7);
        }
        g_signal[i] = static_cast<Real>(static_cast<double>(g_signal[i]) * amplitude + interference);
    }
    PT_CHECK(dsss_demodulate(cfg, std::span<const Real>(g_signal.data(), n),
                             std::span<const Real>(g_chips.data(), chips), 0,
                             g_soft, g_rx) == kNBits);
    double sum = 0;
    for (std::size_t i = 0; i < kNBits; ++i) {
        sum += static_cast<double>(g_soft[i]) * (g_bits[i] ? -1.0 : 1.0);
    }
    const double mean = sum / kNBits;
    double var = 0;
    for (std::size_t i = 0; i < kNBits; ++i) {
        const double v = static_cast<double>(g_soft[i]) * (g_bits[i] ? -1.0 : 1.0) - mean;
        var += v * v;
    }
    return 20.0 * std::log10(std::fabs(mean) / std::sqrt(var / (kNBits - 1)));
}

}  // namespace

PT_TEST(spreading_buys_nothing_against_white_noise_at_fixed_energy_per_bit) {
    // THE CLAIM THIS RELEASE GOT WRONG FIRST TIME, so it is stated plainly.
    //
    // "Spreading a bit over N chips buys 10 log10(N) dB" is true only of a
    // specific comparison, and the obvious way to measure it is not that
    // comparison. Despreading multiplies by a +/-1 sequence; white noise is
    // unchanged in distribution by that, so against a FIXED energy per bit the
    // output SNR cannot improve, whatever N is. Matched filtering already
    // extracts everything a known waveform in white noise has to give.
    //
    // The first version of this test held the chip AMPLITUDE fixed and swept N.
    // Energy per bit is amplitude^2 * N * T_chip, so that sweep quietly
    // multiplied the transmitted energy by N and then reported the resulting
    // 10 log10(N) as processing gain. It was measuring "more energy helps".
    std::printf("       chip rate fixed at 6 kHz; amplitude scaled as 1/sqrt(N) so the\n");
    std::printf("       ENERGY PER BIT is the same for every code length:\n");
    std::printf("       %8s %16s\n", "chips", "output SNR");
    double first = 0, worst_delta = 0;
    for (std::size_t chips : {std::size_t(31), std::size_t(63), std::size_t(127),
                              std::size_t(255), std::size_t(511)}) {
        DsssConfig cfg;
        cfg.sample_rate_hz = 48000;
        cfg.carrier_hz     = 12000;
        cfg.chip_rate_hz   = 6000;
        cfg.chips_per_bit  = chips;
        const double amp = std::sqrt(31.0 / static_cast<double>(chips));
        const double snr = measure_output_snr_db(cfg, chips, amp, 2.5, 0, 4242);
        std::printf("       %8zu %13.2f dB\n", chips, snr);
        if (first == 0) first = snr;
        worst_delta = std::max(worst_delta, std::fabs(snr - first));
    }
    std::printf("       Flat to within %.2f dB over a 16x range of code length.\n", worst_delta);
    std::printf("       Against white noise, spreading is not free SNR and never was.\n");
    PT_CHECK(worst_delta < 1.0);
}

PT_TEST(processing_gain_is_real_against_narrowband_interference) {
    // Where the 10 log10(N) does live: bandwidth expansion at a FIXED data rate.
    // Despreading spreads an interferer that was never spread, so a narrowband
    // one is diluted across the code bandwidth while the signal collapses back
    // to a point.
    //
    // The comparison must therefore hold the DATA rate fixed and let the chip
    // rate -- and so the occupied bandwidth -- grow with N. Holding the chip
    // rate fixed and raising N, as the config's defaults invite, does not spread
    // anything: it trades data rate for integration time.
    std::printf("       data rate held at 100 bps, chip rate = 100*N, so the occupied\n");
    std::printf("       BANDWIDTH grows with N. Energy per bit constant throughout.\n");
    std::printf("       %8s %12s %14s %16s\n", "chips", "chip rate", "white noise", "narrowband tone");

    double tone_first = 0, tone_last = 0, noise_spread = 0, noise_first = 0;
    for (std::size_t chips : {std::size_t(31), std::size_t(63), std::size_t(127),
                              std::size_t(255), std::size_t(511)}) {
        DsssConfig cfg;
        cfg.sample_rate_hz = 192000;
        cfg.carrier_hz     = 40000;
        cfg.chip_rate_hz   = static_cast<Real>(100.0 * static_cast<double>(chips));
        cfg.chips_per_bit  = chips;
        // Samples per bit is fs/chip_rate * N = 1920 for every N, so the
        // integration length -- and with a fixed amplitude the energy per bit --
        // is identical across the sweep. Only the bandwidth changes.
        PT_CHECK(samples_per_bit(cfg) == 1920);

        const double white = measure_output_snr_db(cfg, chips, 1.0, 2.5, 0, 4242);
        // Average over interferer frequencies: a single tone's rejection depends
        // on where it happens to land in one particular code's spectrum, which
        // is a lottery rather than a measurement.
        double acc = 0;
        constexpr int kReps = 24;
        for (int r = 0; r < kReps; ++r) {
            // Sweep only inside the NARROWEST configuration's band (31 chips
            // at 100 bps occupies 40 kHz +/- 1.55 kHz). A tone outside it is
            // rejected equally by every configuration and only dilutes the
            // comparison -- which is what an earlier 38-42 kHz sweep did.
            const double f = 39000.0 + 2000.0 * r / kReps;
            acc += measure_output_snr_db(cfg, chips, 1.0, 2.5, f, 4242);
        }
        const double tone = acc / kReps;
        std::printf("       %8zu %9.0f Hz %11.2f dB %13.2f dB\n",
                    chips, static_cast<double>(cfg.chip_rate_hz), white, tone);
        if (tone_first == 0) { tone_first = tone; noise_first = white; }
        tone_last = tone;
        noise_spread = std::max(noise_spread, std::fabs(white - noise_first));
    }
    std::printf("       White noise: flat to %.2f dB -- unchanged, as it must be.\n", noise_spread);
    std::printf("       Narrowband tone: %+.1f dB from 31 to 511 chips (theory %.1f dB).\n",
                tone_last - tone_first, 10.0 * std::log10(511.0 / 31.0));
    std::printf("       Short of the ideal because this despreader is a bare correlator:\n");
    std::printf("       10 log10(N) assumes a filter matched to the DATA bandwidth after\n");
    std::printf("       despreading, which collects the diluted interferer from a narrow\n");
    std::printf("       slice instead of the whole band. The trend is the real result.\n");
    PT_CHECK(noise_spread < 1.5);
    PT_CHECK(tone_last - tone_first > 7.0);

    // The formula itself is still the formula.
    PT_CHECK_NEAR(static_cast<double>(processing_gain_db(1024)), 30.103, 1e-3);
}

PT_TEST(the_link_budget_is_a_trade_not_a_free_gain) {
    // Quoting processing gain without the throughput it costs is the standard
    // way to oversell a spread-spectrum link.
    DsssConfig cfg;
    cfg.chip_rate_hz = 4000;
    std::printf("       chip rate 4000 Hz:\n");
    std::printf("       %8s %12s %14s\n", "chips", "gain (dB)", "data rate");
    for (std::size_t chips : {std::size_t(31), std::size_t(127), std::size_t(511),
                              std::size_t(2047)}) {
        cfg.chips_per_bit = chips;
        std::printf("       %8zu %12.2f %11.1f bps\n", chips,
                    static_cast<double>(processing_gain_db(chips)),
                    static_cast<double>(data_rate_bps(cfg)));
    }
    cfg.chips_per_bit = 511;
    PT_CHECK_NEAR(static_cast<double>(data_rate_bps(cfg)), 4000.0 / 511.0, 1e-6);
    // 30 dB of gain at 4 kchip/s is under 4 bits per second. That is the honest
    // number, and it is why covert acoustic links carry status rather than data.
    cfg.chips_per_bit = 1024;
    PT_CHECK(data_rate_bps(cfg) < 4);
}

PT_TEST(doppler_is_a_time_scaling_and_it_breaks_long_codes) {
    // Underwater the mobility problem is not a frequency offset, it is a time
    // scaling: v/c with c = 1500 is 6.7e-4 per m/s. Over a long code that
    // accumulates into chip slip, and the coherent sum cancels itself.
    DsssConfig cfg;
    std::printf("       chip slip across one bit, at 1500 m/s sound speed:\n");
    std::printf("       %8s %12s %12s %12s\n", "chips", "1 m/s", "3 m/s", "10 m/s");
    for (std::size_t chips : {std::size_t(31), std::size_t(127), std::size_t(511),
                              std::size_t(2047)}) {
        cfg.chips_per_bit = chips;
        std::printf("       %8zu %12.3f %12.3f %12.3f\n", chips,
                    static_cast<double>(chip_slip(cfg, 1)),
                    static_cast<double>(chip_slip(cfg, 3)),
                    static_cast<double>(chip_slip(cfg, 10)));
    }
    // 2047 chips at 10 m/s slips more than 13 chips -- the code has walked off
    // itself entirely.
    cfg.chips_per_bit = 2047;
    PT_CHECK(chip_slip(cfg, 10) > 13);
    PT_CHECK_NEAR(static_cast<double>(chip_slip(cfg, 1)), 2047.0 / 1500.0, 1e-4);

    // The inverse: how long a code a given speed permits.
    std::printf("       longest code holding slip under a quarter chip:\n");
    for (Real v : {static_cast<Real>(1), static_cast<Real>(3), static_cast<Real>(10)}) {
        std::printf("       %5.0f m/s -> %6zu chips (%.1f dB of gain)\n",
                    static_cast<double>(v), max_chips_per_bit(cfg, v),
                    static_cast<double>(processing_gain_db(max_chips_per_bit(cfg, v))));
    }
    PT_CHECK(max_chips_per_bit(cfg, 1) > max_chips_per_bit(cfg, 10));
    std::printf("       So speed sets a CEILING on processing gain that no amount of\n");
    std::printf("       transmit power can lift. Resynchronise more often, or go slower.\n");
}

PT_TEST(an_hfm_preamble_survives_the_doppler_an_lfm_does_not) {
    // Why the preamble is hyperbolic. An HFM's instantaneous frequency is
    // 1/(a+bt); scale time and it is STILL an HFM, so a scaled copy still
    // correlates strongly with the unscaled replica -- the peak moves rather
    // than collapsing. An LFM has no such invariance.
    constexpr std::size_t kFftSize = 8192;
    const Real fs = 96000;
    auto peak_ratio = [&](PulseType type, Real speed_mps) {
        PulseSpec spec;
        spec.type = type;
        spec.f_start_hz = 18000;
        spec.f_end_hz = 30000;
        spec.duration_s = static_cast<Real>(0.02);
        spec.amplitude = 1;
        spec.taper = Taper::Tukey25;

        static std::array<Complex, kFftSize> replica{};
        const std::size_t rn = render_analytic(spec, fs, replica);
        PT_CHECK(rn > 0);

        static std::array<Real, kFftSize> received{};
        for (Real& v : received) v = 0;
        // Place the pulse INSIDE the buffer, not at index 0. Doppler shifts the
        // matched-filter peak in time, and for an upsweep the shift is NEGATIVE
        // -- v0.2 derived it as dt = -delta * f_end / mu, which here is about
        // -6 samples at 2 m/s. A correlation whose lags start at zero cannot see
        // that peak at all, and reports the pulse as destroyed when it has only
        // moved. This was the second thing wrong with this test.
        constexpr std::size_t kOffset = 512;
        // render_real_doppler takes v/c, NOT 1 + v/c. Passing the latter made
        // the "zero Doppler" baseline a 2x time-compressed pulse -- a hopeless
        // match -- so every other speed measured better than it and the
        // correlation ratios came out above 100%, which a matched filter cannot
        // produce. The absurd number is what exposed the argument convention.
        const std::size_t sn = render_real_doppler(
            spec, fs, speed_mps / static_cast<Real>(1500),
            std::span<Real>(received.data() + kOffset, kFftSize - kOffset));
        PT_CHECK(sn > 0);

        // Pad the received block. A closing target COMPRESSES the waveform, so
        // at positive Doppler it is shorter than the replica -- and
        // correlate_direct refuses a signal shorter than what it correlates
        // against. The padding also gives the peak somewhere to move to, which
        // is the point: an HFM's peak SHIFTS under Doppler rather than
        // collapsing, and a correlation with no spare lags cannot show that.
        PT_CHECK(sn + 2 * kOffset <= kFftSize);
        static std::array<Complex, 2 * kFftSize> corr{};
        const std::size_t cn = correlate_direct(
            std::span<const Real>(received.data(), sn + 2 * kOffset),
            std::span<const Complex>(replica.data(), rn), corr);
        PT_CHECK(cn > 0);
        double peak = 0;
        for (std::size_t i = 0; i < cn; ++i) {
            peak = std::max(peak, static_cast<double>(std::abs(corr[i])));
        }

        // NORMALISE. A raw peak is the wrong figure of merit: Doppler rescales
        // the waveform in time, so the received block has a different ENERGY
        // from the replica and a bare peak moves with speed for that reason
        // alone. The correlation COEFFICIENT removes it -- peak / (||r|| ||s||)
        // measures waveform mismatch and nothing else, and Cauchy-Schwarz
        // bounds it by 1, which makes "above 100%" a usable alarm.
        double er = 0, ep = 0;
        for (std::size_t i = 0; i < sn; ++i) {
            const double v = static_cast<double>(received[kOffset + i]);
            er += v * v;
        }
        for (std::size_t i = 0; i < rn; ++i) {
            ep += static_cast<double>(std::norm(replica[i]));
        }
        return peak / std::sqrt(er * ep);
    };

    std::printf("       12 kHz sweep, 20 ms. Correlation coefficient, as %% of its\n");
    std::printf("       own zero-Doppler value:\n");
    std::printf("       %10s %14s %14s\n", "speed", "LFM", "HFM");
    const double lfm0 = peak_ratio(PulseType::LfmUp, 0);
    const double hfm0 = peak_ratio(PulseType::Hfm, 0);
    double lfm_worst = 1, hfm_worst = 1;
    for (Real v : {static_cast<Real>(2), static_cast<Real>(5),
                   static_cast<Real>(10), static_cast<Real>(20)}) {
        const double l = peak_ratio(PulseType::LfmUp, v) / lfm0;
        const double h = peak_ratio(PulseType::Hfm, v) / hfm0;
        std::printf("       %8.0f m/s %12.1f %% %12.1f %%\n",
                    static_cast<double>(v), 100.0 * l, 100.0 * h);
        lfm_worst = std::min(lfm_worst, l);
        hfm_worst = std::min(hfm_worst, h);
    }
    std::printf("       At 20 m/s the HFM still holds %.0f%% and the LFM has lost %.0f%%.\n",
                100.0 * hfm_worst, 100.0 - 100.0 * lfm_worst);
    std::printf("       The HFM's loss does not GROW with speed, which is the property:\n");
    std::printf("       a time-scaled hyperbolic sweep is still a hyperbolic sweep, so it\n");
    std::printf("       still matches -- the peak moves rather than collapsing. That is\n");
    std::printf("       why the preamble is an HFM and the payload can then be despread\n");
    std::printf("       once the scale it reveals has been corrected.\n");
    PT_CHECK(hfm_worst > 0.95);
    PT_CHECK(lfm_worst < 0.70);
    PT_CHECK(hfm_worst > lfm_worst);
}

PT_TEST(crc32_matches_the_published_check_values) {
    // The standard check values for CRC-32/ISO-HDLC, so this is verified against
    // something published rather than against itself.
    const char* s = "123456789";
    std::array<std::uint8_t, 9> buf{};
    for (std::size_t i = 0; i < 9; ++i) buf[i] = static_cast<std::uint8_t>(s[i]);
    const std::uint32_t got = crc32(buf);
    std::printf("       CRC-32(\"123456789\") = 0x%08X (published 0xCBF43926)\n", got);
    PT_CHECK(got == 0xCBF43926u);
    PT_CHECK(crc32(std::span<const std::uint8_t>()) == 0u);

    // A single flipped bit anywhere must change it.
    std::array<std::uint8_t, 64> data{};
    pt::Rng rng(99);
    for (std::uint8_t& b : data) b = static_cast<std::uint8_t>(rng.uniform01() * 256);
    const std::uint32_t base = crc32(data);
    std::size_t missed = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
        for (int bit = 0; bit < 8; ++bit) {
            data[i] = static_cast<std::uint8_t>(data[i] ^ (1u << bit));
            if (crc32(data) == base) ++missed;
            data[i] = static_cast<std::uint8_t>(data[i] ^ (1u << bit));
        }
    }
    std::printf("       512 single-bit flips, %zu undetected\n", missed);
    PT_CHECK(missed == 0);
}

PT_TEST(reed_solomon_corrects_two_symbol_errors_and_admits_when_it_cannot) {
    pt::Rng rng(2025);
    std::size_t clean = 0, corrected = 0, uncorrectable = 0, wrong = 0;

    for (std::size_t trial = 0; trial < 4000; ++trial) {
        std::array<std::uint8_t, kRsK> data{};
        for (std::uint8_t& d : data) d = static_cast<std::uint8_t>(rng.uniform01() * 16) & 0x0F;
        std::array<std::uint8_t, kRsN> code{};
        PT_CHECK(rs_encode(data, code));

        // Systematic: the data must appear untouched at the front.
        for (std::size_t i = 0; i < kRsK; ++i) PT_CHECK(code[i] == data[i]);

        const std::size_t n_err = trial % 3;   // 0, 1 or 2 errors
        std::array<std::uint8_t, kRsN> rx = code;
        std::size_t hit[2] = {kRsN, kRsN};
        for (std::size_t e = 0; e < n_err; ++e) {
            std::size_t p;
            do { p = static_cast<std::size_t>(rng.uniform01() * kRsN) % kRsN; }
            while (p == hit[0]);
            hit[e] = p;
            std::uint8_t delta;
            do { delta = static_cast<std::uint8_t>(rng.uniform01() * 16) & 0x0F; } while (delta == 0);
            rx[p] = static_cast<std::uint8_t>(rx[p] ^ delta);
        }

        std::size_t fixed = 0;
        const RsResult r = rs_decode(rx, fixed);
        if (r == RsResult::Clean) ++clean;
        else if (r == RsResult::Corrected) ++corrected;
        else ++uncorrectable;

        // Up to two errors must ALWAYS be repaired, exactly.
        if (n_err <= kRsT) {
            PT_CHECK(r == RsResult::Clean || r == RsResult::Corrected);
            for (std::size_t i = 0; i < kRsN; ++i) {
                if (rx[i] != code[i]) ++wrong;
            }
        }
    }
    std::printf("       4000 codewords with 0-2 symbol errors: %zu clean, %zu corrected,\n",
                clean, corrected);
    std::printf("       %zu uncorrectable, %zu symbols left wrong\n", uncorrectable, wrong);
    PT_CHECK(wrong == 0);
    PT_CHECK(uncorrectable == 0);

    // Beyond its power the decoder must say so rather than invent a codeword --
    // but it cannot always, and the reason is the code itself rather than the
    // implementation. RS(15,11) has minimum distance 5. A word with 4 errors
    // sits at distance 4 from the true codeword and can therefore sit at
    // distance 1 or 2 from a DIFFERENT one, which the decoder will then
    // "correct" to, confidently. That is miscorrection, it is a property of the
    // distance, and it is precisely what the CRC is for.
    //
    // So this measures both: how often RS admits defeat, and whether a CRC
    // carried in the frame catches what slips past it.
    std::size_t admitted = 0, miscorrected = 0, caught_by_crc = 0, escaped = 0;
    for (std::size_t trial = 0; trial < 3000; ++trial) {
        std::array<std::uint8_t, kRsK> data{};
        for (std::uint8_t& d : data) d = static_cast<std::uint8_t>(rng.uniform01() * 16) & 0x0F;
        const std::uint32_t frame_crc = crc32(data);

        std::array<std::uint8_t, kRsN> code{};
        rs_encode(data, code);
        std::array<std::uint8_t, kRsN> rx = code;
        for (std::size_t e = 0; e < 4; ++e) {
            const std::size_t p = static_cast<std::size_t>(rng.uniform01() * kRsN) % kRsN;
            std::uint8_t delta;
            do { delta = static_cast<std::uint8_t>(rng.uniform01() * 16) & 0x0F; } while (delta == 0);
            rx[p] = static_cast<std::uint8_t>(rx[p] ^ delta);
        }

        std::size_t fixed = 0;
        const RsResult r = rs_decode(rx, fixed);
        if (r == RsResult::Uncorrectable) { ++admitted; continue; }

        std::array<std::uint8_t, kRsK> out{};
        for (std::size_t i = 0; i < kRsK; ++i) out[i] = rx[i];
        bool same = true;
        for (std::size_t i = 0; i < kRsK; ++i) {
            if (out[i] != data[i]) same = false;
        }
        if (same) continue;                       // RS happened to get it right
        ++miscorrected;
        if (crc32(out) != frame_crc) ++caught_by_crc; else ++escaped;
    }
    std::printf("       3000 codewords with 4 symbol errors (beyond t=2):\n");
    std::printf("         %4zu RS admitted defeat\n", admitted);
    std::printf("         %4zu RS miscorrected to a different valid codeword\n", miscorrected);
    std::printf("              of those, %zu caught by the frame CRC, %zu escaped\n",
                caught_by_crc, escaped);
    std::printf("       Miscorrection is the code's minimum distance (d=5) showing, not a\n");
    std::printf("       defect: 4 errors can land within 2 of a different codeword. The\n");
    std::printf("       CRC is the layer that exists to catch exactly this.\n");

    // The RS decoder must catch the majority on its own...
    PT_CHECK(admitted > miscorrected);
    // ...and the CRC must catch essentially all of the remainder. A 32-bit CRC
    // lets through about 2^-32 of corrupted frames, so zero escapes out of a
    // few hundred is the expected outcome, not a lucky one.
    PT_CHECK(escaped == 0);
}

PT_TEST(a_frame_survives_a_burst_that_a_bit_code_would_not) {
    // Why symbols and not bits. A 4-bit symbol ruined entirely costs the same as
    // one bit flipped, which is the right shape for an acoustic channel where
    // multipath destroys short runs rather than isolated bits.
    std::array<std::uint8_t, kRsK> data{};
    for (std::size_t i = 0; i < kRsK; ++i) data[i] = static_cast<std::uint8_t>(i);
    std::array<std::uint8_t, kRsN> code{};
    PT_CHECK(rs_encode(data, code));

    std::array<std::uint8_t, kRsN> rx = code;
    // An 8-bit burst: two adjacent symbols completely destroyed.
    rx[4] = static_cast<std::uint8_t>(~rx[4] & 0x0F);
    rx[5] = static_cast<std::uint8_t>(~rx[5] & 0x0F);
    std::size_t fixed = 0;
    const RsResult r = rs_decode(rx, fixed);
    std::printf("       8 consecutive bits destroyed -> %s, %zu symbols repaired\n",
                (r == RsResult::Corrected) ? "corrected" : "FAILED", fixed);
    PT_CHECK(r == RsResult::Corrected);
    PT_CHECK(fixed == 2);
    for (std::size_t i = 0; i < kRsN; ++i) PT_CHECK(rx[i] == code[i]);
}
