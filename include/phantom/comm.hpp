// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — spread-spectrum acoustic communication.
//
// The third engine of the original specification, and the last real gap. It is
// the same physics the rest of the library already models, used the other way
// round: where PingAnalyzer pulls a known waveform out of noise to DETECT
// something, this pulls a known waveform out of noise to CARRY something.
//
// Four pieces, in the order a frame meets them:
//
//   PN SEQUENCES   maximal-length LFSR codes. Their two useful properties --
//                  balance, and an autocorrelation of -1/N everywhere off the
//                  peak -- are what make a signal buried under noise
//                  recoverable, and both are verified rather than assumed.
//
//   DSSS           each data bit becomes N chips. What that buys is narrower
//                  than it is usually said to be -- see processing_gain_db.
//
//   SYNCHRONISATION  the hard part underwater, and the reason this module does
//                  not simply reuse the matched filter. See below.
//
//   CODING         CRC-32 to know a frame is wrong, RS(15,11) over GF(16) to
//                  fix it when it is only slightly wrong.
//
// WHY DOPPLER DOMINATES HERE, AND WHY THE PREAMBLE IS AN HFM.
//
// In air, a moving transmitter shifts the carrier and you correct a frequency.
// In water it is not a shift, it is a TIME SCALING of the whole waveform: with
// c = 1500 m/s, 1 m/s of closing speed scales time by 6.7e-4. Over a long PN
// code that accumulates -- a 1023-chip code at 6.7e-4 slips two thirds of a
// chip end to end, and the despreader's coherent sum falls apart.
//
// That is why the preamble is a HYPERBOLIC FM sweep and not a linear one. An
// HFM's instantaneous frequency is 1/(a + bt); scale time and it is still an
// HFM, just with different constants. So a time-scaled HFM still correlates
// against the unscaled replica with a strong peak -- the peak moves in time,
// which costs a range bias you can correct, rather than collapsing, which costs
// the frame. An LFM has no such invariance. The test suite measures both.
#ifndef PHANTOM_COMM_HPP
#define PHANTOM_COMM_HPP

#include "phantom/types.hpp"
#include "phantom/waveform.hpp"

#include <span>

namespace phantom::comm {

// ---------------------------------------------------------------------------
// Maximal-length PN sequences
// ---------------------------------------------------------------------------

// Highest LFSR degree with a tap set in the table. Degree n gives a code of
// 2^n - 1 chips: 31 at degree 5, 32767 at degree 15.
inline constexpr std::uint32_t kMaxLfsrDegree = 15;
inline constexpr std::uint32_t kMinLfsrDegree = 5;

// Length of the code a given degree produces, or 0 if the degree is unsupported.
[[nodiscard]] std::size_t msequence_length(std::uint32_t degree) noexcept;

// Fills `out` with 2^degree - 1 chips of +1/-1.
//
// The tap sets are primitive polynomials, which is what makes the sequence
// maximal length. Rather than trusting that table, the test suite measures the
// period, the balance and the off-peak autocorrelation directly -- a wrong tap
// set produces a short cycle, and a short cycle is exactly what would quietly
// destroy the processing gain.
//
// Returns the number of chips written, 0 on bad degree, zero seed or a span
// too small.
std::size_t generate_msequence(std::uint32_t degree, std::uint32_t seed,
                               std::span<Real> out) noexcept;

// ---------------------------------------------------------------------------
// Direct-sequence spread spectrum
// ---------------------------------------------------------------------------

struct DsssConfig {
    Real sample_rate_hz = 96000;
    Real carrier_hz     = 24000;
    Real chip_rate_hz   = 4000;
    // Chips per data bit. With chip_rate_hz held fixed, raising this divides the
    // data rate and lengthens the integration -- it does NOT spread the signal
    // further, and buys nothing against white noise. See processing_gain_db.
    std::size_t chips_per_bit = 63;
};

// Processing gain in dB: 10 log10(N).
//
// READ THIS BEFORE USING THE NUMBER, because the usual summary of it is wrong
// and this library's first version of the test repeated the mistake.
//
// Against WHITE NOISE at a fixed energy per bit, spreading buys NOTHING.
// Despreading multiplies by a +/-1 sequence, which leaves white noise unchanged
// in distribution; matched filtering already extracts everything a known
// waveform in white noise has to give. Measured: flat to 0.87 dB over a 16x
// range of code length. A sweep that appears to show 10 log10(N) against noise
// is almost always holding the chip amplitude fixed, which multiplies the
// transmitted ENERGY per bit by N -- that is "more energy helps", not gain.
//
// What 10 log10(N) does describe is BANDWIDTH EXPANSION at a fixed data rate,
// and it pays off against:
//
//   NARROWBAND INTERFERENCE -- despreading spreads an interferer that was never
//   spread, diluting it across the code bandwidth. Measured with the data rate
//   held at 100 bps and the chip rate grown with N: +7.1 dB from 31 to 511
//   chips against a tone, while white noise stayed flat.
//
//   COVERTNESS -- the same total power over N times the bandwidth is 10 log10(N)
//   less power spectral density, so the signal sits further under the noise
//   floor of any narrowband listener. This is the "covert" in covert
//   communication, and it is a transmit-side property this library does not
//   measure.
//
// NOTE that raising chips_per_bit while holding chip_rate_hz fixed does NOT
// expand bandwidth -- it trades data rate for integration time. To spread, raise
// the chip rate with the code length.
[[nodiscard]] Real processing_gain_db(std::size_t chips_per_bit) noexcept;

// Data rate in bits per second, chip_rate / chips_per_bit.
[[nodiscard]] Real data_rate_bps(const DsssConfig& cfg) noexcept;

// Samples one modulated bit occupies.
[[nodiscard]] std::size_t samples_per_bit(const DsssConfig& cfg) noexcept;

// BPSK-modulates `bits` (one byte per bit, 0 or 1) onto the spreading code,
// upconverted to the carrier. Returns samples written, 0 if `out` is too small
// or the configuration is degenerate.
std::size_t dsss_modulate(const DsssConfig& cfg,
                          std::span<const std::uint8_t> bits,
                          std::span<const Real> chips,
                          std::span<Real> out) noexcept;

// Coherently despreads. `out_soft` receives the correlator output per bit --
// its sign is the decision and its magnitude is the confidence, which is what a
// soft-decision decoder would want. `out_bits` may be empty if only the soft
// values are wanted.
//
// `phase_rad` is the carrier phase offset. This demodulator does NOT recover
// the phase: it is given it. Carrier recovery is a real subsystem and pretending
// a few lines of arctangent substitute for one would be dishonest -- see the
// limitations in docs/validation.md.
std::size_t dsss_demodulate(const DsssConfig& cfg,
                            std::span<const Real> samples,
                            std::span<const Real> chips,
                            Real phase_rad,
                            std::span<Real> out_soft,
                            std::span<std::uint8_t> out_bits) noexcept;

// Chip slip across one bit for a given closing speed, in chips.
//
// The number that decides whether a code length is usable at a given speed:
// once this approaches 1/2 the coherent sum has cancelled itself and no amount
// of processing gain helps.
[[nodiscard]] Real chip_slip(const DsssConfig& cfg, Real closing_speed_mps,
                             Real sound_speed_mps = static_cast<Real>(1500)) noexcept;

// Longest code whose chip slip stays under `max_slip_chips` at that speed.
[[nodiscard]] std::size_t max_chips_per_bit(const DsssConfig& cfg,
                                            Real closing_speed_mps,
                                            Real max_slip_chips = static_cast<Real>(0.25),
                                            Real sound_speed_mps = static_cast<Real>(1500)) noexcept;

// ---------------------------------------------------------------------------
// CRC-32 (IEEE 802.3, reflected, init/final 0xFFFFFFFF)
// ---------------------------------------------------------------------------

[[nodiscard]] std::uint32_t crc32(std::span<const std::uint8_t> data) noexcept;

// ---------------------------------------------------------------------------
// Reed-Solomon (15, 11) over GF(16)
// ---------------------------------------------------------------------------
//
// 11 data symbols in, 15 out, each symbol 4 bits. Corrects up to 2 symbol
// errors, or detects up to 4. A symbol is 4 bits, so one burst that ruins a
// whole nibble costs the same as one bit flipped -- which is the right shape
// for an acoustic channel, where errors come in bursts from multipath rather
// than singly from thermal noise.
//
// GF(16) is built on x^4 + x + 1. Tables are static and no larger than 16
// entries, so the whole codec allocates nothing.

inline constexpr std::size_t kRsN = 15;   // codeword symbols
inline constexpr std::size_t kRsK = 11;   // data symbols
inline constexpr std::size_t kRsT = 2;    // correctable symbol errors

// Systematic encoding: out[0..10] are the data, out[11..14] the parity.
// `data` must be 11 symbols, each < 16. Returns false otherwise.
bool rs_encode(std::span<const std::uint8_t> data, std::span<std::uint8_t> out) noexcept;

enum class RsResult : std::uint8_t {
    Clean = 0,       // no errors detected
    Corrected,       // errors found and repaired
    Uncorrectable,   // errors detected, too many to repair
    BadInput,
};

// Decodes in place. `errors_corrected` receives the count when the result is
// Corrected.
//
// Uncorrectable is a real and important answer. A decoder that always returns
// something produces plausible wrong data, which is worse than no data -- the
// CRC exists to catch what slips past even this.
RsResult rs_decode(std::span<std::uint8_t> codeword,
                   std::size_t& errors_corrected) noexcept;

}  // namespace phantom::comm

#endif  // PHANTOM_COMM_HPP
