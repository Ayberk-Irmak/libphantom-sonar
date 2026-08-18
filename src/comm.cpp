// SPDX-License-Identifier: Apache-2.0
#include "phantom/comm.hpp"

#include <cmath>

namespace phantom::comm {
namespace {

constexpr Real kZero = static_cast<Real>(0);
constexpr Real kOne  = static_cast<Real>(1);

// Primitive polynomial tap masks, index = degree, for a RIGHT-shifting Fibonacci
// LFSR: for x^n + x^k + 1 the mask carries bit 0 and bit (n-k).
//
// Bit 0 must always be set. The first version of this table used the taps in
// the other convention, which left bit 0 clear -- and a seed of 1 then fed back
// a zero, shifted to the all-zero state on the very first step, and emitted a
// constant from there. The sequence still had the right LENGTH; it just was not
// a sequence. Only the balance and autocorrelation tests caught it.
//
// Degrees below kMinLfsrDegree are zero: a code shorter than 31 chips buys
// under 15 dB and is not worth the framing.
//
// These are a table, and a table can be mistyped. The test suite therefore
// measures the PERIOD of every degree rather than trusting the entry -- a
// non-primitive polynomial produces a short cycle, which would silently destroy
// the processing gain while everything still appeared to run.
constexpr std::uint32_t kTaps[kMaxLfsrDegree + 1] = {
    0, 0, 0, 0, 0,
    0x0005u,   //  5: x^5 + x^3 + 1
    0x0003u,   //  6: x^6 + x^5 + 1
    0x0003u,   //  7: x^7 + x^6 + 1
    0x001Du,   //  8: x^8 + x^6 + x^5 + x^4 + 1
    0x0011u,   //  9: x^9 + x^5 + 1
    0x0009u,   // 10: x^10 + x^7 + 1
    0x0005u,   // 11: x^11 + x^9 + 1
    0x0107u,   // 12: x^12 + x^11 + x^10 + x^4 + 1
    0x0027u,   // 13: x^13 + x^12 + x^11 + x^8 + 1
    0x1007u,   // 14: x^14 + x^13 + x^12 + x^2 + 1
    0x0003u,   // 15: x^15 + x^14 + 1
};

// ---- GF(16), built on x^4 + x + 1 ----------------------------------------
constexpr std::uint8_t kGfExp[32] = {
    1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 1,
    2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 1, 2,
};
constexpr std::uint8_t kGfLog[16] = {
    0, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12,
};

constexpr std::uint8_t gf_mul(std::uint8_t a, std::uint8_t b) noexcept {
    if (a == 0 || b == 0) return 0;
    return kGfExp[(kGfLog[a] + kGfLog[b]) % 15];
}

constexpr std::uint8_t gf_inv(std::uint8_t a) noexcept {
    return (a == 0) ? 0 : kGfExp[(15 - kGfLog[a]) % 15];
}

constexpr std::uint8_t gf_div(std::uint8_t a, std::uint8_t b) noexcept {
    return (b == 0) ? 0 : gf_mul(a, gf_inv(b));
}

constexpr std::uint8_t gf_pow(std::uint8_t a, std::size_t n) noexcept {
    if (a == 0) return 0;
    return kGfExp[(kGfLog[a] * n) % 15];
}

// Generator polynomial (x-a^1)(x-a^2)(x-a^3)(x-a^4), low order first.
// Computed once at namespace scope so it cannot drift from kRsT.
struct Generator {
    std::uint8_t c[kRsN - kRsK + 1];
};

constexpr Generator make_generator() noexcept {
    Generator g{};
    g.c[0] = 1;
    std::size_t deg = 0;
    for (std::size_t i = 1; i <= kRsN - kRsK; ++i) {
        const std::uint8_t root = kGfExp[i];
        // Multiply by (x - root); in GF(2^m) subtraction is XOR.
        for (std::size_t j = deg + 1; j > 0; --j) {
            g.c[j] = static_cast<std::uint8_t>(g.c[j - 1] ^ gf_mul(g.c[j], root));
        }
        g.c[0] = gf_mul(g.c[0], root);
        ++deg;
    }
    return g;
}

constexpr Generator kGen = make_generator();

// ---- CRC-32 --------------------------------------------------------------
constexpr std::uint32_t crc32_entry(std::uint32_t i) noexcept {
    std::uint32_t c = i;
    for (int k = 0; k < 8; ++k) {
        c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
    }
    return c;
}

struct Crc32Table {
    std::uint32_t t[256];
};

constexpr Crc32Table make_crc32_table() noexcept {
    Crc32Table tab{};
    for (std::uint32_t i = 0; i < 256; ++i) tab.t[i] = crc32_entry(i);
    return tab;
}

constexpr Crc32Table kCrcTable = make_crc32_table();

}  // namespace

// ---------------------------------------------------------------------------
// PN sequences
// ---------------------------------------------------------------------------

std::size_t msequence_length(std::uint32_t degree) noexcept {
    if (degree < kMinLfsrDegree || degree > kMaxLfsrDegree) return 0;
    return (std::size_t{1} << degree) - 1;
}

std::size_t generate_msequence(std::uint32_t degree, std::uint32_t seed,
                               std::span<Real> out) noexcept {
    const std::size_t n = msequence_length(degree);
    if (n == 0 || out.size() < n) return 0;
    const std::uint32_t mask = static_cast<std::uint32_t>((std::uint32_t{1} << degree) - 1);
    std::uint32_t state = seed & mask;
    // An all-zero state is the LFSR's fixed point: it would emit zeros forever.
    if (state == 0) return 0;

    const std::uint32_t taps = kTaps[degree];
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint32_t bit = state & 1u;
        out[i] = bit ? kOne : -kOne;
        // Fibonacci LFSR: feed back the parity of the tapped bits.
        std::uint32_t fb = state & taps;
        fb ^= fb >> 16;
        fb ^= fb >> 8;
        fb ^= fb >> 4;
        fb ^= fb >> 2;
        fb ^= fb >> 1;
        state = (state >> 1) | ((fb & 1u) << (degree - 1));
    }
    return n;
}

// ---------------------------------------------------------------------------
// DSSS
// ---------------------------------------------------------------------------

Real processing_gain_db(std::size_t chips_per_bit) noexcept {
    if (chips_per_bit == 0) return kZero;
    return static_cast<Real>(10) * std::log10(static_cast<Real>(chips_per_bit));
}

Real data_rate_bps(const DsssConfig& cfg) noexcept {
    if (cfg.chips_per_bit == 0) return kZero;
    return cfg.chip_rate_hz / static_cast<Real>(cfg.chips_per_bit);
}

std::size_t samples_per_bit(const DsssConfig& cfg) noexcept {
    if (!(cfg.chip_rate_hz > kZero) || cfg.chips_per_bit == 0) return 0;
    const Real spc = cfg.sample_rate_hz / cfg.chip_rate_hz;
    if (!(spc >= kOne)) return 0;
    return static_cast<std::size_t>(spc * static_cast<Real>(cfg.chips_per_bit) + static_cast<Real>(0.5));
}

std::size_t dsss_modulate(const DsssConfig& cfg,
                          std::span<const std::uint8_t> bits,
                          std::span<const Real> chips,
                          std::span<Real> out) noexcept {
    const std::size_t spb = samples_per_bit(cfg);
    if (spb == 0 || chips.size() < cfg.chips_per_bit) return 0;
    const std::size_t total = spb * bits.size();
    if (out.size() < total) return 0;

    const Real spc = cfg.sample_rate_hz / cfg.chip_rate_hz;
    const Real w = static_cast<Real>(2) * kPi * cfg.carrier_hz / cfg.sample_rate_hz;

    for (std::size_t b = 0; b < bits.size(); ++b) {
        // Data is BPSK: bit 0 -> +1, bit 1 -> -1.
        const Real sym = bits[b] ? -kOne : kOne;
        for (std::size_t s = 0; s < spb; ++s) {
            const std::size_t ci = static_cast<std::size_t>(static_cast<Real>(s) / spc);
            const Real chip = chips[(ci < cfg.chips_per_bit) ? ci : cfg.chips_per_bit - 1];
            const std::size_t n = b * spb + s;
            out[n] = sym * chip * std::cos(w * static_cast<Real>(n));
        }
    }
    return total;
}

std::size_t dsss_demodulate(const DsssConfig& cfg,
                            std::span<const Real> samples,
                            std::span<const Real> chips,
                            Real phase_rad,
                            std::span<Real> out_soft,
                            std::span<std::uint8_t> out_bits) noexcept {
    const std::size_t spb = samples_per_bit(cfg);
    if (spb == 0 || chips.size() < cfg.chips_per_bit) return 0;
    const std::size_t n_bits = samples.size() / spb;
    if (out_soft.size() < n_bits) return 0;
    if (!out_bits.empty() && out_bits.size() < n_bits) return 0;

    const Real spc = cfg.sample_rate_hz / cfg.chip_rate_hz;
    const Real w = static_cast<Real>(2) * kPi * cfg.carrier_hz / cfg.sample_rate_hz;

    for (std::size_t b = 0; b < n_bits; ++b) {
        Real acc = kZero;
        for (std::size_t s = 0; s < spb; ++s) {
            const std::size_t ci = static_cast<std::size_t>(static_cast<Real>(s) / spc);
            const Real chip = chips[(ci < cfg.chips_per_bit) ? ci : cfg.chips_per_bit - 1];
            const std::size_t n = b * spb + s;
            // Coherent: multiply by the local carrier and the code, then sum.
            acc += samples[n] * chip * std::cos(w * static_cast<Real>(n) + phase_rad);
        }
        // Normalise by the number of samples so the value is comparable across
        // code lengths; the factor of 2 undoes the cos^2 average.
        const Real soft = static_cast<Real>(2) * acc / static_cast<Real>(spb);
        out_soft[b] = soft;
        if (!out_bits.empty()) out_bits[b] = (soft < kZero) ? std::uint8_t{1} : std::uint8_t{0};
    }
    return n_bits;
}

Real chip_slip(const DsssConfig& cfg, Real closing_speed_mps,
               Real sound_speed_mps) noexcept {
    if (!(sound_speed_mps > kZero)) return kZero;
    const Real scale = closing_speed_mps / sound_speed_mps;
    return std::fabs(scale) * static_cast<Real>(cfg.chips_per_bit);
}

std::size_t max_chips_per_bit(const DsssConfig& cfg, Real closing_speed_mps,
                              Real max_slip_chips, Real sound_speed_mps) noexcept {
    (void)cfg;
    if (!(sound_speed_mps > kZero) || !(max_slip_chips > kZero)) return 0;
    const Real scale = std::fabs(closing_speed_mps) / sound_speed_mps;
    if (!(scale > kZero)) return 0;
    const Real n = max_slip_chips / scale;
    return (n < kOne) ? 0 : static_cast<std::size_t>(n);
}

// ---------------------------------------------------------------------------
// CRC-32
// ---------------------------------------------------------------------------

std::uint32_t crc32(std::span<const std::uint8_t> data) noexcept {
    std::uint32_t c = 0xFFFFFFFFu;
    for (const std::uint8_t byte : data) {
        c = kCrcTable.t[(c ^ byte) & 0xFFu] ^ (c >> 8);
    }
    return c ^ 0xFFFFFFFFu;
}

// ---------------------------------------------------------------------------
// Reed-Solomon (15, 11)
// ---------------------------------------------------------------------------

bool rs_encode(std::span<const std::uint8_t> data, std::span<std::uint8_t> out) noexcept {
    if (data.size() != kRsK || out.size() < kRsN) return false;
    for (const std::uint8_t s : data) {
        if (s > 15) return false;
    }
    constexpr std::size_t kParity = kRsN - kRsK;
    std::uint8_t rem[kParity] = {};

    // Systematic: divide data * x^(n-k) by the generator, keep the remainder.
    for (std::size_t i = 0; i < kRsK; ++i) {
        const std::uint8_t feedback = static_cast<std::uint8_t>(data[i] ^ rem[kParity - 1]);
        for (std::size_t j = kParity - 1; j > 0; --j) {
            rem[j] = static_cast<std::uint8_t>(rem[j - 1] ^ gf_mul(feedback, kGen.c[j]));
        }
        rem[0] = gf_mul(feedback, kGen.c[0]);
    }
    for (std::size_t i = 0; i < kRsK; ++i) out[i] = data[i];
    for (std::size_t i = 0; i < kParity; ++i) out[kRsK + i] = rem[kParity - 1 - i];
    return true;
}

RsResult rs_decode(std::span<std::uint8_t> codeword, std::size_t& errors_corrected) noexcept {
    errors_corrected = 0;
    if (codeword.size() < kRsN) return RsResult::BadInput;
    for (std::size_t i = 0; i < kRsN; ++i) {
        if (codeword[i] > 15) return RsResult::BadInput;
    }

    // Syndromes S_i = r(a^i) for i = 1..4. The codeword polynomial has the data
    // in the high-order positions, so index 0 is the highest power.
    constexpr std::size_t kSyn = kRsN - kRsK;
    std::uint8_t s[kSyn] = {};
    bool any = false;
    for (std::size_t i = 0; i < kSyn; ++i) {
        const std::uint8_t root = kGfExp[i + 1];
        std::uint8_t acc = 0;
        for (std::size_t j = 0; j < kRsN; ++j) {
            acc = static_cast<std::uint8_t>(gf_mul(acc, root) ^ codeword[j]);
        }
        s[i] = acc;
        if (acc != 0) any = true;
    }
    if (!any) return RsResult::Clean;

    // Peterson for t <= 2: solve for the error locator directly. With only two
    // correctable errors the 2x2 system is small enough that a full
    // Berlekamp-Massey would be more code than insight.
    const std::uint8_t det = static_cast<std::uint8_t>(gf_mul(s[0], s[2]) ^ gf_mul(s[1], s[1]));

    std::uint8_t l1 = 0, l2 = 0;
    std::size_t n_errors = 0;
    if (det != 0) {
        // Two errors: [S0 S1; S1 S2] [L2; L1] = [S2; S3].
        l2 = gf_div(static_cast<std::uint8_t>(gf_mul(s[2], s[2]) ^ gf_mul(s[1], s[3])), det);
        l1 = gf_div(static_cast<std::uint8_t>(gf_mul(s[0], s[3]) ^ gf_mul(s[1], s[2])), det);
        n_errors = 2;
    } else if (s[0] != 0) {
        // One error: the locator is degree 1.
        l1 = gf_div(s[1], s[0]);
        l2 = 0;
        n_errors = 1;
        // Consistency: a single error must satisfy S2 = S1^2/S0 and so on.
        if (gf_mul(s[1], l1) != s[2] || gf_mul(s[2], l1) != s[3]) {
            return RsResult::Uncorrectable;
        }
    } else {
        return RsResult::Uncorrectable;
    }

    // Chien search: the locator L(x) = 1 + l1 x + l2 x^2 has roots at the
    // inverses of the error positions.
    std::size_t pos[kRsT] = {0, 0};
    std::size_t found = 0;
    for (std::size_t i = 0; i < kRsN; ++i) {
        const std::uint8_t x = gf_pow(2, (15 - i) % 15);   // a^-i
        const std::uint8_t v = static_cast<std::uint8_t>(
            1u ^ gf_mul(l1, x) ^ gf_mul(l2, gf_mul(x, x)));
        if (v == 0) {
            if (found >= kRsT) return RsResult::Uncorrectable;
            // Position i counted from the HIGH end, matching the syndrome loop.
            pos[found++] = kRsN - 1 - i;
        }
    }
    if (found != n_errors) return RsResult::Uncorrectable;

    // Forney would be the general answer; with at most two errors the magnitudes
    // fall out of the same syndromes directly.
    if (found == 1) {
        const std::uint8_t xi = gf_pow(2, (kRsN - 1 - pos[0]) % 15);
        const std::uint8_t mag = gf_div(s[0], gf_pow(xi, 1));
        codeword[pos[0]] = static_cast<std::uint8_t>(codeword[pos[0]] ^ mag);
    } else {
        const std::uint8_t x1 = gf_pow(2, (kRsN - 1 - pos[0]) % 15);
        const std::uint8_t x2 = gf_pow(2, (kRsN - 1 - pos[1]) % 15);
        // S0 = e1 x1 + e2 x2, S1 = e1 x1^2 + e2 x2^2.
        const std::uint8_t d = static_cast<std::uint8_t>(
            gf_mul(x1, gf_mul(x2, x2)) ^ gf_mul(x2, gf_mul(x1, x1)));
        if (d == 0) return RsResult::Uncorrectable;
        const std::uint8_t e1 = gf_div(
            static_cast<std::uint8_t>(gf_mul(s[0], gf_mul(x2, x2)) ^ gf_mul(s[1], x2)), d);
        const std::uint8_t e2 = gf_div(
            static_cast<std::uint8_t>(gf_mul(s[1], x1) ^ gf_mul(s[0], gf_mul(x1, x1))), d);
        codeword[pos[0]] = static_cast<std::uint8_t>(codeword[pos[0]] ^ e1);
        codeword[pos[1]] = static_cast<std::uint8_t>(codeword[pos[1]] ^ e2);
    }

    // Verify: recompute the syndromes. Peterson can "solve" a system produced by
    // more errors than it can fix and return a confident wrong codeword, and
    // this is the check that turns that into an honest Uncorrectable.
    for (std::size_t i = 0; i < kSyn; ++i) {
        const std::uint8_t root = kGfExp[i + 1];
        std::uint8_t acc = 0;
        for (std::size_t j = 0; j < kRsN; ++j) {
            acc = static_cast<std::uint8_t>(gf_mul(acc, root) ^ codeword[j]);
        }
        if (acc != 0) return RsResult::Uncorrectable;
    }

    errors_corrected = found;
    return RsResult::Corrected;
}

}  // namespace phantom::comm
