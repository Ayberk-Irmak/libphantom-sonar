// SPDX-License-Identifier: Apache-2.0
// FFT verification.
//
// A hand-written transform is exactly the sort of code that "looks right" and
// is subtly wrong -- a swapped twiddle sign or an off-by-one in the bit
// reversal produces output that is still plausible-looking noise. So it is
// checked against a direct O(N^2) DFT, against Parseval, and against transforms
// whose closed form is known.
#include "framework.hpp"

#include "phantom/fft.hpp"

#include <array>
#include <cmath>

using namespace phantom;

namespace {

// Textbook DFT, straight from the definition. Slow on purpose: it shares no
// code with the implementation under test, which is the point.
void direct_dft(std::span<const Complex> in, std::span<Complex> out) {
    const std::size_t n = in.size();
    for (std::size_t k = 0; k < n; ++k) {
        double re = 0.0, im = 0.0;
        for (std::size_t t = 0; t < n; ++t) {
            const double ang = -2.0 * 3.14159265358979323846
                             * static_cast<double>(k) * static_cast<double>(t)
                             / static_cast<double>(n);
            re += static_cast<double>(in[t].real()) * std::cos(ang)
                - static_cast<double>(in[t].imag()) * std::sin(ang);
            im += static_cast<double>(in[t].real()) * std::sin(ang)
                + static_cast<double>(in[t].imag()) * std::cos(ang);
        }
        out[k] = Complex(static_cast<Real>(re), static_cast<Real>(im));
    }
}

double max_rel_error(std::span<const Complex> a, std::span<const Complex> b) {
    double scale = 0.0;
    for (const Complex& c : b) scale = std::max(scale, static_cast<double>(std::abs(c)));
    if (scale <= 0.0) scale = 1.0;
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        worst = std::max(worst, static_cast<double>(std::abs(a[i] - b[i])) / scale);
    }
    return worst;
}

}  // namespace

PT_TEST(fft_matches_direct_dft) {
    constexpr std::size_t kN = 256;
    FftPlan<kN> plan;
    std::array<Complex, kN> x{}, fast{}, slow{};

    pt::Rng rng(12345);
    for (std::size_t i = 0; i < kN; ++i) {
        x[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
    }
    fast = x;
    plan.forward(fast);
    direct_dft(x, slow);

    const double err = max_rel_error(fast, slow);
    std::printf("       max relative error vs direct DFT = %.3g\n", err);
    PT_CHECK(err < pt::tol(1e-13, 1e-5));
}

PT_TEST(fft_round_trip_is_identity) {
    constexpr std::size_t kN = 1024;
    FftPlan<kN> plan;
    std::array<Complex, kN> x{}, y{};

    pt::Rng rng(999);
    for (std::size_t i = 0; i < kN; ++i) {
        x[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
    }
    y = x;
    plan.forward(y);
    plan.inverse(y);

    double worst = 0.0;
    for (std::size_t i = 0; i < kN; ++i) {
        worst = std::max(worst, static_cast<double>(std::abs(y[i] - x[i])));
    }
    std::printf("       round-trip worst absolute error = %.3g\n", worst);
    PT_CHECK(worst < pt::tol(1e-14, 1e-5));
}

PT_TEST(fft_satisfies_parseval) {
    // sum |x[n]|^2 == (1/N) sum |X[k]|^2. Catches any scaling error.
    constexpr std::size_t kN = 512;
    FftPlan<kN> plan;
    std::array<Complex, kN> x{}, X{};

    pt::Rng rng(7);
    double time_energy = 0.0;
    for (std::size_t i = 0; i < kN; ++i) {
        x[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
        time_energy += static_cast<double>(std::norm(x[i]));
    }
    X = x;
    plan.forward(X);

    double freq_energy = 0.0;
    for (const Complex& c : X) freq_energy += static_cast<double>(std::norm(c));
    freq_energy /= static_cast<double>(kN);

    PT_CHECK_REL(freq_energy, time_energy, pt::tol(1e-13, 1e-5));
}

PT_TEST(fft_known_transforms) {
    constexpr std::size_t kN = 64;
    FftPlan<kN> plan;
    std::array<Complex, kN> buf{};

    // A unit impulse at n=0 transforms to a flat spectrum of ones.
    buf.fill(Complex(0, 0));
    buf[0] = Complex(1, 0);
    plan.forward(buf);
    for (std::size_t k = 0; k < kN; ++k) {
        PT_CHECK_NEAR(buf[k].real(), 1.0, pt::tol(1e-14, 1e-6));
        PT_CHECK_NEAR(buf[k].imag(), 0.0, pt::tol(1e-14, 1e-6));
    }

    // A constant transforms to an impulse of height N at k=0.
    buf.fill(Complex(1, 0));
    plan.forward(buf);
    PT_CHECK_NEAR(buf[0].real(), static_cast<double>(kN), pt::tol(1e-12, 1e-3));
    for (std::size_t k = 1; k < kN; ++k) {
        PT_CHECK(static_cast<double>(std::abs(buf[k])) < pt::tol(1e-12, 1e-3));
    }

    // A complex exponential at bin 5 lands entirely in bin 5. Any bit-reversal
    // error scatters it to the wrong bin, which this pins exactly.
    constexpr std::size_t kBin = 5;
    for (std::size_t n = 0; n < kN; ++n) {
        const double ang = 2.0 * 3.14159265358979323846
                         * static_cast<double>(kBin) * static_cast<double>(n)
                         / static_cast<double>(kN);
        buf[n] = Complex(static_cast<Real>(std::cos(ang)), static_cast<Real>(std::sin(ang)));
    }
    plan.forward(buf);
    PT_CHECK_NEAR(std::abs(buf[kBin]), static_cast<double>(kN), pt::tol(1e-11, 1e-2));
    for (std::size_t k = 0; k < kN; ++k) {
        if (k == kBin) continue;
        PT_CHECK(static_cast<double>(std::abs(buf[k])) < pt::tol(1e-11, 1e-2));
    }
}

PT_TEST(fft_is_linear) {
    constexpr std::size_t kN = 128;
    FftPlan<kN> plan;
    std::array<Complex, kN> a{}, b{}, sum{}, fa{}, fb{};

    pt::Rng rng(4242);
    for (std::size_t i = 0; i < kN; ++i) {
        a[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
        b[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
        sum[i] = a[i] + static_cast<Real>(3) * b[i];
    }
    fa = a; fb = b;
    plan.forward(fa);
    plan.forward(fb);
    plan.forward(sum);

    double worst = 0.0;
    for (std::size_t i = 0; i < kN; ++i) {
        const Complex expect = fa[i] + static_cast<Real>(3) * fb[i];
        worst = std::max(worst, static_cast<double>(std::abs(sum[i] - expect)));
    }
    PT_CHECK(worst < pt::tol(1e-12, 1e-3));
}

PT_TEST(fft_rejects_bad_sizes) {
    std::array<Complex, 5> odd{};
    PT_CHECK(!fft_init(odd));           // 10 is not a power of two
    std::array<Complex, 0> empty{};
    PT_CHECK(!fft_init(empty));

    FftPlan<32> plan;
    PT_CHECK(plan.view().valid());
    // Mismatched buffer length is a no-op, not a buffer overrun.
    std::array<Complex, 16> wrong{};
    wrong.fill(Complex(1, 1));
    fft_forward(plan.view(), wrong);
    PT_CHECK_NEAR(wrong[0].real(), 1.0, 1e-12);
}

PT_TEST(fft_sizes_from_2_to_4096) {
    // Every power of two the DSP chain might use, checked by round trip.
    FftPlan<2> p2; FftPlan<8> p8; FftPlan<64> p64;
    FftPlan<512> p512; FftPlan<4096> p4096;

    auto check = [](auto& plan, std::size_t n) {
        static std::array<Complex, 4096> x{}, y{};
        pt::Rng rng(n * 31 + 7);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = Complex(static_cast<Real>(rng.uniform()), static_cast<Real>(rng.uniform()));
        }
        std::span<Complex> sx(x.data(), n), sy(y.data(), n);
        for (std::size_t i = 0; i < n; ++i) sy[i] = sx[i];
        plan.forward(sy);
        plan.inverse(sy);
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            worst = std::max(worst, static_cast<double>(std::abs(sy[i] - sx[i])));
        }
        return worst;
    };

    PT_CHECK(check(p2, 2) < pt::tol(1e-15, 1e-6));
    PT_CHECK(check(p8, 8) < pt::tol(1e-15, 1e-6));
    PT_CHECK(check(p64, 64) < pt::tol(1e-15, 1e-6));
    PT_CHECK(check(p512, 512) < pt::tol(1e-14, 1e-5));
    PT_CHECK(check(p4096, 4096) < pt::tol(1e-14, 1e-4));
}
