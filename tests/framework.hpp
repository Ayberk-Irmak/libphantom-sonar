// SPDX-License-Identifier: Apache-2.0
// A ~90 line test harness. Deliberately hand-rolled: the library advertises
// zero external dependencies, and pulling GoogleTest in just for the test
// binary would make that claim a half-truth.
#ifndef PHANTOM_TEST_FRAMEWORK_HPP
#define PHANTOM_TEST_FRAMEWORK_HPP

#include "phantom/config.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>

namespace pt {

// The algorithms are identical in float and double builds; only the achievable
// precision differs (~1e-7 vs ~1e-16 relative). Tests therefore state BOTH
// tolerances explicitly instead of loosening the double one to match float --
// otherwise a real double-precision regression would slip through unnoticed.
inline constexpr bool kRealIsFloat = sizeof(::phantom::Real) == 4;

constexpr double tol(double double_build, double float_build) noexcept {
    return kRealIsFloat ? float_build : double_build;
}


// Deterministic RNG. The DSP tests are statistical -- estimator variance
// against the Cramer-Rao bound, detection rates -- so they need randomness that
// is identical on every machine and every run, or CI becomes a coin flip.
// A 64-bit LCG is more than enough for Monte Carlo of this kind.
class Rng {
  public:
    explicit Rng(std::uint64_t seed) noexcept : state_(seed * 2862933555777941757ULL + 1) {}

    // Uniform in [0, 1).
    double uniform01() noexcept {
        state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((state_ >> 11) & ((1ULL << 53) - 1))
             / static_cast<double>(1ULL << 53);
    }

    // Uniform in [-1, 1).
    double uniform() noexcept { return uniform01() * 2.0 - 1.0; }

    // Standard normal, Box-Muller. Both values are used, so no sample is wasted.
    double normal() noexcept {
        if (have_spare_) {
            have_spare_ = false;
            return spare_;
        }
        double u1 = uniform01();
        if (u1 < 1e-300) u1 = 1e-300;
        const double u2 = uniform01();
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * 3.14159265358979323846 * u2;
        spare_ = r * std::sin(theta);
        have_spare_ = true;
        return r * std::cos(theta);
    }

  private:
    std::uint64_t state_;
    double spare_ = 0.0;
    bool   have_spare_ = false;
};

using TestFn = void (*)();

struct Case {
    const char* name = nullptr;
    TestFn      fn   = nullptr;
};

inline constexpr std::size_t kMaxCases = 256;
inline Case        g_cases[kMaxCases]{};
inline std::size_t g_case_count = 0;
inline long        g_checks     = 0;
inline long        g_failures   = 0;

struct Registrar {
    Registrar(const char* name, TestFn fn) noexcept {
        if (g_case_count < kMaxCases) g_cases[g_case_count++] = Case{name, fn};
    }
};

inline bool check(bool ok, const char* file, int line, const char* expr) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("    FAIL %s:%d\n         %s\n", file, line, expr);
    }
    return ok;
}

inline bool check_near(double got, double want, double tol,
                       const char* file, int line, const char* expr) {
    ++g_checks;
    const double diff = std::fabs(got - want);
    if (!(diff <= tol)) {
        ++g_failures;
        std::printf("    FAIL %s:%d\n         %s\n"
                    "         got  %.17g\n         want %.17g\n"
                    "         |diff| %.6g  >  tol %.6g\n",
                    file, line, expr, got, want, diff, tol);
        return false;
    }
    return true;
}

inline bool check_rel(double got, double want, double rel_tol,
                      const char* file, int line, const char* expr) {
    const double scale = std::fabs(want) > 1.0 ? std::fabs(want) : 1.0;
    return check_near(got, want, rel_tol * scale, file, line, expr);
}

}  // namespace pt

#define PT_TEST(name)                                                   \
    static void name();                                                 \
    static ::pt::Registrar pt_registrar_##name(#name, &name);           \
    static void name()

#define PT_CHECK(cond) ::pt::check                                      \
    (static_cast<bool>(cond), __FILE__, __LINE__, #cond)

#define PT_CHECK_NEAR(got, want, tol) ::pt::check_near                  \
    (static_cast<double>(got), static_cast<double>(want),               \
     static_cast<double>(tol), __FILE__, __LINE__, #got " ~= " #want)

#define PT_CHECK_REL(got, want, rel) ::pt::check_rel                    \
    (static_cast<double>(got), static_cast<double>(want),               \
     static_cast<double>(rel), __FILE__, __LINE__, #got " ~= " #want)

#endif  // PHANTOM_TEST_FRAMEWORK_HPP
