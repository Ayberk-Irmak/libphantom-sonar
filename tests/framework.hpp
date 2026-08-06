// SPDX-License-Identifier: Apache-2.0
// A ~90 line test harness. Deliberately hand-rolled: the library advertises
// zero external dependencies, and pulling GoogleTest in just for the test
// binary would make that claim a half-truth.
#ifndef PHANTOM_TEST_FRAMEWORK_HPP
#define PHANTOM_TEST_FRAMEWORK_HPP

#include "phantom/config.hpp"

#include <cmath>
#include <cstddef>
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
