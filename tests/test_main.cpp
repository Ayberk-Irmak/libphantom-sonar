// SPDX-License-Identifier: Apache-2.0
#include "framework.hpp"

#include "phantom/config.hpp"

int main() {
    std::printf("libphantom-sonar %s  |  Real = %s  |  %zu test cases\n",
                PHANTOM_VERSION_STRING,
                sizeof(phantom::Real) == 4 ? "float" : "double",
                pt::g_case_count);
    std::printf("--------------------------------------------------------------\n");

    std::size_t failed_cases = 0;
    for (std::size_t i = 0; i < pt::g_case_count; ++i) {
        const long before = pt::g_failures;
        pt::g_cases[i].fn();
        const bool ok = (pt::g_failures == before);
        if (!ok) ++failed_cases;
        std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", pt::g_cases[i].name);
    }

    std::printf("--------------------------------------------------------------\n");
    std::printf("%ld checks, %ld failures, %zu/%zu cases failed\n",
                pt::g_checks, pt::g_failures, failed_cases, pt::g_case_count);
    return pt::g_failures == 0 ? 0 : 1;
}
