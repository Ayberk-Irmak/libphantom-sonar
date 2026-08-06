// SPDX-License-Identifier: Apache-2.0
// libphantom-sonar — compile-time configuration and scalar type selection.
#ifndef PHANTOM_CONFIG_HPP
#define PHANTOM_CONFIG_HPP

#include <cstddef>
#include <cstdint>

#define PHANTOM_VERSION_MAJOR 0
#define PHANTOM_VERSION_MINOR 1
#define PHANTOM_VERSION_PATCH 0
#define PHANTOM_VERSION_STRING "0.1.0"

namespace phantom {

// The library is scalar-type agnostic. `double` is the default because the
// Snell invariant xi = cos(theta)/c must stay stable over tens of thousands of
// arc steps; `float` loses ~7 digits and drifts visibly past ~50 km of range.
// Define PHANTOM_REAL_FLOAT for MCU targets with a single-precision FPU and
// accept the reduced range accuracy documented in docs/validation.md.
#if defined(PHANTOM_REAL_FLOAT)
using Real = float;
#else
using Real = double;
#endif

}  // namespace phantom

#endif  // PHANTOM_CONFIG_HPP
