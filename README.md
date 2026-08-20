# libphantom-sonar

**Real-time ocean acoustics and sonar *simulation* library.**
Zero dependencies, zero heap allocation, C++20.

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![Tests](https://img.shields.io/badge/tests-95846%20checks-brightgreen.svg)
![Validated](https://img.shields.io/badge/validated-Bellhop%20%7C%20UNESCO%20%7C%20ISO%209613--1-brightgreen.svg)
![Portable](https://img.shields.io/badge/portable-x86--64%20%7C%20ARM%20%7C%20RISC--V-blue.svg)

![Munk deep sound channel propagation](data/munk_rays.png)

*61 rays through the canonical Munk deep sound channel, 100 km. Left: the sound
speed profile and its SOFAR axis. Top: ray paths as exact circular arcs. Bottom:
ensonification, with shadow zones in white. Reproduce with
`./build/munk_simulation data && python3 tools/plot_rays.py data`.*

---

## What this is

A library for modelling how sound travels through the ocean, built to run in a
real-time loop on hardware with no operating system and no allocator.

It covers the chain end to end: sound speed and ray tracing, transmission loss
and eigenrays, boundary reflection and reverberation, waveform synthesis and
matched-filter detection, beamforming, multi-target tracking, and spread-spectrum
acoustic communication. It also models **acoustics in air**, because that is the
medium anyone can test in.

It is a simulation and analysis library. It is not deployable firmware, and it
models no specific platform. See [`EXPORT_NOTICE.md`](EXPORT_NOTICE.md).

## What makes it different

**Nothing is verified against itself.** Every claim in this repository is checked
against a closed form, a theoretical bound, or a number published by somebody
else — and where that was not possible, the documentation says so.

| Verified against | Result |
|---|---|
| Bellhop (Acoustics Toolbox), shared `.env` files | turning depths within 0.01 m |
| UNESCO Tech. Paper 44 check value | **1731.9954** vs published 1731.995 |
| UNESCO Tech. Paper 44, 220-value table | 0.0499 m/s worst — the table's own rounding |
| ISO 9613-1:1993 Table 1, 105 values | **0.380%** worst relative error |
| CRC-32 published check value | 0xCBF43926 |
| Bearing estimator vs its Cramér–Rao bound | 0.99–1.13× across 14 dB |
| Tracker NIS vs its χ² distribution | mean 1.951 vs 2.000 |
| FFT vs an independent O(N²) DFT | 2.3e-16 relative |

**The failures are kept, not buried.** Several releases exist because a
measurement was wrong in a way that looked right:

- The sound-speed equation carried **two coefficients from the wrong temperature
  scale** for eleven releases. Mutual agreement between three equations could
  never see it — they agree only to 0.1 m/s and the error was 0.016. Only the
  published check value found it.
- A "processing gain" measurement reproduced theory **exactly** — 12.17 dB
  against 12.17 dB — and was meaningless, because both sides were the same
  tautology. Against white noise at fixed energy per bit, spreading buys nothing.
- A Snell's-invariant test reported **exactly zero** drift, which is suspicious
  rather than reassuring: the tracer derives the angle from the invariant, so the
  test largely measured whether `acos` and `cos` round-trip.
- The C ABI carried **160 kB of `.bss`** — half the RAM of the target part — in a
  buffer that copied data to avoid an assumption that could simply be asserted.

The reasoning for each is in [`docs/validation.md`](docs/validation.md) and
[`CHANGELOG.md`](CHANGELOG.md).

**The ray tracer is analytic, not stepped.** The profile is piecewise *linear* in
depth, and inside a constant-gradient layer a ray is exactly a circular arc — so
each layer crossing is one closed-form step with no integration error and no step
size to tune.

**It allocates nothing.** Proven three ways: a source audit, `nm` over the built
archive, and an ARM cross-build reporting `data == 0` and `bss == 0`. All three
run in CI.

## Quick start

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

./build/phantom_tests           # 95846 checks across 189 cases
./build/phantom_c_tests         # 269 checks, through the C ABI
./build/phantom_bench           # the timing table below

./build/munk_simulation data    # CSV + console report
./build/ping_intercept data     # streaming pulse detection
./build/countermeasure_loop     # intercept -> classify -> reply -> multipath
python3 tools/plot_rays.py data
```

Cross-validate the ray tracer against Bellhop yourself:

```bash
./tools/bellhop_compare/setup_bellhop.sh
python3 tools/bellhop_compare/compare.py \
        --bellhop build/acoustics-toolbox/at/Bellhop/bellhop.exe --check
```

### Using it

```cpp
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

phantom::SoundSpeedProfile<2048> svp;
for (int i = 0; i <= 200; ++i) {
    const phantom::Real z = i * 25.0;
    svp.push(z, phantom::sound_speed::munk(z));
}

phantom::TraceConfig cfg;
cfg.max_range_m = 100000;

std::array<phantom::RayPoint, 4096> path;
const auto result = phantom::trace_ray(svp.view(), 1000, phantom::deg2rad(5), cfg, path);
```

From C, Rust, or an embedded target:

```bash
# C ABI — caller owns the memory, the library allocates nothing
#include "phantom/phantom.h"

# Rust
cd bindings/rust && cargo test

# Cortex-M class part: single precision, no exceptions, no RTTI
cmake -S . -B build-m7 -DCMAKE_TOOLCHAIN_FILE=cmake/cortex-m7.cmake \
      -DPHANTOM_REAL_FLOAT=ON -DPHANTOM_NO_EXCEPTIONS=ON
```

## Measured performance

13th Gen Core i7-13620H, `-O3`, single thread:

| Operation | Time |
|---|---|
| `mackenzie()` | 1.5 ns |
| `chen_millero()` (UNESCO) | 4.0 ns |
| `speed_at()` — 501-layer binary search | 16.6 ns |
| **`trace_ray()` — per arc step** | **52 ns** |
| `trace_ray()` — full 100 km ray (1132 arcs) | 40.8 µs |
| `fft_forward()` — 8192 points | 307 µs |
| `analyze_block()` — 4-template bank, 8192-point block | 2.26 ms |

The ping analyser runs **29× real time on one core** at 96 kHz. Detection latency
is 65 ms and is set by the block length, not the arithmetic.

The budget is quoted **per arc step**, which is the unit that constrains a control
loop. A 100 km trace crosses ~1100 layers and belongs in no sub-microsecond
budget on any hardware — see [`docs/validation.md`](docs/validation.md) for why
the original "< 500 ns per cycle" requirement was restated this way. Re-run
`phantom_bench` on your target before treating any of this as a guarantee.

Library size, ARM 32-bit, `-Os`, float: **64 kB of text, zero static storage.**

## What it does not do

The gaps matter more than the features, so they are listed rather than implied.
The full set is in [`docs/validation.md`](docs/validation.md).

**Nothing here has been measured against a real transducer.** Every number in
this repository is a closed form, a published reference, or a simulation. The
apparatus for a real measurement exists and is self-verified
(`tools/air_bench.py`); the development machine has no working acoustic path, and
that is documented rather than glossed.

- **Ray theory, with its caustics.** Levels near a caustic are flagged, not
  computed; a correct treatment needs a wave solution.
- **Range-independent ocean.** `c` depends on depth alone — no bathymetry,
  fronts or eddies. 2-D range–depth only.
- **Reverberation is a level, not a field**, and does not come from the traced
  eigenrays.
- **Plane-wave, flat-interface reflection.** No beam displacement near the
  critical angle, no sediment layering, no shear.
- **Monostatic only.** Target strength is a dB figure, not a scattering model.
- **No carrier or timing recovery** in the communication engine, and no
  equaliser. The demodulator is *given* the carrier phase.
- **Uniform line arrays only** — no planar geometry, no element directivity.
- **The C ABI is partial.** Beamforming, eigenrays and reverberation are not
  exposed; their C++ interfaces take several spans at once and a C shape has not
  been designed.
- **The Cortex-M7 toolchain file is untested** — this environment's
  `arm-none-eabi` ships no C++ standard library. ARM figures come from
  `arm-linux-gnueabihf`.

## Documentation

| | |
|---|---|
| [`docs/math_spec.md`](docs/math_spec.md) | Every equation, with derivations and the reasoning behind each choice |
| [`docs/validation.md`](docs/validation.md) | What is verified, what is only self-consistent, what is not verified at all |
| [`docs/roadmap.md`](docs/roadmap.md) | Where it came from and where it could go |
| [`docs/hardware.md`](docs/hardware.md) | Building a physical test bench |
| [`CHANGELOG.md`](CHANGELOG.md) | Per-release results, including the corrections |
| [`EXPORT_NOTICE.md`](EXPORT_NOTICE.md) | Scope and classification notes |

## References

- Fofonoff, N. P. & Millard, R. C. (1983). *Algorithms for computation of
  fundamental properties of seawater.* UNESCO Technical Papers in Marine Science
  44.
- Chen, C.-T. & Millero, F. J. (1977). Speed of sound in seawater at high
  pressures. *JASA* 62(5), 1129–1135.
- Wong, G. S. K. & Zhu, S. (1995). Speed of sound in seawater as a function of
  salinity, temperature and pressure. *JASA* 97(3), 1732–1736.
- Mackenzie, K. V. (1981). Nine-term equation for sound speed in the oceans.
  *JASA* 70(3), 807–812.
- Leroy, C. C. & Parthiot, F. (1998). Depth-pressure relationships in the oceans
  and seas. *JASA* 103(3), 1346–1352.
- Munk, W. H. (1974). Sound channel in an exponentially stratified ocean.
  *Deep-Sea Research* 21, 207–217.
- ISO 9613-1:1993. *Attenuation of sound during propagation outdoors — Part 1.*
- Porter, M. B. *The BELLHOP Manual and User's Guide*, Heat, Light and Sound
  Research.
- Locarnini, R. A. et al. (2024). *World Ocean Atlas 2023, Volume 1:
  Temperature.* NOAA Atlas NESDIS 89.

## License

Apache-2.0. See [`LICENSE`](LICENSE).

Bellhop and the Acoustics Toolbox are GPL-3 and are **not** bundled or linked —
`tools/bellhop_compare/` fetches and builds them separately, and the two codes
communicate only through `.env` files and text output.

World Ocean Atlas 2023 data in `data/profiles/` is a US Government work in the
public domain; provenance is recorded in each file's header.
