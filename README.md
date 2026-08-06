# libphantom-sonar

**Real-time ocean acoustics and sonar countermeasure *simulation* library.**
Zero dependencies, zero heap allocation, C++20.

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![Status](https://img.shields.io/badge/status-v0.1%20ray%20tracing-orange.svg)
![Validated](https://img.shields.io/badge/validated-vs%20Bellhop-brightgreen.svg)

![Munk deep sound channel propagation](data/munk_rays.png)

*61 rays through the canonical Munk deep sound channel, 100 km. Left: the sound
speed profile and its SOFAR axis. Top: ray paths as exact circular arcs. Bottom:
ensonification, with shadow zones in white. Reproduce with
`./build/munk_simulation data && python3 tools/plot_rays.py data`.*

---

## What this is

A library for modelling how sound actually travels through the ocean, built to
run in a real-time control loop on hardware that has no operating system and no
allocator.

**v0.1 ships the ocean model and the ray tracer**, cross-validated against
Bellhop. Echo synthesis and covert communication are designed but not implemented — see [the roadmap](docs/roadmap.md).
The repository states what is verified, what is only self-consistent, and what
is not verified at all, in [`docs/validation.md`](docs/validation.md).

It is a simulation and analysis library. It is not deployable countermeasure
firmware, and it models no specific platform. See [`EXPORT_NOTICE.md`](EXPORT_NOTICE.md).

## Why it might interest you

**The ray tracer is analytic, not stepped.** The sound speed profile is stored as
piecewise *linear* in depth, and inside a constant-gradient layer an acoustic ray
is exactly a circular arc. So each layer crossing is one closed-form step with no
integration error and no step size to tune:

```
ξ  = cos θ / c              Snell invariant, constant along the ray
R  = 1 / (ξ |g|)            arc radius,  g = dc/dz
Δr = (sin θ₁ − sin θ₂) / (ξ g)
Δt = ln[ c₂(1 + sin θ₁) / (c₁(1 + sin θ₂)) ] / g
```

That choice is measurable, not aesthetic. Resampling a linear profile from 2
points to 1001 points — the same ocean, described 500× more finely — moves the
answer by **4.6e-12 m over 40 km**. A stepping tracer moves by metres. Full
derivation in [`docs/math_spec.md`](docs/math_spec.md).

**Three sound speed equations, because the usual one is mislabelled.** The
six-term polynomial widely circulated as "the UNESCO equation" is actually
Medwin's simplification, valid only to 1000 m — which is *at or above the deep
sound channel axis*, making it the wrong tool for exactly the problem people
reach for it to solve. This library implements Medwin, Mackenzie (the default,
valid to 8000 m) and Chen & Millero (the real UNESCO algorithm, which takes
pressure rather than depth), and cross-validates them against each other in the
test suite.

**Cross-validated against Bellhop.** The reference ray code of the Acoustics
Toolbox, run on the same `.env` file — `phantom_trace` parses Bellhop's own input
format, so a transcription bug cannot masquerade as agreement.

![Bellhop cross-validation](data/bellhop_comparison.png)

The clean measurement is the **turning depth**: one number per ray, needing no
interpolation, solved for exactly here and approached by Bellhop as its
integrator step shrinks.

| Bellhop step | turning-depth error |
|---|---|
| 500 m | 0.084 m |
| 100 m | 0.0024 m |
| **10 m** | **0.00014 m** |

**0.14 mm over a 101 km path.** The direction is the point: Bellhop converges
*toward* this library's answer, which is what should happen if the analytic arc
solution is the exact limit a stepping integrator approaches.

Full-path RMS deviation is 0.27 m for trapped rays and 0.64 m with surface and
bottom bounces, at Bellhop's default step — but those numbers are a **ceiling,
not a measurement**, and [`docs/validation.md §2`](docs/validation.md) shows why:
the residual halves every time this library's *output sampling* is doubled, so it
is the comparison's chord-versus-arc error rather than a disagreement between the
codes.

**Verified against closed forms, not against baselines.** A recorded baseline
tells you the code did not change. These tell you it is right:

| Property | Result |
|---|---|
| Snell invariant over 29 rays × 100 km | 1.7e-16 relative drift |
| Ray path vs analytic circle, 200 layer crossings | 1.9e-16 relative radial error |
| Turning depth vs `c(z) = c₀/cos θ₀` | < 1e-6 m |
| Range and travel-time budgets | exact to 1e-8 m / 1e-12 s |
| Mackenzie vs Chen-Millero, common validity box | 0.53 m/s max disagreement |

791 checks. Clean under GCC 15 and Clang 21 with `-Wconversion -Wsign-conversion
-Wold-style-cast -Wdouble-promotion -Werror`, and under ASan + UBSan. Passes in
both `double` and `float` builds.

## Measured performance

13th Gen Core i7-13620H, `-O3`, single thread:

| Operation | Time |
|---|---|
| `mackenzie()` | 1.5 ns |
| `chen_millero()` (UNESCO) | 4.0 ns |
| `speed_at()` — 501-layer binary search | 16.6 ns |
| **`trace_ray()` — per arc step** | **52 ns** |
| `trace_ray()` — full 100 km ray (1132 arcs) | 40.8 µs |
| `analyze_sofar()` — 501 points | 543 ns |

The budget is stated **per arc step**, which is the unit that constrains a
control loop. A whole 100 km trace crosses ~1100 layers and does not belong in a
sub-microsecond budget on any hardware — see
[`docs/validation.md §4`](docs/validation.md) for why the original "< 500 ns per
cycle" requirement was restated this way. Re-run `phantom_bench` on your target
before treating any of this as a real-time guarantee.

## Quick start

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

./build/phantom_tests          # 791 checks
./build/phantom_bench          # the table above
./build/munk_simulation data   # CSV + channel analysis
python3 tools/plot_rays.py data

# cross-validate against Bellhop (downloads and builds the Acoustics Toolbox)
./tools/bellhop_compare/setup_bellhop.sh
python3 tools/bellhop_compare/compare.py \
        --bellhop build/acoustics-toolbox/at/Bellhop/bellhop.exe --check
```

```cpp
#include "phantom/channel.hpp"
#include "phantom/profile.hpp"
#include "phantom/ray_tracer.hpp"
#include "phantom/sound_speed.hpp"

using namespace phantom;

// Fixed capacity, no heap. 501 samples of the Munk canonical profile.
SoundSpeedProfile<512> svp;
fill_profile(svp, 0, 5000, 501, [](Real z) { return sound_speed::munk(z); });

// Where is the deep sound channel, and how wide is its trapping cone?
const ChannelInfo ch = analyze_sofar(svp.view());
// -> axis 1300 m, cone +/- 14.381 deg, conjugate depths 0 .. 4800 m

// Trace one ray. The caller owns the buffer; the library never allocates.
std::array<RayPoint, 4096> path;
TraceConfig cfg;
cfg.max_range_m    = 100000;
cfg.bottom_depth_m = 5000;

const TraceResult r = trace_ray(svp.view(), ch.axis_depth_m, deg2rad(10), cfg, path);
// r.status, r.turning_points, r.final_time_s, path[0 .. r.point_count)
```

### Build options

| Option | Default | Effect |
|---|---|---|
| `PHANTOM_BUILD_TESTS` | ON | unit test binary |
| `PHANTOM_BUILD_EXAMPLES` | ON | `munk_simulation` |
| `PHANTOM_BUILD_BENCH` | ON | `phantom_bench` |
| `PHANTOM_BUILD_TOOLS` | ON | `phantom_trace` for the Bellhop comparison |
| `PHANTOM_REAL_FLOAT` | OFF | `Real = float` for single-precision FPUs |
| `PHANTOM_WERROR` | OFF | warnings become errors |
| `PHANTOM_SANITIZE` | OFF | ASan + UBSan on the test binary |

`-ffast-math` is deliberately never enabled: it would let the compiler
reassociate the Snell invariant and break bit-reproducibility across targets,
which is the opposite of what a deterministic library should do.

## What v0.1 does not do

Stated plainly, because the gaps matter more than the features:

- **Geometric acoustics only** — ray paths, not ray intensities. No transmission
  loss, absorption, diffraction or caustic amplitude.
- **Range-independent ocean.** `c` is a function of depth alone; no bathymetry,
  fronts or eddies.
- **Shadow fractions depend on ray density.** At 5761 rays the Munk case is still
  drifting ~0.5 points per doubling. Treat any reported shadow fraction as an
  upper bound and always publish the ray count with it. The convergence table is
  in [`docs/validation.md §6`](docs/validation.md).
- **2-D range–depth only.** No out-of-plane refraction.
- **The Bellhop comparison covers geometry, not amplitude**, and only
  range-independent profiles — because that is all the library models. It is not
  evidence about transmission loss or range-dependent propagation.
- No echo synthesis, no ping analysis, no communication engine — see the
  [roadmap](docs/roadmap.md).

On countermeasures specifically: anti-phase cancellation of a vehicle's acoustic
cross section is **not** achievable in practice with a point transducer against a
broadband ping across all bearings — hull scattering is distributed and
broadband. What v0.3 will implement is narrowband, single-bearing nulling,
documented with its limits.

## Documentation

| Document | Contents |
|---|---|
| [`docs/math_spec.md`](docs/math_spec.md) | every formula, its derivation, its validity range, its citation |
| [`docs/validation.md`](docs/validation.md) | what is verified, what is not, measured numbers |
| [`docs/roadmap.md`](docs/roadmap.md) | v0.2 → v0.6, and what is explicitly out of scope |
| [`docs/hardware.md`](docs/hardware.md) | the hardware-in-the-loop bench: parts, frequencies, latency budget |
| [`tools/bellhop_compare/`](tools/bellhop_compare/README.md) | how the Bellhop cross-validation is set up and what it does not cover |
| [`EXPORT_NOTICE.md`](EXPORT_NOTICE.md) | publication and export-control position |

## References

Jensen, Kuperman, Porter & Schmidt, *Computational Ocean Acoustics*, 2nd ed.,
Springer 2011 · Urick, *Principles of Underwater Sound*, 3rd ed., 1983 ·
Medwin, *JASA* 58(6), 1975 · Mackenzie, *JASA* 70(3), 1981 · Chen & Millero,
*JASA* 62(5), 1977 · Leroy & Parthiot, *JASA* 103(3), 1998 · Munk, *JASA* 55(2),
1974 · Porter, *The BELLHOP Manual and User's Guide*, 2011.

## License

Apache-2.0 — see [`LICENSE`](LICENSE). Apache rather than MIT for the explicit
patent grant, which matters to anyone evaluating this in an industrial context.
