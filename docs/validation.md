# Validation status

What is actually verified, what is only self-consistent, and what is not
verified at all. Stated plainly, because in this domain an unverified claim is
worse than no claim.

Everything below was produced by `phantom_tests`, `munk_simulation` and
`tools/bellhop_compare/` on x86-64, GCC 15.3 and Clang 21.1, `-O3`,
`Real = double` unless noted.

---

## 1. Verified against closed-form results ✅

The tracer is compared against exact analytic solutions, not against a recorded
baseline. A regression baseline only tells you the code did not change; these
tell you it is right. §2 then checks the whole thing against an independent
implementation.

| Property | Method | Result |
|---|---|---|
| Snell invariant `cos θ / c` | held constant along 29 rays × 100 km through a 501-layer Munk profile | **1.7e-16** relative drift |
| Ray path is a circular arc | every point compared against `R = 1/(ξ\|g\|)` about the analytic centre, 200 layer crossings | **1.9e-16** relative radial error |
| Turning depth | compared against `c(z_turn) = c₀ / cos θ₀` | **< 1e-6 m** |
| Layer-refinement invariance | same linear profile sampled at 2 vs 1001 points, traced 40 km | depth differs by **4.6e-12 m** |
| Isovelocity propagation | straight-line geometry and `t = path/c` | exact to **1e-12** |
| Boundary reflection count | 45° sawtooth, analytically 5 surface + 5 bottom hits over 10 km | exact |
| Range budget | `final_range_m` vs configured limit, 25 angle/range combinations | exact to **1e-8 m** |
| Travel-time budget | arc inverted in closed form (§2.5 of `math_spec.md`) | exact to **1e-12 s** |
| SOFAR cone vs tracer | rays at `θ_max ± 0.5°` must respectively stay trapped / escape | consistent |

## 2. Cross-validated against Bellhop ✅

Everything in §1 is internal: closed forms and invariants. Those catch
implementation errors but not a shared misunderstanding of the physics. This
section is the external check.

Reference: **Bellhop**, Acoustics Toolbox `at_2026_7`, built with gfortran 15.3
on x86-64 Linux. Reproduce with `tools/bellhop_compare/`.

Both codes read the **same `.env` file** — `phantom_trace` parses Bellhop's own
input format rather than a transcribed copy, so a transcription bug cannot look
like agreement. Both are configured for **C-linear** SSP interpolation, which is
what libphantom-sonar implements; `phantom_trace` refuses to run on a PCHIP or
spline file, because that would measure the interpolation scheme rather than the
tracer.

### The headline: turning depths

Comparing two polylines means interpolating both, which injects error belonging
to neither code. The turning depth avoids that: one number per ray, solved for
exactly here (`c(z_turn) = 1/ξ`) and approached by Bellhop as its step shrinks.

| Bellhop step | turning-depth error | path RMS |
|---|---|---|
| 500 m | 0.084 m | 0.230 m |
| 200 m | 0.014 m | 0.037 m |
| 100 m | 0.0024 m | 0.012 m |
| 50 m | 0.0017 m | 0.0071 m |
| 20 m | 0.00036 m | 0.0067 m |
| **10 m** | **0.00014 m** | 0.0068 m |

**0.14 mm over a 101 km path** at Bellhop's finest step. The direction matters as
much as the magnitude: Bellhop converges *toward* this library's answer, which is
what should happen if the arc solution is the exact limit of a stepping
integrator.

### Full-path comparison

Canonical Munk profile, source at 1000 m, 101 km, 41 rays, at Bellhop's default
(auto) step size — i.e. what a user gets out of the box:

| Case | rays | boundary hits | mean RMS | worst RMS | max deviation |
|---|---|---|---|---|---|
| Trapped (±13°, refraction only) | 41 | 0 | 0.274 m | 0.687 m | 2.05 m |
| Bouncing (±20°) | 41 | 62 | 0.635 m | 2.61 m | 5.79 m |

### Why those RMS figures are a ceiling, not a measurement

The path-RMS column stops falling at ~0.0068 m no matter how fine Bellhop's step
gets. That floor is **ours**: the tracer emits one point per layer crossing, so
on a 27-point profile a ray is sampled every few kilometres, and the straight
line between those points cuts the corner of each curved arc.

| output sampling | apparent RMS |
|---|---|
| ×8 | 1.055 m |
| ×32 | 0.197 m |
| ×128 | 0.035 m |
| ×512 | 0.0068 m |

Exactly halving per doubling — first order, the signature of chording across the
parabolic nose at a turning point. Refining the sampling does not change the
modelled ocean (`c(z)` is linear within a layer, so inserted points lie on the
same line; `ray_is_invariant_under_layer_refinement` pins that the answer does
not move).

So the true disagreement is **below the smallest value this method can
resolve**. Quoting the RMS without this table would present a limitation of the
measurement as a property of the code.

### What this does and does not cover

- Ray **geometry** only: paths, turning depths, travel times. Not transmission
  loss, beam amplitudes or caustics — v0.1 does not compute those.
- Range-independent profiles only. Bellhop supports range-dependent SSPs and
  bathymetry; this library does not yet.

## 3. Cross-validated between independent equations ✅

Three sound speed equations fitted to different data sets must agree inside
their common validity box (5–25 °C, 34–36 PSU, 0–800 m). A single mistyped
coefficient breaks the agreement immediately.

| Comparison | Max disagreement |
|---|---|
| Mackenzie vs Medwin | **0.295 m/s** |
| Mackenzie vs Chen-Millero (UNESCO) | **0.530 m/s** |

Both are within the published mutual scatter of these formulations. Additional
checks: monotonicity in `T`, `S` and `z` across the full validity range; the
textbook reference point `c(15 °C, 35 PSU, 0 m) ≈ 1507 m/s` reproduced by all
three; depth→pressure ≈ 101 bar at 1000 m.

## 4. Build and runtime hygiene ✅

| Configuration | Result |
|---|---|
| GCC 15.3, `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wold-style-cast -Wcast-qual -Wdouble-promotion -Werror` | clean |
| Clang 21.1, same set | clean |
| AddressSanitizer + UndefinedBehaviorSanitizer | clean, 791 checks |
| `Real = float` (`-DPHANTOM_REAL_FLOAT=ON`) | 791 checks pass |
| Heap allocations during execution | **zero** — no `new`, `malloc`, `std::vector` or `std::string` in `src/` or `include/` |

### float vs double

Same algorithm, different achievable precision. Over a 40 km trace through 1001
layers:

| | `double` | `float` |
|---|---|---|
| Snell drift | 1.7e-16 | 7.1e-08 |
| Layer-refinement depth delta | 4.6e-12 m | 4.0e-03 m |
| Travel time delta | ~1e-13 s | ~1.9e-04 s |

**4 mm of depth error over 40 km** is acceptable for most MCU work. It is not
acceptable if you are differencing travel times for ranging. Use `double` unless
your FPU forces otherwise.

## 5. Performance (measured, not claimed)

13th Gen Intel Core i7-13620H, `-O3`, single thread:

| Operation | Time |
|---|---|
| `mackenzie()` | 1.5 ns |
| `chen_millero()` (UNESCO) | 4.0 ns |
| `speed_at()` — 501-layer binary search | 16.6 ns |
| **`trace_ray()` — per arc step** | **52 ns** (median) |
| `trace_ray()` — full 100 km ray, 1132 arcs | 40.8 µs |
| `analyze_sofar()` — 501 points | 543 ns |

### On the "< 500 ns" requirement

The original specification asked for `< 500 ns per processing cycle` without
defining the cycle. That number is not meaningful for a whole ray trace — a
100 km trace crosses ~1100 layers and cannot fit in a sub-microsecond budget on
any hardware.

The library therefore states its budget **per arc step**, which is the unit that
actually constrains a control loop, and measures **52 ns**. Re-run
`phantom_bench` on your own target before quoting any of this as a real-time
guarantee; these numbers are not a substitute for WCET analysis on the target
silicon.

## 6. Self-consistent but NOT independently verified ⚠️

- **Chen-Millero against the published UNESCO check tables.** The
  implementation is cross-validated against two independent equations (§3),
  which would catch any coefficient typo, but it has not been checked against
  the official UNESCO validation table.
- **Range-dependent propagation.** The Bellhop comparison (§2) covers
  range-independent profiles only, because that is all this library models.
  Nothing here is evidence about range-dependent behaviour.
- **Ray amplitudes.** §2 validates geometry. Transmission loss, beam amplitudes
  and caustics are not computed, so nothing about them is validated.
- **Real measured T/S profiles.** Everything so far uses analytic profiles
  (Munk, linear, isovelocity). Feeding World Ocean Atlas / Argo data is v0.2.

## 7. Known limitations — not bugs, scope ⚠️

- **Shadow fraction depends on ray density.** Gaps between adjacent rays read as
  shadow. Measured convergence for the 100 km Munk case on a 500×250 grid:

  | rays | shadow % |
  |---|---|
  | 61 | 70.1 |
  | 121 | 53.5 |
  | 241 | 37.6 |
  | 481 | 29.2 |
  | 721 | 28.2 |
  | 1441 | 27.5 |
  | 2881 | 26.5 |
  | 5761 | 26.0 |

  Still drifting by ~0.5 points per doubling at 5761 rays. **Treat any reported
  shadow fraction as an upper bound**, and always report the ray count with it.
  A proper treatment needs ray-tube area / caustic handling, not more rays.

- **Range-independent ocean only.** The profile is a function of depth alone. No
  range-dependent bathymetry, no fronts, no eddies.
- **Geometric acoustics only.** No transmission loss, no absorption, no
  amplitude, no diffraction, no caustic corrections. Ray *paths*, not ray
  *intensities*.
- **No bottom interaction model.** Reflection is specular and lossless, or the
  ray is absorbed. Real sediment reflection is grazing-angle and
  frequency-dependent.
- **2-D (range–depth) only.** No azimuthal coupling or out-of-plane refraction.

## 8. Not implemented in v0.1

`PhantomEchoSynthesizer` and `BioMimeticCommEngine` are not in this release. See
`docs/roadmap.md` for what they will require, including the ping-analysis stage
that has to exist before echo synthesis is meaningful.

---

## Reproducing

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/phantom_tests          # 791 checks
./build/phantom_bench          # timing table above
./build/munk_simulation data   # CSV + console report
python3 tools/plot_rays.py data
```

Bellhop cross-validation (§2):

```bash
./tools/bellhop_compare/setup_bellhop.sh
python3 tools/bellhop_compare/compare.py \
        --bellhop build/acoustics-toolbox/at/Bellhop/bellhop.exe --check
```

Convergence sweep from §7:

```bash
for n in 61 121 241 481 721 1441 2881 5761; do
  ./build/munk_simulation /tmp "$n" | grep 'shadow fraction'
done
```
