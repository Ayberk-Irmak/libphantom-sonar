# Bellhop cross-validation

Independent verification of libphantom-sonar's ray tracer against
[Bellhop](https://oalib-acoustics.org/models-and-software/acoustics-toolbox/),
the reference ray/beam code of the Acoustics Toolbox and the de facto standard
in ocean acoustics.

Everything the library's own test suite checks is *internal*: closed forms and
invariants. Those catch implementation errors but cannot catch a shared
misunderstanding of the physics. This directory closes that gap.

## Running it

```bash
./tools/bellhop_compare/setup_bellhop.sh        # downloads + builds Bellhop
cmake -S . -B build && cmake --build build
python3 tools/bellhop_compare/compare.py \
        --bellhop build/acoustics-toolbox/at/Bellhop/bellhop.exe
```

Takes about five seconds once Bellhop is built. Add `--check` to make it exit
non-zero on disagreement (this is what CI runs), and `--json out.json` for a
machine-readable summary.

## Design

**Both codes read the same file.** `phantom_trace` parses Bellhop's own `.env`
format rather than working from a transcribed copy of the profile. A
transcription bug would look exactly like agreement, so the possibility is
removed rather than assumed away.

**The interpolation must match.** Bellhop offers C-linear, N²-linear, PCHIP,
spline and analytic SSP reconstruction. libphantom-sonar is piecewise linear in
`c(z)`, which is Bellhop's `'C'`. Comparing against PCHIP would measure the
difference between two interpolation schemes and say nothing about either
tracer, so `phantom_trace` **refuses to run** on a non-`'C'` file. The canonical
`MunkB_ray.env` shipped with the toolbox uses PCHIP, which is why `make_env.py`
regenerates the same profile with `'CVF'`.

**The beam count is explicit.** With `NBEAMS = 0` Bellhop picks its own beam
count, and the two runs would no longer be launching the same rays.

## What is measured

| Case | Setup | What it isolates |
|---|---|---|
| **A** | 41 rays, ±13° from 1000 m — inside the trapping cone | pure refraction, no boundary contact |
| **B** | 41 rays, ±20° — escapes the channel | surface and sea-floor reflection bookkeeping |
| **C** | Bellhop step from 500 m down to 10 m | Bellhop's discretisation error |
| **D** | libphantom-sonar output sampling from ×8 to ×512 | *this comparison's* own error |

### Case C is the real result

Comparing two polylines requires interpolating both, which injects error that
belongs to neither code. The **turning depth** avoids that entirely: it is one
number per ray, solved for exactly in libphantom-sonar (`c(z_turn) = 1/ξ`) and
approached by Bellhop as its integrator step shrinks.

| Bellhop step | turning-depth error | path RMS |
|---|---|---|
| 500 m | 0.084 m | 0.230 m |
| 200 m | 0.014 m | 0.037 m |
| 100 m | 0.0024 m | 0.012 m |
| 50 m | 0.0017 m | 0.0071 m |
| 20 m | 0.00036 m | 0.0067 m |
| **10 m** | **0.00014 m** | 0.0068 m |

Bellhop converges toward libphantom-sonar's answer, reaching **0.14 mm over a
101 km path**. That is the direction that matters: the analytic arc solution is
the limit a stepping integrator approaches, so this is agreement *and* evidence
that the arc solution is the exact one.

### Case D is why the RMS numbers are honest

The path-RMS column above stops falling at ~0.0068 m. That floor is **ours**,
not Bellhop's. libphantom-sonar emits one point per layer crossing, so on a
27-point profile a ray is sampled every few kilometres, and the straight line
drawn between those points cuts the corner of each curved arc:

| refine | apparent RMS |
|---|---|
| ×8 | 1.055 m |
| ×32 | 0.197 m |
| ×128 | 0.035 m |
| ×512 | 0.0068 m |

Halving exactly with each doubling — first-order, the signature of chording
across the parabolic nose at each turning point. `--refine N` splits each
tabulated layer into `N` sublayers, which **does not change the modelled ocean**
(within a layer `c(z)` is a straight line, so inserting points on that line
reproduces it exactly; `ray_is_invariant_under_layer_refinement` in the test
suite pins that the traced answer does not move). It only raises output density.

So the RMS figures are a **ceiling on the disagreement, not a measurement of
it** — the true difference is below the smallest number the method can resolve.
Reporting them without Case D would be quoting a limitation of the measurement
as if it were a property of the code.

## Files

| File | Role |
|---|---|
| `setup_bellhop.sh` | downloads and builds the Acoustics Toolbox |
| `make_env.py` | writes the `.env` cases from the canonical Munk table |
| `phantom_trace.cpp` | parses a `.env`, traces it with this library, writes CSV |
| `compare.py` | runs both, computes the statistics, draws the figure |

## Licensing

The Acoustics Toolbox is **GPL-3.0**. libphantom-sonar is Apache-2.0 and does
not link against it, vendor any of its code, or redistribute it. `setup_bellhop.sh`
fetches it onto your machine and `compare.py` invokes `bellhop.exe` as a separate
process exchanging text files. No GPL obligation attaches to this library.

## Caveats

- Compares **ray geometry only** — paths, turning depths, travel times. Not
  transmission loss, beam amplitudes or caustics, none of which v0.1 computes.
- Range-independent profiles only. Bellhop supports range-dependent SSPs and
  bathymetry; this library does not yet, so those cases are not exercised.
- Verified against Acoustics Toolbox `at_2026_7` on x86-64 Linux with
  gfortran 15.3. Bellhop is under active development and its defaults change
  between releases; the version is printed in each run's `.prt` file.
