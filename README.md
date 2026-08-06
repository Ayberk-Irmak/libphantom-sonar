# libphantom-sonar

**Real-time ocean acoustics and sonar countermeasure *simulation* library.**
Zero dependencies, zero heap allocation, C++20.

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![Status](https://img.shields.io/badge/status-v0.5%20boundary%20losses-orange.svg)
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

**A ping analyser, because an echo synthesiser with no input is not a module.**
Detecting and characterising the incoming ping is the harder half of a
countermeasure, and the half with the latency budget — so it shipped first. A
bank of matched filters (CW / LFM / HFM), CA-CFAR detection, and a Pulse
Descriptor Word per arrival: time, type, centre frequency, bandwidth, chirp
rate, amplitude, SNR.

```
Filter bank
  [0] CW         12000 -  12000 Hz,  20.0 ms, TB =      0
  [1] LFM-up      8000 -  20000 Hz,  20.0 ms, TB =    240
  ...
Detections
  ToA (s)    type           amp    SNR dB    fc (Hz)   BW (Hz)
  0.20000    LFM-up       0.977      29.5      14000     12000   <- ToA error +1.3 us
  0.54998    CW           0.595      26.3      12000         0   <- ToA error -15.5 us
  0.95000    HFM          0.326      19.0      14000     12000   <- ToA error +3.9 us
```

The arrival-time estimator is checked against the **Cramér-Rao bound** — not
against another implementation, but against the best variance any unbiased
estimator could achieve. It sits at **0.93–1.02× the bound across 20 dB of SNR**.

Getting that comparison right required care: a magnitude detector is bounded by
the RMS bandwidth about the centre frequency, not by `f_rms` about zero. For an
8-20 kHz sweep those differ by 4.16×, so checking an envelope detector against
the coherent bound would make an efficient estimator look four times worse than
it is. Both are computed; [`docs/validation.md §4`](docs/validation.md) shows
which applies and why.

**Underwater Doppler is a time scaling, not a frequency shift.** With `c ≈ 1500
m/s`, `v/c` is ~500× larger than for an airborne radar at the same speed, and
that changes the engineering. Measured peak loss over a 60 ms 8-20 kHz sweep:

| closing speed | CW | LFM | HFM |
|---|---|---|---|
| 5 m/s | — | −5.7 dB | **−0.05 dB** |
| 15 m/s | — | −9.9 dB | **−0.14 dB** |
| 30 m/s | — | −13.1 dB | **−0.30 dB** |

which is why HFM dominates underwater: a time-scaled HFM is a *delayed* HFM.

It also breaks the textbook LFM range-Doppler formula. The narrowband result
uses the sweep's centre frequency; the wideband one uses its **end**:

```
dt = -(v/c) · f_end / mu
```

Verified on upsweeps, downsweeps and narrow sweeps — the narrowband formula is
off by 13 samples on the 8-16 kHz case where the wideband form agrees within 3.
Derivation in [`docs/math_spec.md §6.2`](docs/math_spec.md).

**Boundaries cost what they cost.** Bottom reflection is the Rayleigh
fluid-fluid coefficient, pinned by three independent closed-form limits — the
critical angle from the speed ratio, the impedance ratio at normal incidence,
and unity below critical for a lossless bottom. The critical angle,

```
θ_c = arccos(c_water / c_sediment)     — 24.62° for medium sand
```

is the number that decides which paths survive to long range in shallow water.
Sediment attenuation enters as a complex speed, so `|R| < 1` at *every* angle:
without it, sub-critical rays are trapped forever and shallow-water range is
unbounded, which is not what the ocean does.

The complex square root has two branches and only one is physical. Rather than
reason about sign conventions, the implementation takes the branch that
conserves energy — a 9100-combination sweep over speed, density, attenuation and
angle gives `max |R| = 1.000000000000`.

What that does to a 200 m duct at 3 km, 5 kHz, 8 m/s wind, sand bottom:

| surface | bottom | grazing | TL before | TL now |
|---|---|---|---|---|
| 0 | 0 | — | 71.80 | **71.80** |
| 1 | 1 | 6.83° | 70.67 | 82.85 |
| 2 | 2 | 13.93° | 70.95 | 132.73 |
| 4 | 3 | 24.69° | 71.62 | **199.51** |

Paths that all sat within 1 dB of each other now span 128 dB. The direct path is
unchanged, as it must be.

**Surface loss is capped, and the cap is the honest part.** At 0.5 m seas,
10 kHz and 20° grazing the coherent-scattering formula gives a **700 dB** loss.
That is arithmetic, not physics: the scattered energy goes into a diffuse field
a ray model does not carry, so the path is not 700 dB down — its *specular* part
is. The library reports a capped 30 dB and says why, rather than deleting a path
that is still in the water.

**Eigenrays, and spreading loss that reduces exactly to the closed form.** An
eigenray is a launch angle that reaches a given receiver. The search is exact
rather than interpolated: the tracer already stops on a range budget to machine
precision, so setting that budget to the receiver range makes the ray's final
point *the arrival*.

Geometric spreading then comes from the ray tube — the derivative of arrival
depth with respect to launch angle:

```
TL = -10 log10[ c_rcv cos(θ₀) / (c_src · r · cos(θ_rcv) · |dz/dθ₀|) ]
```

In isovelocity water that must reduce to `20 log10(R)`. It does, to **0.0000 dB**
across nine geometries from 1 to 20 km — with the Jacobian obtained by finite
difference *through the tracer*, so the agreement exercises the search, the
arc-length accumulator and the formula at once.

| range | offset | ray-tube TL | 20 log10(R) | difference |
|---|---|---|---|---|
| 1 km | 500 m | 60.9691 | 60.9691 | 0.0000 |
| 5 km | 2500 m | 74.9485 | 74.9485 | 0.0000 |
| 20 km | 3000 m | 86.1172 | 86.1172 | 0.0000 |

**Multipath the ocean chooses, not the caller.** Feed the eigenrays to the echo
synthesiser and the reply carries the arrival structure of the water it is
transmitted into:

```
Multipath from the ray tracer
      launch     t (ms)      srf      btm    TL (dB)
     -6.498d   1656.190        1        0      77.93
     -5.436d   1655.550        0        0      73.56
      8.952d   1690.677        0        1      72.00
  -> 3 echoes spanning 70.25 ms two-way
```

**Caustics are flagged, not answered.** Where neighbouring rays cross, the tube
collapses and ray theory predicts infinite intensity. That is a failure of the
method, not a property of the ocean, so the library refuses to report a level
rather than returning a number that looks like one.

**Doppler is where underwater differs from radar, and the bank sizing shows it.**
A zero-Doppler template detects a moving target but does not match it. How many
Doppler bins a bank needs is a property of the waveform — derived per family and
checked against measurement. Over ±20 m/s at 1 dB straddling loss, for the same
8-20 kHz 20 ms sweep:

| waveform | bins needed |
|---|---|
| CW | 14 |
| LFM | 5 |
| **HFM** | **2** |

Each bin is another correlation per block and another 128 kB of replica
spectrum, so that ratio is the design argument for hyperbolic sweeps stated as a
cost.

**Echo synthesis, verified by closed loop.** A ping is detected, an echo is
synthesised from the resulting descriptor, and the echo is fed back through the
analyser. A shared sign error between the two halves then shows up as a round
trip that does not close, instead of as two tests that agree and are both wrong.

```
Reply
   range (m)    TS (dB)    v (m/s) delay (ms)
           8         -2          3      10.62
          20         -5          0      26.55
          32         -8         -4      42.48
          44        -11          6      58.40

As received back
      ToA (ms) type            amp    SNR dB    v (m/s)
        12.659 LFM-up        0.525      21.0        4.0   <- ghost at 8 m
        44.472 LFM-up        0.258      17.1      -12.0   <- ghost at 32 m
        60.477 LFM-up        0.227      16.7       12.0   <- ghost at 44 m
```

Delay and Doppler are both **two-way**: `dt = 2dr/c` and
`alpha = (c+v)/(c-v)`, the exact form rather than `1 + 2v/c` (they differ by
0.4% at 30 m/s, which is four samples across a 20 ms pulse). The factor of two
in each is the likeliest bug in an echo path, so the round trip asserts the
reading is closer to the two-way scale than to the ghost's own velocity.

**The two halves cross-check each other.** In the closed-loop example an 8 m/s
ping lands in the 12 m/s bin; the 4 m/s residual produces an arrival-time error
of **+0.0885 ms**, and the wideband coupling formula predicts
`dt = -(v/c)·f_end/mu = +0.0885 ms`. Two results derived independently, agreeing
to four decimals.

**On cancelling the echo rather than faking it.** Every countermeasure
specification asks for anti-phase cancellation. This library quantifies it
instead of implementing it:

| cancellation | timing budget | path budget at 1500 m/s |
|---|---|---|
| −10 dB | 4.2 µs | 6.3 mm |
| −20 dB | 1.3 µs | **2.0 mm** |

Bandwidth barely matters (widening a 12 kHz-centred band from 0 to 12 kHz costs
0.3 dB); timing is everything. Past a quarter period the canceller *adds* up to
6 dB. Add distributed hull scattering — one projector cannot match phase at more
than one bearing — and the conclusion is that **generating false targets is the
achievable countermeasure; cancelling the real one is not.**

**Verified against closed forms, not against baselines.** A recorded baseline
tells you the code did not change. These tell you it is right:

| Property | Result |
|---|---|
| Snell invariant over 29 rays × 100 km | 1.7e-16 relative drift |
| Ray path vs analytic circle, 200 layer crossings | 1.9e-16 relative radial error |
| Turning depth vs `c(z) = c₀/cos θ₀` | < 1e-6 m |
| Range and travel-time budgets | exact to 1e-8 m / 1e-12 s |
| Mackenzie vs Chen-Millero, common validity box | 0.53 m/s max disagreement |
| Ray-tube spreading vs 20 log10(R), 9 geometries | 0.0000 dB |
| Thorp absorption vs the published formula | within 10% (what the fit is worth) |
| Rayleigh \|R\| ≤ 1 over 9100 sediment/angle combinations | max 1.000000000000 |
| Critical angle, normal incidence, sub-critical unity | exact to 1e-10 |
| FFT vs a direct O(N²) DFT sharing no code | 2.3e-16 relative |
| FFT correlation vs direct O(N·L) correlation | 1.1e-15 relative |
| Matched filter peak, width, coherent gain | `A·E/2`, `1/B`, `L/2` — all matched |
| Waveform classification, 4 types at −4.4 dB SNR | 100/100 |

15642 checks. Clean under GCC 15 and Clang 21 with `-Wconversion -Wsign-conversion
-Wold-style-cast -Wdouble-promotion -Werror`, and under ASan + UBSan. Passes in
both `double` and `float` builds. Zero allocation is proven by `nm` over the
built archive, not asserted: the library references only libm and `memset`.

## Measured performance

13th Gen Core i7-13620H, `-O3`, single thread:

| Operation | Time |
|---|---|
| `mackenzie()` | 1.5 ns |
| `chen_millero()` (UNESCO) | 4.0 ns |
| `speed_at()` — 501-layer binary search | 16.6 ns |
| **`trace_ray()` — per arc step** | **52 ns** |
| `trace_ray()` — full 100 km ray (1132 arcs) | 40.8 µs |
| `analyze_sofar()` — 501 points | 596 ns |
| `fft_forward()` — 8192 points | 307 µs |
| `analyze_block()` — 4-template bank, 8192-point block | 2.26 ms |

The ping analyser runs **29× real time on one core** at 96 kHz. Detection
latency is 65 ms and is set by the block length, not the arithmetic — a pulse is
only reported once the block containing it is complete.

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
./build/ping_intercept data    # streaming ping detection, scored against truth
./build/countermeasure_loop    # intercept -> classify -> reply -> multipath
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
- **Cross-template ghosts.** A bank whose waveforms share a band reports one
  arrival more than once — an LFM arrival lights the HFM template ~10.5 dB down
  at a shifted lag. Correct behaviour for a matched filter bank; suppressing it
  needs association logic, which is v0.4.
- **A bank only resolves what it spans.** An echo stretched to 28 ms is detected
  by a 20 ms template and mis-reported in both time and type — shown
  deliberately in `countermeasure_loop`.
- **Doppler is quantised and the range bias is measured, not corrected.**
- **No bearing.** Single channel; bearing needs an array.
- **Ray theory, with its caustics.** Levels near a caustic are flagged rather
  than computed; a correct treatment needs a wave solution through them.
- **No diffuse field.** Surface scattering removes the specular component and
  the loss is capped to acknowledge the energy went somewhere, but where it went
  is not modelled. Reverberation is v0.6.
- **Plane-wave, flat-interface reflection.** No beam displacement near the
  critical angle, no sediment layering, no shear in the bottom.
- **Monostatic only.** Multipath echoes assume the return travels the way the
  ping came; a bistatic geometry needs two path sets.
- **Target strength is still a dB figure**, not the output of a scattering
  model.
- **A Doppler bank is megabytes** — 8 MB at 64 templates and an 8192-point
  transform — and must not go on the stack.
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
