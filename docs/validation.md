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

## 4. Pulse analysis chain — verified against closed forms ✅

Same standard as the ray tracer: every check is against a closed form, an
invariant, or a theoretical bound, never against a recorded output.

### FFT

Hand-written, because "zero dependencies" would otherwise be a half-truth. It is
also exactly the kind of code that looks right while being wrong, so:

| Check | Result |
|---|---|
| vs a direct O(N²) DFT sharing no code, N = 256 | **2.3e-16** max relative error |
| round trip `ifft(fft(x))`, N = 1024 | **4.4e-16** max absolute error |
| Parseval, N = 512 | exact to **1e-15** |
| known transforms (impulse, constant, single bin) | exact |
| linearity, and every size from 2 to 4096 | exact |

### Waveforms

| Check | Result |
|---|---|
| `d(phase)/dt` vs the specified `f(t)`, LFM and HFM | **1.9e-11** relative |
| HFM period linear in time (its defining property) | **1e-9** relative |
| HFM → CW as the sweep vanishes | finite, matches the linear form |
| Doppler is a time *scaling*: 10 m/s shortens a 20 ms pulse | 1920 → 1907 samples, matches `T/(1+v/c)` |

### Matched filter

| Property | Theory | Measured |
|---|---|---|
| FFT correlation vs direct O(N·L) correlation | identical | **1.1e-15** relative |
| Peak lag | the true delay | exact, every delay tested |
| Peak magnitude | `A·E/2` | 480.00 for `E` = 960 |
| Compressed width | `1/B` = 83 µs | 83.3 µs (`B` = 12 kHz, from a 10 ms pulse) |
| Coherent gain | `L/2` = 480 | 480 ± 5% |
| Detection at −13.5 dB input SNR | — | 96/100 within one sample |

### Doppler behaviour of the three waveform families

Measured over a 60 ms 8-20 kHz sweep, peak loss versus closing speed:

| closing speed | CW | LFM | HFM |
|---|---|---|---|
| 2.5 m/s | −9.9 dB | −0.45 dB | — |
| 5 m/s | — | −5.7 dB | **−0.05 dB** |
| 15 m/s | — | −9.9 dB | **−0.14 dB** |
| 30 m/s | — | −13.1 dB | **−0.30 dB** |

The CW figure is at its predicted first null, `v = c/(f0·T)`.

**A correction to the textbook LFM range-Doppler formula.** The narrowband
result `dt = -f_doppler/mu` with the shift taken at the *centre* frequency does
not hold underwater, where Doppler is a time scaling rather than a frequency
shift. The wideband form uses the sweep **end** frequency:

```
dt = -(v/c) · f_end / mu
```

Verified on upsweeps, downsweeps and narrow sweeps, agreeing within 3 samples
while the narrowband formula is off by 13 samples on the 8-16 kHz case
(derivation in `docs/math_spec.md §6.2`).

### Arrival-time estimation vs the Cramér-Rao bound

The strongest check in the library: not against another implementation, but
against the best variance *any* unbiased estimator could achieve.

**Which bound applies matters.** A magnitude (envelope) detector is bounded by
the RMS bandwidth about the centre frequency; only a carrier-phase-coherent
receiver is bounded by `f_rms` about zero. For an 8-20 kHz sweep the two differ
by a factor of 4.16, so comparing an envelope detector against the coherent
bound makes an efficient estimator look four times worse than it is.

| noise σ | measured σ_ToA | envelope CRLB | ratio |
|---|---|---|---|
| 3.0 | 4.30 µs | 4.45 µs | **0.97** |
| 1.0 | 1.48 µs | 1.48 µs | **1.00** |
| 0.3 | 0.41 µs | 0.44 µs | **0.93** |

Flat across 20 dB — the signature of an estimator limited by noise rather than
by a systematic interpolation error. The coherent bound over the same runs is
1.07 µs, i.e. 4.2× tighter and not reachable by this receiver.

The Fisher information is computed two independent ways (analytically from the
instantaneous frequency, and by a spectral derivative of the rendered waveform)
and the two agree to 0.06%.

### Detector behaviour

| Check | Result |
|---|---|
| Classification confusion matrix, 4 waveforms × 25 trials at −4.4 dB | **100/100** correct |
| False alarms, noise only, Pfa 1e-3 → 1e-7 | monotone, < 0.05/block at the tightest |
| Two arrivals, 8 dB apart | both found, amplitude ratio recovered to 0.05 |
| Streaming timestamps across blocks | exact |

### A silent failure mode, reproduced and fixed

CA-CFAR estimates noise from cells either side of a guard band. If the guard is
narrower than the target's response, the target leaks into its own training
cells and raises its own threshold — **the stronger the pulse, the higher the
bar it must clear.** Nothing is reported and the detector looks healthy.

A chirp compresses to ~`fs/B` cells and survives a small guard; a CW does not
compress at all. Measured with a 32-cell guard against a 1920-sample CW
response:

| waveform | detected with a 32-cell guard | with `suggested_cfar_guard()` |
|---|---|---|
| CW | **0/20** | 20/20 |
| LFM | 20/20 | 20/20 |

Three of four waveforms still working is exactly what makes this bug hide.
`analyzer_cfar_guard_must_clear_the_response` pins both halves.

## 5. Doppler bank and echo synthesis — verified ✅

### Doppler bin spacing against measurement

The spacing formula is only useful if it predicts the real loss. Checked
directly against the matched filter: at the `|delta|` the formula reports, the
peak must actually have dropped by about that much.

| waveform | requested | δ | v (m/s) | measured loss |
|---|---|---|---|---|
| LFM | 1 dB | 4.27e-3 | 6.4 | 1.08 dB |
| LFM | 3 dB | 7.39e-3 | 11.1 | 3.15 dB |
| HFM | 1 dB | 0.109 | 163 | 1.55 dB |
| HFM | 3 dB | 0.292 | 438 | 4.12 dB |

Bins needed over ±20 m/s at 1 dB straddling loss, same 8-20 kHz 20 ms sweep:
**CW 14, LFM 5, HFM 2.** That ratio is the design argument for hyperbolic
sweeps, in a number.

### Velocity estimation

A bank over ±12 m/s (3 bins, 12 m/s spacing) recovers radial velocity to within
5.0 m/s, against a half-bin floor of 6.0 m/s — i.e. the estimator is
quantisation-limited, as it should be. The bank also recovers most of the peak a
zero-Doppler template loses: at 10 m/s, −2.54 dB becomes −0.18 dB.

### Echo synthesis — closed loop

Every echo test detects a ping, synthesises an echo from the resulting
descriptor, and feeds it back through the analyser. A shared sign error between
the two halves shows up as a round trip that does not close, rather than as two
tests that agree with each other and are both wrong.

| Quantity | Result |
|---|---|
| Delay, 20 m and 45 m apparent range | 26.656 / 59.990 ms vs 26.667 / 60.000 |
| Target strength, 0 and −6 dB | amplitude 1.000 and 0.501 |
| Two-way Doppler, ghost at 5 m/s | bank reads 10.0 m/s (two-way scale = 10.03) |
| Ghost swarm, 3 targets at 0/−4/−8 dB | all placed exactly; amplitudes 1.000 / 0.631 / 0.398 |
| Extended target, `length_scale` 2 | 1920 → 3840 samples |

The Doppler round trip asserts the reading is closer to the **two-way** scale
than to the ghost's own velocity — the factor of two is the likeliest bug in the
whole echo path, so it is tested for directly.

### Cross-check between the two halves

In `countermeasure_loop`, an 8 m/s ping lands in the 12 m/s bin, and the
resulting 4 m/s residual produces an arrival-time error of **+0.0885 ms**. The
wideband coupling of §4 predicts `dt = -(v/c)·f_end/mu = +0.0885 ms` — the
Doppler bank and the range-Doppler formula, derived independently, agree to four
decimal places.

### Anti-phase cancellation: measured, not claimed

| cancellation | timing budget | path budget |
|---|---|---|
| −6 dB | 6.72 µs | 10.08 mm |
| −10 dB | 4.21 µs | 6.32 mm |
| −20 dB | 1.33 µs | **1.99 mm** |
| −30 dB | 0.42 µs | 0.63 mm |

Bandwidth is nearly irrelevant: with a 2 µs error at 12 kHz, widening the band
from 0 to 12 kHz moves the residual by 0.34 dB. The limit is timing. Past a
quarter period the canceller **adds** up to 6 dB.

Not modelled, and both make it worse: distributed hull scattering (one projector
cannot match phase at more than one bearing) and the latency floor (a canceller
cannot invert a sample it has not received). **Generating false targets is the
achievable countermeasure; cancelling the real one is not.**

## 6. Eigenrays, transmission loss and multipath — verified ✅

### Spreading loss against the closed form

The load-bearing check of the whole module. In isovelocity water with no
boundaries the eigenray is a straight line and the ray-tube formula must reduce
*exactly* to `20 log10(R)`. The Jacobian is obtained by finite difference
through the tracer, so agreement exercises the search, the arc-length
accumulator and the formula at once.

| range | depth offset | ray-tube TL | 20 log10(R) | difference |
|---|---|---|---|---|
| 1 km | 0 / 200 / 500 m | 60.0000 / 60.1703 / 60.9691 | same | **0.0000 dB** |
| 5 km | 0 / 1000 / 2500 m | 73.9794 / 74.1497 / 74.9485 | same | **0.0000 dB** |
| 20 km | 0 / 1000 / 3000 m | 86.0206 / 86.0314 / 86.1172 | same | **0.0000 dB** |

Worst departure across all nine geometries: **0.00000 dB** in double,
0.0009 dB in float.

### Eigenray geometry

Isovelocity, source 1000 m, receiver 1600 m at 5 km:

| quantity | traced | geometric |
|---|---|---|
| launch angle | 6.842773° | 6.842773° |
| path length | 5035.871 m | 5035.871 m (slant) |
| travel time | matched to 1e-7 relative | slant/1500 |
| Jacobian dz/dθ₀ | 5072.000 m/rad | 5072.000 (r/cos²θ₀) |

### Multipath structure

A 200 m duct at 3 km with a ±25° fan yields **14 eigenrays**, ordered in launch
angle, each with distinct bounce counts from 0/0 up to 4/3. Every path is at
least as long as the straight line, every bounced path strictly longer, and
travel times run 1993–2191 ms — a 200 ms spread that is real, not numerical.

### Absorption

Thorp (1967), against the published formula evaluated independently:

| f | computed | published | over 3 km |
|---|---|---|---|
| 1 kHz | 0.069 | 0.069 | 0.21 dB |
| 10 kHz | 1.187 | 1.19 | 3.6 dB |
| 100 kHz | 34.07 | 34.1 | 102 dB |

Monotone across 1–30 kHz in 500 Hz steps. Tolerances are set at ~10%, which is
what Thorp is: a fit at about 4 °C with no temperature or depth dependence.

### Coupling: multipath echoes from traced paths

Eigenrays → `EchoSpec`s → rendered waveform → energy present at every predicted
delay. Delays are referenced to the first arrival (one is exactly zero, none
negative), levels to the strongest path, and the requested target strength lands
on the strongest path exactly.

**The two halves cross-check.** In `countermeasure_loop` an 8 m/s ping falls in
the 12 m/s Doppler bin, and the 4 m/s residual produces an arrival-time error of
+0.0885 ms — exactly what the wideband range-Doppler formula of §4 predicts. CI
asserts both numbers.

### A precision limit worth knowing about

`depth_tolerance_m` **discards** roots it cannot polish, so setting it tighter
than the build can achieve removes paths silently rather than degrading them.
Measured on the 200 m duct at 3 km, 14 paths present:

| tolerance | double | float |
|---|---|---|
| 0.001 m | 14 | **0** |
| 0.010 m | 14 | **3** |
| 0.050 m | 14 | 11 |
| 0.100 m | 14 | 14 |

The default now follows the build (0.01 m double, 0.1 m float). Single precision
cannot resolve the receiver better than ~0.1 m at this range because the
tracer's own depth error there is larger.

### Caustics

Where the ray tube collapses, ray theory predicts infinite intensity. The search
flags such paths and leaves the level at zero rather than reporting a number.
Verified both directly (`spreading_loss_db` with a zero Jacobian returns 0) and
through the search in the Munk channel.

## 7. Boundary reflection — verified ✅

Three independent closed-form limits pin the Rayleigh coefficient. Any one could
be reproduced by a wrong expression; all three together could not.

| limit | check | result |
|---|---|---|
| Critical angle | `arccos(c1/c2)` for 1550 / 1650 / 1800 m/s sediment | 14.593° / 24.620° / 33.557°, exact |
| Normal incidence | `(Z2−Z1)/(Z2+Z1)` over 9 speed × density combinations | exact to 1e-10 |
| Below critical, lossless | `\|R\| = 1` at 5 fractions of the critical angle | exact to 1e-12 |

### Energy conservation, and the branch choice it settles

The complex square root has two branches and only one gives a wave that decays
into the sediment. The implementation picks the one that conserves energy,
because the other yields `|R| > 1`. Swept over 5 speeds × 4 densities ×
5 attenuations × 91 angles = **9100 combinations**:

```
max |R| = 1.000000000000
```

### Attenuation leaks below critical

With a lossless bottom, sub-critical rays are trapped forever and shallow-water
range is unbounded. Measured at half the critical angle:

| attenuation | \|R\| | loss/bounce |
|---|---|---|
| 0.0 dB/λ | 1.000000 | 0.000 dB |
| 0.5 dB/λ | 0.943259 | 0.507 dB |
| 2.0 dB/λ | 0.802190 | 1.915 dB |

Monotone in attenuation. The default sand costs 0.81 dB per bounce there — 8 dB
over ten bounces, which decides whether a path survives.

Bottom loss also grows with grazing angle as it should: 0.19 dB at 2°, 2.04 dB
at the 24.6° critical angle, 9.05 dB at normal incidence. Steep paths die;
grazing paths carry shallow-water range.

### Surface roughness

`|R| = exp(-G²/2)` checked against the closed form across 3 frequencies × 3
angles, agreeing to 1e-9. Flat sea returns exactly 1. Pierson-Moskowitz
`H₁/₃ = 0.0246U²` reproduced; 10 m/s wind gives 0.615 m RMS.

**The cap is deliberate and documented.** At 0.5 m seas, 10 kHz, 20° grazing the
coherent formula gives a **700 dB** loss. That is arithmetic, not physics: the
scattered energy goes into a diffuse field a ray model does not carry, so the
path is not 700 dB down — its *specular* part is. The library reports the capped
30 dB and says why.

### Integration with the eigenray path

Grazing angles at each bounce come from the Snell invariant
(`cos θ = ξ·c` at the boundary), verified against a hand computation: 2 surface
+ 3 bottom bounces at 10° gives 6.7532 dB, matching the sum of the individual
losses exactly.

**The effect this release exists for.** A 200 m duct at 3 km, 5 kHz, 8 m/s wind,
sand bottom:

| surface | bottom | grazing | TL without boundaries | with |
|---|---|---|---|---|
| 0 | 0 | — | 71.80 | **71.80** |
| 1 | 1 | 6.83° | 70.67 | 82.85 |
| 2 | 2 | 13.93° | 70.95 | 132.73 |
| 4 | 3 | 24.69° | 71.62 | **199.51** |

Paths that all sat within 1 dB of each other now span 128 dB. The direct path is
unchanged, as it must be.

## 8. Reverberation — verified ✅

### The range laws

Exact, across three decades of range:

| mechanism | measured decay | theory |
|---|---|---|
| boundary (area ∝ r) | 15.000 / 30.000 / 45.000 / 60.000 dB | 30 log₁₀(r) |
| volume (volume ∝ r²) | 20.000 / 40.000 dB | 20 log₁₀(r) |

### Source level cancels — the result this module exists for

`EL − RL = TS − S_s − 10 log₁₀(A)`, swept over nine combinations of source level
(160 / 200 / 240 dB) and transmission loss (40 / 66 / 90 dB):

| SL | TL | EL | RL | EL−RL |
|---|---|---|---|---|
| 160 | 40 | 90.00 | 84.77 | **5.2288** |
| 200 | 66 | 78.00 | 72.77 | **5.2288** |
| 240 | 90 | 70.00 | 64.77 | **5.2288** |

Identical to four decimals. **Against reverberation, transmit power buys
nothing.** What does: a tenth of the pulse length or a tenth of the beamwidth,
each worth exactly +10.00 dB. Compressing a 20 ms pulse to its 1/B = 83 µs
resolution is worth **+23.8 dB**, which is the time-bandwidth product in dB.

### Lambert and Chapman-Harris

Lambert reproduces `mu + 20 log₁₀(sin θ)` exactly across 5°–90°, returns `mu` at
normal incidence, is monotone over every degree from 2 to 90, and stays finite at
zero grazing.

Chapman-Harris behaves: monotone in wind over 3–40 knots, monotone in grazing
angle over 2°–80°, and its value at exactly 30° matches the wind-and-frequency
term alone — a structural check on the `log₁₀(θ/30)` form. Values fall in the
−24 to −83 dB band.

### Reverberation-limited range

Crossover with ambient noise, verified to actually be a crossover (the two
levels agree to 0.1 dB at the returned range):

| ambient | crossover |
|---|---|
| 40 dB | 24.7 km |
| 60 dB | 5.3 km |
| 80 dB | 1.1 km |

A quiet transmitter into a loud ocean is never reverberation-limited, and the
function returns 0 rather than a number.

### What it does to the detector

A block whose background falls **21.6 dB** end to end. CA-CFAR with a training
window short compared to the decay: **20/20 detections** of a target in the
decayed region, 3 false alarms across 20 empty blocks. A single fixed threshold
at the block mean sits 9.5 dB below the near field and 12.2 dB above the far
field — wrong at both ends, which is precisely the gap CFAR closes.

## 9. Arrays, beamforming and bearing — verified ✅

### Beam pattern closed forms

| property | check | result |
|---|---|---|
| Peak | `B = 1` at the steer angle, and nothing exceeds it over ±90° | exact |
| Nulls | `sin θ − sin θ₀ = m λ/(Nd)` for m = 1…4 | `\|B\|` < 1e-15 |
| Steered nulls | `asin(sin θ₀ + λ/Nd) − θ₀`, at 0° / 30° / 60° | exact, and wider when steered |
| First sidelobe | −12.797 / −13.147 / −13.233 / −13.254 / **−13.260** dB for N = 8…128 | converges to −13.26 |
| Grating lobes | flagged iff a full-height repeat actually exists (`\|B\| > 0.9` away from the beam) | consistent |
| `d ≤ λ/(1+\|sin θ_max\|)` | λ/2 for endfire, λ for broadside-only | exact |

### Array gain, measured not asserted

Monte Carlo through the actual beamformer, 400 trials each:

| N | 10 log₁₀(N) | measured |
|---|---|---|
| 4 | 6.021 | 5.916 |
| 16 | 12.041 | 12.069 |
| 64 | 18.062 | 18.355 |

### Bearing against the Cramér-Rao bound

The spatial twin of §4's arrival-time bound. Both structural properties of
`var ≥ 6/(ρ (kd cos θ)² N(N²−1))` verified:

| N | CRLB | ratio to N/2 (theory 2.828) |
|---|---|---|
| 8 | 1.98991° | — |
| 16 | 0.69939° | 2.8452 |
| 32 | 0.24691° | 2.8326 |
| 64 | 0.08726° | 2.8295 |

and the `1/cos θ` degradation reproduced to 1e-6 at 0° / 30° / 60° / 75°.

**The estimator sits on the bound.** Conventional beamforming with parabolic
peak refinement, 400 trials per point:

| noise σ | measured | CRLB | ratio |
|---|---|---|---|
| 2.0 | 0.56146° | 0.49752° | **1.13** |
| 1.0 | 0.25297° | 0.24876° | **1.02** |
| 0.4 | 0.09866° | 0.09950° | **0.99** |

No factor-of-four surprise this time — unlike v0.2's arrival-time estimator,
where the coherent and envelope bounds differ by 4.16 and picking the wrong one
made an efficient estimator look broken. Here there is no carrier to lose: the
phase gradient across the aperture *is* the signal.

A 32-element half-wave array resolves bearing to **one thirteenth of its own
beamwidth** (0.247° against 3.17°).

### Split-beam

Exact to 1e-16 in the noiseless case across the mainlobe, unambiguous precisely
out to the first null (3.583° for this array), and demonstrably wrapping beyond
it — a source at 8.96° reads as 1.76°. That is the documented limit, not a
defect.

### Joining up with §8

Ten times the elements gives a tenth of the beamwidth gives **exactly 10.00 dB**
of echo-to-reverberation ratio, matching what §8 charges per decade of
ensonified area. The same change buys 10.0 dB of array gain against isotropic
noise. Two independent mechanisms agreeing is a stronger check than either alone.

### A limit the library states rather than hides

Phase steering holds only while the signal stays correlated across the array's
traversal time. For a 32-element half-wave array: 249 Hz of usable bandwidth at
15°, 91 Hz at 45°. **The library's own waveforms are 12 kHz chirps**, so a
phase-steered array cannot handle them off broadside. `narrowband_bandwidth_limit_hz`
returns the number.

## 10. Wideband and adaptive beamforming — verified ✅

### The structural check that caught a sign error

A narrowband tone pushed through the **wideband** path must reproduce
`array_factor` in **closed form**. Worst departure across 21 steer angles:
**5e-3**.

That test earned its place immediately. The first version of the test harness
delayed the element at `+x` instead of advancing it — the mirror of
`synthesize_plane_wave`'s convention — and the beamformer duly steered every
beam to the wrong side. The energy-contrast test only said the contrast was
poor; this one said the pattern was mirrored.

### Wideband steering where phase steering fails

16-element half-wave array, 5 ms LFM over 8–20 kHz, source at 35°:

| | |
|---|---|
| phase-steering bandwidth limit at 35° | 325 Hz |
| signal bandwidth | 12000 Hz — **37× over** |
| beam energy on target | 23.80 dB |
| best off-target beam | 6.79 dB |
| contrast | **17.0 dB** |
| matched filter peak | lag **400**, pulse placed at 400 |

The peak landing exactly on the placement is the part that matters: a
phase-steered beam would smear the chirp and move it.

### A PDW that carries a bearing

The analyser run on each of 17 beams, source at −22° with noise:

| beam | detections | SNR |
|---|---|---|
| −25° | 1 | 37.83 dB |
| **−20°** | 2 | **38.33 dB** |
| −15° | 1 | 24.43 dB |
| 0° | 1 | 14.85 dB |

Strongest beam −20° against a −22° truth, on a 5° scan grid.

### Shading

Peak sidelobe measured from the shaded pattern, against the published window
values:

| window | measured | published | width ×uniform | loss |
|---|---|---|---|---|
| Uniform | −13.26 dB | −13.26 | 1.00 | 0.00 dB |
| Hann | −31.5 dB | −31.5 | 1.64 | −1.76 dB |
| Hamming | −42.7 dB | −42.7 | 1.47 | −1.34 dB |
| Blackman | −58.1 dB | −58.1 | 1.90 | −2.37 dB |

Every window widens the mainlobe and gives up gain: verified as an assertion,
not just observed.

### MVDR resolves what the aperture cannot

16-element half-wave array, conventional resolution limit **7.16°**, two
incoherent sources **4.30°** apart:

| beamformer | peaks found |
|---|---|
| conventional (a^H R a) | **1** |
| MVDR (1 / a^H R⁻¹ a) | **2**, at ±2.145° against a truth of ±2.149° |

And MVDR agrees with conventional on a single source (11.000° vs 10.988° for an
11° truth), so the resolution is not bought with a bias.

**Diagonal loading is verified as necessary, not decorative.** A rank-1
covariance over 16 elements: unloaded, the Cholesky fails and `mvdr_power`
returns 0; with 1% loading it returns a spectrum. Reporting the failure beats
returning noise shaped like an answer.

## 11. Tracking — verified ✅

### Filter consistency: the check that finds what position tests cannot

A Kalman filter whose covariance is wrong still tracks. It just lies about how
well, and then gates correct measurements out or accepts clutter, with nothing
in its output to say so. The normalised innovation squared is chi-square with
2 degrees of freedom when the filter is consistent:

| statistic | measured | expected |
|---|---|---|
| mean NIS (1980 samples) | **1.951** | 2.000 |
| fraction under the 95% gate | **95.0%** | 95% |
| fraction under the 99% gate | **99.0%** | 99% |

### What the filter is worth

| quantity | result |
|---|---|
| RMS position error, raw | 50.97 m |
| RMS position error, filtered | **20.89 m** — 2.44× better |
| velocity (truth 4.00, −7.00) | (4.17, −7.18) |
| closing rate (truth 7.486 m/s) | 7.683 m/s |

Closing rate matters twice over: the Doppler bank of §5 measures it directly,
so the tracker and the bank are two independent routes to the same number.

### Gating

The chi-square gate inverts in closed form (`−2 ln(1−p)`): 5.991 at 95%, 9.210
at 99%, verified to 0.002. On a converged track: NIS **0.268** on target,
**652.7** for a measurement 300 m off, **319.0** for one 20° off.

### Track management

| scenario | result |
|---|---|
| one target, one isolated false alarm, 8 scans | **1 confirmed, 0 tentative** |
| 494 false alarms scattered over 300 scans | **0 confirmed tracks** |
| one missed scan | coasts, then reacquires and re-confirms |
| three missed scans | deleted |
| two well-separated targets | 2 confirmed, each on its own target |

### A roadmap claim, corrected

The v0.8 roadmap said tracking would finally suppress the cross-template ghosts
of §4. **That was wrong**, and the test suite says so rather than the claim
being quietly dropped.

A ghost appears whenever the real arrival does, at a fixed offset set by the
template cross-correlation. It is therefore *exactly* as consistent over time as
the target, moves with it, and forms its own perfectly healthy confirmed track.
Measured: one target plus its ghost gives **two** confirmed tracks, and the test
asserts that outcome so the correction cannot rot.

Tracking kills false alarms — which do not repeat — decisively. Ghosts repeat.
Suppressing them needs the offset and amplitude ratio recognised as a template
artefact, which is a different mechanism. **§12 implements it**; this section is
kept as written because the correction, not just the eventual fix, is the point.

## 12. Fusing what is already measured — verified ✅

Four existing outputs joined to consumers that could already have used them.
Nothing new is measured.

### Range rate: the gate must follow the dimension

Fusing the Doppler bank's closing rate makes the measurement 3-dimensional, so a
gate set for 2 dof is the wrong probability. The 3-dof quantile has no closed
form and is bisected on the exact CDF; it reproduces the textbook values:

```
3-dof gates: 95% at 7.815, 99% at 11.345
```

Consistency survives the fusion. Over 1980 fused updates:

```
mean NIS 2.998 (theoretical 3)
under the 95% gate: 94.7%   under the 99%: 98.9%
```

### What the range rate is worth — and when it stops being worth it

Radial velocity error, the only component a range rate observes:

| scans | position only | with range rate | ratio |
|---|---|---|---|
| 3 | 4.213 m/s | **1.220 m/s** | 3.45× |
| 6 | 2.069 m/s | **0.649 m/s** | 3.19× |
| 20 | 0.365 m/s | 0.368 m/s | 0.99× |

**The gain fades to nothing by twenty scans**, and that is reported rather than
buried. By then the filter's own estimate (0.365 m/s) already beats the
σ = 2 m/s measurement, and a measurement helps exactly as long as it beats the
estimate you have. The value is in the first few scans of a new contact — which
is when a decision usually has to be made.

The test asserts `fused < unfused × 1.05` at every step count, not
`fused < unfused`. The first assertion was written, failed at 20 scans, and was
found to be wrong about the physics rather than the code.

### An unresolved Doppler bin is not a measurement of zero

A bank that fails to resolve a shift reports 0 m/s. That is a different
statement from "the target is stationary", and treating it as a measurement pins
the track's radial velocity to zero:

```
truth closing                8.00 m/s
zero taken as a measurement  5.686 m/s   (28% low)
flag left clear              9.163 m/s
```

Hence the separate `has_range_rate` flag. This is a trap the API could have
walked into silently, so it has a test.

### Ghost recognition, and the check that makes it safe

| case | result |
|---|---|
| target + cross-template ghost | **1** established track |
| different bearing / comparable strength / distant | all kept |
| line-astern formation, **same** waveform | **2 tracks (kept)** |
| same geometry, **different** waveforms | 1 track (suppressed) |

The last two rows differ **only in the waveform label**. Two real targets lit by
one sonar return the same waveform; a ghost is by definition a different
template. Without that check the line-astern formation would be silently
deleted, and deleting a real target is far worse than carrying a ghost.

A note on counting: the test counts Confirmed **or Coasting**. A track that
missed the most recent scan is still a contact, and counting only Confirmed
reports the detector's instantaneous miss rate rather than the number of
targets. The first version of this test counted Confirmed only and failed — on
a run where suppression had correctly not fired and one real track happened to
be coasting.

### Global-cost-ordered association

Two targets crossing at 8 m/s:

```
before crossing: left=1 right=2
after  crossing: left=2 right=1
```

Identities survive. Greedy-by-arrival-order swapped them. The fix is not an
assignment algorithm: it builds every gating pair, sorts by NIS, and assigns
best-first. Still greedy — but over the global cost ordering rather than over
arrival order, which is what this case needs. The gate was not touched; a swap
is not a gating failure, since near the crossing both measurements gate against
both tracks perfectly well.

### Spatial smoothing: MVDR's failure, and its repair

Two **coherent** arrivals at ±6° — which is what §10's multipath produces, one
arrival by two paths:

| | peaks found |
|---|---|
| plain MVDR, 16 elements | 4 (spurious) |
| forward-backward, 10-element subarrays | **2, at −6.00° and +6.00°** |

A single incoherent source is unmoved: 9.00° truth reads 8.987° smoothed.

**The cost is stated, not hidden:** resolution falls from 7.16° (16 elements) to
11.46° (10), and a test asserts that it does. Smoothing spends aperture;
resolving `P` coherent sources needs subarrays longer than `P`, and the library
makes the caller choose.

The first implementation had transposed indices in the backward term and updated
in place over entries it had yet to read. It produced a matrix that still looked
like a covariance and read a 9° source at 3.3°. The invariant that catches this
is **persymmetry**: a correct forward-backward result has it, a transposed one
does not.

## 13. Interacting multiple models — verified ✅

See `docs/math_spec.md` §18 for the algorithm.

### The trade, split by phase

Target turning at 3 °/s between scans 15 and 35, against a single-model EKF
tuned with process noise *between* the IMM's two:

| phase | IMM | single-model EKF | ratio |
|---|---|---|---|
| during the turn | 37.31 m | 36.68 m | **0.98×** |
| after it ends | 25.88 m | 44.53 m | **1.72×** |
| overall | 44.91 m | 50.74 m | 1.13× |

**The gain is in recovery, not in the turn**, and during the manoeuvre the IMM
is marginally worse. That is not the usual description of an IMM, and it is what
the measurement says.

On a straight target it costs nothing: 27.40 m vs 27.88 m, so the win is not
simply more process noise in disguise.

### The manoeuvre detector

Peak `1 − μ_CV` during the turn: **0.991**. Direction recovered correctly
(truth ±3 °/s → ±2.03 / ±2.27 °/s). The estimate saturates at ±ω by
construction and is not a turn-rate measurement.

### Invariants

Over 200 scans of alternating 6 °/s turns: worst `|Σμ − 1|` = 2.22e-16, smallest
model probability 1.4e-2 — no model is ever driven to zero, from which it could
not recover.

The `w → 0` limit converges as 6.946e-N with no jump at the series crossover.
In `float` it floors at 5.0e-4 m, which is the precision limit of differencing
two ~5000 m positions, not a defect; the test is precision-aware.

## 14. Sound speed against published values — verified ✅

The most important result in this release is a **bug that eleven releases of
mutual-agreement testing did not find**.

### The bug

`chen_millero()` held Wong & Zhu's (1995) ITS-90 coefficients for 40 of 42
terms and Chen & Millero's 1977 originals for `A02` and `A03`. It was therefore
neither equation, and it missed the official check value by 0.016 m/s.

Nothing caught it because the only tests were mutual agreement between Medwin,
Mackenzie and Chen-Millero — and those agree only to ~0.1 m/s, which cannot see
a 0.016 m/s error. **Mutual agreement bounds how wrong you can be by how much
your methods differ. That is not the same as being right.**

### The published check value

UNESCO Technical Papers in Marine Science 44 (Fofonoff & Millard 1983), p. 48
and the FORTRAN listing on p. 49:

```
CHECKVALUE: SVEL=1731.995 M/S, S=40 (PSS-78), T=40 DEG C, P=10000 DBAR
```

`chen_millero_1977()` returns **1731.9954**.

### The published table

The same document, p. 50: 220 values over S = 25/30/35/40 PSU, T = 0…40 °C
(IPTS-68), p = 0…10000 dbar. All 220 are in `tests/data/unesco44_table.hpp`.

```
220 published values, worst |delta| = 0.0499 m/s (at S=35 T=0 p=10000 dbar)
the table is printed to 0.1 m/s, so 0.05 is its rounding half-width
```

The agreement is as tight as the published data can express. Note also that a
transcription error would have shown as a large outlier at one cell rather than
a uniform near-miss, which is how the transcription itself was checked.

### Two equations, not one equation and one mistake

| | same T fed to both | after t68 = 1.00024 t90 |
|---|---|---|
| worst |1977 − ITS-90| | 0.0208 m/s | **0.0056 m/s** |

Converting the temperature scale improves the agreement 3.7×, and the 0.0208
figure sits inside Wong & Zhu's own stated "within 0.024 m/s" for their
revision. That is what makes the two versions legitimate rather than a typo.

### Del Grosso: an independent equation, and where a validity box lies

| region | worst |UNESCO − Del Grosso| |
|---|---|
| 0 – 1000 m (most sonar work) | 0.41 m/s |
| 0 – 3000 m | 0.76 m/s |
| 0 – 6000 m | 0.82 m/s |
| 0 – 9810 m (Del Grosso's own limit) | 1.33 m/s |
| full nominal validity box | **3.93 m/s** |

The bottom row is 26 °C water under 1000 bar — a combination no sea contains and
neither equation was fitted to. Quoting a disagreement over a whole (S, T, P)
rectangle measures its corners, not the ocean.

Read the top row against Chen & Millero's quoted standard deviation of
0.19 m/s: two independent equations differ by about twice the uncertainty
either claims. **That is the real accuracy of the speed of sound in seawater**,
and it is why this document refuses to quote ray paths to the metre.

The bounds asserted in these tests are regression bounds set from the measured
values, not predictions from theory. There is no theory that says how far two
empirical fits should drift apart.

## 15. Real ocean profiles — verified ✅

Six WOA23 profiles, fetched from NOAA NCEI by `tools/fetch_woa_profiles.py`.
Every number's provenance is a URL that can be re-fetched.

| site | max depth | c surface | c axis | axis depth |
|---|---|---|---|---|
| aegean | 650 m | 1522.1 | 1514.4 | 175 m |
| black-sea | 2200 m | 1486.9 | 1462.6 | 55 m |
| eq-pacific | 4300 m | 1539.1 | 1484.9 | 950 m |
| levantine | 2400 m | 1530.6 | 1516.1 | 375 m |
| n-atlantic | 4700 m | 1509.0 | 1486.3 | 950 m |
| norwegian | 3200 m | 1477.4 | 1463.7 | 850 m |

**The shipped descriptions are tested against the shipped data.** An earlier
draft called the Norwegian Sea upward-refracting and put the North Atlantic axis
at 1100 m. Neither survived contact with the profiles; both were corrected, and
a test now asserts the structure each description claims.

### Where ignoring salinity stops being safe

Error from assuming 35 PSU instead of the measured value:

| site | worst |
|---|---|
| black-sea | **20.32 m/s** |
| aegean | 4.68 m/s |
| levantine | 4.74 m/s |
| n-atlantic | 0.82 m/s |
| eq-pacific | 0.58 m/s |
| norwegian | 0.19 m/s |

The Black Sea surface is 18.2 PSU. The 20.3 m/s error is the same size as the
whole 24 m/s channel excess measured there — assume 35 PSU and the channel you
trace is not the one that exists.

### A test that was measuring itself

Snell's invariant over real profiles reports **exactly zero** drift, which is
suspicious rather than reassuring: the tracer stores ξ and derives θ from it, so
`cos θ / c` largely measures whether `acos` and `cos` round-trip. It cannot
catch a ray that turns in the wrong *place*.

The turning-depth check can. At θ = 0 the local sound speed must equal
`c_source / cos θ₀`, which is independent of how θ is stored. Over **186 turning
points** across six real profiles the worst error is below **1e-9 m/s**.

The first version of that check reported errors up to 4.3 m/s that were entirely
its own: it detected turns by a sign change across index `i`, and because the
turn itself sits at exactly zero — and zero is not `> 0` — the point *before* the
turn also looked like a reversal.

## 16. Turn-rate estimation and the IMM tracker — verified ✅

### A bracket does not degrade gracefully

Models bracketed at 3 °/s, against a five-state filter that estimates the rate:

| truth | IMM (bracketed) | CTRV (estimated) |
|---|---|---|
| 2 °/s | 1.27 °/s | 1.02 °/s |
| 5 °/s | 2.62 °/s | **4.99 °/s** |
| 8 °/s | 1.82 °/s | **5.64 °/s** |

The 8 °/s row is the interesting one: the IMM reports *less* turn than at 5 °/s.
Once the truth leaves the bracket the ±ω models fit so badly that probability
drifts back to constant velocity and the blend collapses toward zero.

### Zero process noise on omega: confidently wrong

Target straight for 30 scans, then turning at 5.00 °/s:

| `q_w` | reported σ | estimate |
|---|---|---|
| **0** | **0.11 °/s** | **2.08 °/s** |
| 0.005 | 0.91 °/s | 4.97 °/s |
| 0.01 (default) | 1.54 °/s | 4.90 °/s |
| 0.02 | 2.67 °/s | 4.89 °/s |

With no random walk the covariance shrinks, the gain on the fifth state dies,
and the filter reports the **smallest uncertainty about the most wrong answer**.
This is the failure mode a test suite exists for: everything about the output
looks healthy.

The default was 0.02 when the filter was first written and was lowered to 0.01
on this measurement.

### What the fifth state costs

Straight target, 12 runs, against a plain CV EKF: **36.28 m vs 35.74 m** —
a 1.02× penalty. An extra state costs almost nothing when the measurement never
moves it.

### The zero-turn limit, again

Same 0/0 as §13's IMM, but now appearing in the **Jacobian** as well as the
transition, and both must switch to the series at the same threshold. Converges
as 8.322e-N with no jump at the crossover; in `float` it floors at the
precision limit and the test is precision-aware.

### The tracker

`imm_tracker_step()`, gating on the combined estimate:

```
two targets turning in opposite directions, 40 scans -> 2 established tracks
  each detects its own manoeuvre, signs opposite and cancelling
367 false alarms over 200 scans                      -> 0 established tracks
```

M-of-N confirmation survives the swap from a single filter to a mixture.

## 17. Spread-spectrum communication — verified ✅

### PN sequences, measured not trusted

Every degree from 5 to 15: balance exactly +1, periodic autocorrelation exactly
−1 at every non-zero shift (0.0e+00 deviation). Those two properties are what a
maximal sequence has and a short cycle does not.

The first tap table used the wrong convention for a right-shifting LFSR, leaving
bit 0 clear. A seed of 1 fed back zero, the state collapsed to all-zeros on the
first step, and the generator emitted a constant — with the correct *length*.
Only these two measurements caught it.

### Processing gain: the claim this release got wrong first

**Against white noise at fixed energy per bit, spreading buys nothing.** Measured
with amplitude scaled as `1/√N`:

| chips | 31 | 63 | 127 | 255 | 511 |
|---|---|---|---|---|---|
| output SNR | 13.33 | 12.46 | 13.63 | 13.72 | 12.94 dB |

Flat to 0.87 dB over a 16× range.

The first version of this test held the chip **amplitude** fixed while sweeping
N, and reported 12.17 dB against a theoretical 12.17 dB — exact agreement, and
meaningless. Energy per bit is `A²·N·T_chip`, so that sweep multiplied the
transmitted energy by N and reported the result as processing gain. Both sides
of the comparison were the same tautology, which is why they matched perfectly.

This is the third measurement in this project to fail by being compared against
itself, after the sound-speed equations (§14) and the Snell invariant (§15).

**Where the gain is real: bandwidth expansion at fixed data rate.** Data rate
held at 100 bps, chip rate grown with N so the bandwidth grows and the energy
per bit does not:

| chips | chip rate | white noise | narrowband tone |
|---|---|---|---|
| 31 | 3.1 kHz | 21.98 dB | 11.82 dB |
| 511 | 51.1 kHz | 21.28 dB | **18.91 dB** |

White noise flat to 0.75 dB; the tone gains **+7.1 dB**, against an ideal 12.2.
The shortfall is real and explained: a bare correlator collects the diluted
interferer from the whole band, where `10 log10(N)` assumes a filter matched to
the data bandwidth after despreading.

A trap in the module's own API is documented at the point of use: raising
`chips_per_bit` at a fixed chip rate does not spread anything, it trades data
rate for integration time.

### The Doppler ceiling on processing gain

Chip slip at 1500 m/s sound speed, and the longest code holding it under a
quarter chip:

| speed | longest code | gain available |
|---|---|---|
| 1 m/s | 375 chips | 25.7 dB |
| 3 m/s | 125 chips | 21.0 dB |
| 10 m/s | 37 chips | 15.7 dB |

Speed caps the achievable gain, and no transmit power lifts that cap.

### HFM vs LFM under Doppler

Correlation coefficient as a percentage of its own zero-Doppler value, 12 kHz
sweep, 20 ms:

| speed | LFM | HFM |
|---|---|---|
| 10 m/s | 86.0 % | 99.7 % |
| 20 m/s | **54.8 %** | **99.1 %** |

**Three bugs in this one test**, each producing a plausible answer:
`render_real_doppler` takes `v/c` not `1 + v/c` (making the baseline a 2×
compressed pulse, so ratios exceeded 100% — impossible for a matched filter, and
that is what exposed it); raw peaks compared instead of correlation
coefficients, so energy differences masqueraded as mismatch; and a lag window
starting at zero, which missed the *negative* peak shift an upsweep produces.

### Coding, and where the layers hand off

CRC-32 against the published check value: `CRC-32("123456789") = 0xCBF43926`.

RS(15,11) over 4000 codewords with 0–2 symbol errors: **0 left wrong, 0
uncorrectable**. An 8-bit burst destroying two adjacent symbols is corrected
exactly.

With 4 errors — beyond `t = 2` — over 3000 codewords:

```
1836  RS admitted defeat
 992  RS miscorrected to a different valid codeword
        of those, 992 caught by the frame CRC, 0 escaped
```

Miscorrection is the code's minimum distance (`d = 5`) showing, not a defect,
and the CRC is the layer that exists to catch it. The two together let nothing
through in 3000 trials.

## 18. Bindings and portability — verified ✅

### The C ABI

269 checks, compiled by a **C compiler in C11 mode** with
`-Wstrict-prototypes -Wmissing-prototypes -Werror`. Values with published
answers are round-tripped through the boundary rather than merely returned:
the UNESCO check value comes back as 1731.9954, `CRC-32("123456789")` as
0xCBF43926, and Snell's invariant holds across a ray traced through the C API.

### 160 kB of .bss that the size report found

The first C ABI traced rays through a static 8192-point scratch buffer, to avoid
assuming `ph_ray_point` and `RayPoint` share a layout. Cross-compiling and
reporting section sizes showed the cost:

| | before | after |
|---|---|---|
| text | 61216 | 61160 |
| data | 1536 | **0** |
| bss | **163840** | **0** |

Half the RAM of the Cortex-M7 the library is meant to fit on, for a copy that
does nothing — and it made the function non-reentrant. The layout assumption did
not need avoiding, it needed `static_assert`-ing: size, alignment and every
member offset are now checked at build time, so a mismatch breaks the build
rather than scrambling ray paths.

The library now carries **no mutable static storage at all**. CI asserts
`data == 0 && bss == 0` so it cannot return.

### Cross-compilation

| target | word size | Real | result |
|---|---|---|---|
| arm-linux-gnueabihf | 32-bit | float | clean, 0 allocator symbols |
| riscv64-linux-gnu | 64-bit | double | clean, 0 allocator symbols |

Two architectures, two word sizes and both precisions, because the bugs this
catches only appear when `size_t`, pointer width and `Real` change together.

**Not verified: the bare-metal Cortex-M7 build.** `cmake/cortex-m7.cmake` is
provided but this environment's `arm-none-eabi` toolchain ships no C++ standard
library, so it could not be run. The ARM figures above are from
`arm-linux-gnueabihf` — same instruction set, hosted target. The toolchain file
is untested and labelled as such.

### -fno-exceptions -fno-rtti

Builds clean and references no `__cxa_throw`, `_Unwind_*` or
`__gxx_personality` symbols. A no-op that is worth having as a build: if the
library ever starts throwing, this is a link error in CI rather than an unwinder
in a firmware image.

### Rust

8 binding tests pass. They check the **binding**, not the physics — a field in
the wrong order, an enum discriminant off by one, a layout mismatch — by either
round-tripping a published value or comparing against what the library reports
about itself.

Two findings worth recording:

- `ph_profile_speed_at` answers **0 m/s** for a profile with fewer than two
  samples. That is a silent failure a caller can propagate; the safe Rust
  wrapper returns `Err(Error::State)` instead. Found by a binding test.
- Precision mismatch between library and crate is **not** a link error. It
  reinterprets every float crossing the boundary, so `check_precision()` exists
  and every constructor calls it.

## 19. Acoustics in air — verified ✅

Air is here because it is the medium anyone can test in: a speaker and a
microphone are already on the desk. The DSP does not care which medium it is in,
so an air bench exercises the whole chain before any wet hardware is bought.

### Against ISO 9613-1:1993 Table 1

105 published attenuation coefficients at 10 °C, one standard atmosphere:

```
worst relative error 0.380%
worst case at 315 Hz, 90% RH: published 1.30, computed 1.305 dB/km
```

That worst case is the table's own three-figure rounding, not a disagreement.

**Note 5 of the standard, measured.** The table was computed at exact
one-third-octave midband frequencies `1000·10^(k/10)`, not the nominal values in
its own row headings:

| | worst error |
|---|---|
| exact midband | 0.380% |
| nominal | **1.733%** (4.6× worse) |

The transcription was verified the same way as the UNESCO table in §14: the
equations were entered separately from the standard's clause 6.2, and a
transcription slip would appear as one large outlier rather than a uniform
near-miss. The PDF's text layer is OCR-damaged in places, so only unambiguous
cells were included — 105 of them.

### Humidity is not a correction

At 10 °C, 4 kHz: going from 10% to 20% relative humidity **multiplies** the loss
by 1.60. And the sign of the effect reverses with frequency — at 500 Hz damp air
absorbs less, at 8 kHz more. A bench measurement in air that does not record the
humidity means nothing.

### Air is the harsher Doppler environment

`v/c` is 2.91e-3 in air against 6.67e-4 in water, a factor of 4.37. A 511-chip
spreading code slips **1.49 chips per bit at 1 m/s** in air, where in water it
slips 0.34. That makes an air bench a stronger test of §17's comm module than a
water one, not a weaker substitute.

### What the bench cannot do, and one thing it could not do here

`tools/air_bench.py selftest` recovers a simulated direct path at 1.500 ms
(truth 1.500) and an echo at 4.00 ms, so the analysis path is verified without
hardware.

**Real hardware-in-the-loop was attempted on this machine and does not work.**
The audio devices enumerate — an ALC256 with analogue playback and capture — but
there is no working acoustic path. A controlled test settles it: recording while
playing five 50 ms chirps gave *less* energy in the 2–8 kHz band than recording
in silence, by 7.85 dB. Capture is also dominated by a large startup transient.
Whatever the cause, no speaker-to-microphone coupling could be demonstrated, so
no hardware measurement is claimed. The bench tool is written and self-verified;
running it against real hardware is for a machine that has some.

Two limitations hold even on working hardware, and are stated in the tool:

- **Timing and detection only, not level.** Consumer speakers and microphones
  have unknown, uncalibrated, non-flat responses.
- **The measured delay includes the sound card's latency**, typically 10–50 ms
  against 1.5 ms of flight time over half a metre. Only a two-distance
  measurement removes it, because the difference cancels fixed electronic delay.

## 20. Turning a recording into a measurement — verified ✅

### Qualification catches what inspection does not

| channel | excess | clipped | verdict |
|---|---|---|---|
| carrying the probe | +18.92 dB | 0% | **usable** |
| carrying nothing | −0.27 dB | 0% | unusable |
| carrying it, clipping | +31.26 dB | 18.24% | **unusable** |

The third row earns the test: it passes the signal check by a wide margin and is
rejected anyway. A qualifier that only asked "is the probe there" would accept a
recording whose every level is fiction.

This exists because v0.15's dead channel did not look dead — thousands of
distinct sample values, healthy RMS, a matched-filter peak 44 dB above
background, all of it electrical noise. Only the A/B comparison against a silent
recording settled it, and it came out **negative**.

### Two-distance latency cancellation

Truth: c = 343.37 m/s, fixed latency 32.0 ms, distances 0.40 m and 1.60 m.

| | |
|---|---|
| naive single-shot `d1/t1` | 12.1 m/s (28× wrong) |
| two-distance solve | **343.37 m/s**, **32.00 ms** |

Both exact. The solver refuses under-separated distances rather than dividing by
timing noise.

### Impulse response

Three known paths at 200/320/512 samples, gains 1.0/0.5/0.25, recovered at
200.0/320.0/512.0 with −6.03 dB and −12.06 dB against truths of −6.02 and
−12.04.

The arrival picker's separation guard is now a few times the probe resolution
`1/B`. v0.15's tool used 20 ms — seven metres of extra path — and silently
discarded every echo a room produces.

### Still not measured

**No number in this document comes from a real transducer.** §19 records why the
development machine cannot produce one; this section is the apparatus, verified
against channels whose answers are known by construction.

## 21. Build and runtime hygiene ✅

| Configuration | Result |
|---|---|
| GCC 15.3, `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wold-style-cast -Wcast-qual -Wdouble-promotion -Werror` | clean |
| Clang 21.1, same set | clean |
| AddressSanitizer + UndefinedBehaviorSanitizer | clean, 95806 checks |
| `Real = float` (`-DPHANTOM_REAL_FLOAT=ON`) | passes, all four compiler × precision configurations |
| Heap allocations during execution | **zero** — proven by `nm` over the built archive: the library references only libm and `memset` |

### A false verification, and the guard added for it

During v0.10 the four-configuration sweep was run with `-DPHANTOM_USE_FLOAT=ON`.
That option does not exist — the real one is `PHANTOM_REAL_FLOAT` — and CMake
accepts unknown `-D` values silently. Both "float" configurations were ordinary
double builds, and all 41044 checks passed, which is exactly what makes it
dangerous: the failure produces a *clean verification report* for a run that
never used float.

It was caught by reading the `CMakeCache.txt`, not by any test. `CMakeLists.txt`
now rejects that name and eight other near misses with a fatal error. The
genuine float runs were then made, and are the ones reported above; the two
builds are confirmed distinct by the Snell drift, 7.08e-08 against 1.67e-16.

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

## 22. Performance (measured, not claimed)

13th Gen Intel Core i7-13620H, `-O3`, single thread:

| Operation | Time |
|---|---|
| `mackenzie()` | 1.5 ns |
| `chen_millero()` (UNESCO) | 4.0 ns |
| `speed_at()` — 501-layer binary search | 16.6 ns |
| **`trace_ray()` — per arc step** | **52 ns** (median) |
| `trace_ray()` — full 100 km ray, 1132 arcs | 40.8 µs |
| `analyze_sofar()` — 501 points | 596 ns |
| `fft_forward()` — 8192 points | 307 µs |
| `matched_filter_apply()` — one template, 8192-point block | 580 µs |
| `analyze_block()` — 4 templates | 2.26 ms |

At 96 kHz with an 8192-point transform the block advance is 65.3 ms of audio, so
a four-template bank runs about **29× real time on one core**, and end-to-end
detection latency is **65 ms**, dominated by the block length rather than the
arithmetic. Halving the transform halves the latency and halves the longest
detectable pulse.

The FFT is a plain radix-2 with no SIMD and is roughly an order of magnitude off
a vendor-tuned transform. That is the price of the zero-dependency promise; the
correlation is the only hot spot, and it is isolated behind `fft_forward()` /
`fft_inverse()` if you would rather pay a dependency for it.

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

## 23. Self-consistent but NOT independently verified ⚠️

- **Chapman-Harris surface backscatter coefficients.** The formula is written
  out in the header so a reader can check it, and its behaviour is verified
  (monotone in wind and angle, right band, correct structure at 30°). The
  coefficients themselves have NOT been checked against the 1962 paper. Two
  reference values in this repository have already turned out to be misremembered
  rather than wrong in the code — Thorp at 50 kHz and a critical angle — so this
  one is flagged rather than assumed.
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

## 24. Known limitations — not bugs, scope ⚠️

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

## 25. Not implemented in v0.16

`BioMimeticCommEngine` is not in this release. See `docs/roadmap.md`.

Known gaps within what *is* shipped:

- **No verified hardware measurement.** See §19: this machine has no working
  acoustic path, so every number in this document is simulation or published
  reference. Nothing here has been measured against a real transducer.
- **Air sound speed is not Cramer (1993).** No CO2 term, no pressure
  dependence; a few tenths of a m/s over 0-40 C at sea level.
- **The C ABI does not cover everything.** Beamforming, eigenrays and
  reverberation are absent: their C++ interfaces take several spans at once and
  a C shape for them has not been designed. "Stable" means what is there will
  not change shape, not that everything is there.
- **The Cortex-M7 toolchain file is untested.** See §18.
- **No carrier or timing recovery.** The demodulator is GIVEN the carrier
  phase. Carrier recovery is a real subsystem and this module does not pretend
  to one; the HFM preamble reveals the time scale but applying it is the
  caller's job.
- **No equaliser.** Multipath is modelled elsewhere in the library but the
  demodulator does not undo it, so the delay spread a real channel imposes is
  not compensated.
- **CTRV has no tracker of its own.** `imm_tracker_step()` exists; the
  five-state filter is single-target, and multi-target CTRV tracking is not
  implemented.
- **No IMM over CTRV models.** The two approaches are separate; combining them
  (an IMM whose models each estimate their own turn rate) is the standard next
  step and is not done.
- **Profiles are a climatology**, not casts: no eddy, internal wave or diurnal
  layer survives a decadal 1° average.
- **Association is greedy over a global cost ordering**, not a Hungarian
  assignment. It fixes the crossing case (§12) but is not provably optimal.
- **Ghost recognition needs a waveform label.** A caller that reports the same
  label for every detection gets no suppression — correctly, since the label is
  the only thing separating a ghost from a real contact.
- **Spatial smoothing spends aperture.** Resolving coherent sources costs
  resolution, and the subarray length is the caller's choice.
- **Uniform line arrays only**; narrowband element model.
- **Reverberation is a level, not a field**, and does not come from the traced
  paths.
- **Ray theory, with its caustics**; plane-wave flat-interface reflection;
  monostatic multipath; range-independent ocean.

## Reproducing

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/phantom_tests          # 95806 checks
./build/phantom_bench          # timing table above
./build/munk_simulation data   # CSV + console report
./build/ping_intercept data    # streaming ping detection demo
./build/countermeasure_loop    # intercept -> classify -> reply -> multipath
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
