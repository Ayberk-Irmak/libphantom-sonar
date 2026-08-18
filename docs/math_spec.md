# Mathematical specification

Every formula the library implements, with its provenance and its validity
range. If a number in the code does not appear here, that is a bug.

---

## 1. Sound speed in seawater

Three equations are implemented. They exist for different reasons, and picking
the wrong one is the most common modelling error in this domain.

### 1.1 Medwin (1975) — `sound_speed::medwin`

```
c = 1449.2 + 4.6 T − 0.055 T² + 0.00029 T³
  + (1.34 − 0.010 T)(S − 35)
  + 0.016 z
```

`T` in °C, `S` in PSU, `z` in metres. Validity: `0 ≤ T ≤ 35`, `0 ≤ S ≤ 45`,
`0 ≤ z ≤ 1000 m`.

**This equation is frequently mislabelled "the UNESCO equation". It is not.**
It is Medwin's six-term simplification. It is `constexpr`, costs ~1.5 ns, and is
the right choice for shallow water and for baking profiles into flash at compile
time.

Its depth term is linear and its stated validity stops at 1000 m. The deep sound
channel axis sits at 1000–1300 m at most latitudes, so **Medwin must not be used
for SOFAR work.** The pressure-driven rise in sound speed below the axis is
exactly what it does not model.

### 1.2 Mackenzie (1981) — `sound_speed::mackenzie` — **library default**

```
c = 1448.96 + 4.591 T − 5.304e-2 T² + 2.374e-4 T³
  + 1.340 (S − 35)
  + 1.630e-2 D + 1.675e-7 D²
  − 1.025e-2 T (S − 35)
  − 7.139e-13 T D³
```

Validity: `2 ≤ T ≤ 30 °C`, `25 ≤ S ≤ 40 PSU`, `0 ≤ D ≤ 8000 m`. Nine terms,
including a quadratic depth term and a `T·D³` cross term, so it stays honest
through the full water column. Also `constexpr`.

### 1.3 Chen & Millero (1977) — `sound_speed::chen_millero` — the real UNESCO algorithm

```
c(S,T,P) = Cw(T,P) + A(T,P)·S + B(T,P)·S^1.5 + D(T,P)·S²
```

with `Cw` a degree-5×3 polynomial in `(T,P)`, `A` degree-4×3, `B` degree-1×1 and
`D` degree-0×1. Coefficients are in `include/phantom/sound_speed.hpp`.

Validity: `0 ≤ T ≤ 40 °C`, `0 ≤ S ≤ 40 PSU`, `0 ≤ P ≤ 1000 bar`.

**Takes pressure, not depth.** This is the reference and validation path. Not
`constexpr`, because the `S^1.5` term needs `sqrt`, which C++20 does not permit
in constant evaluation.

### 1.4 Depth → pressure — Leroy & Parthiot (1998)

```
h45(Z) = 1.00818e-2 Z + 2.465e-8 Z² − 1.25e-13 Z³ + 2.8e-19 Z⁴     [MPa]
g(φ)   = 9.7803 (1 + 5.3e-3 sin²φ)
k(Z,φ) = (g(φ) − 2e-5 Z) / (9.80612 − 2e-5 Z)
P      = 10 · h45(Z) · k(Z,φ)                                      [bar]
```

Latitude changes the answer by <1 bar at 1000 m (≈0.008 m/s of sound speed).
Small, but free to get right.

### 1.5 Munk (1974) canonical profile — `sound_speed::munk`

```
η    = 2 (z − z_axis) / B
c(z) = c_axis · [1 + ε (η + e^(−η) − 1)]
```

Defaults `z_axis = 1300 m`, `c_axis = 1500 m/s`, `ε = 7.37e-3`, `B = 1300 m`.
This is the standard deep-sound-channel benchmark; Bellhop ships the same case
as `MunkB_ray`, which makes it the natural cross-validation target.

---

## 2. Ray tracing

### 2.1 The invariant

For a horizontally stratified ocean, Snell's law integrates to a constant along
each ray:

```
ξ = cos θ(z) / c(z) = const
```

`θ` is the grazing angle from the horizontal, positive downgoing. Everything
below follows from this one relation, and the test suite verifies it holds to
`< 2e-16` relative over 100 km of propagation.

A ray **turns** where `cos θ = 1`, i.e. where

```
c(z_turn) = 1 / ξ
```

### 2.2 Why constant-gradient arcs, not stepping

The profile is stored as **piecewise linear** in depth, not as a stack of
isovelocity slabs. Inside a layer with `c(z) = c₁ + g·(z − z₁)`, the ray is
**exactly a circle**:

```
radius   R   = 1 / (ξ |g|)
centre   z_c = z₁ − c₁ / g       (the depth where c extrapolates to zero)
```

Because the arc is exact, there is no step size to tune and no integration error
to trade against speed. One arc per layer crossing.

The alternative — isovelocity slabs plus Snell refraction at each interface —
produces stair-step artefacts that only vanish as the layer count grows, i.e.
you pay in both accuracy *and* time. The library's
`ray_is_invariant_under_layer_refinement` test makes the difference measurable:
resampling a linear profile from 2 points to 1001 points moves the answer by
`4.6e-12 m` over 40 km. A stepping tracer would move by metres.

### 2.3 Range increment

From `cos θ = ξ c` and `dz = tan θ · dr`:

```
−sin θ dθ = ξ g dz
dr        = −cos θ dθ / (ξ g)
```

which integrates to

```
Δr = (sin θ₁ − sin θ₂) / (ξ g)
```

### 2.4 Travel time

```
dt = ds / c = −dθ / (g cos θ)
```

integrating with `sec θ + tan θ = (1 + sin θ)/cos θ`:

```
Δt = (1/g) · ln[ c₂ (1 + sin θ₁) / (c₁ (1 + sin θ₂)) ]        (form A)
   = (1/g) · ln[ c₁ (1 − sin θ₂) / (c₂ (1 − sin θ₁)) ]        (form B)
```

The two are algebraically identical via `(1+sin)/cos = cos/(1−sin)`. Form A
loses precision as `θ → −90°`, form B as `θ → +90°`. Within one layer both
angles share a sign, so `sign(sin θ₁ + sin θ₂)` selects the stable branch.
`arc_time()` in `src/ray_tracer.cpp` does exactly this.

### 2.5 Inverting time — exact travel-time budgets

Writing `X(θ) = (1 + sin θ)/cos θ = tan(π/4 + θ/2)`, §2.4 collapses to

```
X₂ = X₁ · e^(−g Δt)
```

so the angle after a given travel time inverts in closed form:

```
θ₂ = 2 · [ arctan( tan(π/4 + θ₁/2) · e^(−g Δt) ) − π/4 ]
```

This is why `TraceConfig::max_time_s` cuts an arc **exactly** rather than
stopping short at the previous layer boundary. Both budgets — range and time —
are honoured to machine precision, and whichever binds first wins.

### 2.6 Degenerate cases the kernel handles explicitly

| Case | Condition | Behaviour |
|---|---|---|
| Isovelocity layer | `\|g\| ≤ 1e-12` | straight-line branch |
| Vertical ray | `ξ ≈ 0` | zero range, depth-only propagation |
| Vertex on a sound speed minimum | turn depth == current depth | axial ray, travels horizontally |
| Horizontal ray in isovelocity water | `sin θ = 0` and `g = 0` | travels horizontally to the range budget |
| Ray launched into a boundary | source on the surface/floor, pointing at it | reflects before the first step |

At a vertex the travel direction is not determined by `sin θ` (it is zero); it is
determined by the gradient — the ray curves toward the **lower** sound speed.
Getting this wrong is the classic way a ray tracer stalls at a layer interface,
so `find_layer()` takes an explicit direction argument and the test suite pins
the behaviour.

---

## 3. Deep sound channel analysis

The axis is the depth of minimum sound speed. A minimum at either endpoint of
the profile is **not** a duct — there is no wall to refract energy back — and
`analyze_sofar()` reports `found = false` for those.

Given an axis speed `c_axis` and a limiting speed
`c_limit = min(c_surface, c_bottom)`, a ray launched on the axis at angle `θ`
turns before escaping when `c_axis / cos θ < c_limit`. Hence the trapping cone:

```
θ_max = arccos( c_axis / c_limit )
```

The **conjugate depths** are where `c(z) = c_limit` above and below the axis.
Energy launched inside the cone never leaves that depth band.

For the default Munk profile this gives an axis at 1300 m, `c_limit = 1548.52
m/s` (surface-limited), a cone of `±14.381°` and conjugate depths of
`0 m … 4800 m` — all reproduced by `examples/munk_simulation.cpp`, and all
cross-checked against the ray tracer in `sofar_trapping_angle_predicts_ray_behaviour`.

---

## 4. Ensonification and shadow zones

A shadow zone is the absence of *every* ray, so it cannot be computed from one
ray. `coverage_mark()` stamps the range–depth cells each ray polyline crosses,
sampling each segment at half-cell resolution so no crossed cell is skipped.

**The reported shadow fraction is a property of your ray fan as much as of the
ocean.** With too few rays, gaps between adjacent rays read as shadow. The
number is only meaningful once it stops moving as you add rays — see
`docs/validation.md` for the convergence data.

This is a purely geometric treatment: no transmission loss, no absorption, no
caustic amplitude. Those are v0.3 work (see `docs/roadmap.md`).

---

## References

- Jensen, Kuperman, Porter, Schmidt, *Computational Ocean Acoustics*, 2nd ed.,
  Springer 2011 — Ch. 3 for ray theory and the constant-gradient arc solution.
- Urick, *Principles of Underwater Sound*, 3rd ed., McGraw-Hill 1983.
- Medwin, "Speed of sound in water: A simple equation for realistic parameters",
  *JASA* 58(6), 1975.
- Mackenzie, "Nine-term equation for sound speed in the oceans", *JASA* 70(3), 1981.
- Chen & Millero, "Speed of sound in seawater at high pressures", *JASA* 62(5), 1977.
- Leroy & Parthiot, "Depth-pressure relationships in the oceans and seas",
  *JASA* 103(3), 1998.
- Munk, "Sound channel in an exponentially stratified ocean, with application to
  SOFAR", *JASA* 55(2), 1974.
- Porter, *The BELLHOP Manual and User's Guide*, HLS Research, 2011.

---

## 5. Active sonar pulse synthesis

Phase is evaluated in closed form at every sample rather than accumulated, so
there is no phase drift to integrate over a long pulse.

### 5.1 CW

```
f(t) = f0
phi(t) = 2 pi f0 t
```

Zero bandwidth, so **no pulse compression and no range resolution**: its matched
filter output is a triangle as wide as the pulse. Excellent Doppler resolution,
which is the trade.

### 5.2 LFM — linear FM

```
f(t)   = f0 + mu t,           mu = (f1 - f0) / T
phi(t) = 2 pi (f0 t + mu t^2 / 2)
```

Range resolution `1/B` regardless of pulse length. The cost is range-Doppler
coupling — see §6.2.

### 5.3 HFM — hyperbolic FM

The *period* is linear in time, hence the name:

```
1/f(t) = 1/f0 + k t,          k = (1/f1 - 1/f0) / T
f(t)   = 1 / (1/f0 + k t)
phi(t) = (2 pi / k) ln(1 + k f0 t)
```

As `k -> 0` this must fall back to the CW form: `ln(1+x) ~ x` gives
`phi -> 2 pi f0 t`, but the coded expression divides by `k`, so the
implementation switches to the linear branch below `|k f0 T| < 1e-9`.

**Why HFM matters underwater.** A time-scaled HFM is, to first order, a *delayed*
HFM — the scaling maps the family onto itself. So the matched filter peak
survives Doppler almost intact. Measured over a 60 ms 8-20 kHz sweep:

| closing speed | LFM peak loss | HFM peak loss |
|---|---|---|
| 5 m/s | −5.7 dB | −0.05 dB |
| 15 m/s | −9.9 dB | −0.14 dB |
| 30 m/s | −13.1 dB | −0.30 dB |

---

## 6. Matched filtering

Correlating a **real** received signal against a **complex** (analytic) replica
yields the complex envelope directly, with no separate basebanding stage:

```
y[n] = sum_k x[n+k] conj(r[k])
```

Implemented as FFT overlap-save. Circular aliasing confines the valid outputs to
lags `[0, M-L]`, which is `M-L+1` usable lags per block and exactly the block
advance.

### 6.1 What the peak means

For `x[k] = A cos(phi_k)` and `r[k] = a_k exp(j phi_k)`:

```
y = sum a_k^2 ( (1 + cos 2phi_k)/2  -  j sin(2phi_k)/2 )  ->  A E / 2
```

the double-frequency terms averaging away. Hence:

- **peak magnitude** `= A E / 2`, which inverts to give the received amplitude
- **coherent gain** `= L / 2` — the factor of two is the half of a real signal's
  energy in the negative-frequency half-plane that an analytic replica discards
- **compressed width** `~ 1/B`, measured at 83 µs for `B = 12 kHz`

### 6.2 Range-Doppler coupling — the wideband form

The textbook narrowband result is `dt = -f_doppler / mu` with the Doppler shift
taken at the **centre** frequency. **That formula is wrong underwater.**

Doppler here is a time *scaling*, not a frequency shift. With `alpha = 1 + v/c`:

```
phi_rx(t) = 2 pi (alpha f0 t + mu alpha^2 t^2 / 2)
```

so the received chirp rate is `mu alpha^2` and the received duration `T/alpha`.
Writing the residual phase against the replica and minimising it over the
overlap — the least-squares linear fit to `t^2` on `[0,T']` has slope `T'` —
gives

```
dt = -(v/c) * f_end / mu,        f_end = f0 + mu T
```

the **end** of the sweep, not its centre. For an 8-16 kHz upsweep the two differ
by 33%, and the error does not vanish at low speed. Verified against the tracer
for upsweeps, downsweeps and narrow sweeps in
`matched_filter_lfm_range_doppler_coupling_matches_theory`.

The derivation assumes the peak survives; above roughly 5 dB of Doppler loss the
peak is too degraded for the expression to mean anything.

---

## 7. Detection and estimation

### 7.1 CA-CFAR

Threshold = `alpha` times the mean power of the training cells either side of a
guard band. For exponentially distributed output power,

```
alpha = N (Pfa^(-1/N) - 1)
```

for `N` training cells.

**The guard band must exceed the width of the response being detected.**
Otherwise the target leaks into its own training cells, raises its own
threshold, and is never detected — the stronger the pulse, the higher the bar it
must clear. A chirp compresses to about `fs/B` cells; a CW does not compress at
all, so its response is as wide as the pulse. `suggested_cfar_guard()` returns
the longest replica length for this reason. See
`analyzer_cfar_guard_must_clear_the_response`, which reproduces the silent
failure and its fix.

### 7.2 Arrival-time estimation and the two Cramér-Rao bounds

For delay estimation in white Gaussian noise of per-sample variance `sigma^2`,
with delay in samples,

```
var(tau) >= sigma^2 / F,     F = sum_n (ds/dn)^2 = (A^2/2) sum_n (2 pi f_n / fs)^2
```

**Which bound applies depends on the receiver**, and confusing the two makes an
efficient estimator look four times worse than it is:

| bound | frequencies measured about | achievable by |
|---|---|---|
| **coherent** | zero — driven by `f_rms` | a receiver tracking absolute carrier phase |
| **envelope** | the centre frequency — driven by RMS *bandwidth* | a magnitude detector |

For an 8-20 kHz sweep these differ by `f_rms / B_rms = 14.4 / 3.46 = 4.16` in
standard deviation. This library is an envelope detector, so the envelope bound
is the relevant one — and carrier-coherent processing is not merely
unimplemented, it needs absolute phase, which Doppler destroys.

Measured across 30 dB of SNR, the estimator sits at **0.93–1.02 ×** the envelope
bound.

---

## References (added for §5-7)

- Van Trees, *Detection, Estimation, and Modulation Theory, Part III*, Wiley 1971
  — ambiguity functions and the delay-estimation bounds.
- Kay, *Fundamentals of Statistical Signal Processing: Estimation Theory*,
  Prentice Hall 1993 — Cramér-Rao bounds, Ch. 3.
- Richards, *Fundamentals of Radar Signal Processing*, 2nd ed., McGraw-Hill 2014
  — pulse compression, CA-CFAR.
- Kroszczynski, "Pulse compression by means of linear-period modulation",
  *Proc. IEEE* 57(7), 1969 — the HFM/Doppler-invariance result.

---

## 8. Doppler banks

A zero-Doppler template detects a moving target but does not *match* it. A bank
replicates each waveform across Doppler bins; how many bins that takes is a
property of the waveform, and the three families differ by orders of magnitude.

Writing `delta = v/c`, the mismatch loss at the midpoint between bins sets the
spacing. Per family:

### 8.1 CW

The matched filter output against a Doppler-shifted CW is `sinc(delta f0 T)`.
Expanding to second order, the `|delta|` giving a loss of `L` dB is

```
delta = sqrt(6 (1 - a)) / (pi f0 T),        a = 10^(-L/20)
```

### 8.2 LFM

A mismatched chirp leaves a residual quadratic phase `B t²` with `B = mu delta`.
The detector absorbs the least-squares *linear* part of that as the range bias
of §6.2; what remains is the deviation of `t²` from its own best-fit line, whose
RMS over `[0,T]` is `T²/sqrt(180)`. So

```
sigma_phi = 2 pi TB delta / sqrt(180)
loss_dB   = 4.343 sigma_phi^2
```

which makes `delta_1dB` almost exactly `1/TB`.

**Validity.** This is a small-mismatch expansion. At 1 dB it predicts 5.49 dB
where 5.73 dB is measured; at 15 m/s on the same waveform it predicts 49 dB
where 9.9 dB is measured, because the true loss saturates while the expression
runs away. It is a bin-spacing tool, not a model of deep mismatch.

### 8.3 HFM

An HFM under time scaling is *exactly* a delayed HFM. Writing
`u = k f0 t` and `v = k f0 tau`, the phase difference against the replica is

```
(2 pi / k) ln[ (1 + alpha u) / (1 + u - v) ]
```

which is constant — a pure phase offset — when `v = (alpha - 1)/alpha`. So there
is no phase residual at all, only the duration mismatch, an amplitude factor of
`(1 - |delta|)`.

The delay that constant corresponds to is large: `tau = (alpha-1)/(alpha k f0)`,
which is **-3.3 ms** at `delta = 0.109` for an 8-20 kHz 20 ms sweep. An HFM under
Doppler moves a long way before it fades.

### 8.4 The consequence

Over ±20 m/s at 1 dB straddling loss, for the same 8-20 kHz 20 ms sweep:

| waveform | bins needed |
|---|---|
| CW | 14 |
| LFM | 5 |
| **HFM** | **2** |

That ratio is the design argument for hyperbolic sweeps stated as a number, and
it is a direct cost: each bin is another correlation per block and another
`FftSize * sizeof(Complex)` of replica spectrum.

---

## 9. Echo synthesis

### 9.1 Delay and Doppler are both two-way

```
dt    = 2 dr / c
alpha = (c + v) / (c - v)  ~  1 + 2 v/c
```

The exact form is used rather than `1 + 2v/c`: at 30 m/s they differ by 0.4%,
which is four samples of delay across a 20 ms pulse — enough for a matched
filter to notice.

The factor of two in each is the single most likely place for a sign or scale
error, so both are verified by a closed loop: an echo synthesised at a stated
range and velocity is fed back through the analyser and must come out with the
range and velocity it went in with.

### 9.2 Anti-phase cancellation

For a residual timing error `tau` between the cancelling pulse and the echo,
the residual power over a flat band of width `B` centred on `f_c` is

```
G^2 = 2 - 2 cos(2 pi f_c tau) sinc(B tau)
```

For any real sonar `f_c >> B`, so the cosine dominates and **bandwidth is
almost irrelevant** — widening a 12 kHz-centred band from 0 to 12 kHz with a
2 µs error moves the residual by 0.3 dB. The limit is timing:

```
tau_max = asin(10^(G/20) / 2) / (pi f_c)
```

| cancellation | timing budget | path budget at 1500 m/s |
|---|---|---|
| −6 dB | 6.7 µs | 10.1 mm |
| −10 dB | 4.2 µs | 6.3 mm |
| −20 dB | 1.3 µs | 2.0 mm |
| −30 dB | 0.42 µs | 0.63 mm |

And past a quarter period the canceller **adds** energy — up to +6 dB.

Two effects this omits, both of which make the real figure worse: hull
scattering is distributed over many wavelengths, so a point projector cannot
match phase at more than one bearing; and a canceller cannot emit the inverse of
a sample it has not yet received, so `tau` has a floor set by acquisition and
processing latency, orders of magnitude above 1.3 µs.

**Conclusion.** Generating false targets is the achievable countermeasure.
Cancelling the real one is not, and the library says so with numbers rather than
implementing something that would not work.

---

## 10. Eigenrays and transmission loss

### 10.1 The search

An eigenray is a launch angle whose ray reaches a given receiver. Because the
tracer honours a range budget exactly (§2.5), setting that budget to the
receiver range makes the ray's final point *the arrival* — there is no
interpolation and none of the chord-versus-arc error that resampling a polyline
would introduce. What remains is a root find on

```
f(theta_0) = z_final(theta_0) - z_receiver
```

A fan brackets the sign changes; bisection polishes each one. Two eigenrays
closer in launch angle than one fan step are seen as one, so the fan resolution
decides completeness, not accuracy.

**A root that cannot be polished to inside `depth_tolerance_m` is discarded.**
That makes the tolerance a completeness knob rather than an accuracy one: set
tighter than the build can achieve and paths vanish silently. Measured on a
200 m duct at 3 km with 14 paths present:

| tolerance | double | float |
|---|---|---|
| 0.001 m | 14 | 0 |
| 0.010 m | 14 | 3 |
| 0.100 m | 14 | 14 |

### 10.2 Arc length

Absorption is quoted per unit *path*, not per unit range, so the tracer
accumulates arc length exactly alongside range and time. On a constant-gradient
segment the ray is a circular arc of radius `R = 1/(xi|g|)`, so

```
ds = R |dtheta| = |theta_1 - theta_2| / (xi |g|)
```

and on an isovelocity segment `ds = |dz| / |sin theta|`. No approximation.

### 10.3 Geometric spreading from the ray tube

Conserve power in the tube between rays launched at `theta_0` and
`theta_0 + dtheta_0`. A point source radiates

```
P cos(theta_0) dtheta_0 / 2
```

into that element. At range `r` the tube crosses an annulus of circumference
`2 pi r`, and its cross-section perpendicular to the ray is `dz cos(theta_rcv)`
— the cosine matters, and leaving it out costs `10 log10(cos theta)`. So

```
I = P cos(theta_0) / (4 pi r cos(theta_rcv) |dz/dtheta_0|)
```

and relative to the 1 m reference `I_1m = P/(4 pi)`,

```
TL = -10 log10[ c_rcv cos(theta_0) / (c_src * r * cos(theta_rcv) * |dz/dtheta_0|) ]
```

**Check.** In isovelocity water a ray is a straight line, so
`dz/dtheta_0 = r/cos^2(theta_0)` and `theta_rcv = theta_0`, giving

```
TL = -10 log10( cos^2(theta_0) / r^2 ) = 20 log10( r / cos theta_0 ) = 20 log10 R
```

exactly spherical spreading. The test suite reproduces this to **0.0000 dB**
across nine geometries from 1 to 20 km, with the Jacobian obtained by finite
difference through the tracer — so the agreement exercises the whole chain, not
just the formula.

### 10.4 Caustics

Where neighbouring rays cross, `dz/dtheta_0 -> 0` and the formula predicts
infinite intensity. That is a **failure of ray theory**, not a property of the
ocean: the correct treatment needs a wave solution through the caustic. The
library flags such paths and refuses to report a level for them, rather than
returning a number that looks like an answer.

### 10.5 Absorption — Thorp (1967)

```
a(f) = 0.11 f^2/(1 + f^2) + 44 f^2/(4100 + f^2) + 2.75e-4 f^2 + 0.003    [dB/km, f in kHz]
```

Boric acid relaxation, magnesium sulphate relaxation, pure-water viscosity, and
a low-frequency floor. A fit at about 4 °C with no temperature or depth
dependence — good to roughly 10%.

| f | a (dB/km) | over 3 km |
|---|---|---|
| 1 kHz | 0.069 | 0.21 dB |
| 10 kHz | 1.19 | 3.6 dB |
| 100 kHz | 34.1 | 102 dB |

At 100 kHz absorption over 3 km exceeds the spreading loss by 30 dB. That single
comparison is why long-range sonar is low-frequency, and why a tank experiment
at 200 kHz can ignore absorption entirely over a metre.

### 10.6 Multipath echoes

For a monostatic geometry the echo returns along the path the ping arrived on,
so each eigenray contributes one arrival at `2t_i` with two-way loss
`2 TL_i`. Delays are referenced to the earliest path, levels to the strongest.

The scale that matters operationally is the multipath spread against the pulse
length: comparable and the paths smear one echo, much larger and they resolve as
separate targets. A 300 m duct at 2.5 km with a ±10° fan spreads arrivals over
70 ms, which a 20 ms pulse resolves; the same geometry with a ±25° fan spreads
them over 400 ms.

---

## 11. Boundary reflection

Until v0.5 every boundary reflected perfectly, which made bounced paths louder
than they are. In shallow water most paths bounce, so this was the largest
remaining overstatement in the transmission loss.

### 11.1 Bottom — the Rayleigh coefficient

A fluid-fluid interface, water (`rho1, c1`) over sediment (`rho2, c2`). With
grazing angles measured from the interface and Snell's law
`cos(th1)/c1 = cos(th2)/c2`:

```
nu = c2/c1,   m = rho2/rho1
sin(th2) = sqrt(1 - nu^2 cos^2(th1))
R = (m nu sin(th1) - sin(th2)) / (m nu sin(th1) + sin(th2))
```

**Three closed-form limits pin this expression.** Any one could be reproduced by
a wrong formula; all three together could not.

| limit | value | why |
|---|---|---|
| normal incidence `th1 = 90°` | `(Z2 - Z1)/(Z2 + Z1)`, `Z = rho c` | the textbook impedance ratio |
| `nu cos(th1) = 1` | `sin(th2) = 0`, so `\|R\| = 1` | the critical angle |
| below critical, lossless | `\|R\| = 1` exactly | total internal reflection |

**The critical angle is the number that matters:**

```
th_c = arccos(c1 / c2)
```

24.62° for medium sand at 1650 m/s. Rays grazing shallower than this are
trapped by the bottom; steeper ones leak into the sediment. It is what decides
which paths survive to long range in shallow water, and it depends on the speed
ratio alone — not on density.

A bottom *slower* than the water has no critical angle at all and leaks at every
angle, which is why mud is acoustically far worse than sand.

### 11.2 Sediment attenuation

Published as `alpha` dB per wavelength, which enters as a complex sediment
speed. From `alpha_dB/lambda = 20 log10(e) * alpha_np * lambda`:

```
delta = alpha / (40 pi log10 e) = alpha / 54.575
c2' = c2 / (1 - i delta)
```

`R` is then complex and `|R| < 1` at **every** angle, including below critical.
That matters: with a lossless bottom, sub-critical rays are trapped forever and
shallow-water range is unbounded, which is not what the ocean does. At half the
critical angle a default sand bottom costs 0.81 dB per bounce — 8 dB over ten.

**Branch choice.** The complex square root has two branches and only one gives a
wave that decays into the sediment. Rather than reason about sign conventions,
the implementation takes the branch that conserves energy: the other yields
`|R| > 1`, which no passive interface can do. A 9100-combination sweep over
speed, density, attenuation and angle confirms `max |R| = 1.000000000000`.

### 11.3 Surface — coherent scattering loss

A flat pressure-release surface reflects perfectly. A rough one scatters the
coherent (specular) component away:

```
|R| = exp(-G^2 / 2),   G = 2 k sigma sin(theta),   k = 2 pi f / c
```

`G` is the Rayleigh roughness parameter and `sigma` the RMS wave height. From
wind, via the Pierson-Moskowitz fully developed sea:

```
H_1/3 = 0.0246 U^2   [U in m/s]      sigma = H_1/3 / 4
```

so 10 m/s of wind gives 0.62 m RMS.

**This is the specular component only.** The scattered energy is not destroyed;
it goes into a diffuse field that a ray model does not represent. At 0.5 m seas,
10 kHz and 20° grazing the formula gives a 700 dB loss, which is arithmetic
rather than physics — the path is not 700 dB down, its *specular* part is. The
library therefore caps the reported surface loss (30 dB by default) and says
why, rather than returning a number that would delete a path that is still in
the water.

### 11.4 Applying it along a path

The grazing angle at each bounce comes from the Snell invariant rather than from
having tracked it, which is exact in a range-independent ocean:

```
cos(theta_boundary) = xi * c(z_boundary)
```

Total boundary loss is then `n_surface * L_s(theta_s) + n_bottom * L_b(theta_b)`.
For a 200 m duct at 3 km, 5 kHz, 8 m/s wind and a sand bottom, this takes paths
that all sat within 1 dB of each other and spreads them over 128 dB — the
four-bounce path is no longer a peer of the direct one.

---

## 12. Reverberation

§11.3 capped the surface scattering loss and said the energy "goes into a
diffuse field the ray model does not carry". This is that field coming back —
and at short range it, not the ambient noise, is what an active sonar competes
against.

### 12.1 Scattering strengths

**Lambert's law**, bottom backscatter:

```
S_s = mu + 10 log10(sin^2 theta) = mu + 20 log10(sin theta)
```

`mu` about −27 dB for many sediments. The `sin^2` is the projected area of a
facet entering twice, once on transmit and once on receive. Angle-only: it
carries no frequency dependence, and it fails near normal incidence where a
smooth bottom gives a specular return far above Lambert.

**Chapman-Harris (1962)**, wind-driven surface backscatter:

```
beta = 158 (v f^(1/3))^(-0.58)                     v in knots, f in Hz
S_s  = 3.3 beta log10(theta/30) - 42.4 log10(beta) + 2.6      theta in degrees
```

The `log10(theta/30)` form means the angular term vanishes at 30°, so the value
there depends only on wind and frequency — a structural property the test suite
checks. **The coefficients have not been verified against the original paper**;
see `docs/validation.md §10`.

### 12.2 Ensonified extent

```
boundary:  A = r phi (c tau / 2)          annulus segment
volume:    V = r^2 Omega (c tau / 2)      shell
```

The `c tau / 2` is the **range resolution**, and that is where pulse compression
enters: a chirp's effective `tau` is `1/B`, not its length. A 20 ms pulse with
12 kHz of bandwidth shrinks the ensonified area by its time-bandwidth product of
240 — **24 dB straight off the reverberation**.

### 12.3 The two range laws

```
RL = SL - 2 TL + S_s + 10 log10(extent)
```

With spherical spreading `TL = 20 log10(r)`:

| mechanism | extent | decay |
|---|---|---|
| boundary | `A ∝ r` | `40 - 10 = ` **30 log10(r)** |
| volume | `V ∝ r²` | `40 - 20 = ` **20 log10(r)** |

Those two exponents are the signature of which mechanism dominates. A measured
decay of 20 rather than 30 says the water column is scattering, not the
boundary.

### 12.4 The result that decides the design

Write the echo and the reverberation with the same source level and the same
two-way loss:

```
EL      = SL - 2 TL + TS
RL      = SL - 2 TL + S_s + 10 log10(A)
EL - RL = TS - S_s - 10 log10(A)
```

**`SL` and `TL` cancel exactly.** In a reverberation-limited geometry, doubling
the transmit power raises the target and the background together and buys
nothing. Verified over nine combinations of source level (160–240 dB) and
transmission loss (40–90 dB): the ratio does not move by a thousandth of a dB.

What *does* help is shrinking `A`: a tenth of the pulse length or a tenth of the
beamwidth each buy exactly 10 dB. That is the whole design argument for pulse
compression, and it is why the analyser's waveforms are chirps rather than
tones.

### 12.5 Where reverberation stops mattering

Reverberation falls as 30 log10(r); ambient noise does not fall at all. Their
crossover is the **reverberation-limited range** — inside it a bigger
transmitter is wasted, outside it the geometry is noise-limited and power helps
again. For a 200 dB source, −30 dB scattering, 0.2 rad beam and 10 ms pulse:

| ambient | crossover |
|---|---|
| 40 dB | 24.7 km |
| 60 dB | 5.3 km |
| 80 dB | 1.1 km |

### 12.6 What this does to CFAR

Reverberation decays by tens of dB across a single processing block. A fixed
threshold cannot be right anywhere: set for the near field it is deaf far out,
set for the far field it is a wall of false alarms close in. Measured on a block
whose background falls 21.6 dB, a single threshold at the block mean sits 9.5 dB
below the near field and 12.2 dB above the far field.

CA-CFAR estimates the background locally and tracks the decay, provided the
training window is short compared to it — 20/20 detections of a target buried in
the decayed region, with 3 false alarms across 20 empty blocks.

---

## 13. Line arrays, beamforming and bearing

Everything before this release was single-channel. A hydrophone knows *when* a
ping arrived and what shape it was, but not where from — and §12 showed that
reverberation scales with the ensonified area, of which the azimuthal half can
only be shrunk with a beam.

Angles are from **broadside**: 0 perpendicular to the array, ±π/2 endfire. That
makes `sin θ` the natural variable, and every formula below is written in it.

### 13.1 The array factor

For a uniform line array of `N` elements at spacing `d`, steered to `θ₀`:

```
psi = k d (sin(theta) - sin(theta_0)),      k = 2 pi / lambda
B   = | sin(N psi / 2) / (N sin(psi / 2)) |
```

Its closed forms are what the test suite checks against, because a beam pattern
that is subtly wrong still looks like a beam:

| property | value |
|---|---|
| peak | `B = 1` at `psi = 0`, exactly |
| nulls | `psi = 2 pi m / N`, i.e. `sin(theta) - sin(theta_0) = m lambda / (N d)` |
| first sidelobe | **−13.26 dB** as `N -> inf` (−12.80 at N=8, −13.26 at N=128) |
| −3 dB width | `0.886 lambda / (N d cos theta_0)` |

The `cos θ₀` in the beamwidth is the projected aperture: a beam steered to 60°
is twice as broad as the same beam at broadside, and the exact null position
follows from `asin(sin θ₀ + λ/(Nd)) − θ₀` rather than the small-angle form.

**Grating lobes** appear where `psi = ±2π`, i.e. `sin θ = sin θ₀ ± λ/d`, which
is real only when that lands inside `[−1, 1]`. Hence

```
d <= lambda / (1 + |sin(theta_max)|)
```

and steering to endfire demands the familiar `d ≤ λ/2`.

### 13.2 Array gain

`AG = 10 log₁₀(N)` against spatially white noise: the signal adds coherently
(amplitude ×N) and the noise does not (×√N), so power SNR improves by N.
Measured by Monte Carlo through the actual beamformer: 5.92 / 12.07 / 18.36 dB
for N = 4 / 16 / 64, against 6.02 / 12.04 / 18.06 predicted.

This assumes the noise is uncorrelated between elements. Isotropic ambient
roughly is at `d ≥ λ/2`; directional interference emphatically is not.

### 13.3 The bearing bound — the spatial twin of §7.2

Exactly the same derivation as the arrival-time bound, with the element index
in place of time:

```
var(theta) >= 6 / (rho (k d cos theta)^2 N (N^2 - 1))
```

where `rho` is the per-element SNR in power. The `N(N²−1)/12` is `Σ n'²` about
the array centre — the spatial counterpart of the waveform's mean-square
bandwidth.

Two consequences read straight off it:

- **`N^(-3/2)`, not `N^(-1/2)`.** More elements buy both more signal *and* more
  aperture. Doubling N improves the bound by `2^1.5 = 2.828`; measured 2.845 /
  2.833 / 2.830 across N = 8 → 64.
- **`1/cos θ`.** A target at endfire is far harder to place than one at
  broadside, because the projected aperture shrinks to nothing.

A 32-element half-wavelength array at 0 dB element SNR is bounded at 0.247°
against a 3.17° beamwidth — **one thirteenth of a beam**. Bearing comes from
phase across the aperture, not from which beam lit up.

The conventional beamformer with parabolic peak refinement measures **1.13 /
1.02 / 0.99 ×** that bound across 14 dB of SNR.

### 13.4 Split-beam

Beamform the two halves separately; the phase between them gives bearing
directly:

```
dphi = k D (sin theta - sin theta_0),     D = (N/2) d
```

`|dphi| = pi` exactly at the first null, so split-beam is unambiguous precisely
out to the mainlobe edge and wraps beyond it. It refines a bearing already known
to within the beam; it cannot find one outside it. Noiseless accuracy is exact
to 1e-16.

### 13.5 The bandwidth limit phase steering carries

A phase-steered array approximates a delay-steered one only while the signal
stays correlated across the array's traversal time `(N−1) d sin θ / c`. Taking a
tenth of its inverse as the usable bandwidth:

| steer | traversal | usable BW |
|---|---|---|
| 0° | 0 | unlimited |
| 15° | 0.40 ms | 249 Hz |
| 45° | 1.10 ms | 91 Hz |

The library's own waveforms are 12 kHz chirps. A phase-steered array **cannot**
handle them off broadside; a real receiver beamforms per frequency bin or steers
with delays. Saying so is cheaper than shipping a beamformer that quietly smears
the beam.

### 13.6 Where this joins §12

Ten times the elements is ten times the aperture is a tenth of the beamwidth is
**exactly 10 dB** of echo-to-reverberation ratio — the same 10 dB per decade
§12.4 charges for the ensonified area. The array also buys 10 dB of gain against
isotropic noise over the same change. Two different mechanisms, and here they
happen to agree.

---

## 14. Wideband and adaptive beamforming

§13 built an array and then measured why it could not be used: phase steering
holds to 91 Hz of bandwidth at 45°, and this library transmits 12 kHz chirps.

### 14.1 Per-bin steering

A phase shift that is wrong across a band is exactly right *within* one bin. So
transform, steer bin by bin, transform back:

```
Y[k] = sum_n  w_n  X_n[k]  exp(-j 2 pi f_k tau_n),     tau_n = x_n sin(theta) / c
```

This is an exact fractional delay, not an interpolation of one. Two details
decide whether it works:

- **Bins above M/2 carry negative frequencies**, `f_k = (k−M)·fs/M`. Steering
  them with the positive value breaks the spectrum's Hermitian symmetry and the
  inverse transform returns noise that still looks like a signal.
- **The forward model's sign must match the steering convention.** The element
  at `+x` sees the wavefront *earlier*, so its time series is further into the
  waveform: `t = (i − at)/fs + x sin θ / c`. Getting this backwards steers every
  beam to the mirror bearing — which looks like a working beamformer until you
  check where it points.

The test suite pins the second one structurally: a narrowband tone pushed
through the wideband path must reproduce `array_factor` in closed form, which a
mirrored convention cannot do.

**Cost.** Transforming every element once costs `N·M` complex values;
re-transforming per beam costs nothing extra in memory but runs `N` FFTs per
beam instead of `N` in total. For 24 elements and 61 beams that is 1464
transforms against 24, so the storage is almost always the right side of the
trade — hence the two-stage API.

### 14.2 Shading

The shaded array factor is the discrete-time Fourier transform of the window, so
its sidelobe level *is* the window's — these are the familiar numbers, not
array-specific ones:

| window | peak sidelobe | mainlobe width | gain given up |
|---|---|---|---|
| Uniform | −13.26 dB | ×1.00 | 0 dB |
| Hann | −31.5 dB | ×1.64 | −1.76 dB |
| Hamming | −42.7 dB | ×1.47 | −1.34 dB |
| Blackman | −58.1 dB | ×1.90 | −2.37 dB |

Shading loss is `10 log₁₀[ (Σw)² / (N Σw²) ]`, always negative: a window trades
gain and resolution for sidelobes, and there is no free direction.

### 14.3 MVDR

Conventional beamforming cannot resolve two sources closer than about one
null-to-peak spacing,

```
delta_theta ~ lambda / (N d cos theta)
```

**whatever the SNR** — the resolution is set by the aperture and nothing else.
MVDR is not bound by it, because it places nulls instead of scanning a fixed
pattern:

```
P(theta) = 1 / (a^H R^-1 a)
```

Solved by complex Cholesky rather than an explicit inverse: `R` is Hermitian
positive definite after loading, so `R = L L^H` and one forward/back
substitution per steering angle.

**Diagonal loading is not optional.** With fewer snapshots than elements `R` is
singular and the factorisation fails outright; even with enough snapshots it is
ill-conditioned whenever a source is strong. The library reports the failure
rather than returning noise shaped like a spectrum. Loading is expressed as a
fraction of the mean diagonal power, and it trades resolution for stability.

Measured on a 16-element half-wave array whose conventional limit is 7.16°:
**two sources 4.30° apart give one conventional peak and two MVDR peaks**, at
±2.145° against a truth of ±2.149°.

MVDR needs the sources to be mutually incoherent. Two coherent arrivals — the
multipath of §10.6, for instance — defeat it, and would need spatial smoothing
that this library does not implement.

---

## 15. Tracking

A Pulse Descriptor Word carries time, type, Doppler and — since §14 — bearing.
Nothing connected them across blocks, which is what turns a list of detections
into a picture.

### 15.1 The model

Constant velocity in Cartesian coordinates, measured in polar. Cartesian because
the dynamics are then **linear**; polar because that is what a sonar produces.
All the nonlinearity lives in the measurement Jacobian.

State `[x, y, vx, vy]` with `y` along broadside, so `bearing = atan2(x, y)` —
the same convention as §13.

```
F = [[1,0,dt,0],[0,1,0,dt],[0,0,1,0],[0,0,0,1]]

Q = q [[dt^4/4, 0, dt^3/2, 0], [0, dt^4/4, 0, dt^3/2],
       [dt^3/2, 0, dt^2,   0], [0, dt^3/2, 0, dt^2  ]],   q = sigma_a^2
```

discrete white-noise acceleration. `sigma_a` is the single knob that sets how
much manoeuvre the filter expects: too small and it lags a turn, too large and
it chases noise.

Measurement `z = [r, theta]` with `r = sqrt(x²+y²)`, `theta = atan2(x, y)`:

```
H = [[ x/r,   y/r,  0, 0],
     [ y/r^2, -x/r^2, 0, 0]]
```

The bearing innovation is wrapped to `[-pi, pi]` before use; without that, a
track sitting near ±180° produces a 2π innovation and diverges on the first
update.

### 15.2 The consistency check that matters

The normalised innovation squared,

```
d^2 = y^T S^-1 y,     S = H P H^T + R
```

is **chi-square with 2 degrees of freedom** when the filter is consistent. Its
mean must be 2 and 95% of samples must fall below 5.991.

This is the test worth running. A filter whose covariance is wrong still tracks
— it just lies about how well, and then gates correct measurements out or
accepts clutter, with nothing in its output to say so. Measured over 1980
samples: **mean 1.951**, with 95.0% under the 95% gate and 99.0% under the 99%.

The chi-square CDF with 2 dof is `1 - exp(-x/2)`, so a gate inverts in closed
form and needs no table:

```
chi2(p) = -2 ln(1 - p)      5.991 at 95%,  9.210 at 99%
```

### 15.3 What tracking is worth

| quantity | result |
|---|---|
| RMS position error, raw measurement | 50.97 m |
| RMS position error, filtered | **20.89 m** (2.44× better) |
| velocity, truth (4.00, −7.00) | estimated (4.17, −7.18) |
| closing rate, truth 7.486 m/s | tracked 7.683 m/s |

Velocity is the interesting one: a single detection cannot know it at all, and
it is what lets a track coast through a missed scan.

### 15.4 What tracking fixes, and what it does not

**False alarms: decisively.** They do not repeat, so M-of-N confirmation
removes them. Measured: **494 false alarms scattered over 300 scans produced
zero confirmed tracks.**

**Cross-template ghosts: not at all.** The roadmap for this release claimed time
consistency would finally suppress the ghosts of §7.6. **That claim was wrong.**
A ghost appears whenever the real arrival does, at a fixed offset set by the
template cross-correlation, so it is exactly as consistent over time as the
target, moves with it, and forms its own perfectly healthy confirmed track.
Nothing about its kinematics is objectionable.

Measured: one target plus its ghost gives **two** confirmed tracks. The test
asserts that outcome so the correction cannot quietly rot.

Suppressing ghosts needs the fixed offset and amplitude ratio to be recognised
as a template artefact — a different mechanism entirely. §16.2 implements it.

### 15.5 Association

*(Superseded by §16.3, which replaced greedy association with a global cost
ordering. Retained because the reasoning is still the reason a global ordering
was needed.)*

Greedy nearest-neighbour on the NIS, gated at `chi2`. Not optimal: a global
assignment does better when two targets cross. But it is
`O(tracks × measurements)` with no allocation, and its failure mode — a swap
during a crossing — is well understood rather than surprising.

`(I − KH)P` is symmetrised after each update. It is symmetric in exact
arithmetic and drifts out of it in floating point; an asymmetric covariance
eventually goes indefinite and the filter diverges with no warning at all.


## 16. Fusing what is already measured

Every quantity in this section was already being produced by an earlier stage
and then discarded at the boundary. Nothing new is measured; four existing
outputs are simply joined to the consumers that could use them.

### 16.1 Range rate as a third measurement

The Doppler bank of §8 estimates a closing rate for every detection, and §15's
filter inferred velocity from position history alone while that estimate sat
unused in the `PulseDescriptor`. Fusing it extends the measurement to
`z = [r, theta, rdot]` with

```
rdot = (x vx + y vy) / r
```

and a third Jacobian row obtained by differentiating it:

```
d(rdot)/dx  = -vx/r + (x vx + y vy) x / r^3
d(rdot)/dy  = -vy/r + (x vx + y vy) y / r^3
d(rdot)/dvx = -x/r
d(rdot)/dvy = -y/r
```

The sign convention is the bank's: **positive closing**, hence the leading
minus. Note that unlike the range and bearing rows, this one depends on the
velocity components — a range rate is the only measurement in the set that
observes the velocity state *directly* rather than through its effect on
successive positions.

**The gate must follow the dimension.** A 2-dof measurement and a 3-dof one at
the same threshold are different gate probabilities, so `TrackerConfig` carries
both and the filter picks by `has_range_rate`. The 3-dof quantile has no
closed form, so it is bisected on the exact CDF

```
F(x; 3) = erf(sqrt(x/2)) - sqrt(2x/pi) exp(-x/2)
```

giving **7.815 at 95%** and **11.345 at 99%** — the textbook values to three
decimals, computed rather than tabulated.

Consistency is preserved: over 1980 fused updates the mean NIS is **2.998**
against a theoretical 3, with 94.7% under the 95% gate and 98.9% under the 99%.

**What it is worth, honestly.** Radial velocity error, the only component a
range rate observes:

| scans | position only | with range rate | ratio |
|---|---|---|---|
| 3 | 4.213 m/s | **1.220 m/s** | 3.45× |
| 6 | 2.069 m/s | **0.649 m/s** | 3.19× |
| 20 | 0.365 m/s | 0.368 m/s | 0.99× |

The gain is large early and vanishes by twenty scans. That is not a defect: by
then the filter's own estimate (0.365 m/s) is already sharper than the σ = 2 m/s
measurement, and a measurement helps exactly as long as it beats the estimate
you already have. Where it matters is the first few scans of a new contact —
which is when a decision usually has to be made.

**An unresolved bin is not a measurement of zero.** A Doppler bank that fails to
resolve a shift reports 0 m/s, which is a completely different statement from
"the target is stationary". Feeding it in pins the track's radial velocity to
zero. Hence the separate `has_range_rate` flag, and the measurement:

```
truth closing                8.00 m/s
zero taken as a measurement  5.686 m/s   (28% low)
flag left clear              9.163 m/s
```

### 16.2 Recognising a cross-template ghost

§15.4 established that time consistency cannot suppress these, because a ghost
is *as consistent as the target*. What distinguishes it is not its kinematics
but its **origin**: it is the same arrival seen through a different matched
filter. So it

* shares the target's bearing and bearing rate (same arrival),
* sits at a fixed range offset (the template cross-correlation lag),
* is weaker (a mismatched filter loses processing gain), and
* carries a **different waveform label**.

A pair is suppressed only when all four hold. **The label check is the whole
safety argument.** Two real targets illuminated by one sonar return the *same*
waveform, so they share a label and are never paired — which is what saves a
line-astern formation, the geometry that is otherwise indistinguishable from a
ghost pair. The test suite demonstrates exactly that:

```
line-astern formation, same waveform -> 2 tracks (kept)
same geometry, different waveforms   -> 1 track (suppressed)
```

The two runs differ *only* in the label. Without that check the first would be
silently deleted, and deleting a real target is a far worse failure than
carrying a ghost.

The bearing tolerance is set from what the array delivers, not from the
geometry. Two tracks on one *true* bearing estimate it independently, so with
1° measurements they routinely differ by 1.7°; the default is 3°. Loosening it
is safe precisely because the label, not the bearing, is doing the separating.

### 16.3 Global-cost-ordered association

Greedy association processes measurements in arrival order and lets the first
one take its best track, which is why two crossing targets swap: at the crossing
the wrong measurement gets first pick.

The fix needs no assignment algorithm. Build every gating `(track, measurement)`
pair with its NIS, sort by cost, and assign best-first, skipping any pair whose
track or measurement is already taken. This is not optimal in the Hungarian
sense — it is still greedy, but greedy over the *global* cost ordering rather
than over arrival order, and that is what the crossing case actually needs. Cost
is `O(TM log TM)` with a fixed-size array and an insertion sort, no allocation.

Measured on two targets crossing at 8 m/s:

```
before crossing: left=1 right=2
after  crossing: left=2 right=1
```

The identities survive. Note also what did *not* need changing: the gate. A
swap is not a gating failure — both measurements gate perfectly well against
both tracks near the crossing, which is exactly why the choice among them has to
be made globally.

### 16.4 Forward-backward spatial smoothing

MVDR (§14) fails on **coherent** sources, and the multipath this library
produces in §10 is precisely that: one arrival reaching the array by two paths,
perfectly correlated. The covariance is then rank-deficient in a way diagonal
loading cannot repair, and the adaptive beamformer nulls the target along with
its own echo.

The remedy averages the covariances of overlapping subarrays, which randomises
the relative phase between the coherent pair and restores rank. Forward-backward
adds the exchange-reversed conjugate, doubling the effective snapshot count for
free on a uniform line array:

```
R_smooth = (1 / 2K) sum_k [ R_k + J conj(R_k) J ]
```

with `J` the exchange matrix and `K = N - L + 1` subarrays of `L` elements.
Elementwise, using the Hermitian symmetry of `R_k`,

```
R_fb[i][j] = R_f[i][j] + R_f[L-1-j][L-1-i]
```

which pairs each entry with exactly one other and gives both the same value.
That the result is **persymmetric** is the check that the indices are the right
way round; getting them transposed produces a matrix that still looks like a
covariance and puts every bearing in the wrong place.

Measured on two coherent arrivals at ±6°:

| | peaks found |
|---|---|
| plain MVDR, 16 elements | 4 (spurious) |
| forward-backward, 10-element subarrays | **2, at −6.00° and +6.00°** |

**The cost is aperture, and it is not small.** An `N`-element array smoothed
with subarrays of `L` has the resolution of an `L`-element array: 7.16° falls
to 11.46° in the case above. Resolving `P` coherent sources needs `L > P`, so
the caller must choose `L` — the library will not choose it for them.


## 17. Real data, and a coefficient that was wrong for eleven releases

### 17.1 Checking against a published value instead of against yourself

Until v0.11 the sound-speed equations were verified by **mutual agreement**:
Medwin, Mackenzie and Chen-Millero must agree inside their common validity box,
so a mistyped coefficient breaks the agreement. That argument is weaker than it
looks. The three equations agree only to about 0.1 m/s, so the check can only
see errors larger than that — and it had been hiding one worth 0.016 m/s since
v0.1.

The primary source settles it. UNESCO Technical Papers in Marine Science 44
(Fofonoff & Millard 1983) publishes both a check value, on p. 48,

```
Check Value: 1731.995 m/s for S = 40, t = 40 C, p = 10000 decibars
```

and a full table on p. 50: **220 values** over S = 25/30/35/40 PSU,
T = 0/10/20/30/40 °C (IPTS-68) and p = 0 … 10000 dbar. Both are now in the test
suite, transcribed in `tests/data/unesco44_table.hpp`.

### 17.2 What the bug was

`chen_millero()` carried Wong & Zhu's (1995) ITS-90 coefficients for 40 of its
42 terms, and Chen & Millero's original 1977 values for the other two:

| | 1977 (IPTS-68) | Wong & Zhu (ITS-90) | what the code had |
|---|---|---|---|
| `A02` | 7.164e-5 | 7.166e-5 | **7.164e-5** |
| `A03` | 2.006e-6 | 2.008e-6 | **2.006e-6** |

So it was neither equation, and its header credited Chen & Millero (1977) for
something that was mostly Wong & Zhu. The fix was not to pick a number but to
implement **both** equations properly:

* `chen_millero_1977()` — IPTS-68, the version the published check value and
  table belong to. It reproduces the check value as **1731.9954** against
  1731.995, and all 220 table entries to a worst error of **0.0499 m/s** — the
  table is printed to 0.1 m/s, so 0.05 is exactly its rounding half-width. The
  agreement is as tight as the published data can express.
* `chen_millero_its90()` — Wong & Zhu (1995), coefficients per NPL. This is the
  library default, because every measurement taken since 1990 is on ITS-90.

That the two are *legitimately different equations* rather than one equation and
one mistake is itself tested: feeding the same number to both differs by
0.0208 m/s, and converting the temperature scale first (`t68 = 1.00024 t90`)
drops that to **0.0056 m/s**. Wong & Zhu state their revision changes the UNESCO
equation by "within 0.024 m/s"; the measured 0.0208 sits inside that.

### 17.3 Del Grosso, and where a validity box lies

Del Grosso (1974), in Wong & Zhu's ITS-90 form, is a genuinely independent
equation — different laboratory, different data, different functional form.
Comparing it to UNESCO says something about the physics in a way that comparing
Medwin to Mackenzie does not.

The result depends almost entirely on **where you compare them**:

| region | worst disagreement |
|---|---|
| 0 – 1000 m (most sonar work) | 0.41 m/s |
| 0 – 3000 m | 0.76 m/s |
| 0 – 6000 m | 0.82 m/s |
| 0 – 9810 m (Del Grosso's limit) | 1.33 m/s |
| full nominal validity box | **3.93 m/s** |

The last row is 26 °C water under 1000 bar. No such sea exists, neither equation
was ever fitted to a sample of one, and they extrapolate apart. A validity box
is a rectangle in (S, T, P) and the ocean is not a rectangle — quoting a
disagreement over the whole box measures the corners, not the water.

Worth noting against Chen & Millero's own quoted standard deviation of
0.19 m/s: in the top kilometre two independent equations differ by about twice
the uncertainty either claims. That is the real accuracy of "the speed of sound
in seawater".

### 17.4 Real profiles

Six World Ocean Atlas 2023 profiles (`data/profiles/`, fetched by
`tools/fetch_woa_profiles.py` from NOAA NCEI over OPeNDAP, one grid cell each):

| site | max depth | c surface | c axis | axis depth |
|---|---|---|---|---|
| aegean | 650 m | 1522.1 | 1514.4 | 175 m |
| black-sea | 2200 m | 1486.9 | 1462.6 | 55 m |
| eq-pacific | 4300 m | 1539.1 | 1484.9 | 950 m |
| levantine | 2400 m | 1530.6 | 1516.1 | 375 m |
| n-atlantic | 4700 m | 1509.0 | 1486.3 | 950 m |
| norwegian | 3200 m | 1477.4 | 1463.7 | 850 m |

These are a **climatology**, not casts: a decadal average over a 1° cell, with no
eddy, internal wave or diurnal layer surviving the averaging. They have the
*shape* of real water, which is what the tests need, and they should not be used
to predict a real range.

The Black Sea earns its place. Its surface is 18.2 PSU — near-fresh — so
assuming the usual 35 PSU makes a **20.3 m/s** error there against at most
4.7 m/s at any other site. For scale, the entire channel excess measured at that
site is 24 m/s: the error is the size of the feature, so the channel you would
trace is not the channel that is there.

### 17.5 What real profiles test that Munk cannot

Munk is analytic and smooth. A real profile is piecewise linear between standard
levels with a kink at every one. Two checks run over all six:

**Snell's invariant** `ξ = cos θ / c` comes out at exactly zero drift — which
should be *suspicious rather than reassuring*. The tracer carries ξ and
reconstructs θ from it as `θ = ±acos(ξc)`, so measuring `cos θ / c` largely
measures whether `acos` and `cos` round-trip. It catches a ray that jumps layers
wrongly; it cannot catch a ray that turns in the wrong *place*.

**Turning depth** can, and is independent of how θ is stored. At a turning point
θ = 0, so the local sound speed must equal `c_source / cos θ₀`. Over 186 turning
points across the six profiles the worst error is **below 1e-9 m/s**: the tracer
cuts the arc exactly at the turn rather than stepping past it and interpolating.

## 18. Interacting multiple models

A constant-velocity filter has two ways to handle a turn, both bad. Quiet
process noise: it lags the manoeuvre, the innovations blow past the gate, and
the track dies on a target that is plainly still there. Loud process noise: it
follows the turn, and spends the other 95% of the time chasing measurement noise.

The IMM refuses the choice. Three models share the four-state Cartesian vector:

```
0   constant velocity,       sigma_a = 0.3 m/s^2
1   coordinated turn at +w,  sigma_a = 2.0 m/s^2
2   coordinated turn at -w,  sigma_a = 2.0 m/s^2
```

A coordinated turn at *known* rate is linear in that same state, which is why no
augmentation is needed:

```
        [ 1  0   sin(wT)/w   -(1-cos(wT))/w ]
F(w,T) = [ 0  1  (1-cos(wT))/w   sin(wT)/w  ]
        [ 0  0    cos(wT)       -sin(wT)    ]
        [ 0  0    sin(wT)        cos(wT)    ]
```

At `w = 0` this is exactly the constant-velocity matrix, but `sin(wT)/w` is 0/0
there, so below `|wT| < 1e-4` the implementation switches to the series. That is
not an approximation to the answer; it is the more accurate way to compute it.
A test walks `w` down through the crossover and requires monotone convergence —
in `double` the difference falls as 6.946e-N with no jump at the threshold.

### 18.1 The four steps

**Mixing.** Before each model predicts, its state is replaced by a blend of all
models' states, weighted by `μ_{i|j} = π_ij μ_i / c_j`. This is what makes an
IMM more than a bank of independent filters: a model that has been idle starts
from somewhere sensible when it takes over. The mixed covariance carries a
spread-of-means term,

```
P0_j = sum_i mu_{i|j} [ P_i + (x_i - x0_j)(x_i - x0_j)^T ]
```

and dropping that outer product makes the filter overconfident exactly when the
models disagree — which is exactly when it should not be.

**Filtering.** Each model runs the ordinary EKF measurement update.

**Re-weighting.** `μ_j ∝ c_j Λ_j` with the Gaussian likelihood

```
L = exp(-d^2/2) / sqrt( (2 pi)^m |S| )
```

The `|S|` term is not decoration. Without it a model could win on residual
simply by being vaguer; `|S|` is what charges it for that. It comes free from
the Cholesky factor already computed for the NIS, as `|S| = (prod L_ii)^2`.

**Combination.** The same mixture formula again, producing the output estimate.

### 18.2 What it is worth, split by phase

A target turning at 3 °/s between scans 15 and 35, against a single-model EKF
tuned with process noise *between* the IMM's two:

| phase | IMM | single-model EKF | ratio |
|---|---|---|---|
| during the turn | 37.31 m | 36.68 m | **0.98×** |
| after it ends | 25.88 m | 44.53 m | **1.72×** |
| overall | 44.91 m | 50.74 m | 1.13× |

**The gain is in recovery, not in the turn.** During the manoeuvre both filters
lag and the single model's larger process noise partly covers for it — the IMM
is marginally *worse*. Afterwards the IMM switches back to a quiet model within
a scan or two while the EKF is still coasting on the velocity the turn left it
with. That asymmetry is the honest description of what an IMM buys, and it is
not the description usually given.

On a straight target it costs nothing: 27.40 m against 27.88 m.

### 18.3 The manoeuvre detector that comes free

`1 − μ_CV` is a manoeuvre probability the filter had to compute anyway for the
mixing. It peaks at **0.991** during the turn, and the probability-weighted turn
rate recovers the direction correctly (truth ±3 °/s → estimate ±2.0…2.3 °/s).

That estimate is *not* a measurement of turn rate. With models at 0 and ±ω it
can only return a value in [−ω, +ω], and a target turning harder saturates it.
Read it as "which way, and roughly how hard".


## 19. Estimating the turn rate, and the IMM in the tracker

### 19.1 Why bracketing is not enough

§18's IMM reports a probability-weighted blend of model turn rates, so its
estimate is confined to `[-ω, +ω]` by construction. Against models bracketed at
3 °/s:

| truth | IMM (bracketed) | CTRV (estimated) |
|---|---|---|
| 2 °/s | 1.27 °/s | 1.02 °/s |
| 5 °/s | 2.62 °/s | **4.99 °/s** |
| 8 °/s | 1.82 °/s | **5.64 °/s** |

Note the 8 °/s row: the IMM reports *less* turn than at 5 °/s. That is not a
bug — the measurement fits the ±3 °/s models so badly that probability drifts
back towards the constant-velocity model, and the blend collapses toward zero.
A bracket does not degrade gracefully once the truth leaves it.

### 19.2 The fifth state

```
x = [ x, y, vx, vy, omega ]
```

with the coordinated-turn transition of §18 applied to the first four, and
`omega` a random walk. The transition is now **nonlinear** — `omega` multiplies
the velocity terms — so the covariance must go through a Jacobian:

```
d(x')/d(omega) = vx d/dw[sin(wT)/w] - vy d/dw[(1-cos wT)/w]
d(y')/d(omega) = vx d/dw[(1-cos wT)/w] + vy d/dw[sin(wT)/w]
d(vx')/d(omega) = -T sin(wT) vx - T cos(wT) vy
d(vy')/d(omega) =  T cos(wT) vx - T sin(wT) vy

d/dw[ sin(wT)/w ]    = (T cos(wT) w - sin(wT)) / w^2      ->  -w T^3/3
d/dw[ (1-cos wT)/w ] = (T sin(wT) w - (1-cos wT)) / w^2   ->   T^2/2
```

with the limits as `w -> 0` on the right. Both the transition and its Jacobian
contain the same 0/0, and they must switch to the series at the **same**
threshold: a state and a covariance that switch at different points describe
different filters. The implementation computes all six quantities in one place
for that reason.

The measurement Jacobian gains a fifth column of **zeros**. Range, bearing and
range rate do not depend on `omega` at all — it is observed only through its
effect on where the target is predicted to be, which is why several scans of
turning are needed before the estimate is worth anything.

### 19.3 The knob that decides whether omega can be learned

The random walk on `omega` is the one parameter with no counterpart in a
constant-velocity filter, and **setting it to zero fails in the way hardest to
notice**. Target flies straight for 30 scans, then turns at 5.00 °/s:

| `q_w` | reported σ | estimate |
|---|---|---|
| 0 | **0.11 °/s** | **2.08 °/s** |
| 0.005 | 0.91 °/s | 4.97 °/s |
| 0.01 (default) | 1.54 °/s | 4.90 °/s |
| 0.02 | 2.67 °/s | 4.89 °/s |

With `q_w = 0` the covariance shrinks monotonically, the Kalman gain on the
fifth state goes to nothing, and the filter can no longer follow a change. It
then reports the *smallest* uncertainty of the four about the *most wrong*
answer. Confident and wrong is a worse failure than uncertain and wrong, because
nothing downstream can tell.

The reported σ settling to a steady state rather than collapsing is therefore
correct, not a defect: `omega` can change at any moment, and a filter that stops
allowing for that has stopped being a tracker.

### 19.4 What the fifth state costs

On a straight target, against a plain CV EKF over 12 runs: **36.28 m against
35.74 m**, a 1.02× penalty. Cheaper than expected — the extra state costs
almost nothing when there is nothing to estimate, because the measurement
simply never moves it.

### 19.5 Choosing between them

They are not ordered.

* **IMM** — reacts faster when a manoeuvre *starts*, because a model that
  already fits is waiting to take over. Cannot report a rate outside its
  bracket. Robust: no model can diverge, only lose probability.
* **CTRV** — measures the rate, with no ceiling. Slower to respond, and carries
  a nonlinearity that an IMM does not.

CTRV *measures* a manoeuvre; an IMM *reacts* to one.

### 19.6 The IMM in the tracker

`imm_tracker_step()` is the IMM counterpart of `tracker_step()`, with the same
global-cost-ordered association and the same M-of-N management. One decision
matters: **gating is on the combined estimate**, not per-model.

A per-model gate would let the worst-fitting model veto a measurement the
mixture accepts happily — and during a manoeuvre that is every model but one.
The mixture's covariance already carries the spread between models (§18.1), so
it widens exactly when the models disagree, which is when the gate should be
generous. A per-model gate throws that property away.

Measured: two targets turning in opposite directions from scan 15 hold two
established tracks through 40 scans, each detecting its own manoeuvre with the
correct sign. And 367 false alarms over 200 scans produce **zero** established
tracks, so M-of-N confirmation survives the swap to a mixture.
