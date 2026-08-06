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
