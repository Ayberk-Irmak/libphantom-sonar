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
