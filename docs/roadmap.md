# Roadmap

The guiding rule: **a finished small thing beats an unfinished big thing.** Each
release below is independently useful and independently publishable. Nothing
ships without the validation section of `docs/validation.md` being updated
first.

---

## v0.1 — Ocean model and ray tracing ✅

- Three sound speed equations (Medwin / Mackenzie / Chen-Millero-UNESCO) plus
  Leroy-Parthiot depth→pressure and the Munk canonical profile
- Fixed-capacity piecewise-linear sound speed profile, zero allocation
- Analytic constant-gradient ray tracer: exact circular arcs, exact turning
  points, exact range **and** travel-time budgets
- SOFAR / duct analysis: axis, trapping cone, conjugate depths, surface duct
- Ensonification grid and shadow zone extraction
- 791 closed-form checks, ASan/UBSan clean, float and double builds
- **Bellhop cross-validation** (`tools/bellhop_compare/`): turning depths agree
  to 0.14 mm over 101 km at Bellhop's finest step, with the residual shown to be
  the comparison's own sampling error rather than a difference between the codes

## v0.2 — Ping analysis ✅

The original specification jumped straight to echo synthesis, but **an echo
synthesiser with no input is not a module.** Detecting and characterising the
incoming ping is the harder half, so it shipped first.

- Radix-2 FFT, zero allocation, verified against a direct DFT to 2.3e-16
- CW / LFM / HFM synthesis with closed-form phase, and true time-scaling Doppler
- FFT overlap-save matched filtering, verified against direct correlation
- `PingAnalyzer`: matched filter bank + CA-CFAR + Pulse Descriptor Words
- ToA estimation at **0.93-1.02x the envelope Cramer-Rao bound** across 20 dB
- 100/100 waveform classification at -4.4 dB input SNR
- 29x real time on one core; 65 ms detection latency, set by the block length

## v0.3 — Doppler bank and echo synthesis ✅

- Doppler bin spacing derived per waveform family and checked against the
  measured loss; over +/- 20 m/s at 1 dB it takes 14 bins for a CW, 5 for an
  LFM and 2 for an HFM
- `PulseBank::add_doppler_bank()`, radial velocity in the PDW
- `EchoSynthesizer`: two-way delay and Doppler, target strength, extended
  targets, ghost swarms -- all verified by closed loop through the analyser
- Anti-phase cancellation quantified rather than implemented: 20 dB needs the
  timing held to 1.3 us, under two millimetres of path
- Cross-check between halves: a Doppler bin error of 4 m/s produces exactly the
  arrival-time bias the wideband coupling formula predicts

## v0.4 — Ray-acoustics coupling ✅

The two halves now meet through the physics, not just the sound speed.

- Eigenray search, exact rather than interpolated: the tracer's range budget is
  set to the receiver range, so the final point *is* the arrival
- Exact arc-length accumulation in the tracer (`R|dtheta|` per segment)
- Ray-tube geometric spreading, reducing to `20 log10(R)` in isovelocity water
  to **0.0000 dB** across nine geometries from 1 to 20 km
- Thorp absorption; at 100 kHz over 3 km it exceeds spreading by 30 dB
- Caustics flagged, not reported as levels
- Multipath echo synthesis: arrival structure from the traced paths
- A measured precision limit: `depth_tolerance_m` discards roots, so a float
  build silently loses 11 of 14 paths at the double-precision default. The
  default now follows the build.

## v0.5 — Boundary losses ✅

The largest remaining overstatement in the transmission loss: every boundary
used to reflect perfectly, and in shallow water most paths bounce.

- Rayleigh bottom reflection, pinned by three closed-form limits: the critical
  angle from the speed ratio, the impedance ratio at normal incidence, and
  unity below critical for a lossless bottom
- Sediment attenuation as a complex speed, so sub-critical rays leak; without
  it shallow-water range is unbounded
- Energy conservation over 9100 combinations settles the complex-sqrt branch
- Rough-surface coherent loss, with the cap the diffuse field demands
- Grazing angles at each bounce from the Snell invariant, exact
- A 200 m duct at 3 km: paths that sat within 1 dB now span 128 dB

## v0.6 — Reverberation ✅

Where the energy v0.5 capped away actually goes, and what an active sonar is
really fighting at short range.

- Lambert bottom backscatter and Chapman-Harris surface backscatter
- Ensonified area and volume, and the two decay laws that identify the
  mechanism: 30 log10(r) for boundary, 20 log10(r) for volume
- The echo-to-reverberation ratio, in which source level and transmission loss
  cancel exactly -- verified over nine SL/TL combinations
- Reverberation-limited range: where a bigger transmitter starts helping again
- A reverberation envelope generator, and a demonstration that CA-CFAR tracks a
  21.6 dB decay across a block where no fixed threshold can

## v0.7 — Arrays and bearing ✅

The last purely structural gap: everything before this was single-channel.

- Line array geometry, the array factor and its closed forms -- peak, nulls,
  the -13.26 dB first sidelobe, grating-lobe spacing
- Array gain, measured through the beamformer rather than asserted
- The bearing Cramer-Rao bound, the spatial twin of v0.2's arrival-time bound:
  N^(-3/2) scaling and 1/cos(theta) degradation both verified, and the
  estimator lands at 0.99-1.13x of it
- Split-beam bearing, exact across the mainlobe and wrapping beyond it
- The bandwidth limit phase steering carries -- 91 Hz at 45 degrees for an
  array whose own waveforms are 12 kHz wide

## v0.8 — Beams meet the detector ✅

v0.7 built an array and measured why it could not be used with this library's
own waveforms. This closes that, and joins the array to the analyser.

- Wideband per-bin beamforming: exact fractional delay, verified against the
  closed-form array factor to 5e-3
- A 12 kHz chirp steered 37x beyond the phase-steering limit, with the matched
  filter peak landing exactly where the pulse was placed
- The analyser run per beam, so a Pulse Descriptor Word carries a bearing
- Shading, with the sidelobe/beamwidth/gain trade measured against the
  published window values
- MVDR: two sources 4.30 degrees apart resolved where the conventional
  aperture limit is 7.16, by Cholesky rather than an explicit inverse
- Diagonal loading verified as necessary rather than decorative

## v0.9 — Tracking ✅

- An EKF over a constant-velocity target: Cartesian dynamics, polar
  measurements, with the nonlinearity confined to the measurement Jacobian
- Filter consistency verified where it counts -- NIS mean 1.951 against a
  2-dof expectation of 2.000, with 95.0% and 99.0% under their gates. A filter
  that tracks well and reports the wrong covariance passes every position test
  and fails this one.
- Chi-square gating, greedy association, M-of-N confirmation and coasting
- 494 false alarms over 300 scans produce zero confirmed tracks
- **A roadmap claim corrected**: v0.8 said tracking would suppress the
  cross-template ghosts of v0.2. It does not, and the test suite now asserts
  the honest outcome. Ghosts repeat; false alarms do not, and only the latter
  are what time consistency removes.

## v0.10 — Fusing what is already measured ✅

Several quantities were measured by one subsystem and discarded at the boundary
of the next. Joining them is cheaper than any new physics and worth more.

- [x] **The Doppler bank's closing rate feeds the tracker** as a third
      measurement, with the chi-square gate following the dimension (3-dof
      quantiles bisected on the exact CDF: 7.815 / 11.345). Consistency holds,
      mean NIS 2.998 against a theoretical 3. Radial-velocity error improves
      **3.45×** at 3 scans, 3.19× at 6 — and **0.99× at 20**, which is reported
      rather than buried: by then the filter's own estimate already beats the
      measurement. An unresolved bin is explicitly *not* a measurement of zero.
- [x] **Ghost suppression by template-artefact recognition**, delivering what
      v0.8 wrongly expected of tracking. The distinguishing feature is origin,
      not kinematics: shared bearing, fixed offset, weaker, and a **different
      waveform label**. The label check is the safety argument — a line-astern
      formation with the same waveform is kept, the same geometry with different
      waveforms is suppressed, and the two runs differ only in that field.
- [x] **Global-cost-ordered association.** Every gating pair is built, sorted by
      NIS and assigned best-first. Crossing targets keep their identities where
      greedy-by-arrival-order swapped them. Not Hungarian-optimal, and said so.
- [x] **Forward-backward spatial smoothing**, so MVDR survives the coherent
      multipath v0.4 produces. Two coherent arrivals at ±6° go from 4 spurious
      peaks to 2 correct ones. The cost is aperture — resolution falls from
      7.16° to 11.46° — and a test asserts that it does.
- [ ] *Deferred:* an interacting-multiple-model filter, so a manoeuvre is
      modelled rather than absorbed as process noise. Moved to v0.11; it is new
      estimation machinery rather than a fusion of what already exists, which is
      what this release was about.

## v0.11 — Real data and manoeuvre ✅

- [x] **An interacting-multiple-model filter.** Three models sharing one
      four-state vector: constant velocity, and coordinated turns at +/-omega.
      Against a single-model EKF on a 3 deg/s turn the gain is **1.72x AFTER
      the manoeuvre and 0.98x during it** -- the win is in recovery, not in the
      turn, which is not the usual claim and is what the measurement says.
      Costs nothing on a straight target (27.40 m vs 27.88 m). Ships a free
      manoeuvre detector: peak P(manoeuvre) 0.991, direction recovered.
      NOT yet wired into tracker_step(); that is v0.12.
- [x] **Chen-Millero against the published UNESCO check table.** This found a
      real bug: the equation was Wong & Zhu's ITS-90 coefficients with two
      terms (A02, A03) left at their 1977 values, so it was neither equation and
      missed the official check value by 0.016 m/s. Eleven releases of
      mutual-agreement testing could not see it, because the equations only
      agree to 0.1 m/s. Now both versions are implemented properly and verified
      against the primary source: check value **1731.9954** vs published
      1731.995, and all **220 table values** within 0.0499 m/s -- the table's
      own rounding half-width.
- [x] **Del Grosso (1974)** added as a genuinely independent fourth equation.
      Its disagreement with UNESCO is 0.41 m/s in the top kilometre and 3.93 m/s
      over the full nominal validity box -- the difference being that the box
      contains 26 C water under 1000 bar, which no sea does.
- [x] **Real T/S profiles.** Six WOA23 sites (Levantine, Aegean, Black Sea,
      North Atlantic, Norwegian, equatorial Pacific), fetched over OPeNDAP by
      `tools/fetch_woa_profiles.py` so every number's provenance is a URL. The
      Black Sea earns its place: at 18.2 PSU surface salinity, assuming the
      usual 35 PSU is a 20 m/s error -- the size of the whole channel feature.
- [x] **A turning-depth check** that is independent of how the tracer stores
      its state, replacing a Snell-invariant check that was largely measuring
      whether `acos` and `cos` round-trip. 186 turning points, all exact.

## v0.12 — IMM in the tracker, and the turn rate measured ✅

- [x] **The IMM wired into a full tracker.** `imm_tracker_step()` carries the
      same global-cost-ordered association and M-of-N management as
      `tracker_step()`, gating on the COMBINED estimate rather than per model --
      a per-model gate would let the worst-fitting model veto a measurement the
      mixture accepts, which during a manoeuvre is every model but one. Two
      targets turning in opposite directions hold two tracks through 40 scans
      with the correct manoeuvre signs; 367 false alarms produce zero tracks.
- [x] **The turn rate estimated rather than bracketed**, as a fifth state with a
      nonlinear transition. Against models bracketed at 3 deg/s and a truth of
      5 deg/s: **4.99 vs 2.62**. At 8 deg/s the IMM reports *less* turn than at
      5, because once the truth leaves the bracket the models fit so badly that
      probability drifts back to constant velocity -- a bracket does not degrade
      gracefully, which is the finding that justifies the fifth state.
- [x] **The trap found and documented**: with zero process noise on omega the
      covariance shrinks, the gain dies, and the filter reports the SMALLEST
      uncertainty of any setting about the MOST wrong answer (0.11 deg/s sigma,
      2.08 deg/s estimate against a truth of 5.00). The default was lowered from
      0.02 to 0.01 on the strength of the sweep.
- [x] **Cost of the fifth state measured**: 1.02x on a straight target, which is
      cheaper than expected -- an extra state costs almost nothing when the
      measurement never moves it.

## v0.13 — Spread-spectrum acoustic communication ✅

The third engine of the original specification, and the last real gap.

- [x] **DSSS, with processing gain stated correctly rather than flatteringly.**
      The first measurement here was WRONG in the project's recurring way: it
      held chip amplitude fixed while sweeping N, which multiplies energy per
      bit by N, and reported the resulting 12.17 dB against a theoretical 12.17
      as a triumph. Both sides were the same tautology. Corrected: against white
      noise at fixed energy per bit spreading buys NOTHING (flat to 0.87 dB over
      16x), and the real gain is bandwidth expansion against narrowband
      interference (+7.1 dB from 31 to 511 chips at a fixed data rate).
      The link budget is published as a TRADE -- 30 dB at 4 kchip/s is under
      4 bits per second, which is why a covert link carries status and not data.
- [x] **PN sequences whose maximality is verified**, degree 5 to 15: balance
      exactly +1 and autocorrelation exactly -1 at every shift. The first tap
      table was in the wrong convention and produced a constant with the right
      LENGTH; only these two measurements caught it.
- [x] **Doppler as a time scaling**, with the ceiling it puts on processing
      gain: 10 m/s permits 37 chips and 15.7 dB, and no transmit power lifts
      that. HFM vs LFM preambles under Doppler measured -- at 20 m/s the HFM
      holds 99.1% of its correlation and the LFM 54.8%.
- [x] **CRC-32 verified against its published check value**, and RS(15,11) over
      GF(16) correcting 2 symbol errors with 0 failures in 4000 trials. Beyond
      t=2 the two layers were measured TOGETHER: RS miscorrects 992 of 3000
      four-error words to a valid-but-wrong codeword -- the minimum distance
      showing, not a defect -- and the CRC catches all 992.
- [ ] *Deferred:* carrier and chip-timing recovery, and an equaliser. The
      demodulator is given the carrier phase. A few lines of arctangent would
      not be a recovery loop and the header says so rather than implying one.
- [ ] *Deferred:* bio-mimetic waveform shaping. It needs an ambient-noise and
      biological-transient model the library does not have, and shaping against
      a spectrum you have not measured is decoration.

## v0.14 — Bindings and portability ✅ *(current)*

- [x] **Stable C ABI** (`include/phantom/phantom.h`), hand-written, no
      generator. Caller-owned storage throughout -- the library allocates
      nothing and the ABI does not quietly change that. Tested by a C compiler
      in C11 mode, not by C++ in C mode: 269 checks.
- [x] **Rust `phantom-sonar-sys` + safe `phantom-sonar` crate.** 8 tests. The
      safe layer earns its place: a `Profile` owns its storage so the C
      contract becomes the borrow checker's problem, and `speed_at` returns an
      error where C answers 0 m/s for a one-sample profile -- a silent failure a
      binding test found.
- [x] **Cross-compilation verified** on ARM 32-bit/float and RISC-V
      64-bit/double, with size reporting. This found a real defect: the C ABI
      was carrying **160 kB of .bss** in a ray-trace scratch buffer -- half the
      RAM of the target part -- for a copy that did nothing. Replaced with
      static_asserts on the layout; the library now has data == 0 and bss == 0,
      asserted in CI.
- [x] **`-fno-exceptions -fno-rtti` build mode**, verified to reference no
      unwinder symbols.
- [ ] *Not done:* the bare-metal Cortex-M7 build. The toolchain file is
      provided but this environment's arm-none-eabi has no C++ standard
      library, so it is UNTESTED and labelled so rather than claimed.

## v0.15 — Hardware in the loop

See `docs/hardware.md` for the bench design and the bill of materials. Bench and
tank only — a real acoustic transmitter in open water is a regulatory and
marine-mammal issue, not a weekend experiment.

- [ ] Red Pitaya STEMlab acquisition path, measured end-to-end latency
- [ ] Tank measurements at 40–200 kHz, ping detection latency histogram
- [ ] Measured vs simulated arrival times, published as a figure

---

## Explicitly out of scope

- Anything requiring classified or export-controlled parameters. The library
  models physics that is in the open literature; it does not encode any
  platform's signature, any real system's waveform library, or any national
  capability.
- Deployable countermeasure firmware. This is a simulation and analysis library.
