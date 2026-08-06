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

## v0.4 — Ray-acoustics coupling ✅ *(current)*

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

## v0.5 — Wave-acoustics and boundaries

The gaps v0.4 documented rather than closed.

- [ ] Bottom loss: grazing-angle and frequency-dependent reflection, so a
      bottom-bounced path is not reported louder than it is.
- [ ] Surface scattering and bubble loss, both of which matter in sea state.
- [ ] Beam displacement and a caustic treatment that returns a level instead of
      a flag.
- [ ] Bistatic geometries: `echoes_from_eigenrays` assumes the echo returns the
      way the ping came.
- [ ] Range-dependent bathymetry (piecewise-linear bottom).
- [ ] Cross-template ghost suppression: association logic across detections.
- [ ] Bearing estimation, which needs an array rather than the single channel
      assumed so far.

## v0.6 — Real data

- [x] **Bellhop cross-validation.** Done in v0.1.1. Both codes read the same
      `.env`; Bellhop converges toward the analytic arc solution as its step
      shrinks. See `docs/validation.md §2`.
- [ ] **Chen-Millero against the published UNESCO check table**, as a unit test
      with the official values rather than mutual agreement between equations.
- [ ] **Real T/S profiles.** World Ocean Atlas (WOA23) and Argo float ingest;
      ship two or three real regional profiles in `data/` with provenance.

## v0.7 — Covert acoustic communication

- [ ] DSSS modulator/demodulator with configurable chip rate and explicit
      processing gain `10·log10(N)`
- [ ] **Doppler-tolerant synchronisation.** In water this dominates: `v/c` with
      `c ≈ 1500 m/s` means 1 m/s of closing speed is 6.7e-4 — and it is a
      *time-scaling* of the waveform, not merely a frequency shift. HFM sweeps
      for the preamble, not LFM.
- [ ] Multipath-tolerant framing; CRC-32 plus RS(15,11) over GF(16) with static
      tables
- [ ] Bio-mimetic waveform shaping (spectral masking against ambient noise and
      biological transients)

## v0.8 — Bindings and portability

- [ ] Stable C ABI (`phantom.h`, hand-written, no generator)
- [ ] Rust `phantom-sonar-sys` + safe wrapper crate
- [ ] Cortex-M7 / RISC-V cross-compilation in CI with size and WCET reporting
- [ ] `-fno-exceptions -fno-rtti` build mode

## v0.9 — Hardware in the loop

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
