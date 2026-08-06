# Roadmap

The guiding rule: **a finished small thing beats an unfinished big thing.** Each
release below is independently useful and independently publishable. Nothing
ships without the validation section of `docs/validation.md` being updated
first.

---

## v0.1 — Ocean model and ray tracing ✅ *(current)*

- Three sound speed equations (Medwin / Mackenzie / Chen-Millero-UNESCO) plus
  Leroy-Parthiot depth→pressure and the Munk canonical profile
- Fixed-capacity piecewise-linear sound speed profile, zero allocation
- Analytic constant-gradient ray tracer: exact circular arcs, exact turning
  points, exact range **and** travel-time budgets
- SOFAR / duct analysis: axis, trapping cone, conjugate depths, surface duct
- Ensonification grid and shadow zone extraction
- 791 closed-form checks, ASan/UBSan clean, float and double builds

## v0.2 — External validation and real data

The credibility release. Nothing new gets built until the existing code is
checked against something that was not written here.

- [ ] **Bellhop cross-validation.** Run the Munk `MunkB_ray` case through the
      Acoustics Toolbox, overlay ray paths, publish RMS path deviation and the
      comparison figure. This is the highest-value single task in the project.
- [ ] **Chen-Millero against the published UNESCO check table**, as a unit test
      with the official values rather than mutual agreement between equations.
- [ ] **Real T/S profiles.** World Ocean Atlas (WOA23) and Argo float ingest;
      ship two or three real regional profiles in `data/` with provenance.
- [ ] Eigenray search — launch angles that connect a source to a given receiver.
      Needed by everything in v0.3.
- [ ] Range-dependent bathymetry (piecewise-linear bottom).

## v0.3 — Ping analysis and echo synthesis

Note the ordering. The original specification jumped straight to echo synthesis,
but **an echo synthesiser with no input is not a module.** Detecting and
characterising the incoming ping is the hard part and it comes first.

- [ ] `PingAnalyzer` — replica-correlator / matched-filter bank producing a
      Pulse Descriptor Word: type (CW / LFM / HFM), centre frequency, bandwidth,
      pulse width, chirp rate, time of arrival, estimated bearing. Detection
      latency is the binding constraint, target < 1 ms.
- [ ] Transmission loss along traced paths: spherical/cylindrical spreading plus
      Thorp absorption, so target strength numbers mean something.
- [ ] `EchoSynthesizer` — delay/Doppler/target-strength shaped echo generation,
      `Δt = 2Δd/c`; multi-target pulse trains.
- [ ] Ray-tube area and caustic handling, which is what a physically defensible
      shadow zone actually requires.

**Honest scope note.** Anti-phase cancellation of a vehicle's acoustic cross
section is *not* achievable in practice with a point transducer against a
broadband ping across all bearings — hull scattering is distributed and
broadband. What can be simulated is narrowband, single-bearing nulling, and that
is what will be implemented, documented with its limits. Overclaiming here is
how a project loses the audience it is trying to reach.

## v0.4 — Covert acoustic communication

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

## v0.5 — Bindings and portability

- [ ] Stable C ABI (`phantom.h`, hand-written, no generator)
- [ ] Rust `phantom-sonar-sys` + safe wrapper crate
- [ ] Cortex-M7 / RISC-V cross-compilation in CI with size and WCET reporting
- [ ] `-fno-exceptions -fno-rtti` build mode

## v0.6 — Hardware in the loop

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
