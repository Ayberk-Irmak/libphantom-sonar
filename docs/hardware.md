# Hardware-in-the-loop bench

The plan for turning `libphantom-sonar` from a simulation into something with
measurements behind it. Written to be executable: specific parts, specific
frequencies, specific reasons.

Prices are order-of-magnitude as of 2026 and will drift. Nothing here needs a
research vessel.

---

## 0. The one design decision that matters

**Pick the operating frequency for the tank, not for realism.**

A tank is an acoustically terrible room. Standing waves and reverberation
dominate unless the tank is many wavelengths across. In water at `c ≈ 1500 m/s`:

| frequency | wavelength | 1 m tank is… |
|---|---|---|
| 5 kHz | 30 cm | ~3 wavelengths — hopeless |
| 20 kHz | 7.5 cm | ~13 wavelengths — marginal |
| **40 kHz** | **3.75 cm** | **~27 wavelengths — workable** |
| 200 kHz | 7.5 mm | ~130 wavelengths — comfortable |

So: **work at 40–200 kHz**, document the scaling to the 1–20 kHz band that real
systems use, and never quietly imply the scaled result is the real one. Scaling
frequency up by 10× and range down by 10× keeps the wavelength-to-geometry ratio
honest; absorption does not scale the same way, and that has to be stated.

---

## 1. Bill of materials

### Stage A — signal path only, no water (~€150)

Prove the DSP before adding acoustics. Every latency bug is easier to find here.

| Item | Purpose | ~Cost |
|---|---|---|
| Any Linux PC you already own | host | — |
| USB audio interface, 192 kHz (Focusrite Scarlett 2i2 or similar) | analogue loopback at up to ~90 kHz | €130 |
| Coax, BNC/TRS adapters, 50 Ω terminators | loopback rig | €20 |

Deliverable: ping generated → played out → captured → detected, with a measured
detection-latency histogram. No water involved.

### Stage B — water tank (~€700 on top of A)

| Item | Purpose | ~Cost |
|---|---|---|
| Aquarian Audio H1a or H2a hydrophone (with preamp) | calibrated receiver, flat to ~100 kHz | €250–400 |
| Piezo ceramic disc transducers, 40 kHz, 2–4 pcs | projector; cheap and replaceable | €20 |
| Small class-D or lab amplifier, 10–20 W | driving the projector | €60 |
| Glass or acrylic tank, ≥ 80 × 40 × 40 cm | the ocean | €120 |
| Acoustic absorber (closed-cell foam, anechoic wedges) | kills tank reverberation | €80 |
| Thermometer / conductivity probe | so `mackenzie(T, S, z)` has real inputs | €60 |
| Positioning rig (aluminium extrusion, clamps) | repeatable geometry — matters more than it sounds | €80 |

Deliverable: measured two-way travel time vs. hydrophone separation, compared
against `trace_ray()` with the *measured* T and S. This is the figure that turns
the repo from "well-tested code" into "validated against reality".

### Stage C — real-time embedded (~€500 on top of B)

| Item | Purpose | ~Cost |
|---|---|---|
| **Red Pitaya STEMlab 125-14** | Zynq-7010 FPGA + Linux + 2×ADC / 2×DAC at 125 MS/s | €400–500 |
| Alternative: STM32H7 Nucleo + external ADC | if you want bare-metal WCET numbers | €80 |
| Alternative: Raspberry Pi 5 | if you only need the Linux-side pipeline | €80 |

The Red Pitaya is the right single box: FPGA for the front end, Linux for the
library, both ADC and DAC on the same clock. That last part is what makes an
honest end-to-end latency measurement possible.

Deliverable: ping arrives → detected → synthetic echo transmitted, with the full
loop latency measured on a scope, not estimated.

---

## 2. Latency budget

This is the whole engineering problem. A synthetic echo that arrives late is a
second target, not a decoy.

| Stage | Budget | Notes |
|---|---|---|
| Acoustic propagation, 10 m round trip | 13.3 ms | this is the *window*, not a cost |
| ADC + buffering | < 200 µs | keep buffers small; this dominates on USB audio |
| Ping detection (matched filter) | **< 1 ms** | the binding constraint |
| Ray trace / delay computation | < 50 µs | measured at 40.8 µs for a 100 km ray; a tank case is far shorter |
| Waveform synthesis | < 200 µs | |
| DAC + amplifier | < 200 µs | |

The 13.3 ms window for a 10 m offset is generous. **Detection latency is where
the project succeeds or fails**, which is another way of saying `PingAnalyzer`
(v0.3) is the real work, not the echo synthesiser.

USB audio interfaces will not get you under ~5 ms round trip. That is fine for
Stage A/B (measuring propagation) and not fine for Stage C (closing the loop) —
hence the Red Pitaya.

---

## 3. Measurement protocol

1. **Calibrate the geometry first.** Measure hydrophone separation with a ruler,
   to the millimetre. At 1500 m/s, 1 mm is 0.67 µs; if your timing resolution is
   better than that, your ruler is now the error term.
2. **Measure T and S every run.** Tank water warms during a session. A 1 °C
   change moves sound speed by ~4.6 m/s, which is 0.3% — larger than most of the
   effects being measured.
3. **Characterise the tank before trusting it.** Record the impulse response;
   find the direct arrival and the first surface/wall reflections. Anything
   arriving after the first reflection is tank, not physics. Gate accordingly.
4. **Publish the residual, not the agreement.** "Simulated and measured travel
   times agree" means nothing. "Mean residual 4 µs, σ = 11 µs over 40
   separations from 0.1 to 0.7 m" means something.

---

## 4. Regulatory and ethical constraints

- **Tank and bench only.** Do not transmit into open water. Underwater acoustic
  transmission at power is regulated in most jurisdictions and is a genuine
  marine-mammal harassment issue — the same literature that makes bio-mimetic
  signalling interesting is the literature on why it disturbs animals.
- Bio-mimetic waveform work should stay at simulation and tank scale. Playing
  synthesised marine mammal vocalisations into the sea is not a demo, it is an
  incident.
- Check `EXPORT_NOTICE.md` before publishing measurement data from any hardware
  configuration that resembles a fielded system.

---

## 5. What to build first

If you want one deliverable that makes the repository credible, it is **not** the
full countermeasure loop. It is Stage B, step 4: a scatter plot of measured
versus simulated travel time across 40 hydrophone separations, with the residual
quoted, taken in a 60 cm tank with a €400 hydrophone.

That plot is small, cheap, honest, and almost nobody publishing acoustic
simulation code on GitHub has one.
