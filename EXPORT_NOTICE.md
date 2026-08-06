# Export control and publication notice

## What this project is

`libphantom-sonar` implements ocean acoustic propagation physics that is
published in the open scientific literature — the sound speed equations of
Medwin, Mackenzie and Chen & Millero, Snell's law in a stratified medium, and
the standard constant-gradient ray solution found in Jensen et al.,
*Computational Ocean Acoustics*. Every formula in the code is cited in
`docs/math_spec.md` with its published source.

It is a **simulation and analysis library**. It contains no platform signature
data, no measured target strengths, no waveform library from any fielded system,
and no parameters derived from any non-public source.

## What this project is not

- Not deployable countermeasure firmware
- Not a model of any specific vehicle, sensor or national capability
- Not validated for, or intended for, operational use

## For contributors and users

Underwater acoustic detection and countermeasure technologies appear in national
export control schedules, including lists derived from the Wassenaar
Arrangement. Publicly available and fundamental research generally falls outside
those controls in most jurisdictions, and this project is intended to stay
firmly inside that space.

That intention is not a legal determination, and this file is not legal advice.
Specifically:

1. **If you work in the defence sector, get written publication approval before
   contributing.** Most defence employment contracts contain intellectual
   property and publication clauses that apply to personal projects in an
   overlapping technical field. Ask first; a retroactive problem is much worse
   than a delayed merge.
2. **If you intend to contribute anything derived from non-public sources —
   measured data, system parameters, internal reports — do not.** Contributions
   must be traceable to published literature or to your own bench measurements
   of hobbyist-grade hardware.
3. **Check with your organisation's export control function** if you are unsure
   whether a specific contribution changes the character of the project.

## Contribution policy

Pull requests that add capability derived from classified, controlled, or
proprietary sources will be rejected regardless of technical quality. If a
contribution cannot cite a public reference or an open measurement, it does not
belong here.
