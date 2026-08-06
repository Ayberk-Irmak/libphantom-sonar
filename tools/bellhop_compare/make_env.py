#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Generate Bellhop .env files for the cross-validation cases.

The sound speed table is the canonical Munk profile shipped with the Acoustics
Toolbox as `at/tests/Munk/MunkB_ray.env`, reproduced verbatim so the comparison
runs against the reference community benchmark rather than a profile invented
here.

Every generated case uses **C-linear** interpolation ('C'). libphantom-sonar
stores the profile as piecewise linear in c(z), which is exactly what Bellhop's
'C' option does. Running against PCHIP or spline would compare interpolation
schemes, not ray tracers, and would produce a disagreement that says nothing
about either code's correctness.
"""

import argparse
from pathlib import Path

# at/tests/Munk/MunkB_ray.env, verbatim.
MUNK_SSP = [
    (0.0, 1548.52), (200.0, 1530.29), (250.0, 1526.69), (400.0, 1517.78),
    (600.0, 1509.49), (800.0, 1504.30), (1000.0, 1501.38), (1200.0, 1500.14),
    (1400.0, 1500.12), (1600.0, 1501.02), (1800.0, 1502.57), (2000.0, 1504.62),
    (2200.0, 1507.02), (2400.0, 1509.69), (2600.0, 1512.55), (2800.0, 1515.56),
    (3000.0, 1518.67), (3200.0, 1521.85), (3400.0, 1525.10), (3600.0, 1528.38),
    (3800.0, 1531.70), (4000.0, 1535.04), (4200.0, 1538.39), (4400.0, 1541.76),
    (4600.0, 1545.14), (4800.0, 1548.52), (5000.0, 1551.91),
]

BOTTOM_DEPTH = 5000.0


def write_env(path, title, source_depth_m, angle_deg, n_beams,
              range_km=100.0, step_m=0.0, ssp=MUNK_SSP):
    """Write one ray-trace case.

    step_m = 0 lets Bellhop choose its own step size, which is what a user gets
    by default and therefore what the headline comparison should use.
    """
    lines = [
        f"'{title}'\t\t! TITLE",
        "50.0\t\t\t! FREQ (Hz)",
        "1\t\t\t! NMEDIA",
        # C-linear interpolation, Vacuum above the surface, dB/m/kHz attenuation.
        "'CVF'\t\t\t! SSPOPT - C-linear to match libphantom-sonar",
        f"{len(ssp)}  0.0  {BOTTOM_DEPTH:.1f}\t\t! NPTS, SIGMA, bottom depth (m)",
    ]
    for z, c in ssp:
        lines.append(f"  {z:8.1f}  {c:8.2f}  /")
    lines += [
        "'A' 0.0\t\t\t! bottom BC: acousto-elastic halfspace",
        f" {BOTTOM_DEPTH:.1f}  1600.00 0.0 1.0 /",
        "1\t\t\t! NSD",
        f"{source_depth_m:.1f} /\t\t! SD (m)",
        "51\t\t\t! NRD",
        f"0.0 {BOTTOM_DEPTH:.1f} /\t\t! RD (m)",
        "1001\t\t\t! NR",
        f"0.0  {range_km:.1f} /\t\t! R (km)",
        "'R'\t\t\t! run type: ray paths",
        # An explicit beam count matters: with 0 Bellhop picks its own, and the
        # two runs would no longer be launching the same rays.
        f"{n_beams}\t\t\t! NBEAMS (explicit, so both codes launch the same rays)",
        f"{-angle_deg:.4f} {angle_deg:.4f} /\t\t! ALPHA1,2 (degrees)",
        f"{step_m:.1f}  {BOTTOM_DEPTH + 500.0:.1f}  {range_km + 1.0:.1f}"
        "\t\t! STEP (m), ZBOX (m), RBOX (km)",
    ]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir", nargs="?", default="bellhop_run",
                    help="directory to write .env files into")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    written = []

    # Case A -- pure refraction. Source on the 1000 m isobath, launch angles
    # inside the trapping cone (+/- 14.19 deg from this depth), so no ray ever
    # touches a boundary. This isolates the refraction kernel: any disagreement
    # is the arc solution versus Bellhop's stepping, with nothing else mixed in.
    written.append(write_env(
        out / "cmp_trapped.env", "Munk C-linear, trapped rays",
        source_depth_m=1000.0, angle_deg=13.0, n_beams=41))

    # Case B -- with boundary interactions. Steeper launch angles escape the
    # channel and bounce off the surface and the sea floor, which exercises the
    # reflection bookkeeping as well.
    written.append(write_env(
        out / "cmp_bounce.env", "Munk C-linear, surface and bottom bounces",
        source_depth_m=1000.0, angle_deg=20.0, n_beams=41))

    # Case C -- a single ray at several Bellhop step sizes. The point is not
    # just that the two agree, but that Bellhop's answer *converges toward* the
    # analytic arc solution as its step shrinks. That is the direction of the
    # error, and it is the part worth publishing.
    for step in (500.0, 200.0, 100.0, 50.0, 20.0, 10.0):
        written.append(write_env(
            out / f"cmp_step_{int(step):04d}.env",
            f"Munk C-linear, step {step:g} m",
            source_depth_m=1000.0, angle_deg=10.0, n_beams=11, step_m=step))

    for p in written:
        print(p)


if __name__ == "__main__":
    main()
