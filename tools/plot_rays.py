#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Render the CSV written by the munk_simulation example.

This is tooling, not library code -- the zero-dependency promise applies to
libphantom-sonar itself, not to the plotting helpers around it.

    ./build/munk_simulation data
    python3 tools/plot_rays.py data
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_profile(path):
    depth, speed = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            depth.append(float(row["depth_m"]))
            speed.append(float(row["munk_mps"]))
    return np.array(depth), np.array(speed)


def read_rays(path):
    rays = defaultdict(lambda: ([], []))
    launch = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rid = int(row["ray_id"])
            launch[rid] = float(row["launch_deg"])
            r, z = rays[rid]
            r.append(float(row["range_m"]) / 1000.0)
            z.append(float(row["depth_m"]))
    return rays, launch


def read_coverage(path):
    r, z, hits = [], [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            r.append(float(row["range_m"]) / 1000.0)
            z.append(float(row["depth_m"]))
            hits.append(int(row["hits"]))
    ru = np.unique(r)
    zu = np.unique(z)
    grid = np.array(hits, dtype=float).reshape(len(zu), len(ru))
    return ru, zu, grid


def main(argv):
    data = Path(argv[1] if len(argv) > 1 else "data")
    out = data / "munk_rays.png"

    depth, speed = read_profile(data / "profile.csv")
    rays, launch = read_rays(data / "rays.csv")
    ru, zu, cov = read_coverage(data / "coverage.csv")

    axis_depth = float(depth[int(np.argmin(speed))])

    fig = plt.figure(figsize=(16, 7.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[1, 1])

    # --- sound speed profile ------------------------------------------------
    ax0 = fig.add_subplot(gs[:, 0])
    ax0.plot(speed, depth, lw=1.8, color="#1f4e79")
    ax0.axhline(axis_depth, color="#c0392b", ls="--", lw=1.0)
    ax0.annotate(
        f"SOFAR axis\n{axis_depth:.0f} m",
        xy=(speed.min(), axis_depth),
        xytext=(6, -34),
        textcoords="offset points",
        fontsize=9,
        color="#c0392b",
    )
    ax0.invert_yaxis()
    ax0.set_xlabel("sound speed  [m/s]")
    ax0.set_ylabel("depth  [m]")
    ax0.set_title("Munk profile", fontsize=11)
    ax0.grid(alpha=0.25)

    # --- ray fan ------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 1])
    angles = np.array([launch[k] for k in sorted(rays)])
    norm = plt.Normalize(angles.min(), angles.max())
    cmap = plt.get_cmap("coolwarm")
    for rid in sorted(rays):
        r, z = rays[rid]
        ax1.plot(r, z, lw=0.6, color=cmap(norm(launch[rid])), alpha=0.9)
    ax1.axhline(axis_depth, color="#c0392b", ls="--", lw=1.0, zorder=5)
    ax1.invert_yaxis()
    ax1.set_ylabel("depth  [m]")
    ax1.set_title(
        f"{len(rays)} rays launched from the channel axis, "
        "analytic constant-gradient arcs",
        fontsize=11,
    )
    ax1.grid(alpha=0.2)
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax1,
        label="launch angle [deg]",
        pad=0.01,
    )

    # --- ensonification / shadow zones --------------------------------------
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)
    shaded = np.ma.masked_where(cov > 0, np.ones_like(cov))
    ax2.pcolormesh(
        ru, zu, np.log10(cov + 1.0), cmap="viridis", shading="auto", rasterized=True
    )
    ax2.pcolormesh(ru, zu, shaded, cmap="Greys_r", vmin=0, vmax=1, shading="auto")
    ax2.invert_yaxis()
    ax2.set_xlabel("range  [km]")
    ax2.set_ylabel("depth  [m]")
    shadow_pct = 100.0 * float((cov == 0).sum()) / cov.size
    ax2.set_title(
        f"ensonification (log ray-hit count); white = shadow zone, "
        f"{shadow_pct:.1f}% of cells",
        fontsize=11,
    )

    fig.suptitle(
        "libphantom-sonar - deep sound channel propagation, 100 km",
        fontsize=13,
        fontweight="bold",
    )
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
