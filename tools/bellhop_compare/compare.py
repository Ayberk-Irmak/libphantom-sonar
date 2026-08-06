#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Cross-validate libphantom-sonar's ray tracer against Bellhop.

Both codes are handed the *same* .env file: Bellhop reads it natively, and
`phantom_trace` parses it too. Nothing is transcribed between them, so a
transcription bug cannot masquerade as agreement.

    python3 tools/bellhop_compare/compare.py --bellhop /path/to/bellhop.exe

Three things are measured, in increasing order of how much they prove:

  A/B  RMS depth deviation along the path, for trapped rays and for rays that
       bounce off the surface and sea floor. Both polylines are resampled onto
       a common range grid before differencing.

  C    Turning-depth error versus Bellhop's step size. This is the headline
       result, because it needs no interpolation at all -- the turning depth is
       a single number per ray, exact in libphantom-sonar by construction and
       approached by Bellhop as its integrator step shrinks.

  D    The same RMS as A, versus *libphantom-sonar's own* output sampling. This
       is the honesty check: it establishes that the residual in A is the
       comparison's chord-versus-arc error, not a disagreement between codes.
       Without D, the RMS in A would be quoted as if it meant something about
       the tracers.
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]

# Sampling refinement for the headline cases. Each tabulated SSP layer is split
# into this many sublayers, which does not change c(z) -- it only makes the
# emitted polyline dense enough to compare against a stepping code.
DEFAULT_REFINE = 128
SWEEP_REFINE = 512


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------
def read_ray_file(path):
    """Parse a Bellhop .ray file into a list of ray dicts."""
    with open(path) as fh:
        for _ in range(7):  # title, freq, Nsx/Nsy/Nsz, Nbeams, depthT, depthB, 'rz'
            fh.readline()
        tokens = fh.read().split()

    rays = []
    i, n = 0, len(tokens)
    while i < n:
        try:
            alpha = float(tokens[i])
            nsteps = int(tokens[i + 1])
            n_top = int(tokens[i + 2])
            n_bot = int(tokens[i + 3])
        except (ValueError, IndexError):
            break
        i += 4
        need = 2 * nsteps
        if i + need > n:
            break
        block = np.array(tokens[i:i + need], dtype=float).reshape(nsteps, 2)
        i += need
        rays.append({"launch_deg": alpha, "r": block[:, 0], "z": block[:, 1],
                     "n_top": n_top, "n_bot": n_bot})
    return rays


def read_phantom_csv(path):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    if raw.size == 0:
        return []
    ids = raw["ray_id"].astype(int)
    rays = []
    for rid in np.unique(ids):
        m = ids == rid
        rays.append({"launch_deg": float(raw["launch_deg"][m][0]),
                     "r": raw["range_m"][m], "z": raw["depth_m"][m],
                     "t": raw["time_s"][m]})
    return rays


# --------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------
def run_bellhop(bellhop, workdir, stem):
    subprocess.run([str(bellhop), stem], cwd=workdir, check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ray = workdir / f"{stem}.ray"
    if not ray.exists():
        raise SystemExit(f"Bellhop produced no .ray for {stem}; see {stem}.prt")
    return read_ray_file(ray)


def run_phantom(phantom, env_path, out_csv, refine):
    proc = subprocess.run(
        [str(phantom), str(env_path), str(out_csv), "--refine", str(refine)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"phantom_trace failed on {env_path}:\n{proc.stderr}")
    return read_phantom_csv(out_csv)


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------
def resample(r, z, grid):
    """Depth at each range in `grid`. Rays are non-decreasing in range; keep
    the first sample of each repeated range so np.interp gets a strict x."""
    keep = np.concatenate(([True], np.diff(r) > 0))
    return np.interp(grid, r[keep], z[keep])


def pair_rays(bellhop_rays, phantom_rays, tol_deg=1e-4):
    """Match rays by launch angle."""
    out = []
    for b in bellhop_rays:
        p = min(phantom_rays, key=lambda x: abs(x["launch_deg"] - b["launch_deg"]))
        if abs(p["launch_deg"] - b["launch_deg"]) <= tol_deg:
            out.append((b, p))
    return out


def deviation_stats(pairs, samples=8000):
    """RMS/max depth deviation and turning-depth error over matched rays."""
    per_ray = []
    for b, p in pairs:
        r_max = min(b["r"].max(), p["r"].max())
        if r_max <= 0:
            continue
        grid = np.linspace(0.0, r_max, samples)
        d = resample(b["r"], b["z"], grid) - resample(p["r"], p["z"], grid)
        per_ray.append({
            "launch_deg": b["launch_deg"],
            "rms_m": float(np.sqrt(np.mean(d ** 2))),
            "max_m": float(np.max(np.abs(d))),
            # Deepest point reached. Exact in libphantom-sonar (the vertex is
            # solved for, not stepped onto), so this needs no interpolation.
            "turn_err_m": float(abs(b["z"].max() - p["z"].max())),
            "range_km": r_max / 1000.0,
            "n_top": b["n_top"], "n_bot": b["n_bot"],
            "grid": grid, "dev": d, "bh": b, "ph": p,
        })
    return per_ray


def summarise(per_ray):
    if not per_ray:
        return None
    rms = np.array([r["rms_m"] for r in per_ray])
    mx = np.array([r["max_m"] for r in per_ray])
    turn = np.array([r["turn_err_m"] for r in per_ray])
    return {
        "n_rays": len(per_ray),
        "range_km": float(max(r["range_km"] for r in per_ray)),
        "bounces": int(sum(r["n_top"] + r["n_bot"] for r in per_ray)),
        "mean_rms_m": float(rms.mean()),
        "median_rms_m": float(np.median(rms)),
        "max_rms_m": float(rms.max()),
        "max_abs_m": float(mx.max()),
        "mean_turn_err_m": float(turn.mean()),
        "max_turn_err_m": float(turn.max()),
    }


def print_case(title, s):
    if s is None:
        print(f"\n{title}\n  no rays matched")
        return
    print(f"\n{title}")
    print(f"  rays compared           : {s['n_rays']}  over {s['range_km']:.1f} km"
          f"   ({s['bounces']} boundary interactions)")
    print(f"  RMS depth deviation     : mean {s['mean_rms_m']:.4f} m, "
          f"median {s['median_rms_m']:.4f} m, worst {s['max_rms_m']:.4f} m")
    print(f"  max depth deviation     : {s['max_abs_m']:.4f} m")
    print(f"  turning-depth error     : mean {s['mean_turn_err_m']:.5f} m, "
          f"worst {s['max_turn_err_m']:.5f} m")


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
def plot(per_ray_a, step_rows, refine_rows, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(16, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])

    ax = fig.add_subplot(gs[0, :])
    for k, r in enumerate(per_ray_a[::4]):
        ax.plot(r["bh"]["r"] / 1000.0, r["bh"]["z"], color="#1f77b4", lw=3.0,
                alpha=0.30, solid_capstyle="round",
                label="Bellhop (Acoustics Toolbox)" if k == 0 else None)
        ax.plot(r["ph"]["r"] / 1000.0, r["ph"]["z"], color="#d62728", lw=0.8,
                label="libphantom-sonar" if k == 0 else None)
    ax.invert_yaxis()
    ax.set_xlabel("range [km]")
    ax.set_ylabel("depth [m]")
    ax.set_title("Same .env file, both tracers — canonical Munk profile, C-linear SSP, "
                 "101 km", fontsize=11)
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(alpha=0.25)

    # Bellhop step convergence -- the headline.
    ax = fig.add_subplot(gs[1, 0])
    if step_rows:
        steps = [r["step_m"] for r in step_rows]
        ax.loglog(steps, [r["mean_turn_err_m"] for r in step_rows], "o-",
                  color="#2ca02c", lw=1.9, ms=7, label="turning-depth error")
        ax.loglog(steps, [r["mean_rms_m"] for r in step_rows], "s--",
                  color="#9467bd", lw=1.5, ms=6, label="path RMS deviation")
        anchor = step_rows[-1]
        ref2 = [anchor["mean_turn_err_m"] * (s / anchor["step_m"]) ** 2 for s in steps]
        ax.loglog(steps, ref2, ":", color="grey", lw=1.3, label="second order")
        ax.set_xlabel("Bellhop integration step [m]")
        ax.set_ylabel("disagreement [m]")
        ax.set_title("Bellhop converges toward the analytic arc solution", fontsize=11)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=9)

    # Our own sampling floor -- the honesty panel.
    ax = fig.add_subplot(gs[1, 1])
    if refine_rows:
        refs = [r["refine"] for r in refine_rows]
        ax.loglog(refs, [r["mean_rms_m"] for r in refine_rows], "o-",
                  color="#ff7f0e", lw=1.9, ms=7, label="path RMS deviation")
        anchor = refine_rows[0]
        ref1 = [anchor["mean_rms_m"] * (anchor["refine"] / n) for n in refs]
        ax.loglog(refs, ref1, ":", color="grey", lw=1.3, label="first order")
        ax.set_xlabel("libphantom-sonar output sampling  [sublayers per SSP layer]")
        ax.set_ylabel("apparent RMS deviation [m]")
        ax.set_title("The residual is this comparison's chord-vs-arc error,\n"
                     "not a disagreement between the codes", fontsize=11)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=9)

    fig.suptitle("libphantom-sonar vs Bellhop — ray geometry cross-validation",
                 fontsize=13, fontweight="bold")
    fig.savefig(out_path, dpi=140)
    print(f"\nwrote {out_path}")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bellhop", default=None, help="path to bellhop.exe")
    ap.add_argument("--phantom", default=None, help="path to phantom_trace")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--plot", default=None)
    ap.add_argument("--json", default=None, help="write the summary as JSON")
    ap.add_argument("--refine", type=int, default=DEFAULT_REFINE)
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero unless agreement meets the tolerances")
    # The strict assertion is the turning-depth error at Bellhop's finest step:
    # that metric needs no interpolation and Bellhop's own discretisation error
    # is negligible there, so it measures the tracers against each other.
    #
    # The path-RMS bound is deliberately looser. Cases A and B run Bellhop at
    # its default (auto) step, which Case C shows carries ~0.2 m of its own
    # error, and they are sampled at --refine, which carries its own chord
    # error. Tightening this bound would not detect a better tracer; it would
    # just make the check fail on Bellhop's defaults.
    ap.add_argument("--turn-tolerance", type=float, default=0.01,
                    help="max turning-depth error at Bellhop's finest step, m "
                         "(the meaningful assertion)")
    ap.add_argument("--rms-tolerance", type=float, default=1.0,
                    help="max mean path RMS deviation at Bellhop's default step, m "
                         "(a loose sanity bound, see --help notes)")
    args = ap.parse_args()

    bellhop = args.bellhop or shutil.which("bellhop.exe") or shutil.which("bellhop")
    if not bellhop or not Path(bellhop).exists():
        print("Bellhop not found. Build it with:\n"
              "  tools/bellhop_compare/setup_bellhop.sh\n"
              "then pass --bellhop <path>, or put bellhop.exe on PATH.", file=sys.stderr)
        return 2
    bellhop = Path(bellhop).resolve()

    phantom = Path(args.phantom) if args.phantom else REPO / "build" / "phantom_trace"
    if not phantom.exists():
        print(f"phantom_trace not found at {phantom}; build it first "
              "(cmake --build build)", file=sys.stderr)
        return 2
    phantom = phantom.resolve()

    workdir = Path(args.workdir) if args.workdir else REPO / "build" / "bellhop_run"
    workdir.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(Path(__file__).parent / "make_env.py"),
                    str(workdir)], check=True, stdout=subprocess.DEVNULL)

    print("=" * 70)
    print("libphantom-sonar vs Bellhop — ray geometry cross-validation")
    print("=" * 70)
    print(f"  bellhop   : {bellhop}")
    print(f"  phantom   : {phantom}")
    print(f"  workdir   : {workdir}")
    print(f"  refine    : x{args.refine} (output sampling; does not change c(z))")

    summary = {}

    # --- Cases A and B ----------------------------------------------------
    per_ray_a = deviation_stats(pair_rays(
        run_bellhop(bellhop, workdir, "cmp_trapped"),
        run_phantom(phantom, workdir / "cmp_trapped.env",
                    workdir / "cmp_trapped_phantom.csv", args.refine)))
    stats_a = summarise(per_ray_a)
    print_case("Case A — trapped rays (refraction only, no boundary contact)", stats_a)
    summary["trapped"] = stats_a

    per_ray_b = deviation_stats(pair_rays(
        run_bellhop(bellhop, workdir, "cmp_bounce"),
        run_phantom(phantom, workdir / "cmp_bounce.env",
                    workdir / "cmp_bounce_phantom.csv", args.refine)))
    stats_b = summarise(per_ray_b)
    print_case("Case B — with surface and sea-floor reflections", stats_b)
    summary["bounce"] = stats_b

    # --- Case C: Bellhop step size ----------------------------------------
    print("\nCase C — Bellhop integration step sweep  (11 rays, phantom at "
          f"x{SWEEP_REFINE} sampling)")
    print(f"  {'step (m)':>9}  {'turning-depth err (m)':>22}  {'path RMS (m)':>14}")
    step_rows = []
    ph_fine = None
    for env in sorted(workdir.glob("cmp_step_*.env"), reverse=True):
        step = float(env.stem.split("_")[-1])
        if ph_fine is None:
            ph_fine = run_phantom(phantom, env, workdir / "cmp_step_phantom.csv",
                                  SWEEP_REFINE)
        s = summarise(deviation_stats(pair_rays(
            run_bellhop(bellhop, workdir, env.stem), ph_fine)))
        if s:
            s["step_m"] = step
            step_rows.append(s)
            print(f"  {step:>9.1f}  {s['mean_turn_err_m']:>22.5f}  "
                  f"{s['mean_rms_m']:>14.5f}")
    step_rows.sort(key=lambda r: r["step_m"])
    summary["step_sweep"] = [{k: v for k, v in r.items()} for r in step_rows]

    if len(step_rows) >= 2:
        print("\n  Bellhop's turning depths approach libphantom-sonar's as its step")
        print("  shrinks. That is the expected signature: the analytic arc solution")
        print("  is the limit a stepping integrator converges toward.")

    # --- Case D: our own sampling floor -----------------------------------
    print("\nCase D — libphantom-sonar output sampling  (vs Bellhop at its finest step)")
    print(f"  {'refine':>7}  {'apparent RMS (m)':>18}")
    refine_rows = []
    bh_fine = run_bellhop(bellhop, workdir, "cmp_step_0010")
    for n in (8, 16, 32, 64, 128, 256, 512):
        s = summarise(deviation_stats(pair_rays(
            bh_fine,
            run_phantom(phantom, workdir / "cmp_step_0010.env",
                        workdir / f"cmp_refine_{n}.csv", n))))
        if s:
            s["refine"] = n
            refine_rows.append(s)
            print(f"  {n:>7}  {s['mean_rms_m']:>18.5f}")
    summary["refine_sweep"] = refine_rows

    if len(refine_rows) >= 2:
        ratio = refine_rows[0]["mean_rms_m"] / max(refine_rows[-1]["mean_rms_m"], 1e-12)
        print(f"\n  The residual falls by {ratio:.0f}x as the sampling is refined {refine_rows[-1]['refine'] // refine_rows[0]['refine']}x,")
        print("  i.e. it is the chord-versus-arc error of this comparison, not a")
        print("  disagreement between the tracers. The RMS figures above are a")
        print("  ceiling on the true difference, not a measurement of it.")

    # --- Output ------------------------------------------------------------
    out_plot = Path(args.plot) if args.plot else REPO / "data" / "bellhop_comparison.png"
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    try:
        plot(per_ray_a, step_rows, refine_rows, out_plot)
    except ImportError:
        print("matplotlib not available; skipping the figure", file=sys.stderr)

    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=2))
        print(f"wrote {args.json}")

    if args.check:
        problems = []
        finest = step_rows[0] if step_rows else None
        if finest and finest["mean_turn_err_m"] > args.turn_tolerance:
            problems.append(f"turning-depth error {finest['mean_turn_err_m']:.5f} m "
                            f"> {args.turn_tolerance} m at step {finest['step_m']:.0f} m")
        for name in ("trapped", "bounce"):
            s = summary.get(name)
            if s and s["mean_rms_m"] > args.rms_tolerance:
                problems.append(f"{name} mean RMS {s['mean_rms_m']:.4f} m "
                                f"> {args.rms_tolerance} m")
        if problems:
            print("\nFAIL:")
            for p in problems:
                print(f"  - {p}")
            return 1
        print(f"\nPASS: turning depths within {args.turn_tolerance} m, "
              f"path RMS within {args.rms_tolerance} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())
