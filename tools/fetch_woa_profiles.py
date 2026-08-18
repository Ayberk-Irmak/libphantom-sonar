#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Fetch real temperature/salinity profiles from the World Ocean Atlas 2023.

Every profile shipped in data/profiles/ was produced by this script. It talks to
NOAA NCEI's OPeNDAP server and pulls ONE grid cell rather than the 100+ MB
global field, so re-running it is cheap and the provenance of every number in
the repository is a URL you can fetch yourself.

    python3 tools/fetch_woa_profiles.py --out data/profiles

WOA23 is a US Government work and is in the public domain. Cite it as:

    Locarnini, R. A. et al. (2024). World Ocean Atlas 2023, Volume 1:
    Temperature. NOAA Atlas NESDIS 89.
    Reagan, J. R. et al. (2024). World Ocean Atlas 2023, Volume 2:
    Salinity. NOAA Atlas NESDIS 90.

A NOTE ON WHAT THIS DATA IS. WOA23 is a CLIMATOLOGY -- a decadal average over
all observations in a 1-degree cell. It is not a cast, not a forecast, and not
what a sonar would see on any particular day. Real profiles have interleaving,
eddies, internal waves and a diurnal surface layer, none of which survive
averaging. Use these to test that the library handles the SHAPE of real water;
do not use them to predict a real range.
"""

import argparse
import pathlib
import sys
import urllib.request

BASE = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA"
# decav = the 1955-2022 "decadal average" climatology; t00/s00 = annual mean.
URLS = {
    "t": BASE + "/temperature/netcdf/decav/1.00/woa23_decav_t00_01.nc",
    "s": BASE + "/salinity/netcdf/decav/1.00/woa23_decav_s00_01.nc",
}

# Chosen so the set spans the sound-speed structures the library must handle,
# not so it spans the globe. The comment on each is what it is FOR.
SITES = [
    # Every note below was written AFTER looking at the data it describes. An
    # earlier draft called the Norwegian Sea upward-refracting and put the North
    # Atlantic axis at 1100 m; neither survived contact with the profiles, and
    # shipping a description the shipped data contradicts is worse than shipping
    # no description. Axis depths quoted here are what the fetched profile gives.
    ("levantine",   34.5,  32.5,
     "Eastern Mediterranean. The saltiest water here (39.3 PSU at the surface), "
     "which lifts sound speed everywhere and makes the salinity term matter "
     "instead of being a small correction. Channel axis near 375 m."),
    ("aegean",      38.5,  25.5,
     "Aegean. Only 650 m deep in this cell, strongly stratified and bounded on "
     "every side -- the case where a range-INDEPENDENT ocean is least defensible, "
     "included so that limitation has something real to fail against. Axis 175 m."),
    ("black-sea",   43.5,  34.5,
     "Black Sea. Near-fresh at the surface (18.2 PSU) under a permanent "
     "halocline, so the profile is driven by SALINITY rather than temperature -- "
     "the reverse of every other site here. Axis 55 m, very shallow."),
    ("n-atlantic",  45.5, -40.5,
     "North Atlantic. The textbook deep sound channel, axis near 950 m -- the "
     "closest real analogue to the Munk profile the tests use."),
    ("norwegian",   70.5,   5.5,
     "Norwegian Sea. Cold throughout and weakly stratified: temperature falls "
     "monotonically from 6.6 C to below zero, so the near-surface water refracts "
     "DOWNWARD and a surface source has a shadow zone under it. Axis 850 m. "
     "(The seasonal surface duct real ships see is averaged away by a "
     "climatology -- which is itself worth seeing.)"),
    ("eq-pacific",   0.5, -150.5,
     "Equatorial Pacific. Warm-topped with a thick thermocline over a deep "
     "channel at 950 m -- the long-range propagation case."),
]



def cell_index(lat_deg, lon_deg):
    """WOA23 1-degree grid: lat[0] = -89.5, lon[0] = -179.5."""
    return int(round(lat_deg + 89.5)), int(round(lon_deg + 179.5))


def fetch_column(var, lat_i, lon_i, timeout):
    url = "{}.ascii?{}_an[0][0:101][{}][{}]".format(URLS[var], var, lat_i, lon_i)
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        text = fh.read().decode("utf-8", "replace")

    values, depths, lat, lon = [], [], None, None
    section = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("{}_an.{}_an".format(var, var)):
            section = "values"
            continue
        if line.startswith("{}_an.depth".format(var)):
            section = "depth"
            continue
        if line.startswith("{}_an.lat".format(var)):
            section = "lat"
            continue
        if line.startswith("{}_an.lon".format(var)):
            section = "lon"
            continue
        if line.startswith("{}_an.".format(var)):
            section = None
            continue
        if not line:
            continue
        if section == "values" and line.startswith("["):
            # "[0][7][0], 18.3" -- the value follows the last comma.
            values.append(float(line.rsplit(",", 1)[1]))
        elif section == "depth":
            depths.extend(float(x) for x in line.split(",") if x.strip())
        elif section == "lat" and lat is None:
            lat = float(line.split(",")[0])
        elif section == "lon" and lon is None:
            lon = float(line.split(",")[0])
    return depths, values, lat, lon


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/profiles", help="output directory")
    ap.add_argument("--timeout", type=float, default=120.0)
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for name, lat, lon, note in SITES:
        lat_i, lon_i = cell_index(lat, lon)
        print("{:<12} lat {:>6.1f} lon {:>7.1f}  (grid {}, {})".format(
            name, lat, lon, lat_i, lon_i), end="", flush=True)
        try:
            zt, t, glat, glon = fetch_column("t", lat_i, lon_i, args.timeout)
            zs, s, _, _ = fetch_column("s", lat_i, lon_i, args.timeout)
        except Exception as exc:                                 # noqa: BLE001
            print("  FAILED: {}".format(exc))
            continue

        # WOA marks land and below-bottom levels with a huge fill value. Keep
        # only levels where BOTH variables are present: a sound speed needs both,
        # and the two fields do not always stop at the same level.
        rows = []
        for i in range(min(len(zt), len(t), len(s))):
            if abs(t[i]) > 1e30 or abs(s[i]) > 1e30:
                break
            if zs[i] != zt[i]:
                raise SystemExit("depth axes differ between t and s")
            rows.append((zt[i], t[i], s[i]))
        if len(rows) < 2:
            print("  no water (land cell?)")
            continue

        path = out / "{}.csv".format(name)
        with path.open("w") as fh:
            fh.write("# World Ocean Atlas 2023, decav annual mean (t00/s00), 1-degree grid\n")
            fh.write("# site: {}\n".format(name))
            fh.write("# {}\n".format(note))
            fh.write("# grid cell centre: lat {:.1f}, lon {:.1f}\n".format(glat, glon))
            fh.write("# source: {}\n".format(URLS["t"]))
            fh.write("#         {}\n".format(URLS["s"]))
            fh.write("# fetched by tools/fetch_woa_profiles.py\n")
            fh.write("# WOA23 is a US Government work, public domain.\n")
            fh.write("# NOTE: a climatology, not a cast. See the script docstring.\n")
            fh.write("depth_m,temperature_c,salinity_psu\n")
            for z, tv, sv in rows:
                fh.write("{:.1f},{:.4f},{:.4f}\n".format(z, tv, sv))
        print("  {} levels to {:.0f} m -> {}".format(len(rows), rows[-1][0], path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
