#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Turn the fetched WOA23 CSVs into a header the tests can compile in.

The library allocates nothing and reads no files, so the test suite cannot go
looking for data at run time without giving up both properties. Embedding the
profiles keeps the tests hermetic: they measure the same water on every machine,
with no network, no working directory and no I/O error path to handle.

    python3 tools/make_profile_header.py data/profiles tests/data/woa_profiles.hpp
"""
import pathlib
import sys


def main():
    src = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "data/profiles")
    dst = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "tests/data/woa_profiles.hpp")
    files = sorted(src.glob("*.csv"))
    if not files:
        raise SystemExit("no CSVs in {}".format(src))

    out = ["// SPDX-License-Identifier: Apache-2.0",
           "// Real temperature/salinity profiles, embedded for the test suite.",
           "//",
           "// World Ocean Atlas 2023 (decav annual mean, 1-degree grid), a US Government",
           "// work in the public domain. Fetched by tools/fetch_woa_profiles.py, converted",
           "// by tools/make_profile_header.py. Do not edit by hand.",
           "//",
           "// A climatology, not a cast: a decadal average over a 1-degree cell. It has the",
           "// SHAPE of real water, which is what these tests need, but no eddy, no internal",
           "// wave and no diurnal surface layer survives the averaging.",
           "#ifndef PHANTOM_TEST_WOA_PROFILES_HPP",
           "#define PHANTOM_TEST_WOA_PROFILES_HPP",
           "",
           "#include <cstddef>",
           "",
           "namespace woa23 {",
           "",
           "struct Level { double depth_m; double temperature_c; double salinity_psu; };",
           "",
           "struct Site {",
           "    const char* name;",
           "    const char* note;",
           "    double latitude_deg;",
           "    double longitude_deg;",
           "    const Level* levels;",
           "    std::size_t level_count;",
           "};",
           ""]

    sites = []
    for f in files:
        name = f.stem
        note, lat, lon = "", 0.0, 0.0
        rows = []
        for line in f.read_text().splitlines():
            if line.startswith("# site:"):
                continue
            if line.startswith("# grid cell centre:"):
                body = line.split(":", 1)[1]
                lat = float(body.split("lat")[1].split(",")[0])
                lon = float(body.split("lon")[1])
                continue
            if line.startswith("#"):
                if not note and "Atlas" not in line and "source" not in line \
                        and "fetched" not in line and "Government" not in line \
                        and "NOTE" not in line and "degree grid" not in line:
                    note = line[1:].strip()
                continue
            if line.startswith("depth"):
                continue
            rows.append(tuple(float(x) for x in line.split(",")))
        ident = name.replace("-", "_")
        out.append("// {}".format(note))
        out.append("inline constexpr Level k_{}[] = {{".format(ident))
        for z, t, s in rows:
            out.append("    {{{:8.1f}, {:9.4f}, {:9.4f}}},".format(z, t, s))
        out.append("};")
        out.append("")
        sites.append((name, note, lat, lon, ident, len(rows)))

    out.append("inline constexpr Site kSites[] = {")
    for name, note, lat, lon, ident, n in sites:
        out.append('    {{"{}",'.format(name))
        out.append('     "{}",'.format(note.replace('"', '\\"')))
        out.append("     {:.1f}, {:.1f}, k_{}, {}}},".format(lat, lon, ident, n))
    out.append("};")
    out.append("")
    out.append("inline constexpr std::size_t kSiteCount = sizeof(kSites) / sizeof(kSites[0]);")
    out.append("")
    out.append("}  // namespace woa23")
    out.append("")
    out.append("#endif  // PHANTOM_TEST_WOA_PROFILES_HPP")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(out) + "\n")
    print("{} sites, {} levels -> {}".format(len(sites), sum(s[5] for s in sites), dst))


if __name__ == "__main__":
    main()
