# data/

`munk_rays.png` is committed — it is the README figure.

The CSV files it is built from (`rays.csv`, `coverage.csv`, `profile.csv`) are
several megabytes and are regenerated on demand rather than tracked:

```bash
./build/munk_simulation data
python3 tools/plot_rays.py data
```

Real measured T/S profiles (World Ocean Atlas, Argo) will land here in v0.2,
each with a provenance note recording the source, region and date.
