# Contributing

## Before you open a PR

Read [`EXPORT_NOTICE.md`](EXPORT_NOTICE.md). Short version: everything here must
be traceable to published literature or to your own measurements on hobbyist
hardware. If you work in the defence sector, get written publication approval
from your employer first — personal projects in an overlapping technical field
are commonly covered by IP and publication clauses.

## Non-negotiables

These are what the project is *for*. A change that breaks one of them will be
rejected however good it otherwise is.

1. **Zero dependencies.** Standard C++ library only. No Boost, no Eigen, no
   GoogleTest. The test harness is 90 lines in `tests/framework.hpp` precisely
   so the claim stays true. Python is allowed in `tools/` — that is tooling, not
   library code.
2. **Zero allocation in `include/` and `src/`.** No `new`, no `malloc`, no
   `std::vector`, no `std::string`. Callers own all buffers. CI enforces this
   twice: `tools/audit_no_alloc.py` on the sources, and `nm` over the built
   archive, which is the check that actually cannot be fooled.
3. **No `-ffast-math`**, ever. It permits reassociation of the Snell invariant
   and breaks bit-reproducibility across targets.
4. **No virtual dispatch in hot paths.** Templates and free functions.
5. **New physics needs a closed-form test.** See below.

## Testing standard

The bar is not "the output looks right" or "it matches what it printed
yesterday". A regression baseline only proves the code did not change.

Every physics change must be verified against something independent:

- an **analytic solution** (isovelocity → straight lines; constant gradient →
  exact circular arcs; known turning depths),
- an **invariant** that must hold regardless of implementation (`cos θ / c`
  constant; refining a linear profile must not move the answer),
- or a **published reference value** with its citation.

Look at `tests/test_ray_tracer.cpp` for the pattern. Each test states in a
comment what closed form it checks and why that form is the right one.

Tolerances must be precision-aware: use `pt::tol(double_tol, float_tol)` so the
`float` build does not force the `double` build's tolerance to be loosened —
otherwise a genuine double-precision regression slips through unnoticed.

## Documentation standard

- Any new formula goes in `docs/math_spec.md` with its derivation, validity
  range and citation.
- Any new claim about accuracy or speed goes in `docs/validation.md` **with the
  measurement that supports it**, and in the right section: verified, or only
  self-consistent. Moving something into the "verified" column requires an
  external reference, not internal agreement.
- If you find a limitation, document it rather than quietly narrowing the test.
  The honest limitations section is the most valuable part of this repository.

## Before pushing

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPHANTOM_WERROR=ON
cmake --build build && ./build/phantom_tests

# float precision
cmake -S . -B build-float -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DPHANTOM_REAL_FLOAT=ON -DPHANTOM_WERROR=ON
cmake --build build-float && ./build-float/phantom_tests

# sanitizers
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPHANTOM_SANITIZE=ON
cmake --build build-asan && ./build-asan/phantom_tests

# allocation audit
python3 tools/audit_no_alloc.py
```

Clang is worth running too — it catches implicit float widening that GCC's
`-Wdouble-promotion` misses.

## Style

Follow the surrounding code: 4 spaces, 100 columns, `snake_case` functions,
`PascalCase` types, trailing `_` on private members. Comments explain *why* a
choice was made, especially where the obvious approach is wrong — the branch
selection in `arc_time()` is the model to imitate.
