# GS mixing configuration isolation — 2026-11-16

## Scope and change

Internal retained GS now captures method, initial mixrate/alpha_mb/beta_p,
Pulay memory, automatic adjustment switch/ratio and initial no-mixing steps in
`s_gs_run_config%mixing`. The state copies that value before ordinary or DC
initialization. `s_mixing%input` supplies method dispatch, Pulay memory, adaptive
halving and resets; mutable rates/residual/history are separately owned.
No mixing formulas, reductions, update ordering or fixed 21-slot storage changed.
Legacy initialization captures globals; separate DFT-MD reset gating remains global.

## Local evidence

Apple Silicon/macOS, GNU Fortran, OpenBLAS, static Libxc, CPU serial HSE ON;
-O2 -fallow-argument-mismatch, without unlimited-line-length flags. Build jobs=2;
GS OMP/BLAS=1. No network download, GitHub Actions or physical convergence campaign.

- Red: prior metric build rejects probe because the mixing configuration field
  does not exist (`/tmp/salmonpy-gs-mixing-red-20261116/build.log`).
- Build: `/tmp/salmonpy-gs-mixing-build-20261116.json`; source/object forbidden
  eigensolver checks passed. Reference:
  `/tmp/salmonpy-gs-metric-reference-20261116.json` (original pinned fast HSE).
- 13 final lifecycle cases: simple/Broyden/Pulay each with pre-SCF disabled and
  consecutive counts 1/2/3, plus simple_potential/norm_pot without pre-SCF.
  Density methods enable auto adjustment (ratio=0 to exercise halving); simple
  also uses two initial no-mixing iterations. A coefficients are .35/.65/.55
  and Pulay memory=2; B coefficients .175/.325/.275 and memory=1.
- Each A/B output, eigenvalues, k data and **every restart payload** exactly
  matched a separate original solver run under its respective input. Printed
  per-iteration residual labels and values matched in sequence.
- Legacy globals are poisoned after capture. Invalid method, nonfinite rate and
  out-of-capacity memory are rejected before activation. A caller's config is
  changed after A's run. Retained A's arrays, mixing configuration, adaptive
  coefficients and previous residual remain byte-identical before/after B's
  run and close. A separate two-coefficient check verifies halving/reset leaves
  the other state unchanged.
- Final case records: `/tmp/salmonpy-gs-mixing-final-{simple,broyden,pulay,simple_potential}-20261116/result.json`.
  Both reference hash sets and residual traces are retained. Temporary density/wavefunction payloads are removed after comparison.

## Broader checks

- Source Python suite: 435 tests, 1 skip (`/tmp/salmonpy-gs-mixing-python-20261116.log`).
- Explicit sdist-to-wheel build succeeded; installed wheel suite: 435 tests,
  2 skips (`/tmp/salmonpy-gs-mixing-final-installed-tests-20261116.log`). Existing local
  dependency environment reused; this does not test online dependency resolution.
- Fast-HSE TDDFT interleaving regression: all 8 cases matched original output,
  2/3 k points, Taylor4/PTCN/factor modes, different dt/OpenMP/capacity,
  12 steps (`/tmp/salmonpy-gs-mixing-tddft-20261116`).
- Human guides and embedded AI guide copies synchronized.

## Limits

A tiny Si diagnostic fixture (8³, 2 k points, loose stopping threshold) checks
migration and state ownership, not converged optical or material properties.
Mixed methods within one process, spinor simple_dm, DC, MPI, geometry optimization,
DFT-MD, GPU and other platforms were not newly executed in this increment.
Other material/input/parallel/timer globals remain; no general simultaneous or
arbitrary different-material GS API is claimed. Settings are internal Fortran;
no new Python public entry point is exposed.
