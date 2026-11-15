# Explicit real-space grid initialization — 2026-11-15

The explicit fast-HSE profile splits input capture (grid_config_from_input) from
configuration-only grid bounds/halo/index generation (init_grid_whole_config).
The legacy wrapper captures existing input and preserves native termination on
failure. Cell/k/pseudo/XC/HSE input extraction is NOT complete; full TDDFT retains
its common-grid/k/XC restriction. Coordinates and MPI maps are unchanged and
remain the responsibility of subsequent initialization routines.

The explicit entry validates geometry, periodicity, finite values, consistency,
integer conversion/halo bounds and total grid-point integer range before committing
state. It builds trial index arrays and uses move_alloc only on success. Input
errors return status1; handled allocation failure returns status2. Integer limits
are not a RAM budget. No simulated allocation failure was used in the test.

## Local verification

- Red probe failed on missing config type/entry in the previous build:
  `/tmp/salmonpy-grid-config-red-20261115/build.log`.
- Final native probe includes a provenance-attributed frozen original routine
  from fast-PR c9e8ce431f43eefd4a88998134adb4796f8ec837. This avoids comparing only
  two wrappers of the same implementation. Even/odd periodic/isolated grids,
  point-count/spacing mode, independent shapes, old-global poisoning, Maxwell
  consistency and invalid-call retention pass. Independent expected bounds,
  halo widths and the original near-epsilon behavior are checked (length8-5e-11,
  spacing1: periodic count8; isolated count9). This preserves original rounding,
  not a claim of physical usefulness for that boundary configuration.
- Full-context regression:8 cases, OMP2, different dt/OpenMP/capacities; separate
  reference, reverse order and close/reopen match bit-for-bit. Grid remains common.
- Legacy native-driver regression:32 cases vs immutable fast-PR reference;
  formatted current and energy differences are zero.
- Source and all application-object eigensolver audits pass before linking.
  Inventory hashes match fresh sources; no new executable native line exceeds132.
- Python source suite305 passed. Final sdist-to-wheel and new installed-wheel
  tests pass (305 tests, one optional plotting skip). Installed preparation
  passes audits; every replacement byte and augmented probe/tool match source.
- Independent read-only review found no implementation blocker. Its test-quality
  concern was addressed with the frozen original and independent boundary checks.
- English/Japanese/AI guides are synchronized. No Actions launched. Changes are
  uncommitted at report creation.

## Limits

Local synthetic Apple Silicon tests only. The low-level isolated-grid check adds
no isolated TDDFT capability. No heterogeneous full solver, MPI/GPU/concurrent
calls, AVX/SVE, force/microAc validation or physical-convergence claim.

## Receipts

- `/tmp/salmonpy-grid-config-frozen-probe-20261115/result.json`; SHA256 `73e9c582ef8b590c9e70644ad15f9a2ee08ba9fcf5ad9e32bf1ffe99cc8c7ff3`.
- `/tmp/salmonpy-grid-config-context-20261115/result.json`; SHA256 `44088174847a5cd12970769b2448f4a120269881ee9e227be88c13a47af7faa8`.
- `/tmp/salmonpy-grid-config-default-20261115/result.json`; SHA256 `176cbbb0173b2a680b6823934ef2c06c6b9735dbcdb97532bd45c469d887ae8c`.
- `/tmp/salmonpy-grid-config-build-20261115.json`; SHA256 `a268b4bae5b4be118b238f7658be13417f2a6e8e743aa3b50a31c270692c52b1`.

Logs: `/tmp/salmonpy-grid-config-python.log`,
`/tmp/salmonpy-grid-config-final-packaging.log`,
`/tmp/salmonpy-grid-config-installed-tests.log`,
`/tmp/salmonpy-grid-config-wheel-prepare.log`.
