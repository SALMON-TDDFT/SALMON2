# Structure/grid preparation validation — 2026-11-15

Implemented optional ASE StructureInput, immutable GridSelection and
maximum-width grid resolution for positive axis-aligned cells. Existing ASE DFT/
TDDFT entry points accept grid_spacing plus explicit grid_spacing_unit, or the
original grid counts. No native numerical code or integration order was changed.

- Nine new tests: detached arrays/nested metadata, CIF/extxyz frame and SHA256,
  preservation of expanded cell, ceiling/alignment/Bohr conversion, input and
  geometry rejection, mixed/partial occupancy rejection at loader and DFT boundary,
  and generated DFT input text equality with explicit equivalent counts.
- Initial module tests failed before implementation. Review regressions reproduced
  shallow ASE info sharing and silently accepted mixed occupancy before fixes.
- Full source suite: 337 tests, pass, 1 skip (case-sensitive filesystem), 8.410 s.
- Local source distribution → final wheel built with two jobs. Fresh venv isolated
  installed suite: 337 tests, pass, 2 skips (filesystem + optional plotting), 8.301 s.
- Installed preparation example ran with socket audit events forbidden. A supplied
  four-atom supercell (10×6×7 Å) resolved at max 0.4 Å to 25×15×18 points,
  actual widths 0.4, 0.4, 0.3888888888888889 Å. Original file was deleted, bundle
  moved, hashes checked and recipe re-resolved identically from the transferred
  path. Material ID film-Si, atom count and cell remained intact.
- Installed structure/ASE modules, example, paired guide and AI JSON match source
  bytes. Guide synchronization and git diff whitespace checks pass.
- ASE emits a NumPy 2.5 shape deprecation warning during its own Atoms operations;
  tests complete successfully. No full solver or production physics run occurred.

Artifacts outside checkout:
`/tmp/salmonpy-structure-source-tests-20261115.log`,
`/tmp/salmonpy-structure-v2-build-20261115.log`,
`/tmp/salmonpy-structure-installed-tests-20261115.log`,
`/tmp/salmonpy-structure-offline-smoke-20261115/validation.json`.

The explicit-count DFT path retains existing geometry behavior; new width
resolution rejects nonorthogonal/rotated cells. Periodicity is checked at the
ASE DFT boundary and in the example. CIF symmetry expansion belongs to parsing,
not implicit supercell construction. Source hashes cover the whole file. Direct
ASE width calls retain adopted counts only; save GridSelection.as_dict when
request/rounding provenance is required. Preparation example saves that report.

No resource/convergence guarantee, MPI decomposition, automatic potential choice,
automatic backend registration or scientific input completeness is claimed.
All tests are local Apple Silicon/macOS; no Actions were run. New work is uncommitted.
