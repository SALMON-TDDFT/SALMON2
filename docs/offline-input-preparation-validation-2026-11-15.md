# Offline input preparation validation — 2026-11-15

Implemented `salmonpy.preparation.prepare_inputs` and `verify_inputs`, plus the
verification CLI and executable `examples.prepare_inputs`. Preparation reuses
MaterialRegistry/resolve_materials and copies an explicit local asset mapping.
Initial scope is configuration/file integrity, not native execution or automatic
scientific dependency resolution. ASE and grid-spacing milestones remain open.

Evidence (local Apple Silicon/macOS, Python 3.14):

- Initial nine focused tests failed because the module was missing, then passed.
- Review found Unicode-normalization and implicit manifest-name collisions.
  A Unicode regression reproduced callback execution before rejection and passed
  after comparing normalized names. Original filename spelling is retained.
- Final focused suite: 11 tests, one case-sensitive-filesystem test skipped on
  this case-insensitive volume. Linux filesystem execution was not performed.
- Full source suite: 328 tests, pass, 1 skip, 8.348 s.
- Source distribution → wheel built successfully with two build jobs. Fresh venv
  installed the wheel, then isolated-import suite: 328 tests, pass, 2 skips,
  8.346 s (case-sensitive filesystem and optional plotting dependency).
- Installed sample ran with Python socket audit events forbidden. Original
  coefficient file was deleted and bundle moved. Manifest hash verification and
  recipe resolution from the transferred location matched resolved.json.
- CLI verified the relocated bundle. Installed module, example, paired guides
  and AI JSON were byte-equal to source. No solver or network acquisition ran.
- AI manual synchronization and git whitespace checks passed. No Actions started.

Artifacts outside checkout:

- `/tmp/salmonpy-preparation-v2-source-tests-20261115.log`
- `/tmp/salmonpy-preparation-v2-build-20261115.log`
- `/tmp/salmonpy-preparation-v2-installed-tests-20261115.log`
- `/tmp/salmonpy-preparation-offline-smoke-20261115/validation.json`

The first pre-review wheel is not the final validated artifact. No claim is made
about production structure interpretation, physical convergence, actual HPC
installation, MPI/OS portability of runtimes, or automatic dependency completeness.
User callbacks remain trusted code; file verification does not call them.
