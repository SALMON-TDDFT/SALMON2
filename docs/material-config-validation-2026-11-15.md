# Material configuration catalog validation — 2026-11-15

Scope: `salmonpy.config`, schema `salmonpy.materials/1`, explicit custom resolver
registry, immutable resolved definitions and structurally validated initial-state
references. Configuration only; no native numerical changes in this milestone.
Earlier uncommitted GS/SCF/native changes are separate work.

- Added 12 focused tests: exact IDs, duplicate JSON keys, unknown fields and
  references before callbacks, source/prepare exclusivity, finite/JSON-safe values,
  cycles/Unicode/large integers, immutable copies, backend registration, resolver
  errors/outputs, Python/JSON parity, file base directory and definition hashes.
- First run failed because the module was absent. Review identified oversized
  Python integer serialization as an inconsistent error path; a reproducing test
  failed before the fix and passed afterward.
- Source suite: 317 tests, all pass, 8.277 s.
- Local source distribution → wheel: successful, Python 3.14 / Apple Silicon /
  macOS 26 wheel. Fresh venv installed that wheel; isolated import test suite:
  317 tests, pass, one optional plotting dependency skip, 8.352 s.
- Installed module, paired guides, example and AI JSON are byte-equal to source.
  Installed example output matches source example.
- AI guide synchronization and `git diff --check` pass.
- Read-only peer review found no important functional issue within stated scope.

Artifacts are outside checkout under `/tmp/salmonpy-material-config-*20261115*`:
source-v2 log `salmonpy-material-config-v2-source-tests-20261115.log`,
wheel build log `salmonpy-material-config-v2-build-20261115.log`,
installed suite log `salmonpy-material-config-installed-tests-20261115.log`.
The first suite run caught a missing required AI task metadata field; corrected
before the final build. The initial pre-fix wheel is not the validated artifact.

No Actions were started. No new MPI or solver-physics validation was performed.
No external structure/GS data were loaded. ASE/file adapters, grid-spacing
resolution, built-in parameter schemas, execution capabilities and ID propagation
through point states/output/restart remain checklist M2–M6. Initial-state recipe
contents are not yet certified compatible. Hashes cover resolved parameter values,
not automatically referenced files or physical accuracy.
