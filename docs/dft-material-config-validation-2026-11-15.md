# DFT material configuration and TDDFT bridge validation — 2026-11-15

Added explicit salmonpy.dft/dft configuration resolver and identity-carrying
ConfiguredGroundState. Fixed scientific parameters remain separate from runtime,
SCF stopping controls, threads, timeout and warm-start execution arguments.
Existing DFT/GS/TDDFT native code and time evolution are unchanged.

- Six dedicated tests cover equivalent original DFT input text, grid-width report,
  source hashes and alteration rejection, paths/units/unknown keys, mocked GS→RT
  identity forwarding and cleanup, wrapped warm start, and UPF element mismatch.
- Module support was absent in the initial failing test run. Reviewer-identified
  warm-start forwarding and mislabeled UPF regressions failed before fixes.
  Existing warm-start scientific compatibility checks remain in DFT.run.
- Full source suite: 350 tests, pass, one filesystem skip, 8.394 s.
- Final source distribution→wheel built locally (two jobs). New venv isolated
  installed suite: 350 tests, pass, two skips (filesystem/optional plotting),
  8.391 s. No GitHub Actions.
- Source and installed configured_si examples used the existing audited local
  runtime, executable SHA256
  29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed.
  No full runtime compilation or network acquisition occurred.
- Final installed example ran with Python socket audit events forbidden; staged
  structure, pseudopotential, license and provenance were verified after execution.
  Si starter energy and currents had zero difference from packaged expected values;
  comparator still uses its declared tolerances. Energy: -31.74653905752293 Ha.
- Two TDDFT states received opposite drives in interleaved order; subsequent
  same-initial-state replay was bit-identical. Material ID Si-starter and definition
  digest survived model→GS→initial→states. Grid resolved to 16³. Each process used
  one thread, at most two workers; elapsed 11.905 s in final installed run.
  Peak memory was not measured or hard-limited.
- Installed module, paired guide, sample and AI JSON match source. Manual sync
  and git whitespace checks pass. Early sample syntax error was fixed before
  the numerical runs; pre-final wheels are not the final validation artifacts.

Artifacts outside checkout:
`/tmp/salmonpy-dft-config-source-tests-20261115.log`,
`/tmp/salmonpy-dft-config-final-build-20261115.log`,
`/tmp/salmonpy-dft-config-final-installed-tests-20261115.log`,
`/tmp/salmonpy-dft-config-final-smoke-20261115/validation.json`.
Generated native restart payloads from the three smoke runs were removed after
verification; manifests retain their hashes. Those diagnostic GS directories are
not reusable restart bundles after cleanup. Current/report/input files remain.

Scope: local Apple Silicon/macOS, short PZ Si connection regression with SCF
convergence checked by existing API. Not grid/k/time-converged optics, new Libxc
force validation, or cross-platform/MPI evidence. KY element and PP/XC provenance
are not automatically certified. Generic saved-GS/initial_state_id selection,
point/run IDs, MPI assignment and automatic disk ID propagation remain open.
