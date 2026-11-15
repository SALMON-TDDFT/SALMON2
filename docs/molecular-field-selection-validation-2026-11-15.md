# Molecular sampled-drive selection validation

The pure `select_molecular_field` API implements explicit provided_drive and
own-source vacuum subtraction. No FDTD/TDDFT update or native code changed.
Provenance is caller-declared; matching identifiers do not certify separate
solver runs. Subtraction also removes vacuum radiation reaction. This is not
validation of two-way molecular dynamics, emission or a local-field model.

Four regression tests verify own-only cancellation, retention of additive
external/other-molecule fields, mismatch refusal for IDs/time/position/units/
sampling/source histories and footprints, immutable detached metadata and
numerical overflow. Artificial additive fields are used, not native propagation.
Independent read-only review found no required corrections.

Local Apple Silicon/macOS, Python 3.14. Source test log:
`/tmp/salmonpy-molecular-field-tests.log`.
Distribution built explicitly from sdist into
`/tmp/salmonpy-molecular-field-dist-20261115`; build log:
`/tmp/salmonpy-molecular-field-build.log`.
No Actions, remote publication or physical convergence sweep was run.

Results: source 435 tests/one skip; installed wheel 435 tests/two skips, no
failures. The new example ran from installed Python with `-I`, source imports
excluded. Output: `/tmp/salmonpy-molecular-field-example.json`; installed log:
`/tmp/salmonpy-molecular-field-installed-tests.log`. A fresh environment reused
existing dependencies offline via a test-only path; dependency resolution was
not exercised. No new native/MPI runs are claimed for this pure selector.
