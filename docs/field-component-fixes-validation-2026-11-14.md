# Field-component fixes: local validation, 2026-11-14

Base: `salmonpy-field-components`, commit `0575a770`.
Environment: local macOS arm64, Python 3.14, GNU Fortran and Apple Accelerate.
This is a local engineering regression, not a new physical-convergence study.

Before fixes:
- Built an unrepaired wheel with pip wheel. The bundled dylib SHA256 differed
  from its .sha256 sidecar after CMake installation.
- New endpoint test failed for all three uniform field models: 11 samples could
  not initialize a 10-step state.
- New boundary test showed incident waves accepted with PEC/periodic left edges.

After fixes:
- CMake computes the installed native checksum after install-name/RPATH changes.
  Build-tree hashing and source/object audits remain; repair still rehashes.
- Uniform fields evaluate only the consumed endpoints 0..steps.
- Explicit incident drives require a left ABC boundary; non-ABC runs without an
  incident drive remain accepted.
- Built sdist then wheel using `python -m build --no-isolation`.
- Unrepaired wheel native hash passes. Installed outside checkout: 135 tests pass.
- delocate-repaired wheel native hash passes. Separate install: 135 tests pass.
- Existing wheel RECORD/hash repair unit test passes; manual tests pass.

CI now checks and installs the unrepaired sdist-built wheel on Linux before
uploading it, in addition to existing repaired wheel tests. This CI change has
not been pushed or run remotely at the time of this record. Windows/Linux local
builds were not performed. Existing external native MPI/TDDFT computations were
not rerun: no native physics kernel or runtime patch changed.

The proposed shared coupling schedule is documented separately in
`plans/2026-11-14-unified-coupling-design.md`; it is not implemented by these fixes.
