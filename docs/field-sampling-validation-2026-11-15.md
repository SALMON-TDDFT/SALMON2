# Explicit field sampling validation — 2026-11-15

SampledField owns immutable times, positions and values. sample() requires an
explicit spatial convex stencil and exact/linear temporal policy. sample_packet
uses the separate clock of each required port and returns FieldPacket. Only point
values are supported. Missing history, extrapolation, wrong position or units,
and missing ports reject. Native solver stepping and current feedback are unchanged.

Tests use manufactured affine fields: spatial/temporal exactness, separate A/E
clocks, immutable arrays, missing histories, no implicit H, and invalid stencils.
Independent review found a scalar position tolerance allowed a large x translation
to hide off-plane y errors. Componentwise tolerance and regression fix this case.

Apple Silicon macOS/Python 3.14:
- Source suite before final regression: 415 tests, one skip.
- Final focused sampling tests: four passed.
- Final sdist-to-wheel build/install: successful.
- Final installed full suite: 416 tests, two skips (optional plot and filesystem).
- Installed example field_sampling: manufactured A/E values and clocks match.

Artifacts: /tmp/salmonpy-field-sampling-final-dist-20261115/.
Logs: /tmp/salmonpy-field-sampling-final-build.log and
/tmp/salmonpy-field-sampling-final-installed-tests.log.
No Actions, native rebuilds, MPI changes, or physical convergence sweeps.

This does not complete native history extraction, the H sampling accessor,
clock mapping of returned j, or a self-consistent packet coupling driver.
