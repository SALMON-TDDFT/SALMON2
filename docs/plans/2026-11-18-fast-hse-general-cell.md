# Fast-PR general-cell initialization

The rotated-cell full-k path initializes, while reduced-k scalar FFT-orbit reuse
fails ACE Hermiticity. Disabling orbit reuse diagnostically makes the identical
full/reduced calculation agree. Scalar orbit gauge, k lookup/rotation and wrap
phase still use componentwise Cartesian box lengths; the spinor path already
uses full cell/Cartesian rotations.

Implement complete replacement exx_k_orbits.f90 in the explicit-routing source
stage, preserving pinned upstream files. Use full cell transformations only for
non-axis-aligned cells; retain the exact old arithmetic on positive diagonal
cells. Export a versioned native capability from the corrected core. Admit
nonsingular general cells following native column norms and grid priorities;
require that core capability for these cells before consuming native lifetime.
Do not reorder primitive axes, alter wfn records, project fields or change dt.

Regression: old core fails rotated reduced-k initialization. New core must match
full-k for rotated/right/left-handed and xy-sheared cells, internal MPI1/2,
independent material groups, direct own save/resume and existing diagonal traces.
Retain native symmetry metric/grid/atom compatibility checks and current signed
permutation scope; general-cell support does not imply arbitrary space groups.
Run source/installed tests and sdist-to-wheel after complete-source layout change.
No Actions, no ordinary-HSE engine, no mandatory physical-convergence campaign.

## Actual-field follow-up

The host interface is Cartesian A[3], not fixed z. Only attached-host zero-impulse
response initialization skips the hardwired z check; the actual vector is checked
before every step using owned Cartesian symmetry. Native batch behavior remains
unchanged. Add an explicit field-core capability for old libraries and x/C2x,
y/C2y,circular-xy/mirror-z full/reduced MPI/save-resume regressions. Incompatible
retained symmetry still rejects before updating; use suitable full-k GS data
instead of projecting the field or changing symmetry mid-run.

The sheared-cell first-step refusal was measured as2.855e-18 between fractional
A and reconstructed reciprocal B symmetry tables. Match their existing init
tolerance1e-10 in the same owner; keep B/global owned-input identity exact.

Nonbinary diagonal length7.3/grid12 also triggers the old admission bug. Compare
that case to old optimized kernels with only admission corrected; label that
library comparison-only and do not distribute it as another supported engine.
