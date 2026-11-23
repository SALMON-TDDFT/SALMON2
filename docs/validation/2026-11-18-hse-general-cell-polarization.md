# General-cell fast-PR HSE and actual-field polarization

The previous rotated-cell reduced-k route failed ACE Hermiticity while the full-k
route completed. Diagnostic selection of reference exchange within the same fast-PR
engine restored agreement. Scalar FFT-orbit reuse used aligned Cartesian lengths
for gauges, k lookup/rotations and wrap phases. The complete replacement now uses
cell matrices/Cartesian rotations in general cells, following the existing spinor
formulas. Off-diagonal entries select the branch; positive diagonal cells retain
old arithmetic even when n*(length/n) differs from length at roundoff.

The first sheared-cell step was separately refused because immutable fractional A
and reciprocal B symmetry tables differed by2.855e-18. Their owner comparison now
uses the same1e-10 rule as map initialization; exact B/global owned-input identity
checks remain. No tolerance for current/energy agreement was increased.

The host API receives Cartesian A[3]. Attached-host, zero-impulse response starts
at A=0 without a z-direction assumption; actual vectors are checked before every
step. Ordinary batch initialization keeps its prior checks. General-cell symmetry
and non-z-preserving groups require explicit corrected native-core capabilities.
Old missing-capability libraries reject before acquiring a lifetime. No field
projection, primitive-axis reorder, automatic GS conversion or new HSE engine.

Final explicit local build: /private/tmp/salmonpy-general-cell-release-build/build.json
(prepared source /private/tmp/salmonpy-general-cell-release-source). OpenBLAS and
O2 match the preceding reference. Fixed binaries/pinned snapshots were not edited.
Full files and before/after hash manifests retain source provenance/licensing.

Verification (OMP/BLAS1; synthetic8³, full8/reduced4-k,32 steps):

- Rotated orthogonal, left-handed, xy-sheared and z-permuted cells; each compares
  full/reduced at internal k ranks1 and2.
- x/C2x, y/C2y and circular-xy/mirror-z fields; also internal k ranks1 and2.
- All14 comparisons pass. Maximum current difference1.734723475976807e-18;
  maximum energy difference1.3322676295501878e-15 (atomic units).
- Each family runs two independent MPI2 material groups with different screening.
  Within-group traces agree; different material groups have distinct traces.
- Own fresh-process7→32 direct resume matches continuous traces/final checkpoint
  bytes, with zero accepted-step replay. Incompatible field at step4 is rejected
  before update; saved state bytes stay unchanged and the next valid step succeeds.
- Original diagonal-cell trajectory bytes equal the previous library. For
  length7.3/grid12, the old library itself refuses admission; new results match
  old optimized kernels with only that admission check corrected. This auxiliary
  reference is comparison-only, not a supported or distributed alternative engine.
- Normal mixed FDTD/LD/SBE/reduced-HSE, PBE, ordinary Taylor/PT-CN and seven-rank
  protected mixed tests pass against previous references.
- 587 source tests pass, one optional skip. Package/build evidence includes the
  source-distribution-to-wheel/native hash check and isolated installed tests.
- No GitHub Actions launched. Molecular work remains paused. This is numerical
  connection/initialization validation, not converged real-material physics.

Evidence is retained outside the repository in developer outputs directory
salmonpy-hse-general-cell-20261118. Tools/test_fast_hse_general_cell_mpi.py records
all commands and per-case differences; its optional roundoff-reference-build is
explicitly labeled. Same binary/model/MPI decomposition remains required on resume;
checkpoint array/format layout is unchanged. Arbitrary space groups, symmetry
learned/PT-CN, unrestricted state handles and checkpoint repartitioning remain open.
