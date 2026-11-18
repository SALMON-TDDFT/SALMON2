# Fast-PR input assets admission

Next checklist subset: pre-native structural validation of generated k meshes,
retained stars/weights, referenced PP assets and the supported legacy GS records.
No input normalization, native formula change, download or compile on import.
Existing collective create_state error agreement remains the boundary.

- Match native shifted uniform-mesh order and representative weights. Reject
  incompatible symmetry/shift meshes, explicit file_kw (native sets kmesh=-1),
  unsupported even gamma-centred endpoint grids, empty k ownership and overflows.
- Resolve exact case-preserving bundle-local pseudo filenames for every species.
  Check readable/nonempty/native-recognized format; plain scalar KY numeric
  header/table structure. Other format-content validation remains partial.
- Supported GS layout: native-endian4-byte Fortran record markers and integer/
  logical fields, float64 occupations, single global real64/complex128 wfn.bin.
  Verify metadata bounds, effective k/state/grid sizes, occupations/electron sum,
  file lengths/markers and finite occupation values using bounded streaming reads.
  Occupation range/electron-sum differences warn without clipping or renormalization;
  startup does not rescan the entire orbital file. Explicit deep diagnostics remain separate.
- Do not infer cell/XC/pseudopotential/k-order identity from legacy headers that
  lack those fields. This is structural admission, not physical-GS provenance.

Unit failures first; malformed one-rank assets in two independent MPI2 groups,
repair/retry, normal traces and own direct resume. Preserve normal non-symmetric
Taylor/PT-CN and mixed samples. Update JP/EN/AI/checklist, source and installed
suites, sdist-to-wheel. Local only; no Actions.

## 2026-11-18 warning-first review

Resource budgets are advice, not execution limits. Layer rank planning now warns
while retaining all points and reporting required ranks; actual MPI size must
still match. Structural impossibility/corruption and confirmed native STOP
conditions remain errors. Do not expand unsupported-profile guards merely because
a case is untested. Pseudopotential parser coverage remains conservative and needs
further comparison against native-accepted syntax; do not claim this audit proves
that all existing profile restrictions are minimal. SCF non-convergence currently
remains a failed-result contract, not an advisory precision estimate.
