# Same-model reduced-k fast-PR HSE sharing

## Implementation and boundary

NativeFastHSEModel and normal ResponseSession now accept explicitly opted-in
reduced-k symmetry states of one immutable material/cell/k mesh/retained group.
Python requires shared_symmetry_protocol1 before attaching or consuming the native
lifetime. First handle initialization binds the model symmetry option; later
allocations inherit it. Existing settings and exact symmetry-table checks remain.
The actual A/c of every step must be invariant under the retained operations.

The first real two-context run exposed repeated allocation of density symmetry
maps. A complete-file HSE override now releases/rebuilds construction-only maps
synchronously. The arithmetic is unchanged. Propagated density/wavefunctions,
exchange maps/kernels/history, output, transactions and staged restores remain
context-owned. Child communicators remain attachment-owned until final detach;
no arbitrary subgroup retirement, concurrent entry or heterogeneous-model claim.
Original bundled sources and provenance are retained; the override and manifest
hashes identify the new preparation. Fixed binaries were not modified.

## Local verification

Mac64GiB/18 cores, macOS26 Apple Silicon, Python3.14, Open MPI; OMP/BLAS1. Explicit
fresh preparation/build from verified bundled sources; mandatory source/object
PZHEEV prohibition audits pass. No GitHub Actions dispatched.

- Numerical matrix: Taylor strict, Taylor learned, PT-CN off, PT-CN learned ×
  internal MPI1/2 × full/reduced-k. Three solver groups per run: one group shares
  two independently driven states; the other two use the legacy independent
  initializer.32 accepted steps, different nonzero z drives, direct saves at8/9.
  Shared/independent results agree exactly. Both fresh restores continue exactly
  without replay and do not alter the sibling. Symmetry-breaking x drive rejects
  without changing either original. After closing all states, a fresh model state
  still starts from original GS; a surviving lease works after model.close.
- Full/reduced-k comparison keeps the existing tolerances: current rtol2e-10,
  atol2e-12; energy rtol2e-12, atol2e-11. Maximum measured energy difference across
  the matrix is2.220446049250313e-15 Hartree. No thresholds relaxed.
- Ordinary FDTD/LD/SBE/two-HSE mixed examples, learned Taylor and learned PT-CN,
  128 steps saved at63: separate6 ranks vs shared4 ranks have exactly equal
  A/c, E and j. Both HSE energy streams also match continuous execution after
  direct restart;65 new steps execute,0 replay. Each point also matches the
  independent full-k serial reference within existing tolerances.
- Late restore candidate failure on the last internal rank preserves original
  field/LD, SBE, both electronic binary states, Python clocks/arrays/flags/object
  identities and output streams. Retry in the same session succeeds and matches
  continuous propagation. PT-CN with two field ranks (5 total) matches one field
  rank exactly.
- Four old full-k/PBE/Taylor/PT-CN/combined protected regressions,32 steps, match
  previous baselines with the new runtime. Existing default behavior retained.
- Source and isolated installed-wheel suites:603 tests each, one external-condition
  skip each. Source distribution→wheel/native hashes pass. Installed wheel can
  prepare verified bundled sources with the new override. Its normal reduced-k
  mixed save/resume/late-failure retry, standalone Python example, and installed
  MPI2 shared/independent/idle-recreation/direct-resume probe all pass.
- Paired guides and AI JSON synchronized. Read-only review found no blocker.

## Important negative cases

The ordinary localized 3D Taylor mixture requested256 steps. After145 steps its
local field ceased to preserve c2z and the next protected step rejected it with
original_preserved=True. This is an actual-field restriction, not a sharing
failure; it also occurs with separate groups. The non-symmetric PBE combined
fixture emits incompatible transverse components much earlier, also safely
rejected. A previously used alternate PBE fixture did not cure this restriction.
Do not project away components, loosen relative checks or silently shorten user
runs. The128-step LD/SBE/HSE campaign is explicitly a limited compatible window;
use a genuinely compatible field/geometry/group or prepared full-k input beyond
it. Full-k combined PBE/HSE remains covered by the positive legacy regression.

The original shared-symmetry Python veto, density-map repeated allocation, initial
wheel packaging path typo and installed probe missing its explicit symmetry flag
are retained in developer logs alongside successful corrected runs. These are
synthetic ownership/interface regressions, not converged material benchmarks or
cross-platform certification. Evidence retained outside Git under developer
outputs salmonpy-shared-symmetry-20261118; raw checkpoints are not portable or
canonical archives.
