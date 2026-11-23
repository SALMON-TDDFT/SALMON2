# Internal fast-PR symmetry TDDFT connection

Goal: Connect validated core symmetry/checkpoint support to complete TDDFT state
through an explicit internal rt_open(..., allow_symmetry=.true.) opt-in. Existing
Python worker/native-C/FDTD entry calls omit it and remain rejected. Single
internal k rank, scalar periodic CPU/full-q, factor-off Taylor/ACE only.

Architecture: retain existing initialization and propagation arithmetic. Bound
symmetry profile before initialization; owner map checks actual A before any
backup/propagation. A broken field little group returns a distinct local status34
without mutating anything; other admission failures preserve existing semantics.
Use existing staged full-context + exchange V2 checkpoint path, with fresh-owned
operator preparation and no accepted-step replay. Symmetry configuration remains
process-common/immutable; no heterogeneous same-process contexts are promised.

1. Add synthetic full8/reduced4-k fixtures (8^3 grid, symmetry-preserving orbitals,
   origin atom, C2z or mirror-x, z drive). No SCF/physical-convergence claim.
2. Fail-first explicit native rt_open probe; continuous/save/resume at both odd
   and even splits with exact API/internal traces, staged rejection and forbidden
   A update preservation/retry; compare reduced/full fast-PR traces within tolerance.
3. Implement opt-in guard and distinct field status. Preserve public callers'
   default false, bound unsupported factor/PT-CN/MPI2 before initialization.
4. New external adapter library/probe; run full/ACE, C2z/mirror-x, OMP1/2,
   independent MPI1/2 processes with COMM_SELF; negative and non-reduced regression.
5. Update JP/EN/AI/checklist, source/installed suites, sdist-wheel/native hashes;
   persist evidence, skip-ci commit/push. Next is native Python opt-in + protected
   mixed integration and reduced internal-k MPI; do not claim these now.

Symmetry archive copies canonicalize overlap-cache dirty flags to true, because
overlap caches are not serialized. Live flags and non-reduced V1 bytes stay unchanged.
