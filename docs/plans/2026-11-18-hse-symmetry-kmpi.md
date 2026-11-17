# Reduced-k internal MPI implementation plan

Extend only the explicit fast-PR scalar periodic/full-q/factor-off Taylor profile
from one to two internal k ranks. Existing exchange star/orbit distribution is
retained; remove admission restrictions only after direct numerical tests. No
conventional-HSE fallback or field projection. Default symmetry denial remains.

Add a separate native capability query for symmetry k-MPI; old binaries retain
one-rank admission. Validate nonempty contiguous representative ownership and
rank-local phase coverage in checkpoint support. Archives remain bound to the
same MPI rank count/ownership; repartitioning a checkpoint is not supported.

Acceptance: original phase0 full8/reduced4 k, one/two internal ranks; compare
current and energy over32 steps, OMP1/2. Direct save7/16 and fresh resume must
match continuous traces and checkpoint bytes without replay. Retry rejection
must preserve per-rank state. Check two independent groups and normal non-symmetry
regressions. Update Japanese/English/AI docs and checklist, package checks.
Local testing only; do not run GitHub Actions.
