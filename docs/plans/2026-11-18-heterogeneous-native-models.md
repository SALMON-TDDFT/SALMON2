# Heterogeneous standard TDDFT model registry and activation

Approved next stage. Keep legacy single-model ABI/defaults; add explicit registry
ABI and Python registry owning one MPI attachment. Separate negative generation
model tokens from positive state tokens. Register owns parser/workdir, dt/capacity/
initial A, post-initialization input view and optimization/XC scratch caches; spawn
selects the requested model. All handle step/info/checkpoint/close paths select
that state's owner. No copies of wavefunctions or Libxc/MPI owners during switching.

1. Red public entry test and real two-model/independent MPI test.
2. Native registry, model/state lifetimes and explicit activation. Park/restore
   optimization modulo tables and XC scratch via move_alloc; reset PP format flags
   per species. Keep process timers/MPI attachment common and synchronous.
3. Python open_registry/open_model/new_state connection, one work copy/model,
   explicit close with live child leases, per-model dt/capacity/decomposition.
4. Compare interleaved PBE/TB-mBJ/PZ/Libxc-PZ/fixed TB-mBJ with independent old
   controls, MPI1/orbital2/spatial2 and mixed splits, opposite order, direct resume,
   stale tokens, sibling close and model reuse. Reuse existing GS only.
5. Paired documentation/AI/checklist, source/wheel and old-default regression.
   All new sources/objects audited before use; no Actions; canonical calc2026
   selected results only, no GS collection or portable-checkpoint claims.

Supported profile stays full-k periodic unpolarized fixed-ion middlepoint standard
semilocal TDDFT. No HSE engine changes, no thread concurrency, no general parser/
allocation/I/O/MPI crash recovery. Native registration still uses synchronous
compatibility views; model-specific derived numerical caches are explicit owners.

## Completed verification

Eight real heterogeneous MPI/OMP cases pass, including16³/8³ cells/atom/k counts,
five XC choices, dt0.003/0.004, opposite registration order and mixed internal
splits. Interleaved charge currents and direct restore match independent previous
runtime controls exactly; replay count0. Additional source/installed capacity
cases reach B56 with A capacity48. Fifteen legacy API comparisons retain exact
RT and energy columns. Source and installed wheel621 tests each, one skip.
Installed two-rank example and six-rank independent comparison pass.

Review found no remaining blockers after sorted cleanup with attempt-all/retry,
collective queued-model selection and coordinate-block profile parsing fixes.
Prototype link/profile failures are retained as logs. Native source/object audits
and wheel hashes pass. No Actions, new GS calculations or GS collection. Selected
evidence: canonical calc2026/salmonpy-heterogeneous-models-20261118. Local writes
verified; cloud synchronization unverified. GNU15 Apple Silicon only, Fujitsu and
other platforms remain unmeasured. Automatic ordinary-layout heterogeneous group
coalescing, full GS contexts and specialized SOC/+U/DC resources remain separate.
