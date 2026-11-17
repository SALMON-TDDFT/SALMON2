# Fast-PR symmetry checkpoint implementation plan

**Goal:** Extend internal same-binary numerical checkpoints to the already admitted
scalar CPU/full-q/internal-k1 fast-PR symmetry exchange core, preserving all current
non-reduced archive bytes and public symmetric TDDFT/FDTD rejection.

**Architecture:** V1 remains unchanged for non-reduced owners. V2 carries the
value-only symmetry map and cache alongside existing kernel identity/numeric history.
Fresh independently initialized owners supply trusted static layout; read exact
static scalar/layout/array identity before loading bounded cache arrays. Stage
numeric/map values before committing; do not serialize, copy, or replace FFT/C/MPI
handles. Capture/install carry map/cache values for rollback without replay.

**Tech stack:** Fortran same-binary stream IO; existing checkpoint adapter/native
fast-PR probes; Python source and installed-wheel tests; local MPI only.

1. Add native fail-first probe using real independent C2z/mirror-x core owners.
   Save action/history/cache, perturb/advance, validate-only and staged read,
   install then compare exact current/future actions and energy. Keep other owner
   unchanged. Reject different map/partition, truncated archive and missing tail
   without numerical or map mutation. Exercise absent/present phase caches.
2. Extend numeric payload with allocatable value-only map; V2 read/write guarded
   by initialized static map/full phase coverage; compare all static fields and
   array layouts against fresh map; bound cached dimensions by trusted metadata.
   Preserve V1 bytes and public symmetry guards. Add malformed cache tests.
3. Fresh external library: recompile checkpoint/context/fast-MPI against matching
   modules; run new probe full/ACE × star/orbits × OMP1/2 × MPI1/2 independent
   COMM_SELF processes, with distinct parameters/schedules. Existing exchange,
   mixed failure/retry/direct restore and old-vs-new V1 archive bytes unchanged.
4. Update paired guides/AI JSON/checklist with exact internal scope. Run source
   and installed545 tests, sdist-to-wheel/native hashes and installed MPI smoke.
   Persist evidence outside tmp; commit/push with skip-ci, no Actions.

Full symmetric TDDFT/FDTD, learned/PT-CN symmetry and reduced internal-k MPI2
remain separate work. Initialization may construct operators but never replays
accepted propagation. Raw streams require caller-authenticated archive/config.
