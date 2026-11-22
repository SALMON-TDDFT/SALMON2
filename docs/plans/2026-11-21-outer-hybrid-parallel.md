# Persistent outer q OpenMP team

Goal: retain MPI k distribution and single-thread worker BLAS, but move dense
screening's OMP parallel region outside all bounded q batches for one target k.
Avoid repeated fork/join and per-worker allocation; share original-H refinement
parts once per target. Keep bounded source exchange buffers and physics/checks.

SALMON initializes MPI_THREAD_FUNNELED. MPI calls inside the new parallel region
are confined to omp master (the original MPI-initializing thread), separated
from worker activity by barriers. No single region may perform MPI, because an
arbitrary worker could be chosen. Dense H preparation and optical/grid callbacks
stay outside this region: their x/hx/target_info are shared and not reentrant.
Matrix-free/SOC callback paths remain unchanged until thread-private adapters are
available. Existing Mac sequential-OpenBLAS vendor guard is retained.

1. Add a persistent-worker test of static/pump q responses across several batches,
   changed G counts and occupations, versus the current C-ABI call and independent
   oracle. Observe unavailable worker API before implementation.
2. Extract the unchanged one-q solver into a module worker with explicit local
   context/scratch lifetimes. Keep C ABI2 as a wrapper. Parent BLAS policy and
   original-H slices are target-local, worker BLAS state/scratch are thread-local.
3. Add driver dense helper: prepare context outside OMP; one parallel team per
   target k; primary-thread MPI/packing/collective rejection; omp do q vertex and
   response work; primary-thread chi accumulation and bounded buffer recycling.
   All threads cross the same barriers on success, empty-k ranks and rejection.
4. Audited native tests and numeric/refinement regressions (including4MPI12OMP),
   uneven MPI distributions and collective invalid-input rejection. Full audited
   SALMON build and actual Si comparisons, with frozen validated references.
5. Confirm region-count reduction and thread teams; compare whole W/phase timings
   at matched resources. No claimed acceleration without measurement. Save only
   derived diagnostics/records in canonical calc2026, no GS payload collection.

## Implementation and review outcomes

- Kept one dense q team per target and thread-private reusable scratch. Added
  primary-thread MPI capability wrapper in real and dummy communication modules.
- Initial multi-MPI test caught a termination-flag race: a next-iteration master
  update could overtake the previous bottom flag read. The final loop reads done
  only after the top master/barrier handshake. The omp-do end barrier protects
  old data; nonmaster threads touch no shared state before the next top barrier.
  Redundant initial/bottom barriers were subsequently removed and tests repeated.
- Review caught unallocated optional ASSOCIATE selectors. Associations now exist
  only inside their allocation branches; recursive worker/helper procedures make
  local lifetimes explicit. The initial binaries/sources/stack evidence are
  preserved as unapproved test evidence, never executed or reused again.
- Checked4MPI×12OMP and uneven5MPI×2OMP pass. Collective rejection covers bad
  vertices detected inside the team, nonfinite inputs, missing FUNNELED support
  and only rank0 requesting above its OMP thread limit.9MPI with8k points and
  zero box-pair ownership also matches the old complete W. Dummy communication
  compiled/probed with the project's legacy compatibility flags.
- Independent90-digit original-H refinement still passes4MPI×12OMP; no numerical
  thresholds or physical parameters changed. C ABI2 preserved.
- Final same-backend native1MPI×12OMP interleaved3-trial medians: static W
  0.0668→0.0583s (1.146×); integer-occupation pump/refinement W1.9727→1.9317s
  (1.021×). Regions64 versus192batches. These are small-model local measurements,
  not production-size Si or Fugaku scaling results.
