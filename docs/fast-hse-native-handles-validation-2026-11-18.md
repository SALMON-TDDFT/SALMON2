# Fast-PR native C handles — 2026-11-18

The new protocol1 C entry uses generation-checked opaque int64 tokens for multiple
independent states of one fixed model in one attached MPI group. Existing kernels
and legacy singleton functions remain. The adapter selects only the state pointer,
accepted A and failure flag during one synchronous main-thread call, then returns
those values to the owning slot. It never copies whole owning contexts. Registry
resizing preserves target ownership; a retired token does not select its successor.

## Scientific and ownership checks

Six cases: learned Taylor and learned PT-CN, each internal MPI1/OpenMP1,
MPI2/OpenMP1 and MPI2/OpenMP2. Each case runs five separate jobs: isolated A,
isolated B, mixed handles, preceding legacy singleton A and B. Every job has two
independent model groups with different HSE omega. Total30 MPI invocations.

A uses dt=.001/capacity32 and z-polarization; B uses dt=.002/capacity34 and rotating
xy polarization. Current/total-energy trajectories, eight local history diagnostics
and the entire serialized HSE numeric payload agree exactly with independent
and old-singleton controls. Mixed jobs additionally exercise:

- simultaneous per-handle protection/rollback and staged restore commit/cancel;
- a recoverable deep numeric overflow without changing accepted state/output;
- close A with B's restore candidate alive, then continue B;
- reusing A's slot without revalidating A's old token;
- direct loading step9 into fresh independent states, then23/24 new calls through
  A32/B33, with every future current/energy equal to the continuous control;
- expanding to eight simultaneously live slots while a transaction and another
  restore candidate remain live, then independent cleanup;
- invalid0/negative/large tokens, different valid tokens across MPI ranks, and
  mismatched positive initial dt/capacity followed by corrected initialization.

All handles close before host detach. Host MPI remains initialized after detach.
The C example compiles with the matching host MPI and runs without Python.

Reproduce a case explicitly (two material groups):

```sh
mpiexec -n 4 python tools/test_fast_hse_handles_mpi.py \
  BUILD FULL_K_FIXTURE NEW_OUTPUT --group-size 2 --mode mixed
```

Use separate fresh output directories and `--mode single-a`/`single-b` controls.
`--legacy` is an explicit isolated reference path only; never a fallback for mixed
handles. Fixture inputs are synthetic connection tests, not converged material GS.

## Archive diagnosis

An initial comparison of entire independently generated checkpoint files failed,
although propagated observables/history agreed. A read-only decoder located all
examined differences in existing inactive periodic fixed-ion workspace:
`tpsi%update_zwf_overlap`, unwritten `rt%Qp_e(:,:,1:)`, projection-only `rt%E_old`,
isolated multipole `poisson%npole_total`/`n_multipole_xyz`, and disabled +U `energy%E_U`.
Restore deliberately invalidates overlap caches. The examined propagated orbital,
density/potential and HSE payload regions matched. No unknown byte offsets are
masked and no checkpoint format is changed. Verification compares the complete
HSE payload and actual future propagation from restored state instead of imposing
whole-file canonicality on old workspace. Raw archives and decoding evidence stay
in private local outputs; they are not published fixtures.

## Remaining scope

Python NativeFastHSEBackend retains one native lifetime. Arbitrary heterogeneous
models, multiple shared symmetry initializations, threads, shared I/O/timer
extraction and per-context MPI subgroup retirement are not completed. Currently
subgroups are attachment-owned until detach. Stable token identity does not imply
portable/canonical archives or unlimited create/close cycling. Local validation
is not a scaling, cross-platform or physical convergence certificate. Actions
are not used for this campaign.
