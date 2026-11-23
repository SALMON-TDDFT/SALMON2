# Reduced-k fast-PR HSE internal MPI2 (2026-11-18)

Scope: scalar periodic CPU/full-q, factor-off Taylor/ACE, explicit symmetry,
nonempty k-only groups of one/two ranks, one native lifetime per process. Fast-PR
only. Full-action remains test-only. Synthetic8³/full8-vs-reduced4 C2z input with
the original phase0; not a physical-convergence certificate. macOS Apple Silicon,
GNU Fortran15, OpenMPI5, Libxc5.2.3; local execution only, no Actions.

The previous small-gradient fix was committed/pushed as e6020b90 before this work.
Git/ref-log writes and Libxc reads briefly waited for OneDrive; synchronization
restored access. Failed initialization/comparison logs are retained separately.

## Implementation

- New `sp_fast_hse_symmetry_kmpi_protocol()==1` is required for group2 in addition
  to the existing symmetry protocol. Missing/unknown capability rejects before
  native lifetime acquisition. Default init continues to reject symmetry.
- Existing fast-PR star/orbit distribution is retained. Each rank owns contiguous
  representatives and the corresponding full-star phase span.
- Checkpoint maps validate rank-local cache first/count/offsets. Read allocation
  sizes derive from the independently initialized expected map; cached ownership
  and kernel phase span must match. Existing archive headers pin MPI rank and
  partition. No archive repartitioning is supported.
- MAX rather than SUM consensus preserves state-preserving status34 over two
  ranks. Rank-local rejection/backup paths remain collective.
- Ordinary mixed sample accepts --hse-ranks2 with explicit symmetry/reference.

## Checks

1. A new unit test fails before capability support and passes afterward. Old
   libraries refuse group2 before native initialization (deny-old-k2).
2. Original phase0 local FDTD A/c, dt=.003,32 steps, two independent groups of
   two internal ranks: full8/reduced4 max current difference4.34e-19. Energy
   difference4.44e-16 with OMP1 and8.89e-16 with OMP2. Existing tolerances remain
   rtol2e-10/atol2e-12 for current; no field projection or fixture perturbation.
3. Two independent two-rank groups with omega=.11/.12: direct save7→32 fresh
   resume reproduces current/energy traces and complete final checkpoint bytes
   exactly. No accepted-step replay. Forbidden fields leave before/after state
   bytes and Python metadata unchanged before retry.
4. Four-rank FDTD/LD/SBE/HSE mixed continuous32 and fresh16→32 resume pass.
   Four-rank failure/retry and five-rank protected restore preserve original
   field/LD/electronic states; three future retry steps match exactly and replay0.
5. Existing PBE, fast-PR Taylor, PT-CN and combined non-symmetric five/seven-rank
   regressions pass against former binaries with unchanged tolerances.
6. Source Python suite553 tests, one skip. Final installed suite and
   source-distribution-to-wheel/native-hash verification are retained with the
   external evidence. No other platform was executed.

## Reproduction

Prepare NEW_SOURCE with tools/prepare_hse_rt_state.py --explicit-routing; build a
new matching host-MPI library explicitly. Preserve the former library for the
negative test. Both groups must use the same new build and equivalent full/reduced
initial states; all ranks in a group participate in every collective call.

```sh
mpiexec -n 4 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD FULL_BUNDLE REDUCED_BUNDLE NEW_COMPARE --group-size 2 \
  --drive ORIGINAL_LOCAL_FIELD_NPZ --point-index 1
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD REDUCED_BUNDLE NEW_SAVE --mode save --group-size 2 --checkpoint CHECKPOINT
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD REDUCED_BUNDLE NEW_RESUME --mode resume --group-size 2 --checkpoint CHECKPOINT
```

Maintain the same MPI decomposition on resume. Initial operator construction is
still required; temporal history is loaded rather than replayed. Learned/PT-CN
symmetry, larger/empty/real-space/orbital partitions, broader initialization and
converged real-material benchmarks remain outside this verified scope.

Local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-hse-symmetry-kmpi-20261118/`.
