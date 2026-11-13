# Validation for the HSE/SOC proposal

## Submitted regression

Only the existing small Si HSE system is used: eight Si atoms,12^3 real-space
points,4^3 k grid,32 electrons, MPI4/OMP2. Case420 produces a converged GS;
case421 reads that GS and propagates64 impulse-response steps. Both follow the
ordinary `preparation` / `verification` / `create_mpi_test()` convention.

The GS verifier requires normal termination, convergence below the input1e-8
tolerance, finite eigenvalues and restart data. The RT verifier checks normal
termination, expected samples, finite current/energy/response, electron count,
energy drift and native dielectric/conductivity output. It is not a converged
optical spectrum and is not a performance acceptance test.

```sh
python3 tools/check_forbidden_eigensolvers.py --source-root src &&
cmake -S . -B build -DUSE_HSE=ON -DUSE_MPI=ON -DCMAKE_BUILD_TYPE=Release &&
cmake --build build --parallel 4 &&
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 ctest --test-dir build -R '42[01]_bulk_Si_hse' --output-on-failure -j1
```

## Evidence and boundaries

The 2026-11-13 preparation deliberately excludes developmental giant-Gamma
routes. Its fresh verification is recorded in
[the reduced-tree scope report](hse-pr-scope-20261113.md#verification-on-the-reduced-tree).
No new GS data was collected; the Si trajectory results below describe the
earlier proposal and do not certify this reduced tree.

Before submission cleanup, GNU15 CPU serial/MPI builds, application-object audits,
40 focused checks, six automatic-layout starts and the64-step Si RT passed.
Those extra development probes have been archived, rather than included in the
public regression suite. The reduced-tree Release build and full fresh Si GS/RT passed for this preparation.

Historical Fugaku tests compare different immutable snapshots; they do not
certify this complete merged tree. No broad runtime speedup is claimed.
Fujitsu/A64FX and GPU execution of this candidate remain unverified. Source and
application-object guards remain active; the public test reduction does not
remove scientific acceptance checks from the implementation.

## Reduced-tree result (2026-11-06)

All six standard preparation/run/verification stages passed. GS converged in98
iterations with density residual9.0503234e-9. GS run262.39s, RT run295.97s;
these are single-run timings, not a speedup claim. [CTest record](validation/si-hse-20261106.txt)
normalizes only the build-directory pathname. The exported source configured
successfully with MPI/HSE enabled. Source and all linked application-object
subtrees passed the prohibition checker. No additional regression system is
included in the proposal.

## OpenACC compile guard follow-up

The merged upstream OpenACC fixes (including PR1292) are retained. A subsequent
source audit found a new CPU-only MGGA action with an unguarded rank-eight
`flux` declaration. It now follows upstream's `add_xc_tau_operator` pattern:
OpenACC compiles a CPU-only rejection branch, with no operation for an inactive
payload; the original declarations and action remain in the CPU branch.

GNU preprocessing of111 changed Fortran files found no rank-eight-or-higher
assumed/deferred-shape declarations in the OpenACC view. A focused check failed
before the guard and passed afterward; CPU preprocessed nonblank lines remain
identical. The OpenACC view of this module passed GNU syntax checking against
the existing interfaces, and the CPU Release rebuild and object audit passed.
These are local compile-structure checks, not NVHPC/Miyabi compilation or GPU
execution. No additional public regression system was introduced. The fresh Si
result above predates this guard-only change; its CPU action is unchanged.
