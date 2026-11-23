# Existing GS → in-memory standard TDDFT validation

Apple Silicon/macOS, GNU/OpenMPI, OpenBLAS/Libxc. Baseline d824b4c6.

The new model owns input/orbitals/occupations, uses the old RT restart-read boundary,
and follows the same remaining initialization/update operations. Existing GS import
runs initialization only: no SCF and no new GS restart bundle. Imported convergence
is unknown; uncomputed forces are rejected. GS/RT timers initialize once per host
attachment, and input/cache views are parked during interleaved GS loads.

Five cases passed exactly: PZ MPI1, orbital-MPI2 and spatial-MPI2; Libxc-PBE
orbital-MPI2; TB-mBJ spatial-MPI2. Each compares 16 driven current steps to the
existing-file path on the same runtime. Oppositely driven siblings do not change
the primary trace. A fresh independent state matches the trace; direct save8→16
restore matches without replay. GS closes before later state creation. An existing
GS import between RT steps preserves the RT owner. This does not validate new
SCF execution interleaved with RT, nor establish physical material convergence.

Analytic real/complex arrays on MPI2 verify empty-band orbital redistribution
(6 source bands split3+3 → 4 occupied target bands split2+2), without SCF/files.
Source637 and installed-wheel637 tests pass, one skip each. The installed MPI2
`salmonpy.examples.memory_gs` sample passes. Native source/all application-object
prohibition checks and wheel native hash verification pass. No Actions dispatched.

Reproduce against an existing matching input/GS bundle:

```sh
mpiexec -n 4 python tools/test_memory_gs_mpi.py BUILD EXISTING_GS NEW_OUTPUT --size 2 --split spatial
mpiexec -n 2 python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_OUTPUT
```

No new SCF/reference GS data was generated. A proposed fresh-reference operation
was rejected by automatic approval review under the user's no-new-GS instruction;
validation instead used existing local fixtures and synthetic array inputs.
Selected diagnostic logs/hashes, runtime and wheel are stored in canonical
`calc2026/salmonpy-memory-gs-tddft-20261118`; no GS restart/orbital/density bundles
were collected. Local OneDrive files verified, cloud synchronization unverified.

Remaining: different subgroup/repartitioning, retained-GS checkpoint, HSE/SOC/+U/DC
and threaded/concurrent ownership; Linux/Intel Mac runtime verification.
