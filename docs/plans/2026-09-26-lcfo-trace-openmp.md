# LCFO trace and OpenMP optimization

User approved eliminating duplicated diagnostics and parallelizing masks.
Use the exact retained-ACE half trace -dv/2 sum_n f_n |sum_r F_r^H C_r|².
Reduce the overlap, never reconstruct global ACE action for this diagnostic.
Keep diagnostics current on every refresh (no stale energy or schedule change).
Reorder spherical masking by WF column and add OpenMP over independent WFs.
Compute loss per WF, then sum in deterministic order outside parallel regions.
Validate complex unequal MPI partitions/nonunit dv/fractional occupation against
full ACE action. Test spherical masks with OMP1/2/4. Native regression and
Si128 R9/N4 before-after (same spherical code), one numerical job at a time.
Keep BLAS1; do not assume an oversubscribed MPI16/OMP2 setup is faster.
Resume the preserved long process after validation/benchmark.
