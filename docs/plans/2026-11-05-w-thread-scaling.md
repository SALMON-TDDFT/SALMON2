# W MPI / numerical-library thread scaling

Approved scope: compare 8 MPI x 1, 4 MPI x 2, 2 MPI x 4 under the same nominal eight-core budget; verify numerical consistency and select useful settings. Preserve physics and remote production jobs.

1. Identify actual NumPy/SciPy BLAS/LAPACK providers. Current wheels use Apple Accelerate, so do not label native library threading as proven OpenMP or infer utilization from OMP_NUM_THREADS alone.
2. Benchmark fresh processes on the existing 6^3 synthetic fixture and a 12^3-grid synthetic operator fixture. Use repeated interleaved runs, correct provider-specific environment before imports, no concurrent workloads, and record CPU/wall utilization and numerical file differences. Keep all failed artifacts.
3. Use actual provider query/control APIs where available to verify the selected threading model; distinguish requested upper bound from observed cores. Assess large tridiagonal reductions and RHS products separately if needed.
4. Document best observed configuration and reproducible launch settings, with limitations for Si production/Fugaku and non-matrix-free scaling. Store reports under canonical calc2026; no GS collection or executable replacement.
