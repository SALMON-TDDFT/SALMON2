# Taylor+ACE BLAS OpenMP Implementation Plan

Goal: support internally threaded vendor BLAS, including 12 threads.
Architecture: all HSE BLAS calls outside application OpenMP regions; retain elementwise OpenMP. Honor vendor environment instead of forcing OpenBLAS to one thread. No additional orbital replicas.

User approved this approach. Execute inline in TDCDFT.

1. Extend production-size ACE oracle to OMP12.
2. Remove forced BLAS thread setter and outer GEMM parallelism.
3. Build and run Taylor+ACE MPI1 x OMP1/2/4/12, two repeats, orbital parity and RSS.
4. Document measured local results; Fugaku performance remains unmeasured.

Completed: native build, 4 native tests, real-MPI test with OMP12, and eight
Taylor+ACE benchmark runs passed. Maximum final-orbital relative error 4.18e-12.
MPI1 OMP1/2/4/12 mean seconds per step: 5.0462/3.6481/3.1933/4.6873.
Fugaku build/performance remains untested; no machine-specific speed claim.
