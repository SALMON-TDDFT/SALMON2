# HSEX review fixes and CPU scaling

Goal: fix the three reviewed defects, then verify numerical equivalence and measured CPU speedup without altering remote production runs.

Approved design: retain the original translation-space exchange kernel for phase derivatives; select the insulating head only after checking occupations across k and a positive gap; separately normalize macroscopic charge-channel convergence. Reuse existing checkout and occupied-only RT.

1. Add regressions for even-grid finite-A derivative, integer electron/hole pockets, and charge-channel errors masked by local fields. Observe failures before fixing.
2. Correct each path and build a corrected serial-loop baseline in a new scratch build, with mandatory source/object audits and SHA recording.
3. Parallelize independent local k loops (MPI stays outside OMP); use ZGEMM for dense products, per-thread scratch and serial BLAS under OMP. Preserve scalar IEEE checks.
4. Build optimized executable separately, audit and record hash. Run unit tests, serial/MPI/OMP checks, finite-difference derivative and same-input RT timings. Compare corrected baseline, optimized OMP1 and optimized OMP4; multiple repeats, report full RT and targeted costs honestly.
5. Update manual and record results under canonical calc2026; no GS collection, no changes to queued jobs, no unrelated commits.
