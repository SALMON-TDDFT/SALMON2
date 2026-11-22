# Scalar Sternheimer screening and profiling implementation plan

**Goal:** Add an opt-in matrix-free static scalar-insulator W solver and detailed timing without changing the running baseline or default numerical path.

**Architecture:** Keep `hsex_w_solver='auto'` on the existing dense scalar and matrix-free SOC paths. `sternheimer` enables the real-space callback for scalar insulating GS screening only; fractional scalar updates stay on the existing checked path. Normalize scalar occupations to one spin for the existing Sternheimer band solver and restore the factor two on its response. Both routes retain the optical-limit head/wings, local fields, q mesh, finite checks and acceptance thresholds.

**Tech stack:** Fortran, existing SALMON hpsi/MPI wrappers, LAPACK/BLAS, Python offline regression harness only.

1. Add a synthetic scalar W regression with known dense/full-spectrum reference. Cover q=k and q<k, optical head, replicated/distributed layouts and invalid input rejection. Demonstrate rejection before extending the driver.
2. Extend the matrix-free driver validation and spin normalization; retain positive-gap/eigenpair requirements and reject fractional scalar populations. Add opt-in input, broadcast, variables-log and manual entries. Share the existing real-space Hamiltonian callback between scalar and SOC paths.
3. Add monotonic per-rank phase accumulation and rank-maximum timing logs: H construction, reduction/factorization, optical response, source exchange, vertices, solves, and matrix redistribution/checks. Do not change existing three aggregate timers.
4. Run source/native-source and all-application-object prohibition audits before every new build/link. Run synthetic tests, scalar IEEE checks, existing SOC regressions, a GNU full build, and actual Si comparisons using fresh scratch. Preserve the running 0cda655f executable and all current jobs.
5. Compare native scalar W, optical dielectric tensor, and identical RT currents against the dense reference at unchanged tolerances. Record phase timings and matched-resource wall time; do not claim speedup from differently loaded machines. If scalar hpsi/eigenpair checks reject the physical fixture, diagnose rather than lower acceptance criteria.

Deferred until evidence supports them: distributed box-matrix redistribution/reciprocity exchanges, phase/previous-solution caches, block/preconditioned solves and repeated W-update scheduling. These are separate changes with their own numerical and capacity checks.

Results and records go under canonical calc2026; no GS payload collection. Production W control/solve/I/O remains Fortran. Source defaults and current jobs stay unchanged.

## Measured outcome and scope adjustment

Synthetic scalar W (including arbitrary unused-empty-state rotations), q=k
and q<k, 1/uneven5MPI and invalid-population/NaN/Infinity tests passed.
The scalar implementation projects occupied states only, as the dense
static reference does. Requiring unused high empty eigenpairs caused an
initial physical-Si rejection; failure evidence is preserved in scratch.

Actual Nr8^3, Nk=Nq=2^3, boxes2^3, cutoff1Ha, 1MPI1OMP: the full matrix-free
route agrees with dense W to 1.23e-13 but does not accelerate the whole small
case. Optical work improves ~4.9x while individual finite-q CG solves slow
down. Therefore add opt-in `hybrid`: dense finite-q BLAS + real-space optical
head/wings. Auto remains unchanged. Interleaved three-run TB-mBJ median W
time improves16.0854->9.2133s (1.746x), response3.198x. Hybrid W difference
2.92e-16;16-step RT-current difference7.30e-19au. LDA single-pair W time
16.3432->9.3965s(1.739x), W error4.86e-16. Machine was also running two
baseline jobs: these are matched-resource, contended-host measurements.

Default dense and fractional pump-reference regressions and SOC regressions
pass. An existing CMake directory-scope bug left sibling test link-checker
paths empty. Promote the fixed mandatory checker path to INTERNAL cache;
previously failing guarded probes now build with actual PZHEEV checks, and
a generated-guard regression CTest passes. No prohibition checks disabled.

Nr12^3,Nk=Nq=4^3,boxes6^3,cutoff3Ha matched8MPI1OMP dense/hybrid follow-up
is waiting on the existing spectrum pipeline COMPLETE marker, avoiding
extra concurrent production jobs. Its executable is a frozen audited copy.
Fujitsu build, large-box capacity and production-scale speedup remain unverified.
