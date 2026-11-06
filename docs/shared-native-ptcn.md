# Shared native PT-CN

`src/rt/ptcn_core.f90` supplies the functional-independent callback solver.
`ptcn_native.f90` supplies native density, Hartree, XC and Hamiltonian callbacks.
`exx_ptcn_core` and `exx_ptcn` remain compatibility interfaces.

Build with `USE_HSE=ON` for `propagator='hse_ptcn'`; use `USE_SCALAPACK=ON`
for the distributed matrix/occupation-block backend and `ptcn_mgga23`.
FFTW is required by either configuration and is discovered by CMake.
Both propagators require `yn_predictor_corrector='n'` in `&propagation` and
`gram_schmidt_interval=0` in `&tgrid`. Existing functional-specific input guards
remain authoritative; sharing the core does not certify new HSE SOC, thermal,
Wannier or laser configurations.

The native adapter uses FFT finite-difference kinetic preconditioning, coupled
MPI-reduced Anderson coefficients, and equal-occupation-block projection.
The non-Anderson compatibility solver retains its previous tighter inner stop.
HSE freezes only exact exchange during inner iterations and refreshes the full
physical exchange before final acceptance. The optional learned-inner mode is described below.
No full-column Lowdin rotation is applied to fractionally occupied states.

Final gates: normalized fresh residual <1e-10, Gram error <=1e-8, electron-number
error <=1e-7, and finite values. This implementation requires a full local spatial
domain. The submitted regression is the small Si GS/RT system; see [validation](PR_VALIDATION.md).
The 2026-10-28 local validation is saved under canonical calc2026/ptcn-shared-hse-20261028;
short tests do not establish long-time accuracy or Fujitsu portability.


## Learned ACE inside CN, exact endpoint acceptance

For the certified scalar integer-occupation CPU dense-k native HSE route:

```fortran
&functional
 xc='hse06'
 exx_factor_mode='learned'
 exx_factor_exact_interval=1
 exx_factor_n_teach=2
/
&propagation
 propagator='hse_ptcn'
 yn_predictor_corrector='n'
/
&tgrid
 gram_schmidt_interval=0
/
```

The first CN right-hand-side action is exact. Once enough accepted teachers are
available, a three-frame learned factor history predicts the endpoint ACE for
the first inner solve only. A fresh physical exchange action checks convergence;
failure retains the newly rebuilt exact ACE for correction, without reinstalling
the prediction. Prediction installation explicitly invalidates the physical-source
cache, including unchanged-wavefunction cases. Only after residual, Gram and
electron-number gates pass are exact endpoint factors aligned and committed.

CN requires `exx_factor_exact_interval=1`: every step has an exact endpoint
check and teacher. This does not use Taylor's skip-exact-step schedule. With two
teachers, prediction starts at step five if the learned system is ready; otherwise
ordinary exact ACE is used. Constant dt is required. RT restart is rejected until
history persistence is separately certified. SOC, thermal occupations, Gamma-only,
Wannier/localized support, spatial/orbital distribution and GPU learned CN remain
rejected. Existing Taylor factor-history behavior is unchanged.

`HSE_PT_CN_EXACT` reports actual native full exchange/ACE rebuilds after source-cache
checks for each entire step; `HSE_PT_CN ... builds=` counts fresh callback actions,
which can include exact-source cache hits. Do not equate the two. The benchmark
script alternates ordinary/learned runs, includes warmup, and reports measured
external and RT-step times. A short weak-impulse benchmark cannot certify long
laser trajectories or performance on Fugaku.


## Restored HSE endpoint Lowdin (2026-10-28)

The original Fugaku orbital route performs Lowdin retraction after the nonlinear
CN solve. That call was missing from the first shared adapter and is now restored
for the certified scalar integer-occupation HSE path. The MGGA partial-occupation
path is unchanged. Serial and orbital-distributed Lowdin use the existing helpers.
The raw Gram/electron acceptance guard is retained, then Lowdin is applied and
Gram/electron guards are checked again. The exact H/ACE is rebuilt at the actual
orthonormal output before collecting an exact teacher.

`HSE_PT_CN_ORTH` separately reports pre/post Gram, pre/post residual of the SAME
CN equation, retraction norm and Lowdin time. Raw nonlinear-root residual <1e-10
remains the original solver gate. Post-retraction residual is an independent
finite diagnostic: geometric retraction changes the root and can exceed1e-10.
Do not describe a case as satisfying an additional strict post-residual gate if
the reported post-residual exceeds that value. Existing thresholds are not relaxed.
Core `builds=` counts CN fresh callbacks including the starting RHS; it excludes
the separately reported one post-Lowdin fresh action. Actual ACE rebuild counters
include the post-Lowdin reconstruction after physical-source cache checks.

The restored-vs-original comparison is in canonical calc2026/hse-lowdin-25-50as-20261028.
Old no-Lowdin tests and source/binary snapshots remain valid records of that
previous variant; they are not the restored original orthogonalized algorithm.

## Occupied-only scalar HSE endpoint certification

Zero-temperature scalar RT initialization sets `system%no=nelec/2` for a
closed shell, independently of the GS `nstate`. GS empty states are not evolved.
The standard HSE adapter collectively requires every propagated occupation to
be finite and equal to 2 before applying the original full-column endpoint
Lowdin. It does not add empty-state propagation or a mixed-occupation rotation.
Unequal, empty, fractional and nonfinite propagated occupations are rejected.
The MGGA occupation-block machinery and existing distribution guards remain.
The projected fixed-point CN experiment is separate and is not selected here.

The local occupied-only guard and shared regressions passed. Independent
12.5/25/50as comparisons are in progress; Fujitsu portability and full-pulse
accuracy of this exact version remain unverified. Existing jobs and fixed
binaries are unchanged.
