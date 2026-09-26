# Si64 chain locality and compiler isolation plan

Goal: rearrange the existing Si64 crystal to 8x1x1 conventional cells and DC
fragments, measure axial WF tails against exchange errors, and isolate the user's
GNU Fortran 15/AArch64 loop-vectorization report before interpreting SCF failure.
User confirmed the 64-atom rearrangement. Keep a=10.26 bohr, h=.64125 bohr,
128x16x16 total grid, 16^3 core, Gamma and 8 MPI ranks. Start with x buffers
8/16/24 (no y/z buffer), making fragment periods 2a/3a/4a. These avoid the
noncommensurate 1.75a fragment boundary in the former six-grid buffer. This does
not by itself prove buffer convergence.

1. Keep BLAS/OpenMP threads=1 and external BLAS linkage fixed. Compare O3 against
   O3 -fno-tree-loop-vectorize with identical inputs. Verify actual build flags.
   Distinguish loop vectorization from SLP and external BLAS implementation.
2. Add opt-in Gamma fragment Psi/Hpsi diagnostics and an independent residual/
   orthogonality analysis. First demonstrate missing-output failure on a small
   HSE case, then build both variants and verify parsed dimensions and residuals.
3. Run buffer-8 paired pilots, compare early and final density residuals,
   eigenstate residuals, overlap matrices, density and eigenvalue differences.
   Do not infer a compiler bug merely from differing late SCF trajectories.
4. Converge buffer 8/16/24 using the justified build and identical physical
   controls. Keep iteration limits and failures explicit. Separately localize
   saved states and verify density/exchange invariance.
5. Measure periodic x-center, axial tail fractions and norm-containment radii
   for Phi and occupation-weighted Q. Sweep axial support radii and measure
   frozen-density exchange/target-action errors without feeding an inconsistent
   operator into SCF. There is no universal acceptable spread threshold.
6. Preserve inputs, hashes, reports, tests and unresolved limitations. Review
   implementation, run relevant regressions and commit locally.

## Execution constraint after user feedback

Run only one simulation at a time. Four overlapping 8-rank jobs were excessive;
three were stopped, retaining only the bounded PZ damping diagnostic. Per-rank
CPU usage rose from roughly 50% to 98% after removing contention. This is not a
measured HSE speedup. Keep BLAS threads fixed, and compare MPI/OpenMP settings
serially using iteration wall time before scaling the test matrix. Do not
interpret manually terminated jobs as numerical failures.

The user's GNU15/AArch64 report matches docs/hse-platforms.md: Netlib 3.12.1
ZLARF1L, with an existing fallback-only -fno-tree-loop-vectorize workaround.
The actual current library is Homebrew OpenBLAS. The repository's reproducer
hse_lapack_eigenvectors passes in both tested SALMON builds with relative
residual 2.0512586062782591e-15. Turning off SALMON loop vectorization is not
equivalent to rebuilding that previously faulty Netlib routine.

## User-directed full-cell baseline

The user requested an undivided 8x1x1 full calculation first, then gradually
smaller support. The next run uses yn_dc=n, the actual Si64 82.08x10.26x10.26
cell/grid128x16x16, Gamma, 128 occupied orbitals (zero-temperature closed-shell),
full-support Wannier HSE and ACE. This differs from the fractional-occupation DC
pilots; it is not a controlled same-temperature DC comparison. One MPI rank
(current backend requires full grid/orbitals), 8 OpenMP threads, 8 BLAS threads
(non-nested). No concurrent simulations. Preserve every 10th SCF checkpoint and
use normal time-limit exit so ordinary full-cell restart is available. Establish
SCF and eigenstate residuals before calling any support cutoff acceptable.

## Corrected meaning of full calculation

User clarified that the 8 DC fragments must remain, requiring 8 MPI ranks.
The previous yn_dc=n / MPI1 run was a mistaken interpretation of full and is
not the requested baseline. It has finished and must not be resumed for this
request. Full means no spatial support truncation and no pair pruning within
the 8-fragment DC calculation. Use one job, MPI8; keep total OpenMP/BLAS
thread count bounded, with BLAS1. Establish this reference before shrinking
WF support.

## Implemented verification checkpoint

- Independent Gamma/single-k eigen-residual export and analysis; H4 occupied
  RMS residual about 3e-8 Ha, orthogonality about 1e-15 at converged density.
- Axial support sweep uses one common Fock operator for consistently masked
  source factors, reports changed-density self-exchange separately from the
  original-density expectation, and never certifies SCF from these tests alone.
- OpenMP target-column FFTs: worker-private plans/buffers, serial planning and
  destruction, fixed per-target accumulation order. Review found no blocker;
  an empty-target guard was added after review.
- 15 unit tests pass, including independent Bloch exchange, serial/2-worker
  equality, source selection, support endpoints, Hermiticity and offline
  localization. HSE422 and complex LCFO130 CTest prep/run/verify: 6/6 pass with
  OMP2/BLAS1. Serial compilation remains exercised by the offline-localizer test.
- Current eight-rank DC run uses roughly 1.6-1.7 CPU cores per rank; this is CPU
  activity, not a controlled runtime speedup measurement. No concurrent jobs.
- Full-support DC SCF remains unconfirmed; the observed mixed-density residual
  is not an unmixed fixed-point residual and does not replace eigenstate checks.
