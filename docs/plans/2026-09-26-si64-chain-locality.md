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

## DC fragment boundary correction (2026-09-26)

The first MPI8 full-support baseline stopped normally at iteration 545 after
600.868 s without convergence (mixed density residual 1.5019237e-3). Final-current-H
occupied RMS eigen residuals were 0.604--2.347 Ha despite orthogonality within
2.721e-13. These are not frozen-H inner-solver diagnostics.

Inspection found a separate, reproducible geometry defect: equivalent fragments
had atom counts [16,16,18,14,16,16,16,16], giving initial electron counts 72 and
56 in fragments 3 and 4 instead of 64. The direct binary floating-point boundary
comparison reproduced those counts, while decimal arithmetic gave all 16.

Use a shared per-axis roundoff allowance 32*epsilon*total-cell-length on both
faces of each half-open fragment box, including lower and excluding upper
boundary atoms without modifying coordinates. It is a numerical boundary
convention, not a physical buffer expansion. The actual MPI8 one-iteration
regression in testsuites/dc_fragment_boundary/check.py failed on the old binary
with the observed 18/14 atom counts and passed after recompilation. Complex
LCFO130, HSE422 and LAPACK eigenvector CTest checks: 7/7 passed serially.

A single unchanged-input baseline is rerunning in dc-full-boundary-fixed with
MPI8/OMP2/BLAS1, full support and no pair pruning. Geometry correctness is proven;
SCF convergence and a valid WF cutoff remain to be established separately.

## Literature-guided next convergence investigation

- Lin, Adaptively Compressed Exchange Operator, JCTC 12, 2242--2249
  (2016), DOI 10.1021/acs.jctc.6b00092.
  https://math.berkeley.edu/~linlin/publications/ACE.pdf
  Distinguish a fixed-ACE linear eigensolver from nested density SCF with
  fixed exchange followed by an outer exchange update.
- Hu, Lin, Yang, Projected Commutator DIIS Method for Accelerating Hybrid
  Functional Electronic Structure Calculations, JCTC 13, 5458--5467
  (2017), DOI 10.1021/acs.jctc.7b00892.
  https://math.berkeley.edu/~linlin/publications/PCDIIS.pdf
  Gauge-invariant projected density-matrix/commutator mixing; compatible
  with ACE, demonstrated for HSE06 silicon. Adaptation to DC with global
  chemical potential and fractional occupations needs separate derivation
  and validation; the projector formulas cannot be copied unmodified.
- Kudin, Scuseria, Cances, A black-box self-consistent field convergence
  algorithm: One step closer, JCP 116, 8255--8261 (2002),
  DOI 10.1063/1.1470195. EDIIS followed by DIIS is a broader robustness
  reference, not direct validation of this DC implementation.

Actual current call path: CG -> hpsi -> hse_add_action -> hse_ace_apply;
refresh_wannier builds ACE outside CG after density/occupation updates.
There is no explicit converged inner density SCF at fixed ACE. Proposed
next diagnostic: measure fixed-H residual immediately after solve_orbitals,
then after occupation, local-potential and exchange updates. Do not attribute
all nonconvergence to a particular update until those stages are measured.
No new mixing scheme has been implemented in this boundary-fix change.

## Boundary-fixed baseline outcome

Completed normally by native time shutdown after 601.339 s, iteration 573.
Final mixed density residual: 9.8743525e-4 (target 1e-7), NOT converged.
Final-current-H occupied RMS eigen residual across fragments: 0.449--1.358 Ha;
maximum orthogonality error: 1.077e-13. This does not establish a CG inner
failure because the Hamiltonian is refreshed after the solve. The corrected
geometry alone does not resolve SCF instability. No WF cutoff sweep certified.

The expanded one-iteration MPI regression also verified all 16 periodic atom
sites are distinct and equal across the eight fragments. It passed after the
baseline finished; no concurrent simulations were run. Original input, output,
executable/source/input SHA256 provenance and per-fragment binary diagnostics
are retained under work/si64-chain/dc-full-boundary-fixed (workspace root).


## Locality measurements on the converged reference

Offline localization retained all38 positive-occupation states (no occupation
cutoff). A1e-6 gradient target stalled near1e-5 in line search; it is NOT
claimed reached. At a declared2e-5 target all8 fragments converged and exchange
was invariant to rounding. A1e-4 comparison gives the same spread to about1e-9.

Six Q factors have circular x-center reliability of order1e-7, versus about.91
for the other32. Cropping all38 gives arbitrary-center sensitivity up to1.06
meV/atom when changing localization tolerance. The diagnostic now optionally
retains factors below reliability.1 in full support. A regression with a
uniform factor and a localized factor failed before the option and passes
afterward. With these six factors protected, tolerance sensitivity is below
5.1e-9meV/atom and all8 fragments agree closely.16 unit tests pass.

Worst whole-fragment exchange errors at axial halfwidths9,8,7bohr are2.3833,
6.4016,20.6859meV/atom; corresponding discarded norm fractions are.0004311,
.0008825,.0029261. These are fixed-density full-periodic-16-atom-fragment
exchange differences, NOT DC-core total-energy errors or truncated-SCF results.
No certified support cutoff or production pair-pruning speedup is claimed.
Data: samples/dc_hse/si64-chain/locality-results.json. Full snapshots and raw
sweeps: workspace work/si64-chain/localized. No production solver changes.
