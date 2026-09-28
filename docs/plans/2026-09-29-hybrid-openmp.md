# Hybrid local-loop OpenMP implementation plan

Goal: extend the user-approved collapse/private/reduction approach from the DC
core exchange integral to independent local hybrid grid loops.

Architecture: retain MPI/FFT/BLAS ownership and numerical algorithms. Parallelize
only disjoint grid writes or scalar reductions, with explicit data scoping and
static scheduling. Keep the contiguous coordinate innermost. No large per-thread
arrays, input changes, or nested FFT/BLAS threading.

Tech stack: existing Fortran/OpenMP, MPI, FFTW, Libxc and ScaLAPACK CPU build.

1. Audit EXX native/spatial/orbital/ACE and rVV10 loops. Select native pack/unpack,
   action addition, spatial reciprocal multiplier, and rVV10 local spline/input
   and potential/output loops. Existing rVV10 reciprocal convolution is already
   parallel. Leave MPI/FFT source-pair loops, prefix-compaction loops and BLAS
   matrix products unchanged.
2. Replace pack/unpack reshapes and action's running grid index with direct
   indices, collapse four outer loops, retain contiguous x. Add collapsed spatial
   kernel loop with direct z-pencil index. Parallelize independent rVV10 points;
   spline_basis uses only local scratch and immutable tables.
3. Build MPI/HSE/Libxc/ScaLAPACK and serial HSE. Compare short DC GS for four
   functionals against the pre-change result, and native RT at OMP1/2 for three
   global hybrids including multi-k and Gamma adaptive source ACE. Check charge,
   convergence, currents/energies, and numerical tolerances. Keep total test
   concurrency bounded because spectra are running in other processes.
4. Record verification and limitations. No scaling improvement claimed without
   a controlled benchmark; this change does not alter asymptotic complexity.

## Results

Implemented six additional regions (pack, unpack, action addition, reciprocal
multiplier, rVV10 input and output), plus the previously requested core integral.
The reciprocal multiplier keeps z contiguous, while native grid loops keep x
contiguous. Shared arrays have disjoint writes; only the core scalar uses a
reduction. Small nq-sized spline temporaries remain thread-local.

Validation on 2026-09-29:
- MPI/HSE/Libxc/ScaLAPACK build and non-MPI HSE build passed.
- All three changed Fortran modules compiled without OpenMP enabled.
- `test_conventional_rt`: 10 tests passed, including MPI2 x OMP1/2 comparisons
  for DC HSE06/PBE0/PBEh40/rVV10 and Gamma adaptive/full-k native RT.
- DC comparison requires converged SCF, charge within 1e-10 electrons, and
  final exchange within 1e-9 Ha. RT current/energy comparison uses atol=1e-8,
  rtol=1e-7, with electron-number conservation checked independently.
- Independent static review found no required changes; diff whitespace checks
  passed. OpenACC/GPU configurations were not run (these native hybrid routes
  retain their existing CPU-only restrictions).

No controlled timing benchmark was run. Small loops can be dominated by thread
startup, and pack/unpack speed is ultimately memory-bandwidth limited. This is
local threading, not a change to pair counts, MPI communication, or N-scaling.
