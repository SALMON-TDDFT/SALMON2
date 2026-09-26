# DC-HSE with Wannier rotations and ACE

Build with `USE_HSE=ON` and MPI; see `docs/hse-build.md`. Dependencies are FFTW,
Libxc and BLAS/LAPACK. The H4 example is in
`testsuites/422_H_dcdft_hse/inputfile` (copy `testsuites/pseudo/H_rps.dat` beside
it). It uses two fragments with buffers, a 1x2x1 k mesh, four MPI ranks,
fractional occupations, HSE fragment-SCF and complex LCFO. Run with one OpenMP
and BLAS thread for a small reproducible test:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 /path/to/salmon < inputfile > outputfile
```

The numbered case verifies SCF convergence and finite complex LCFO eigenvalues.
The additional comparison programs below require Python with NumPy and `mpiexec`;
each invocation needs a new output directory. Every MPI process has a 180-second
subprocess timeout. They do not consume stored reference wavefunctions.

```sh
python3 samples/dc_hse/validate.py /path/to/salmon /tmp/dc-hse-comparison
python3 samples/dc_hse/validate_rt.py /path/to/salmon /tmp/dc-hse-rt-comparison
python3 -m unittest discover -s testsuites/unit_hse_wannier
```

`validate.py` compares k parallelism, one-fragment/no-buffer DC against ordinary
HSE with the same backend, checks Gamma, and runs a 32x8x8 total grid with
24x8x8 fragment+buffer grids. It saves energies, LCFO eigenvalues and elapsed
smoke-test times in `results.json`. Those times are not scaling benchmarks.
`validate_rt.py` produces a fresh occupied-only GS and compares four Taylor4+ACE
steps with spread-minimization intervals 1 and 10. This verifies in-memory U
reuse and a restart after step 2; it does not validate long-time stability or
DC-to-global RT conversion.

Compiled unit tests compare rectangular/shifted/shuffled-grid exchange with an
independent reciprocal Bloch sum, fractional and zero occupations, full complex
rotations, arbitrary target Hermiticity, ACE interpolation, polar transport,
analytic spread gradients and continuous phases across ±pi. Set `FC`,
`FFTW_ROOT`, `OPENBLAS_ROOT` for your compiler/installations; the tests default
to Homebrew locations. Unit compilation requires the indicated libraries.

Read `docs/inputs/hse.md` for controls and limitations. This version retains full
spatial support and all exchange pairs. It is a correctness baseline for later
controlled locality approximations, not a demonstrated linear-scaling HSE
implementation. Some coarse-mesh/extra-state fixtures do not converge the MLWF
gradient within 200 iterations; exact full-support exchange remains valid and
this is printed explicitly.

## Si64 locality diagnostic

`si64/` contains the historical DG Si64 geometry/grid/buffer adapted to HSE06,
with file hashes and its 400-state cost implications documented. The full HSE
baseline must be converged before comparing physical SCF results.

`pair_screening.py` analyzes frozen **Gamma** density factors exported as NPZ.
It screens symmetric pairs by an L2 overlap bound using max(V_G)=pi/omega^2,
keeps all self pairs, and compares exact and screened column actions. Reported
energy is spin-summed fragment-periodic short-range HSE exchange; it is not the
DC core-weighted total energy. The script calculates full pair potentials for
verification, so its runtime is not a production screening speedup benchmark.
Its metric gate is only necessary for ACE, not an SCF correctness certificate.
No screening switch has been enabled in the Fortran SCF/LCFO path by this script.


## Snapshot and post-SCF localization

Set `SALMON_HSE_WANNIER_SNAPSHOT=1` for a native final density-factor snapshot.
For DC it is written separately under `data_dcdft/fragments/NNNNNN/`. A normal
time-limit shutdown can write an **unconverged** snapshot; successful process
exit does not establish SCF convergence. Preserve the input and output log.

```sh
python3 samples/dc_hse/read_snapshot.py hse_wannier_snapshot.bin factors.npz
```

The converter refuses an unconverged SCF reference unless explicitly given
`--allow-unconverged`. Such data may diagnose pair locality at a frozen density,
but must not be presented as converged HSE observables. The converter currently
accepts Gamma only. `highest_source_occupation` is not the occupation of the
highest originally retained eigenstate: exactly empty bands have been removed.

The optional `localize_snapshot.f90` program rebuilds a Gamma gauge offline,
without rerunning SCF. Compile it together with `hse_wannier_gauge.f90` and
`hse_wannier.f90`, linking FFTW3 and BLAS/LAPACK; then run:

```sh
localize_snapshot input.bin localized.bin 0 3000 1e-6
```

The last three arguments are occupation cutoff, maximum localization iterations,
and gradient tolerance. Cutoff 0 removes only exact zeros and preserves the
density. A positive cutoff is an additional approximation; discarded occupation
is printed and a changed density cannot inherit the SCF-converged flag. Check
localization status and original/relocalized exchange before screening. A unitary
gauge that has not met the localization criterion still gives exact full-support
exchange but must not be called a converged MLWF set.

Binary version 1 uses a native-endian marker and 14 int32 entries: marker,
version, grid(3), k mesh(3), source count, refresh count, localization iterations,
localization status, SCF iteration, SCF-converged flag. These are followed by 9
real64 entries: spacing(3), omega, spread, gradient, minimum transport singular
value, periodic exchange energy, SCF residual; then occupation, U, Phi and Q in
Fortran order. The reader detects endian order. The format is diagnostic, not a
restart format: it omits original band indices/retained-state count, actual k
vectors, and the SCF criterion/threshold. Preserve the generating input/log and
do not use reduced U rows as original band indices. Transport-only refreshes
record localization status 2 and spread/gradient -1 (not evaluated).

The multi-budget pair diagnostic shares pair FFTs across budgets. It still
computes the full pair reference, so its wall time does not measure production
screening performance. Empty-state source reduction and polar-only refreshes
are implemented acceleration steps; spatial pair pruning remains diagnostic.

## CPU execution and full-support reference

Run one simulation at a time. Eight DC fragments require eight MPI ranks with
the current one-rank-per-fragment Gamma layout. MPI8 x OpenMP2 / BLAS1 bounds
the active computational threads at 16 on the development machine. This is not
a measured optimal configuration for other hardware. OpenMP parallelism now
distributes independent target-column FFTs within each fragment. FFTW plans and
buffers are private to each worker; planning/destruction is serial, following
https://www.fftw.org/fftw3_doc/Thread-safety.html . Each target retains the same
source accumulation order. The serial build remains supported.

`SALMON_HSE_EIGEN_DIAGNOSTIC=1` exports the final Psi and refreshed-Hamiltonian
Hpsi on a single-k, full-grid/full-orbital layout before LCFO. The diagnostic
reader in locality.py computes eigen-residuals and orthogonality independently.
The v1 export does not contain actual k coordinates or convergence flags; keep
the input and matching Wannier snapshot. Its single-k guard does not itself
prove the point is Gamma. Use the documented Gamma inputs for locality studies.

`locality.py snapshot.bin --radii ... --natom 64 --output report.json` requires
both SCF and localization convergence by default. It applies an axial x mask to
every source factor consistently, without renormalizing. It reports discarded
norm, self-exchange of the changed density, exchange expectation on the original
density, and Fock-action error on the original factors. A common Hermitian Fock
operator is preserved, but this remains a frozen-density test, not a variational
SCF/force/RT accuracy certificate. The minimum circular-center reliability is
reported because an almost uniform axial density has no well-defined center.
A support radius of at least half the periodic x length is the uncut reference.


Axial support diagnostics can keep factors with undefined circular centers
in full support: `locality.py ... --min-center-reliability 0.1`. The threshold
is a dimensionless magnitude of the normalized first circular density moment,
not a WF amplitude cutoff. The same mask defines the exchange operator on all
targets. `full_support_factor_count` reports protected factors;
`max_tail_fraction` is the geometric tail before protection, whereas
`max_truncated_factor_tail_fraction` excludes protected factors and
`discarded_norm_fraction` measures the actual discarded total norm.
No support-sweep result certifies truncated SCF, forces or propagation.
