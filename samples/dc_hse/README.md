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
