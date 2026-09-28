# Hybrid distribution kernel measurements

This harness measures synthetic production-kernel workloads, **not** water SCF,
MD trajectories, full SALMON runtime, or a production-scale memory bound.
It uses GNU Fortran/MPI build objects on macOS/Linux and the existing
`unit_lcfo_rt/peak_rss.c`. It does not modify SALMON production code.

## Run

Use an existing GNU MPI build with HSE and ScaLAPACK enabled. The usual configure
switches are `--enable-mpi --enable-scalapack`; this harness is not validated
with the Fujitsu compiler.

```sh
python3 testsuites/benchmark_hybrid_distribution/run.py \
  --build /absolute/path/to/build --output /absolute/path/to/smoke.json --smoke

python3 testsuites/benchmark_hybrid_distribution/run.py \
  --build /absolute/path/to/build --output /absolute/path/to/scaling.json \
  --grid 32 --states 16 --dense-dim 1024 --repeat 3

python3 testsuites/benchmark_hybrid_distribution/run.py \
  --build /absolute/path/to/build --output /absolute/path/to/more-states.json \
  --grid 32 --states 64 --repeat 3 --phase exchange
```

Run measurements sequentially on an otherwise idle allocation. Set any required
process-placement flags in `--mpiexec`, e.g. an MPI launcher command and its
options. The harness forces OpenMP and common BLAS thread-count variables to one.
If vendor ScaLAPACK and FFTW flags cannot be recovered from the build cache,
provide `--link-flags` with compatible FFTW/ScaLAPACK/BLAS/LAPACK libraries.
For a standard system installation an example is
`--link-flags='-lfftw3 -lscalapack -llapack -lblas'`.
`--phase dense` skips exchange, and `--phase exchange` skips the dense solver.
Sizes and timeout are user controls, not implementation memory ceilings.

## Workloads and checks

Exchange state counts must be even; the near-cubic factorization must have between
two and `grid` modes on each axis. These layouts require an even grid size.

Exchange input is a deterministic orthonormal tensor Fourier wave-packet set with pair rotations generated directly
on owned grid points and orbital columns; there is no replicated reference
wavefunction. The grid spacing is 1 bohr, occupations are two, and MLWF refresh
uses at most three minimizer iterations. Iteration count, status, spread and
gradient are recorded. A nonconverged unitary gauge can still yield the correct
full-support exchange operator; its timing does not represent converged water
localization. Screened exchange uses omega=0.11; Coulomb uses the existing
half-cell cutoff. No hybrid mixing, semilocal functional, or rVV10 is evaluated.
ACE action must reproduce its exact source action, and the unscaled exchange
trace must agree between decompositions. The same streamed-column algorithm is
used even for one orbital group, so its one-rank reference is not the legacy
serial native path.

Dense input is a deterministic complex Hermitian matrix. **Every rank generates
all input blocks** in a temporary N-by-32 buffer, retaining only its local
block-cyclic entries. This input-generation timing is not a benchmark of the
production fragment/halo assembly. Solver time includes full diagonalization,
orthogonality and residual checks; its baseline is one-rank ScaLAPACK, not
LAPACK or CheFSI. All eigenvalues are requested. The trace is compared across
layouts. These checks complement, rather than replace, the independent
`unit_hse_ace/validate_exchange.py` and `unit_lcfo_scalapack/validate.py` oracles.

## Interpretation

Every repetition starts fresh MPI processes. Phase times include cold first-call
work such as FFT plan initialization and validation collectives. Report the
maximum rank time per phase, its median, and its range over repeats. These are
not steady-state Taylor-step timings. Different rank layouts may incur different
communication, FFT buffer, and library workspace costs.

RSS is the lifetime high-water resident set of each MPI process, in bytes.
The baseline is recorded after MPI initialization. Neither subtracting that
baseline nor adding peaks isolates a specific allocation: **sum of rank peaks
is not simultaneous node RSS**. Shared pages, MPI/BLAS/FFTW runtime state and
per-rank overhead contribute. Per-array analytical storage and process peak RSS
must not be conflated.

JSON contains every rank's raw measurements, phase summaries, workload settings,
compiler/MPI/build metadata, and SHA-256 hashes of linked production objects and
probe sources, with their source snapshots. External library resolution remains environment-dependent.
A completed file ends with `complete: true`; interrupted runs retain only the
completed records. Archive the JSON together with the source version and the
actual build configuration.
