# Native HSE06 on the TDCDFT branch

This experimental CPU implementation integrates HSE06, full screened exchange,
ACE and Taylor4 predictor/corrector into SALMON. It is not an upstream release.
The initial supported configuration is fixed-ion, unpolarized periodic systems,
fully occupied spin pairs, a cubic real-space grid and a uniform cubic k mesh (full or supported symmetry-reduced mesh). Only k-point MPI decomposition is supported. NLCC, spin-orbit, DFT+U,
spatial/orbital MPI decomposition, finite-temperature occupations and ionic
dynamics are rejected. The RT path supports a transverse impulse. On full k meshes, linear Acos2
laser pulses with an optional impulse probe are supported; pulse RT restart
is rejected. Symmetry-reduced impulse RT requires a field-preserving z group.

HSE RT automatically selects Taylor4 + ACE with predictor/corrector when
`propagator` is omitted. The entire `&propagation` namelist may be omitted.
The fourth-order exponential polynomial averages the local and exchange
operators between the initial and predicted endpoint. Its self-consistent
predictor/corrector has generally second-order global time accuracy. Assess
time-step convergence, electron count and overlaps; orbitals are not renormalized.

See the [implementation and efficiency note](../../docs/hse-implementation-notes.md)
for the current Taylor+ACE baseline, algorithm details, and measured tradeoffs.

## Build

HSE is enabled by default in a normal SALMON CPU CMake build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON
cmake --build build -j 8
```

Compatible installed libraries are used automatically. Missing FFTW, Libxc or
BLAS/LAPACK are downloaded and built in the build directory. No dependency paths
or compiler workaround flags are normally needed. MPI itself must be installed
when `USE_MPI=ON`; omit that option for serial execution.
See [dependency details](../../docs/hse-build.md).

## Si example

From the repository root, prepare a separate calculation directory:

```sh
mkdir -p calculations/my_native_hse/gs calculations/my_native_hse/rt
cp testsuites/pseudo/Si_rps.dat calculations/my_native_hse/Si_rps.dat
cp samples/hse_native/gs.inp calculations/my_native_hse/gs/inputfile
cp samples/hse_native/rt.inp calculations/my_native_hse/rt/inputfile
```

Run in `gs`, then in `rt`, using the absolute path to the built executable:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 8 /path/to/salmon < inputfile > run.log
```

Both inputs use Si8, 12³ grid points and 4³ k points. Change `nproc_k` and the
MPI process count together. `xc='hse06'` selects 25% short-range exchange with
omega=0.11 bohr^-1 by default, and the matching Libxc semilocal remainder.
RT uses Taylor4 + ACE, dt=0.16 au and a z-directed impulse of 10^-4 au.
See [HSE input](../../docs/inputs/hse.md) for configurable screening and units.
The 100-step RT example is a stability pilot, not a resolved optical spectrum.

For SCF initialized from *only* imported wavefunctions, set
`yn_restart='y'`, `read_gs_restart_data='wfn'`,
`yn_reset_step_restart='y'`, and `nscf_init_no_diagonal=0`.
Resetting the iteration initializes density mixing history; initial subspace
diagonalization also permits a noncanonical occupied gauge. Do not use these
reset settings for genuine RT continuation.

For RT continuation, retain the original input physics, set `yn_restart='y'`,
point `directory_read_data` to an accepted `checkpoint_rt_XXXXXX/`, and increase
`nt` to the desired total step count. Keep `hse_restart.bin` with the checkpoint.
The code verifies dt, impulse, geometry, k points, occupations and exact
pseudopotential bytes, then rebuilds ACE and FFT caches. Resetting the time or
freezing the functional is rejected. The propagator does not renormalize wavefunctions.

## Validation and timing

Standalone numerical parity tests:

```sh
SALMON_TEST_MPI=1 python3 -m unittest discover -s samples/hse_mlwf_reference -p 'test_native*.py'
```

The native tests use a Fortran compiler, FFTW, OpenBLAS and Libxc, with paths
configurable through `FC`, `FFTW_ROOT`, `OPENBLAS_ROOT`, and `LIBXC_ROOT`.
Optional profiling uses `SALMON_HSE_PROFILE=1`. The full sampled exchange
kernel is retained; these speedups do not rely on a new MLWF support cutoff.
Taylor4 + ACE at dt=0.32 au was unstable on the tested Si grid. Use a smaller step and check time-step convergence for your system.

## What ACE does

At an exchange refresh, the full screened Fock operator is applied to the
current occupied orbitals U to obtain W=K[U]U. ACE constructs low-rank factors
from the occupied metric -U^dagger W (with grid weights). Repeated actions on
trial vectors then use two BLAS products, -B(B^dagger X), instead of another
full exchange calculation. This reproduces K[U] on the construction subspace;
it is approximate for general trial directions. It is neither a local
potential nor a Wannier localization/support-cutoff approximation.

Taylor keeps the initial ACE for its predictor and averages initial/predicted
ACE operators for its corrector, then refreshes at the accepted state. Only exactly identical
source arrays bypass refresh through the cache. See the
[ACE construction and update schedule](../../docs/hse-implementation-notes.md#exchange-and-ace).

## Distributed exchange memory

Source, target, action, phase and ACE factors now remain on their owning k
ranks. Only density-matrix row tiles are transposed for the exchange FFT.
Taylor snapshots are freed after use. This reduces memory but adds communication;
see the [implementation note](../../docs/hse-implementation-notes.md).

### Internally threaded BLAS

HSE calls BLAS outside application OpenMP regions. Link a threaded BLAS and
configure its runtime (for the tested OpenMP OpenBLAS: `OMP_NUM_THREADS=12`,
`OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=12`). HSE does not override the vendor
thread count. These are local OpenBLAS settings, not a certified Fugaku job script.
MPI/FFTW calls remain outside OpenMP; no per-thread orbital replicas are allocated.
Measure the MPI/OpenMP balance on the target machine.

## Input and restart smoke tests

Use an existing converged GS matching the Si example:

```sh
python3 samples/hse_native/check_input_smoke.py --exe /path/to/salmon \
  --gs /path/to/data_for_restart --work /tmp/new-hse-smoke-directory
```

Add `--off-exe /path/to/serial-hse-disabled-salmon` for the PZ odd-restart and
stale-metadata regression checks. Custom-omega and PZ runs reuse orbitals only
for wiring checks, not for physical spectra. Each run advances one or two steps.
