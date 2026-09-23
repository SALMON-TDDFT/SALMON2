# Native HSE06 on the TDCDFT branch

This experimental CPU implementation integrates HSE06, full screened exchange,
ACE and self-consistent PT-CN into SALMON. It is not an upstream release.
The initial supported configuration is fixed-ion, unpolarized periodic systems,
fully occupied spin pairs, a cubic real-space grid and a full uniform cubic
k mesh. Only k-point MPI decomposition is supported. NLCC, spin-orbit, DFT+U,
spatial/orbital MPI decomposition, finite-temperature occupations and ionic
dynamics are rejected. The RT path currently supports a transverse impulse;
laser pulses and pump-probe HSE are not yet enabled.

## Build

Enable `-DUSE_HSE=ON -DUSE_MPI=ON` in a normal SALMON CPU CMake build.
Supply `-DFFTW_INSTALLDIR=/path/to/fftw` and
`-DLIBXC_INSTALLDIR=/path/to/libxc` if these libraries are not on the search path.
FFTW's `fftw3.f03` and the Libxc C library are required; Libxc's compiler-specific
Fortran modules are not required. Use a BLAS/LAPACK library such as OpenBLAS.
The default `USE_HSE=OFF` build has no new dependency.

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
omega=0.11 bohr^-1, and the matching Libxc semilocal remainder. RT selects
`propagator='hse_ptcn'`, dt=0.32 au and a z-directed impulse of 10^-4 au.
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
freezing the functional is rejected. Accepted steps must satisfy the fresh
full-exchange endpoint residual (<10^-10), electron count and overlap gates;
the propagator does not renormalize accepted wavefunctions.

## Validation and timing

See `docs/results/si-hse-native/README.md`. Standalone numerical parity tests:

```sh
python3 -m unittest discover -s samples/hse_mlwf_reference
```

The native tests use a Fortran compiler, FFTW, OpenBLAS and Libxc, with paths
configurable through `FC`, `FFTW_ROOT`, `OPENBLAS_ROOT`, and `LIBXC_ROOT`.
The runtime prints `HSE_PT_CN` for residual/iteration/whole-step timing and
`HSE_TIMING` for exchange, ACE, exchange communication, local Hamiltonian and
preconditioner costs. The full sampled exchange kernel is retained: these
speedups do not rely on a new MLWF support cutoff.
