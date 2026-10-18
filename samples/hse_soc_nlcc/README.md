# Native HSE06 with SOC and fixed-core NLCC

This developer fixture validates a two-atom periodic Si cell using an unmodified,
fully relativistic OpenMX Si_PBE19.vps with nonzero core density. It is a deliberately
coarse 8x8x8 real-space Gamma calculation, not a converged Si or WSe2 optical spectrum.
The source is a development branch based on 0f1bddf3c; it does not replace any Fugaku
production executable.

## Supported path

- CPU, periodic fixed nuclei, no DC; HSE06 with `spin='noncollinear'` and `yn_spinorbit='y'`.
- Full uniformly weighted native k meshes, k-only MPI, `nproc_ob=1`, `nproc_rgrid=1,1,1`.
- `nstate=nelec`: each occupied spinor has occupation 1 (not the scalar spin-pair occupation 2).
- GS and Taylor4 + ACE or direct-reference Taylor4 RT, including the existing SOC Hamiltonian/current.
- NLCC adds the fixed, unpolarized core density and its gradient to semilocal XC only.
  Exact exchange uses valence spinors only. XC energy and both density/gradient derivatives
  use the same total density. Ionic motion/relaxation with hybrid NLCC is rejected.
- Reduced q meshes, crystal symmetry, learned ACE, fractional occupations, finite spatial
  exchange support, GPU exchange, DC, and spinor checkpoint continuation are outside this path.

The spinor exchange forms all four spin-density blocks, including the two off-diagonal
blocks, through the existing distributed density-tile kernel. ACE overlaps sum both spin
components. This correctness-first implementation performs four scalar density-block passes;
large WSe2 dense-k performance and memory usage have not been established.

The noncollinear local XC matrix is time-centered in the HSE predictor/corrector, alongside
Hartree and Fock. Its density eigenbasis has a defined identity rotation at zero magnetization.
The SOC ionic nonlocal energy is evaluated directly from the pseudopotential, so it must not
receive the scalar path's subtraction of twice the Fock energy.

## Reproduce

Use an independent local build with MPI, HSE, FFTW and the Libxc C library enabled. From the
source root, obtain the original potential from the URL in `pseudo.json`, verify its SHA256,
and run (the output path must not already exist):

```bash
python3 samples/hse_soc_nlcc/validate.py \
  --salmon /absolute/path/to/new-build/salmon \
  --pseudo /absolute/path/to/Si_PBE19.vps \
  --output /absolute/path/to/new-validation-directory --threads 2
```

This runs scalar and spinor NLCC GS to density residual <1e-7, then scalar RT and spinor
ACE/direct/zero-field/half-time-step RT. Each weak-field trajectory covers 2 atomic time units
(100 steps at dt=0.02, or 200 at dt=0.01). A larger impulse and dt=0.08/0.04/0.02 test the
corrector's time order. The stronger-field test is a numerical stress test, not a linear spectrum.
Logs, exact inputs, binary/potential hashes and `results.json` are retained.

The independent compiled probes are:

```bash
python3 -m unittest discover -s samples/hse_mlwf_reference -p 'test_native_spin_nlcc.py'
```

Set `RUN_MPI_TESTS=1` to include two-rank k distribution. Local default library paths match
`test_native_exchange.py`; override `FC`, `MPIFC`, `FFTW_ROOT`, `OPENBLAS_ROOT`, `LIBXC_ROOT`,
and `MPIEXEC` as needed. Unit probes cover frozen-core Libxc agreement, energy finite
differences, equal-spin/unpolarized limits and vacuum, off-diagonal spinor Fock action against
an explicit grid sum, SU(2) covariance, closed-shell exchange, ACE interpolation and k distribution.

## Interpretation

GS convergence is checked explicitly; process exit 0 and restart files alone are insufficient.
The short RT runs verify the implementation, not optical energy resolution, mesh convergence,
or accuracy of a particular material/pseudopotential pair. Si SOC is weak; exchange algebra is
also tested with generic spinors having strong off-diagonal spin density. The discrete ionic
operators used by the legacy scalar and j-dependent pseudopotential implementations are not
certified as interchangeable; this work does not claim equality of their total energies after
manually zeroing SO projectors. No pseudopotential is modified in the reproducible suite.

Before WSe2 production, select fully relativistic norm-conserving W/Se potentials with recorded
valence counts, reconverge GS, and establish real-space, vacuum, k-point, time-step and duration
convergence for the intended spectrum. Use a new source/build area on Fugaku. Preserve the old
binary, job 52445014, and its inputs and locks.
