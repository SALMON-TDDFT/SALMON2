# Native HSE06 / ACE / MPI validation (2026-09-24)

Implemented in the TDCDFT branch: native Fortran screened exchange, ACE,
Libxc HSE06 semilocal remainder, SCF Hamiltonian/energy integration and a
self-consistent PT-CN real-time solver. This is experimental branch support,
not upstream acceptance. The full sampled exchange kernel is retained.

## Matched short-run performance

Si8, 32 electrons, 16 occupied orbitals, cubic 12³ grid, 4³ k mesh, dt=0.32 au,
z impulse=10^-4 au; same converged initial orbitals and endpoint tolerances.
Apple M5 Pro, 18 cores, 64 GB; gfortran 15.2, Open MPI 5.0.9, FFTW, Libxc 7,
OpenBLAS. OpenMP and BLAS threads fixed to one per MPI process.

| Implementation | MPI ranks | Median seconds/step, first 3 steps |
|---|---:|---:|
| Python reference | 8 | 13.802 |
| Native Fortran | 1 | 14.577 |
| Native Fortran | 4 | 4.628 |
| Native Fortran | 8 | 2.947 |

Native MPI8 is **4.68× faster** than matched Python MPI8 and **4.95× faster**
than native MPI1. These are three-step pilots, with the independent Python8
production reference running in the background, not isolated-node scaling
measurements. Propagation timing excludes checkpoint I/O. No k convergence
scan was performed. Earlier Python timings at later trajectory steps are not
used as the denominator.

Native MPI8 first-step costs (seconds): full exchange 1.750; ACE build 0.013;
ACE applications 0.138; exchange gather/reduction 0.076; remaining local
Hamiltonian/potential work 0.594; FFT preconditioner 0.081. Other orchestration
costs account for the remainder. Timings are local rank-zero wall times;
communication includes the large exchange data collectives, not every solver
reduction. Logical PT-CN refresh counts include cache hits; there are three
actual full exchange builds per step in this pilot.

Full source orbitals are replicated during exchange refresh; row blocks have
unique MPI owners and their actions are reduced. ACE factors and repeated local
Hamiltonian applications stay on the owning k rank. A resident-memory snapshot
during the separate MPI8 stability run was 295–298 MiB per rank (sum 2.31 GiB).
This is a snapshot, not peak RSS or unique physical memory; shared pages can be
counted repeatedly. The full source arrays still scale with k count, so this
implementation does not establish asymptotically local scaling.

## Numerical evidence

- Fixed Si reference: native Hpsi relative error 8.90e-11; total energy
  difference 6.39e-14 Ha; total XC energy difference 3.73e-14 Ha.
  Total energy -31.31474671493668 Ha; screened exchange -1.65293882461363 Ha.
- The Libxc semilocal derivative can change by a few 10^-10 under one-ulp
  density perturbations. The cross-language potential/action test uses 2e-9;
  the PT-CN fresh full-exchange residual gate remains 1e-10.
- Native MPI1/4/8 after three steps: wavefunction relative differences from
  Python below 5.7e-12; current absolute differences below 3.7e-15 au.
- Continuous three steps versus two steps plus a fresh-process restart:
  wavefunctions are bitwise identical. ACE/FFT caches are rebuilt.
- The 100-step native MPI8 impulse pilot (0.774 fs) completes with maximum
  fresh full-exchange residual 9.47e-11, final electron-count error 5.79e-11
  and maximum orbital-overlap error 7.40e-11. The post-impulse energy range is
  2.60e-12 Ha at ten-step sampling. The initial impulse deposits 1.51e-7 Ha;
  that physical energy increase is excluded from the drift measure.
- Changed dt, impulse or pseudopotential bytes, missing/truncated HSE metadata and
  frozen-functional input are rejected before an accepted RT step.
- Native SCF started from the independently converged reference returns to the
  same fixed point in five iterations (density criterion 7.70e-10).
  For a wavefunction-only import, reset SCF iteration/mixing history and use
  initial subspace diagonalization; see the sample README.
- Native SCF from SALMON's own initial orbitals converges after 87 updates
  (density criterion 7.37e-9). An independent full-exchange evaluation gives
  -31.31474671493524 Ha, only 1.43e-12 Ha above the Python reference, with
  projected residual 9.82e-7 Ha. The imported-state SCF endpoint has projected
  residual 9.52e-8 Ha. These checks use the unmixed final orbital density.
- 78 Python/numerical tests pass, including standalone Fortran full-kernel,
  ACE, semilocal and nonlinear PT-CN parity. Existing Si/TDCDFT CTest suite:
  six tests pass. Both HSE-enabled and HSE-disabled builds succeed.

Native E_xc includes semilocal and exact exchange once. The nonlocal-ion energy
inferred from orbital expectations subtracts twice the exchange energy, so
adding exchange to Hpsi does not double count it. The GGA potential uses
SALMON's negative-divergence convention. Accepted PT-CN steps are not rescaled;
electron-count and overlap checks remain active.

Machine-readable results: `operator_parity.json`, `step_parity.json`,
`restart_parity.json`, `restart_guards.json`, `scf_validation.json`,
`stability.json`, and `timings.json`.
Raw development fixtures are under ignored `calculations/si_hse_native/`.
Reproducible inputs and build instructions are in `samples/hse_native/`.

Initial support is k-only MPI, CPU, cubic uniform grids/k meshes, fixed ions,
unpolarized fully occupied spin pairs, no NLCC, and transverse impulse RT.
Laser HSE, broader geometry/spin/process layouts and physical support cutoffs
require separate validation. These short native tests do not establish a
resolved exciton spectrum; the existing long Python HSE/TDCDFT comparison
remains the independent reference.
