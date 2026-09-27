# PBEh(40)+rVV10 validation — 2026-09-27

Base `4048d67f`; CPU Apple arm64, GNU Fortran 15, FFTW, Libxc, OpenBLAS. MPI-enabled build used for water with one rank and one OpenMP/BLAS thread. `USE_HSE=OFF`, non-MPI build also compiles. Existing compiler warnings in legacy sources remain.

## Kernel and integration checks

- 10 compiled MLWF tests: inherited tests plus unscreened Gamma, multi-k/fractional occupations, explicit Coulomb radius and too-large-radius rejection. Independent Bloch-sum comparison, unitary-gauge invariance, Hermiticity, negative exchange metric and occupied-space ACE reproduction are included.
- 2 compiled functional tests: 0.6 PBE X + PBE C versus independent Libxc calls; rVV10 Fourier kernel versus independent radial quadrature, rho/sigma finite differences, uniform-density cancellation, vacuum, periodic translation and q-mesh convergence at 16/32/64/128 channels.
- 6 executable rejection tests: unconverged BOMD, legacy snapshot, DC+rVV10, restart, invalid rVV10 channels, retained HSE MD rejection.
- 4 inherited native exchange/semilocal/PT-CN/ACE tests pass.
- Existing 422 DC-HSE fixture on 4 MPI ranks completes; final density change 5.1200636e-9, below 1e-8.

## Water force and short NVE results

One water molecule, periodic 5 Å cube, 24^3 grid, Gamma, 4 occupied states, supplied H/O test pseudopotentials, b=5.3, C=0.0093, 32 spline channels, automatic exchange radius 2.5 Å. SCF density threshold 1e-10. These parameters are for consistency checks, not converged liquid-water predictions.

| Check | Result |
|---|---:|
| H x-force, analytic | -0.82023734 eV/Å |
| H x-force, central energy difference (±0.001 Å) | -0.82024237 eV/Å |
| H absolute discrepancy | 5.03e-6 eV/Å |
| O x-force, analytic | 3.01790690 eV/Å |
| O x-force, central energy difference (±0.001 Å) | 3.01750304 eV/Å |
| O absolute discrepancy | 4.04e-4 eV/Å |
| NVE max energy deviation, dt=0.1 fs, 4 steps | 1.827e-5 eV |
| NVE max energy deviation, dt=0.05 fs, 8 steps | 3.326e-6 eV |

Both NVE paths cover only **0.4 fs** with identical deterministic initial velocities at 300 K. Halving dt reduces the maximum deviation by about 5.5 times. This does not validate long trajectories, thermodynamics or diffusivity. All SCFs in the eight accepted calculation logs converged. Machine-readable results and log hashes accompany this file. Long calculation output stays in the external work directories.

## Why projector forces changed

Before the PBEh force correction, the 24^3 H force differed from its energy derivative by 0.5396 eV/Å. Omitting rVV10 gave essentially the same 0.5402 eV/Å discrepancy. The inherited force formula applies finite-difference gradients to orbitals; at finite grid spacing it does not exactly differentiate the interpolated projector used in the Hamiltonian. Differentiating that radial cubic and its solid spherical harmonic directly reduced the discrepancy to the value above. This change is restricted to PBEh.

An exploratory 32^3 job crossed a binary rebuild and is not included as validation evidence; no grid-convergence claim is made. A separate initial 24^3 MD job demonstrated the old larger force error; it is not included among the accepted trajectories.

The inherited BOMD path also advanced with deliberately unconverged `nscf=1`. PBEh now rejects that case before accepting a force/ionic step; the executable regression test checks the diagnostic.

## Remaining work

Total-density DC+rVV10, DC-MD forces, restart metadata, real-space distribution, controlled pair pruning, large-water-cell scaling, and physical cell/radius/grid/time convergence remain unimplemented or unvalidated as detailed in the [input documentation](../../inputs/pbeh40-rvv10.md). The inherited MLWF/ACE path is retained; no linear-scaling claim is made.

## Static DC total-density rVV10

`dc-static.json` records the 16x8x8 hydrogen-cell integration check from
`testsuites/unit_pbeh_rvv10/dc_hydrogen.inp`. The two-fragment buffers cover the
full cell, so this checks communication and energy accounting, not convergence
of a genuinely truncated fragment approximation. One-fragment DC agrees with
conventional PBEh40+rVV10 within 2e-6 eV; two/four-rank calculations also agree
within that tolerance. The regression additionally checks a nonzero energy
change when rVV10 is removed. The periodic wrapper's complete density derivative
(including density gradients) passes the 2e-8 atomic-unit tolerance.

Reproduce with an HSE-enabled executable and a local MPI launcher:

```
SALMON_TEST_EXE=/absolute/path/salmon SALMON_TEST_MPIEXEC=/absolute/path/mpiexec \
  python3 -m unittest discover -s testsuites/unit_pbeh_rvv10
```

This validates static DC only. Root-only nonlocal FFT performance, realistic
water fragment/buffer convergence, DC forces, and checkpoint restart remain open.

## Shared EXX controls and spherical source support

The canonical controls are now `exx_mlwf_interval/maxiter/tolerance`, with legacy
HSE input aliases, and `exx_mlwf_radius` in input length units. See
[the input contract and approximation limits](../../inputs/exx-mlwf.md).
Standalone probes cover HSE and PBEh kernels, zero/full limits, periodic masks,
protected ambiguous centers, Hermiticity, compact/dense action agreement and
loss of overlap during gauge transport. Input tests cover alias equality,
conflicts, units and rejected modes, as well as a converged small conventional
finite-radius HSE/PBEh SCF.

Finite-radius DC is not certified by the previous full-support numbers. A new
check exposed MPI-dependent initial seeds and unconverged MLWF gauges. The
seeds now use physical k indices for finite support. The same first DC update
agrees across two/four ranks for Gaussian and random starts. A truncated SCF
result is rejected if the last MLWF minimization did not converge, including
a failed gauge transport that invalidates earlier convergence. This prevents a
small density residual from being mistaken for localization convergence.
