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

## Exact local-box exchange FFT

`local-fft.json` compares `exx_local_fft='auto'` and `'off'` at the **same**
finite source radius (4 bohr). This is an implementation-parity check, not a
physical cutoff-convergence result. The native four-H fixture uses a 32x16x16
bohr cell/grid and 1x2x1 k mesh. Both runs converge at step 30; the last exchange
refresh uses 54,000 local pair-grid points instead of 262,144 global points.
The regression `test_local_fft_native_scf` reconstructs both inputs and checks
the final energy difference below 1e-8 eV and actual local-path use.

The fixed-source benchmark on a 32x24x20 grid compares the whole PBEh Wannier
exchange action, including gather/scatter, with four compact sources and eight
targets (one identically zero). The relative action difference is 2.94e-16.
28 pairs use 20,412 padded-grid points instead of 430,080; one-thread warm
exchange actions took approximately 1.24 ms locally versus 5.92 ms globally in
this run. These are five-action averages excluding kernel/plan setup, not a
production-water or end-to-end SCF scaling result. Local pairs currently execute
serially, and dense targets may still require every source-target pair.

Reproduce the standalone correctness and bounded benchmark checks:

```
python3 -m unittest discover -s testsuites/unit_hse_wannier -p test_local_fft.py
python3 -m unittest discover -s testsuites/unit_hse_wannier -p test_local_fft_benchmark.py
```

The standalone direct sum includes complex kernels, wrapped support, moving box
origins, full-grid fallback then local cache reuse, singleton support and tiny
nonzero densities. Operator checks cover HSE/PBEh, translated multi-k sources,
Hermiticity and extended-source fallback. This engine does not activate PBEh
LCFO optical propagation or DC-MD; those adapters and force consistency remain
separate requirements.


## DC-LCFO response integration

`lcfo-response.json` records the H4 impulse fixture for PBEh40 and PBEh40+rVV10,
with x/y/z spatial partitions and 2 versus 4 ranks (orbital groups). At 0.08 au,
0.02/0.01-au step refinement gives current differences around 1e-15 au and
post-impulse energy widths around 5e-14 Hartree. Integrated electron count is
4 within 2e-13. Grid-spacing input parity and per-fragment functional/run-ID
mismatch rejection are also covered by `test_pbeh_response.py`.
These tiny, short tests establish plumbing consistency, not water spectra or
production scaling. Those measurements used the original root FFT backend.
The subsequent distributed FFT tests below add direct single-grid potential comparisons.


## Distributed rVV10 convolution

The native FFTE adapter compares local energy, vrho and vsigma against serial
FFTW, then assembles the complete discrete potential for comparison with
`rvv10_periodic`. Tests cover a rectangular 16x12x8 grid, axis and combined
layouts, 2/4 ranks, 1/2 OpenMP threads, and q grids 8/16/32. They also verify
collective rejection of a density invalid on one rank, unsupported-grid
fallback, and preservation of Poisson FFT tables for a different grid size.

Native DC→LCFO response is compared with the committed root-backend trajectory
fixture; gradient/divergence use the production halo routines. The exact
measured bounds and validation status are recorded in `distributed-fft.json`.
The FFT channel storage is divided by Py*Pz, with replication across Px.
DC still assembles the scalar potential globally for fragment mapping. No
large-water timing or end-to-end scaling claim is made by these small tests.

## Cached FFTW pencils

`fftw-pencils.json` records three trials per grid for the optional
`rvv10_fft='fftw'` backend. Six measured FFTW plans per batch size are reused
across forward/inverse calls. Direct tests include 2/3/4/5 channels (including
a 4+1 tail), changing spatial layouts, comparison with FFTE and inverse
normalization. All 54 distributed functional/potential cases pass; the maximum
full-potential difference from the serial reference is 1.18e-13. Native SCF
backend selection, 25 unit tests, y/z DC-LCFO trajectories, and HSE ON/OFF builds
pass. Independent review found no critical or important defects.

Median warm forward/inverse times for 32 channels, in milliseconds:

| Grid | MPI ranks | FFTE | FFTW |
|---|---:|---:|---:|
| 32x24x16 | 2 | 5.15 | 5.59 |
| 32x24x16 | 4 | 3.61 | 3.38 |
| 64x48x32 | 2 | 47.24 | 79.52 |
| 64x48x32 | 4 | 30.52 | 43.29 |

The current FFTW packing dominates on the larger grid. Whole-functional
medians for that grid are 0.639/0.631 s (FFTE/FFTW, 2 ranks) and
0.312/0.318 s (4 ranks), with appreciable run-to-run variation. There is no
consistent whole-functional advantage; FFTE remains the default. These local
OMP=1 measurements compare current adapters, including FFTE's per-transform
table initialization, and do not establish production scaling. Setup and
component timings, all trials, ranges and log hashes are retained in the JSON.

## Retained Z spectra

`transposed-spectrum.json` records the subsequent FFTW layout optimization.
Forward spectra stay in z,x,y memory order for the rVV10 kernel; inverse FFTs
run z,y,x and return the original real-space layout. Each pair now requires
four redistributions instead of eight. The same-process, same-plan comparison
uses three trials, 32 channels and OMP=1. Median warm pair times (ms):

| Grid | MPI ranks | Old FFTW X | New FFTW Z | Reduction |
|---|---:|---:|---:|---:|
| 32x24x16 | 2 | 6.25 | 4.41 | 29.5% |
| 32x24x16 | 4 | 3.32 | 2.30 | 30.8% |
| 64x48x32 | 2 | 78.27 | 53.37 | 31.8% |
| 64x48x32 | 4 | 40.61 | 27.76 | 31.6% |

Complete-functional medians on 64x48x32 are 0.612/0.603 s
(FFTE/new FFTW, 2 ranks) and 0.317/0.311 s (4 ranks). These are small
differences relative to trial variation; FFTE stays the default. This does
not measure old-versus-new FFTW whole-functional speedup. The generic FFTW
X-to-X API remains available as a reference; `rvv10_fft='fftw'` automatically
uses retained Z spectra in SCF/DC/LCFO.

All 54 direct functional/potential cases pass, including full-potential
differences below 1.20e-13. Individual Z-layout Fourier coefficients agree
with assembled FFTE spectra; 2/3/4/5-channel inverse normalization and plan
reuse pass. The native 25-test suite and y/z DC-LCFO response fixtures pass.
Independent review found no critical or important defects. Phase maxima may
come from different ranks and must not be added as an exact breakdown.
These remain local, small-system measurements, not production-water scaling.

## DC-MD force gate (not yet MD support)

The opt-in static `yn_dc_force_diagnostic` prints the explicit nuclear derivative
at fixed fragment orbitals and occupations. Global atom IDs/images are preserved;
total-grid electrostatics and core-weighted nonlocal projector terms are
assembled with the same rank ownership as DC energy. It does not write these
uncertified values into the ordinary ionic force field.

`dc-md-force-audit.json` records the initial energy-only reference;
`dc-force-response-audit.json` adds the explicit-force comparison;
`dc-force-thermal-audit.json` also compares the diagnostic core-weighted E-TS.
At 300 K, buffer2/1 leave approximately0.05/0.08 eV/angstrom force discrepancies
after converged SCF and displacement-step refinement. The entropy term changes
but does not remove them. These are electronic-response residuals, not solely
orbital response. The full-cell buffer limit agrees within~8e-6 eV/angstrom.

Independent review found no explicit-force formula or reduction defect, but
identified the need to isolate the asymmetric projector derivative and thermal
response. The production projector contraction is now tested against frozen
complex s/p/d-channel energy finite differences with asymmetric, full and empty
core masks; the thermal audit measures both E and E-TS. This synthetic test does
not replace a native truncated-oxygen force certification. Periodic wrapping
checks do not certify moving atom-list updates or fragment-face crossings.

The force gate is not passed for truncated DC. No DC-MD integrator is enabled.
The next response formulation and outstanding validation are described in
[the response plan](../../plans/2026-09-27-dc-force-response.md).

The final force-gate validation has32 passing tests (31 regressions plus the
separately run24x24x24 water one-fragment limit); HSE-enabled MPI and HSE-disabled
serial builds pass. The water limit covers the actual H/O radial projectors,
while the asymmetric contraction is isolated by the synthetic frozen-state
test. Earlier16x16x16 water settings did not meet the requested SCF tolerance
and were rejected; those runs are not force certifications. See
`dc-force-validation.json` for measured water parity and log hashes.
