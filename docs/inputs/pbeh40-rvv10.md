# PBEh(40)+rVV10 through MLWF and ACE

Experimental implementation on `pbeh40-rvv10-water-md`, based on `dc-hse-mlwf-ace` (`4048d67f`).

## Model and inputs

```
&functional
  xc='pbeh40_rvv10'        ! or 'pbeh40' to omit nonlocal correlation
  pbeh_coulomb_radius=0    ! input length; 0 selects half the shortest BvK supercell side
  rvv10_b=5.3
  rvv10_c=0.0093
  rvv10_nq=32
  exx_mlwf_interval=5
  exx_mlwf_maxiter=100
  exx_mlwf_tolerance=1d-7
  exx_mlwf_radius=0        ! input length; 0 retains full orbital support
  exx_local_fft='auto'    ! exact local convolution when compact support is cheaper
/
&parallel
  nproc_rgrid=1,1,1
  nproc_ob=1
  nproc_k=1
/
```

The functional is 0.4 Fock exchange + 0.6 PBE exchange + PBE correlation + rVV10 nonlocal correlation. It is not HSE06 with a modified mixing fraction. `hse_omega` is not used for the PBEh exchange kernel. MLWF is enabled automatically for both PBEh functional names; the inherited occupation weighting, gauge transport, localization updates, full-support pair-density FFT and ACE remain in use.

The short-range damping parameter `b=5.3` is the water-hybrid choice; `C=0.0093` is dimensionless. Check the parameter convention against the target study. For example, see [Evolution of Aqueous Electron with Varying Temperature](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/62ab9df7f5524a36fb1528e8/original/evolution-of-aqueous-electron-with-varying-temperature.pdf).

### Exchange boundary convention

The unscreened Coulomb interaction is spherically truncated at R:

`v(G) = 8*pi*sin(|G|*R/2)^2 / |G|^2`, `v(0)=2*pi*R^2`.

R must not exceed half the shortest side of the Born–von Karman supercell (`num_rgrid*num_kgrid*hgs`). The default is that upper bound, printed in bohr. Positive input radii are converted from the selected input length unit. The kernel includes Coulomb within R without HSE screening. A finite R is still an approximation to the infinite periodic global hybrid: converge **cell/k mesh and radius together**. Never interpret this cutoff as evidence for converged full-range exchange in a small cell. No arbitrary zeroing of the G=0 exchange term occurs.

MLWF gauge localization by itself retains full support. `exx_mlwf_radius>0`
selects the shared [EXX spherical source approximation](exx-mlwf.md) for static
DFT. This radius truncates occupation-weighted Wannier sources; it is distinct
from `pbeh_coulomb_radius`, which changes the interaction kernel. At radius 0,
localization iteration-limit messages do not invalidate full-support exchange.
Finite-radius results depend on localization and must be converged against the
support radius. ACE reproduces the chosen (possibly masked) operator on its
reference subspace; it does not remove the truncation error.

### rVV10 evaluation

In Hartree atomic units:

- `omega0 = sqrt(4*pi*n/3 + C*sigma^2/n^4)`, `sigma=|grad n|^2`.
- `kappa = (3*pi/2)*b*(n/(9*pi))^(1/6)`, `q=omega0/kappa`, `theta=n/kappa^(3/2)`.
- `Phi(q,q',r)=-3/[2*(1+q*r^2)*(1+q'*r^2)*(2+(q+q')*r^2)]`.
- `E_nl = 1/2 integral theta(r)*Phi*theta(r') dr dr' + beta integral n dr`, `beta=(3/b^2)^(3/4)/32`.

The implementation uses natural cubic spline channels, a logarithmic q grid from 1e-4 to 0.5, and the usual 12-term smooth q saturation. Increase `rvv10_nq` (8–128) to check interpolation convergence. Values below the lower q bound are clamped with zero q derivative. Points at density <=1e-18 bohr^-3 contribute only the beta term. The energy and potential differentiate this same regularized discrete functional.

The rational kernel's analytic three-dimensional Fourier transform is used, including its G=0 value and equal-q limit. No radial table or image cutoff is needed. Periodic convolution costs O(nq^2*G+nq*G*log G), not a G-by-G pair matrix.

Compatible grids now use the native Poisson FFTE layout: full x lines and distributed y/z pencils. Density and sigma are collected only within the x communicator. Channel storage per rank is O(nq*G/(Py*Pz)); the work is replicated across Px. **Use y and/or z decomposition to reduce per-rank FFT memory.** An x-only decomposition retains full channel storage on each rank. FFT tables are private to rVV10 and do not alter the Poisson solver's saved tables.

The automatic path requires uniform block decomposition, 2/3/5-smooth grid dimensions from 2 through 4096, and FFTE transpose divisibility (Nx divisible by Py, Ny divisible by Pz, in addition to each dimension's own process count). Incompatible grids retain the root FFTW reference path. The chosen backend and x/y/z process counts are logged. No additional namelist switch is required.

DC evaluates the mixed **total** density with native halo gradients and divergence; the global nonlocal energy is summed once. Its scalar potential is still assembled globally for fragment/buffer mapping. Thus DC still holds full scalar grids, while the nq-channel convolution is distributed. This does not by itself establish whole-program weak scaling.

`vrho` and `vsigma` are added before SALMON's GGA divergence. Density gradients and the corresponding negative divergence therefore use the same finite-difference operator. The formula follows [Sabatini, Gorni and de Gironcoli, PRB 87, 041108(R) (2013)](https://doi.org/10.1103/PhysRevB.87.041108).

## Fixed-cell water dynamics

Keep `exx_mlwf_radius=0` and set `theory='dft_md'` for Born–Oppenheimer MD. A fresh calculation performs its initial SCF; restart input is currently rejected. The supplied [water_md.inp](../../samples/pbeh40_rvv10/water_md.inp) is a **small consistency fixture**, not an equilibrated liquid or a converged production setting. Copy `H_rps.dat` and `O_rps.dat` from `testsuites/pseudo` to its run directory. For scientific water simulations, choose and validate appropriate pseudopotentials, grid, supercell, exchange cutoff, SCF tolerance and timestep.

PBEh forces differentiate the same cubic radial projector and solid spherical harmonics used in the Hamiltonian. This avoids the finite-grid inconsistency of moving the nonlocal projector derivative onto a finite-difference orbital gradient. All other functionals retain their existing force path. At full support (`exx_mlwf_radius=0`), no explicit ionic derivative of exchange/rVV10 is needed for fixed-cell, fixed-grid, fully self-consistent orbitals without NLCC. A PBEh BOMD run stops if an SCF does not converge, before accepting an ionic step.

## Supported and rejected combinations

- Periodic, orthorhombic, unpolarized conventional DFT and fixed-cell BOMD, CPU, k-only MPI; uniform full k meshes.
- Static `pbeh40` and `pbeh40_rvv10` use the inherited DC MLWF+ACE exchange path; DC convergence against fragment/buffer size is still needed.
- **Not supported:** DC MD, conventional RT/Ehrenfest, ionic optimization, spin polarization, NLCC, OpenACC, variable-cell stress/NPT, restarting PBEh checkpoints, or legacy HSE Wannier snapshot export. Projector angular momentum above f is rejected by the PBEh force routine.
- DC+rVV10 convolves the **total density**. The initial fragment orbital preparation omits this term until the first regular total-density SCF update. DC-MD remains disabled: this static integration does not establish variational forces for truncated fragments.
- No production scaling, long liquid trajectory, diffusivity, RDF, density, or exchange-cutoff convergence claim follows from the bounded tests below.

## DC-LCFO electronic response

A fixed-nuclei, Gamma-point, occupied-only response can now start from the
complex LCFO files of a static DC calculation with the same `xc`. Set
`theory='tddft_response'`, `yn_dc='n'`, `yn_conventional_from_dcdft='y'`,
and `yn_hse_lcfo_rt='y'`. The historical LCFO flag name is retained.
Omit electronic temperature and use `nstate=nelec/2`. Each spatial rank must
match one saved fragment core; orbital MPI groups are supported. The default
propagator is `hse_taylor4` with predictor/corrector. Use an impulse field.
This is electronic response in a fixed LCFO subspace, with no initial SCF
reconvergence; it does not enable DC-MD.

Both the projected exchange and its continuity diagnostic use 40% unscreened
spherical-cutoff Coulomb exchange. Automatic Coulomb radius is half the shortest
**fragment plus buffer** side. Transported MLWF sources and coefficient-space
ACE use the existing LCFO algorithms. This first PBEh RT route requires
`exx_mlwf_radius=0` in the generating GS and `hse_lcfo_wf_radius=0` in RT.
A finite-source-radius RT model needs separate convergence/energy validation.

For rVV10, the density gradient and potential divergence use the same halo
exchange and finite-difference stencils as the spatially decomposed Laplacian.
The convolution uses the same native FFTE pencil adapter as SCF/DC on compatible
grids. Derivatives return to each owned core and enter the existing halo
divergence. Orbital groups do not multiply nonlocal energy. The reference
root FFT remains a compatibility fallback, as described above.

New complex LCFO exports include a `functional.txt` per fragment, bound to the
binary run ID. Reconstruction checks the functional, relevant exchange/rVV10
parameters and source radius for every fragment. PBEh requires these files;
legacy HSE data without them retain their previous route. Radius input values
are compared after unit conversion (automatic zero and an explicit equivalent
radius are deliberately distinct input settings). RT restart/checkpoint output
remains unsupported.

Run the integration fixture with:

```
python3 testsuites/unit_lcfo_rt/test_pbeh_response.py \
  --binary /absolute/path/to/salmon --pseudo /absolute/path/to/H_rps.dat
```

Direct distributed energy, density derivatives and full-potential comparisons
against serial FFTW are available with:

```
python3 testsuites/unit_pbeh_rvv10/test_distributed.py --build /absolute/path/to/build
```

This uses the configured MPI/HSE build objects and tests 2/4 ranks, individual
and combined spatial axes, q grids of 8/16/32, OpenMP 1/2, collective rejection,
fallback and preservation of independent Poisson FFT tables.

Use `--xc pbeh40` for exchange-only hybrid coverage and `--axis y` or `--axis z`
to exercise the corresponding domain boundaries. The fixture checks time-step
refinement, orbital decomposition and functional metadata rejection, not a
converged water absorption spectrum.

## Build and verification

Build with `USE_HSE=ON` (FFTW, Libxc and BLAS/LAPACK). Existing HSE dependency discovery is reused; `USE_LIBXC=ON` is not additionally required for this C-ABI semilocal evaluator. Example local configuration:

```
cmake -S . -B build -DUSE_HSE=ON -DUSE_MPI=ON \
  -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 6
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s testsuites/unit_hse_wannier -p test_wannier.py -v
SALMON_TEST_EXE=/absolute/path/to/build/salmon \
  python3 -m unittest discover -s testsuites/unit_pbeh_rvv10 -v
python3 samples/pbeh40_rvv10/validate.py build/salmon /fresh/force-H --grid 24
python3 samples/pbeh40_rvv10/validate.py build/salmon /fresh/force-O --grid 24 --atom O
python3 samples/pbeh40_rvv10/validate.py build/salmon /fresh/nve --grid 24 --md
```

The standalone kernel tests currently use the same macOS Homebrew toolchain layout as existing native HSE tests. Run directories must be fresh. Results and limitations are recorded in [validation results](../results/pbeh40-rvv10/README.md).
