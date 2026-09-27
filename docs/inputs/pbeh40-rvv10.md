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
  hse_mlwf_interval=5
  hse_mlwf_maxiter=100
  hse_mlwf_tolerance=1d-7
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

MLWF gauge localization does not itself truncate orbital support here. No new distance/pair pruning is enabled. Exact occupied-space reproduction by ACE is checked independently of localization convergence. Localization iteration-limit messages therefore do not invalidate full-support exchange, but do indicate that locality is not fully optimized.

### rVV10 evaluation

In Hartree atomic units:

- `omega0 = sqrt(4*pi*n/3 + C*sigma^2/n^4)`, `sigma=|grad n|^2`.
- `kappa = (3*pi/2)*b*(n/(9*pi))^(1/6)`, `q=omega0/kappa`, `theta=n/kappa^(3/2)`.
- `Phi(q,q',r)=-3/[2*(1+q*r^2)*(1+q'*r^2)*(2+(q+q')*r^2)]`.
- `E_nl = 1/2 integral theta(r)*Phi*theta(r') dr dr' + beta integral n dr`, `beta=(3/b^2)^(3/4)/32`.

The implementation uses natural cubic spline channels, a logarithmic q grid from 1e-4 to 0.5, and the usual 12-term smooth q saturation. Increase `rvv10_nq` (8–128) to check interpolation convergence. Values below the lower q bound are clamped with zero q derivative. Points at density <=1e-18 bohr^-3 contribute only the beta term. The energy and potential differentiate this same regularized discrete functional.

The rational kernel's analytic three-dimensional Fourier transform is used, including its G=0 value and equal-q limit. No radial table or image cutoff is needed. Periodic convolution costs O(nq^2*G+nq*G*log G), with O(nq*G) arrays, not a G-by-G pair matrix. The conventional implementation uses a full grid on each k rank. In DC, owned portions of the mixed total density are gathered, the total-communicator root evaluates rVV10, and its potential is broadcast and mapped to every fragment (including buffers). The global nonlocal energy is added once, after semilocal core accumulation. This is not a distributed real-space rVV10 solver: the root uses O(nq*G) storage and all ranks hold full scalar grids.

`vrho` and `vsigma` are added before SALMON's GGA divergence. Density gradients and the corresponding negative divergence therefore use the same finite-difference operator. The formula follows [Sabatini, Gorni and de Gironcoli, PRB 87, 041108(R) (2013)](https://doi.org/10.1103/PhysRevB.87.041108).

## Fixed-cell water dynamics

Set `theory='dft_md'` for Born–Oppenheimer MD. A fresh calculation performs its initial SCF; restart input is currently rejected. The supplied [water_md.inp](../../samples/pbeh40_rvv10/water_md.inp) is a **small consistency fixture**, not an equilibrated liquid or a converged production setting. Copy `H_rps.dat` and `O_rps.dat` from `testsuites/pseudo` to its run directory. For scientific water simulations, choose and validate appropriate pseudopotentials, grid, supercell, exchange cutoff, SCF tolerance and timestep.

PBEh forces differentiate the same cubic radial projector and solid spherical harmonics used in the Hamiltonian. This avoids the finite-grid inconsistency of moving the nonlocal projector derivative onto a finite-difference orbital gradient. All other functionals retain their existing force path. No explicit ionic derivative of exchange/rVV10 is needed for fixed-cell, fixed-grid, fully self-consistent orbitals without NLCC. A PBEh BOMD run stops if an SCF does not converge, before accepting an ionic step.

## Supported and rejected combinations

- Periodic, orthorhombic, unpolarized conventional DFT and fixed-cell BOMD, CPU, k-only MPI; uniform full k meshes.
- Static `pbeh40` and `pbeh40_rvv10` use the inherited DC MLWF+ACE exchange path; DC convergence against fragment/buffer size is still needed.
- **Not supported:** DC MD, LCFO projection/RT, Ehrenfest/RT, ionic optimization, spin polarization, NLCC, OpenACC, variable-cell stress/NPT, restarting PBEh checkpoints, or legacy HSE Wannier snapshot export. Projector angular momentum above f is rejected by the PBEh force routine.
- DC+rVV10 convolves the **total density**. The initial fragment orbital preparation omits this term until the first regular total-density SCF update. DC-MD remains disabled: this static integration does not establish variational forces for truncated fragments.
- No production scaling, long liquid trajectory, diffusivity, RDF, density, or exchange-cutoff convergence claim follows from the bounded tests below.

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
