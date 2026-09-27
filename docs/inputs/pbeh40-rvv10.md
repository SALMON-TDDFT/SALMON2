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

The automatic path requires uniform block decomposition, 2/3/5-smooth grid dimensions from 2 through 4096, and FFTE transpose divisibility (Nx divisible by Py, Ny divisible by Pz, in addition to each dimension's own process count). Incompatible grids retain the root FFTW reference path. The chosen backend and x/y/z process counts are logged. The optional `rvv10_fft='ffte'|'fftw'` selects the compatible-grid backend; the default is `ffte`. Both choices preserve this fallback contract.

DC evaluates the mixed **total** density with native halo gradients and divergence; the global nonlocal energy is summed once. Its scalar potential is still assembled globally for fragment/buffer mapping. Thus DC still holds full scalar grids, while the nq-channel convolution is distributed. This does not by itself establish whole-program weak scaling.

`vrho` and `vsigma` are added before SALMON's GGA divergence. Density gradients and the corresponding negative divergence therefore use the same finite-difference operator. The formula follows [Sabatini, Gorni and de Gironcoli, PRB 87, 041108(R) (2013)](https://doi.org/10.1103/PhysRevB.87.041108).

## Fixed-cell water dynamics

Keep `exx_mlwf_radius=0` and set `theory='dft_md'` for Born–Oppenheimer MD. A fresh calculation performs its initial SCF; restart input is currently rejected. The supplied [water_md.inp](../../samples/pbeh40_rvv10/water_md.inp) is a **small consistency fixture**, not an equilibrated liquid or a converged production setting. Copy `H_rps.dat` and `O_rps.dat` from `testsuites/pseudo` to its run directory. For scientific water simulations, choose and validate appropriate pseudopotentials, grid, supercell, exchange cutoff, SCF tolerance and timestep.

PBEh forces differentiate the same cubic radial projector and solid spherical harmonics used in the Hamiltonian. This avoids the finite-grid inconsistency of moving the nonlocal projector derivative onto a finite-difference orbital gradient. All other functionals retain their existing force path. At full support (`exx_mlwf_radius=0`), no explicit ionic derivative of exchange/rVV10 is needed for fixed-cell, fixed-grid, fully self-consistent orbitals without NLCC. A PBEh BOMD run stops if an SCF does not converge, before accepting an ionic step.

## Supported and rejected combinations

- Periodic, orthorhombic, unpolarized conventional DFT and fixed-cell BOMD, CPU, k-only MPI; uniform full k meshes.
- Static `pbeh40` and `pbeh40_rvv10` use the inherited DC MLWF+ACE exchange path; DC convergence against fragment/buffer size is still needed.
- **Not supported:** direct truncated-fragment MD, unrestricted conventional RT/Ehrenfest, ionic optimization, spin polarization, NLCC, OpenACC, variable-cell stress/NPT, restarting PBEh checkpoints, or legacy HSE Wannier snapshot export. Projector angular momentum above f is rejected by the PBEh force routine.
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
The convolution uses the same native pencil adapter as SCF/DC on compatible
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

## Reusable FFTW pencil backend

Set `rvv10_fft='fftw'` in `&functional` to select cached FFTW transforms for
SCF, static DC and the supported LCFO response route. FFTW uses measured
`plan_many` plans, serial FFTs per MPI rank, and batches of at most four channels.
Plans and packing maps are reused until the grid, process coordinates or batch
size changes. Calls must enter serially per rank; axis communicators must be
coordinate ordered and collective dimensions/channel counts must agree.
The forward transform now retains Z pencils with contiguous z lines and
local shape `(Nz,Nx/Py,Ny/Pz)` in z,x,y order. The rVV10 kernel uses the
corresponding physical wavevectors; inverse FFTs run in z,y,x order and
return the original real-space layout. This halves redistribution calls
from eight to four per forward/inverse pair. The generic transform still
supports its original X-to-X layout for reference comparisons.
No FFTW MPI or FFTW threads library is required. The backend does not change
the functional metadata or enable DC-MD.

FFTE remains the default. Benchmark the selected adapter: transform speed
and complete-functional speed can differ, and packing/communication still matter. Measure the complete functional on the
target machine before selecting a backend. Setup, warm forward/inverse pairs,
local FFT, communication, packing and full rVV10 timings are separated by:

```
python3 testsuites/unit_pbeh_rvv10/test_distributed.py --build /path/to/build --benchmark --scale 1
python3 testsuites/unit_pbeh_rvv10/test_distributed.py --build /path/to/build --benchmark --scale 2
```

The grids are 32x24x16 and 64x48x32, with 32 channels, 2/4 MPI ranks and one
OpenMP thread. FFTE timings include its current per-transform private-table
initialization; these compare implemented paths, not isolated library speed.
Use `--fft-backend fftw` with the LCFO integration fixture to check that route.

## Static DC force diagnostics (direct fragment MD remains disabled)

For static PBEh DC calculations, `yn_dc_force_diagnostic='y'` in `&dc` prints
a separately labelled frozen-orbital nuclear derivative in Hartree/bohr. It
requires a converged SCF, full MLWF support, and automatic fragment atom lists.
It sums total-grid electrostatics and the two derivatives of each
core-restricted nonlocal projector product, with global atom identities.
It does not populate the regular ionic force or advance any ion.

The diagnostic also reports total electron residual, core-weighted occupation
`TS`, and `E-TS`. These are diagnostic quantities; stationarity of `E-TS` for
truncated fragments has not been established. Fragment orbital and occupation
response terms are absent from the printed force. Periodic-image bookkeeping
is present, but moving atom-list updates and boundary-crossing MD are not.

Run the reproducible small H4 audit with:

```
python3 testsuites/unit_pbeh_rvv10/audit_dc_force.py \
  --binary /path/to/salmon --mpiexec mpiexec --output /path/to/results.json
```

At 300 K electronic temperature, full-cell buffers agree with energy finite
differences. Truncated buffers leave about 0.05–0.08 eV/angstrom disagreement
even after the tested core-weighted entropy term is subtracted. Therefore
`theory='dft_md', yn_dc='y'` is still rejected. This diagnostic is preparation
for an energy-consistent DC force, not production MD support.

For positive `temperature_k`, PBEh DC occupations now enforce the core-weighted
electron count with a bounded chemical-potential solve. Insufficient weighted
state capacity stops with a request to increase `nstate_frag`; failure is no
longer silently accepted. The finite-temperature response formulation keeps
electronic temperature fixed and differentiates E-TS at fixed total electron
count. Occupation/entropy derivative kernels are available internally, but the
coupled orbital/density force correction and ionic integration are still pending.
The existing zero-temperature and non-PBEh occupation paths are unchanged.

## DC initial state to real-space Ehrenfest

A bounded native route now uses `theory='tddft_response'`, `yn_md='y'`,
`yn_dc='n'`, `yn_conventional_from_dcdft='y'`, and `yn_hse_lcfo_rt='n'`.
The existing reader reconstructs DC-LCFO initial states onto the real-space
mesh and checks functional/run metadata. Thereafter the wavefunction is the
ordinary `spsi%zwf` mesh array; it is never projected back to the LCFO subspace.
MLWF/ACE accelerates the time-dependent exchange action. Initial integer
occupations stay fixed; no electronic-temperature fitting or SCF occurs in RT.

Use `ensemble='NVE'`, `step_update_ps=1`, `out_rt_energy_step=1`, full EXX support,
occupied-only states, the default hybrid Taylor4 predictor/corrector, and
`ae_shape1='impulse'` or `ae_shape1='Acos2'` without a second field. For Acos2,
use `theory='tddft_pulse'`, positive `omega1` and `tw1`, nonnegative `t1_start`,
and linear transverse polarization. Supply initial velocities normally.
Checkpoint continuation/output remain rejected in this route. Initial force uses the post-impulse nonlocal phases. The initial
energy output retains SALMON's pre-impulse electronic reference; evaluate
post-excitation conservation using Eall+Tion after the kick. E_work in the RT
file is ionic mechanical work, not laser work.

Ionic positions are advanced by the existing Verlet steps. During electronic
propagation, local/nonlocal pseudopotentials use midpoint positions; endpoint
pseudopotentials are rebuilt before energy and force evaluation. There are no
moving-basis terms in this fixed real-space mesh representation.

See `samples/pbeh40_rvv10/water_ehrenfest_gs.inp` and
`water_ehrenfest_rt.inp` for a verified small water setup. Native exchange now supports the Gamma y/z spatial layout described below.
Orbital distribution, finite EXX support and long trajectories are not certified.
This is DC preparation followed by total-grid Ehrenfest, not propagation of
independent truncated-fragment forces. The old direct DC-BOMD guard remains.

For Acos2, propagation uses midpoint A and nuclear positions; energy, current
and forces use endpoint A. The ionic electric field is the centered difference
`E(t)=-(A(t+dt)-A(t-dt))/(2*dt)`, including initialization. Ionic current uses
completed Verlet velocities. Independently integrate
`volume*(Jion-Jmatter) dot E` to compare external work against Eall+Tion.
No electronic temperature is assigned during excitation. The H4 pulse fixture
shows second-order energy/work and final-current convergence; this does not
establish long-time or liquid-water accuracy.


### Spatial ACE algebra (development stage)

`hse_ace_build` and `hse_ace_apply` accept an optional `sum_grid` callback,
which sums a complex matrix over the spatial communicator in place. Each rank
stores only its local grid rows of factors, sources and targets. Construction
reduces the occupied-state metric before factorization; application reduces
factor/target overlaps before the local action. No full-grid array is gathered.
Peers must use matching orbital/k dimensions, grid volume and call order;
local row counts may differ and may be zero. Average ACE combines local factor
columns and retains the same spatial partition. Memory is O(local_grid *
ACE_rank * local_k), with small replicated metric/overlap matrices; all occupied
columns remain present locally, so this does not provide orbital decomposition.

This API is verified independently with MPI1/2/4, complex states, non-unit grid
volume, unequal/empty partitions, midpoint operators, zero exchange and bad
local data. It is also used by the Gamma spatial mesh RT route below; other layouts
retain their existing admission restrictions.


### Gamma spatial mesh Ehrenfest RT

For DC-initialized PBEh40/PBEh40+rVV10 native mesh RT, set `nproc_k=1`,
`nproc_ob=1`, and for example `nproc_rgrid=1,2,2` with four MPI ranks.
The path requires an unshifted Gamma point, orthorhombic periodic cell,
fully occupied spin pairs, `exx_mlwf_radius=0`, and the default Taylor4+ACE.
The x dimension is complete on each pencil; y/z are distributed. Grid sizes
must satisfy the FFTW pencil divisibility constraints. Static DC fragments,
projected LCFO RT and multi-k calculations retain their prior layout support.

MLWF localization reduces six band-overlap matrices; polar transport reduces
the temporal overlap matrix. Each rank retains only local grid rows of the
localized source and previous frame. An identity initial gauge replaces the
serial pivoted seed; at full support this changes localization history but not
the exchange operator. The usual exx_mlwf_interval/maxiter/tolerance apply.
The full-support pair convolution uses the existing distributed FFTW backend,
with up to four target columns per batch and the same spherical Coulomb cutoff
as the serial path. ACE construction and action reduce only small band matrices.
No full wavefunction-grid gather occurs during these exchange RT operations.
`exx_local_fft` does not select compact source convolutions on this path.

Initial DC-to-mesh reconstruction now streams bounded destination-grid chunks
as described below. Validated fragment basis and coefficient records remain resident.
All occupied columns remain local to each spatial rank, and small band matrices
are replicated. The verified MPI1/2/4 H4 pulse and MPI1/2 water smoke cases are
correctness checks, not evidence for giant-system timing or long-time accuracy.
Use `samples/pbeh40_rvv10/h4_ehrenfest_pulse_spatial.inp` after its matching GS.


### Bounded DC-to-mesh reconstruction

Complex DC initialization validates all input metadata and fragment core
coverage before reconstructing orbitals. Coverage and complex contributions
are processed in flat destination-grid chunks of at most65,536 points. Each
chunk is reduced only to its owning spatial/orbital rank. No global coverage
array or global single-orbital reduction buffer is allocated by reconstruction.
The complex scratch consists of two buffers (at most2 MiB combined with the
tested 16-byte complex representation). Two integer coverage buffers are freed
before these complex buffers are allocated. Small index/metadata arrays and
resident fragment basis/coefficient records are additional memory.

The diagnostic `DC_LCFO_TILE scratch_points/global_points` reports buffer
capacity in grid points and the global grid size. This is not total process
memory. File metadata/provenance validation and projected LCFO configuration
remain unchanged. Wrapped grid maps are preserved, and missing or duplicated
core coverage is rejected collectively. Each destination/chunk is handled in
sequence; fragment intersections are rescanned across chunks, so this memory
change is not a claim of faster initialization. Streaming fragment input and
more efficient sparse routing remain separate work.
