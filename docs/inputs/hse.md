# Native HSE input

HSE is enabled by default in ordinary CPU CMake builds. Compatible installed
Libxc, FFTW and BLAS/LAPACK are selected automatically; missing dependencies
are downloaded and built locally. See [build instructions](../hse-build.md).
The ordinary system, pseudopotential, real-space/k-space grid, SCF and parallel
namelists remain necessary. The HSE-specific functional selection for both GS
and RT is:

```fortran
&functional
  xc = 'hse06'
  hse_omega = 0.11d0
/
```

`hse_omega` is the positive finite screening parameter in **bohr^-1** for
`unit_system='au'` or `'a.u.'`, and **Angstrom^-1** for `'A_eV_fs'`.
Internally, SALMON converts it to bohr^-1 before exchange or Libxc evaluation:
`omega_internal = omega_input * 0.52917721067` in A_eV_fs.
The omitted default is always standard HSE06 (0.11 bohr^-1), equivalent to
approximately 0.2078698738 Angstrom^-1. For example:

```fortran
&units
  unit_system = 'A_eV_fs'
/
&functional
  xc = 'hse06'
  hse_omega = 0.2078698738003611d0
/
```

All other dimensionful inputs must also follow the selected unit system. The SR kernel is
`erfc(hse_omega * |r-r'|) / |r-r'|`; increasing omega shortens its range.
It is distinct from any numerical spatial truncation radius.
The short-range exact-exchange fraction remains 0.25. Omega=0 is unsupported.
Changing omega defines a modified HSE functional rather than standard HSE06.
Both the Fock kernel and Libxc's HF/PBE screening parameters use the input value.

For GS use `&calculation theory='dft' /`. For Taylor4+ACE RT use:

```fortran
&calculation
  theory = 'tddft_response'
/
&control
  directory_read_data = 'path/to/gs/'
/
&tgrid
  dt = 0.16d0
  nt = 6000
/
&emfield
  trans_longi = 'tr'
  ae_shape1 = 'impulse'
  e_impulse = 1d-4
  epdir_re1 = 0d0, 0d0, 1d0
/
```

The time step and length above are examples from the Si tests, not universal
accuracy settings. HSE RT defaults to Taylor4 + ACE with predictor/corrector;
the entire `&propagation` namelist may be omitted. No ACE switch is necessary.
Explicitly incompatible polynomial orders or disabling predictor/corrector are
rejected. Other functionals retain their existing propagation defaults. Current native restrictions include
periodic, unpolarized, fixed-ion systems, occupied-only states and fixed occupations.
Do not add `xname` or `cname` to HSE.

Use the same omega for the converged GS and subsequent RT. Changing it requires
reconverging the GS. RT checkpoint metadata records omega and rejects a changed
value on restart; GS wavefunction loading alone does not certify functional
consistency. The input value is broadcast to MPI ranks and written to variables.log.

Validation: the native semilocal regression covers default 0.11 and custom 0.06,
0.20 against direct Libxc calls (both screening parameters set), and rejects zero,
negative, NaN and infinity. The full native exchange test file passes four tests.

MPI8 one-step input smoke tests also pass: omitted and explicit 0.11 produce
identical current/energy files, 0.20 changes the result, and zero is rejected
before propagation. These tests reuse a default-omega GS only to check wiring;
the changed-omega trajectory is not a physical production result.

Unit conversion validation (MPI8, one step): A_eV_fs with omitted omega and
with explicit 0.2078698738003611 Angstrom^-1 reproduces the atomic-unit total
energies within 1e-10 Ha. Input 0.2 Angstrom^-1 is logged internally as
0.105835442134 bohr^-1. Default and explicit-value tests both passed.

Default-propagator validation: MPI8 one-step Si calculations with omitted
`&propagation` and explicit Taylor4/predictor-corrector give identical current
and energy files. Explicitly disabling predictor/corrector is rejected. The
non-HSE PZ default remains middlepoint without predictor/corrector. Developer
propagation paths remain available for historical validation.

## Fragment-periodic DC-HSE and the Wannier backend

This branch adds `yn_dc='y'` with `xc='hse06'`. Each fragment **including its
buffer** has periodic screened exchange and its own ACE. Hartree/density assembly
retains the existing DC implementation. Complex orbitals are used at Gamma too.
See the runnable [small example and validation](../../samples/dc_hse/README.md).

The following `&functional` inputs select/control the new backend:

| Variable | Default | Meaning |
|---|---|---|
| `yn_hse_wannier` | `'n'` | Select the rectangular full-mesh Wannier backend. Automatically enabled for DC-HSE. |
| `exx_mlwf_interval` | `10` | Positive number of exchange refreshes between spread minimizations. Polar transport/reconstruction is performed on every changed source. |
| `exx_mlwf_maxiter` | `200` | Positive maximum number of unitary spread-descent iterations per minimization. |
| `exx_mlwf_tolerance` | `1d-6` | Positive finite spread-gradient norm tolerance, always in atomic units (bohr squared), independent of `unit_system`. |
| `exx_mlwf_radius` | `0` | Periodic spherical source support, in input length units; 0 is full support. Positive values automatically enable the Wannier backend and support static DFT only. |

`exx_local_fft='auto'` selects a reduced padded FFT for compact sources when its
volume is smaller; `'off'` retains the global reference path. Both use the same
discrete global interaction kernel.

The three old `hse_mlwf_*` control names remain input aliases. See
[shared EXX controls and radius limitations](exx-mlwf.md).

The refresh count includes predictor and corrected states in RT, and is **not**
a count of physical time steps. Unchanged orbitals *and* occupations reuse the
cached exchange. Occupation changes invalidate the cache. `status=0` in
`HSE_WANNIER` diagnostics means the gradient criterion was met; `status=1` means
it was not. `status=2` denotes transport only (spread and gradient are -1, not
evaluated). A finite iteration limit does not guarantee an MLWF minimum.
At the default `exx_mlwf_radius=0`, unconverged localization retains full-support exact exchange with an explicit
message. The first refresh has no previous overlap (reported as zero).

For DC, provide a nonnegative electronic `temperature` or `temperature_k` and
enough `nstate_frag` to hold the occupied fragment space. Fractional occupations
and extra states are supported. States with positive occupation at any k define
the localization frame `Phi=Psi U`. The density factors used for exchange are
`Q=Psi sqrt(f/2) U`, **not** `Phi` with artificially equal occupations. Q need not
be orthonormal and is not itself an occupied-only MLWF basis. This preserves
`Psi (f/2) Psi†` exactly for unitary U. No disentanglement is implemented.

Requirements: orthogonal periodic cells, full uniformly weighted standard
rectangular k meshes, CPU, unpolarized fixed ions, `nproc_ob=1`, and
`nproc_rgrid=1,1,1` within each fragment. Fragments and their k points can run in
parallel. Existing DC restrictions on split-axis k points and buffer directions
still apply. Symmetry-reduced/custom k meshes, GPU, orbital/grid decomposition,
DFT+U, spin-orbit and microscopic vector potentials are unsupported.

Ordinary (non-DC) GS may explicitly select this backend, including fractional
occupations and extra states. Its RT path currently requires occupied-only fixed
occupations and default Taylor4+ACE. Polar temporal transport uses the previous
localized orbitals: `U=polar(Psi_current† Phi_previous dv)`, followed by occasional
spread minimization. Accepted link phases are continued through ±pi during the
line search. U/history are held in memory; they are reinitialized after restart,
which leaves the full-support exchange operator invariant. No new checkpoint
format or DC-to-global-RT projection is supplied.

DC energy uses core-weighted exchange: the full core expectation is removed from
the inferred ionic nonlocal term and half is added to XC once. LCFO applies the
full retained-fragment exchange to its basis, rather than extrapolating the ACE
operator outside its construction space. Buffer-size convergence remains a
required physical convergence study; small smoke tests do not establish it.

**Default performance boundary (`exx_mlwf_radius=0`):** this is the exact full-support Wannier
baseline. It transforms a fragment's full k mesh to its Born–von Karman
supercell and evaluates all source/target pairs without spatial or distance
cutoffs. Each fragment's k root performs full exchange; ACE applies locally on
all k owners. DC and ACE are active, but pair pruning, local Poisson boxes,
distributed Wannier-pair scheduling and large-system speedup are not yet
implemented or demonstrated. Source-only truncation from the reference scripts
is deliberately not used as a variational SCF approximation.


Exactly empty bands (zero occupation at every k) are omitted from the exchange
source. No positive occupation, however small, is discarded. The band union is
common to all k points. Changes to this set reset the localization gauge/history;
all retained states remain targets of the full action and ACE construction.
This is an exact density-rank reduction, not spatial pair screening. For a zero
density, one zero-weight source is retained to represent the zero operator.

For diagnostic export set the environment variable
`SALMON_HSE_WANNIER_SNAPSHOT=1`. A collective final refresh writes
`hse_wannier_snapshot.bin` in each fragment directory before LCFO processing.
See `samples/dc_hse/README.md` for the binary version, conversion, and provenance
limitations. The export records SCF and localization convergence separately.

## Spatial mesh RT from DC initialization

HSE now shares the PBEh Gamma spatial MLWF/FFTW/ACE implementation for
`theory='tddft_response'` and `theory='tddft_pulse'` with
`yn_dc='n'`, `yn_conventional_from_dcdft='y'`, `yn_hse_lcfo_rt='n'`.
DC initializes the mesh orbitals; time propagation acts directly on those
orbitals, without an LCFO projection. `yn_hse_wannier` is enabled automatically
for spatial layouts; use `yn_hse_wannier='y'` for the serial reference.

Use `nproc_ob=1`, `nproc_k=1`, `num_kgrid=1,1,1`, and y/z spatial layouts such
as `nproc_rgrid=1,2,1` or `1,2,2`. The x grid must divide by Py; y by both
Py and Pz; z by Pz. Cells must be orthogonal, the k point unshifted Gamma,
occupations fixed at two per orbital, and `exx_mlwf_radius=0` (full support).
The existing `exx_mlwf_interval/maxiter/tolerance` controls apply.

The HSE reciprocal kernel is `4*pi*(1-exp(-G^2/(4*omega^2)))/G^2`, with
zero mode `pi/omega^2`; the PBEh Coulomb cutoff does not affect this kernel.
The HSE mixing fraction remains 0.25. Distributed FFT and ACE store only
local grid rows. Supported fields are impulse and Acos2 without a second
field, using Taylor4+ACE. Checkpoints and snapshots are not supported here.

This migration covers fixed-ion RT and the conventional Gamma SCF scope below.
HSE MD remains disabled until its forces are validated. Legacy multi-k, finite-support and projected routes remain available under their
existing restrictions; they will be removed only after their replacements are
implemented and verified.

Validation: `testsuites/unit_hse_ace/validate_exchange.py` compares screened
and Coulomb exchange against serial Wannier on 1/2/4 ranks;
`testsuites/unit_pbeh_rvv10/test_hse_spatial.py` compares HSE DC-initialized
impulse and pulse histories against the serial route.

## Conventional spatial HSE SCF

For `xc='hse06'`, `theory='dft'`, `yn_dc='n'`, y/z spatial layouts now use
the same screened-exchange, MLWF and ACE implementation. Wannier activation
is automatic when the spatial process count exceeds one. The Gamma,
orthogonal-cell, FFT divisibility and full-support restrictions above apply.
Use occupied-only states (`2*nstate=nelec`), fixed occupations (omit electronic
temperature), fixed ions, and `yn_hse_lcfo_rt='n'`.

During this migration use `write_gs_restart_data='no'` and
`yn_self_checkpoint='n'`. Restart, checkpoints, Wannier snapshots and the
full-grid eigen/solver diagnostic exporters are not supported by this route.
Final SCF energy and eigenvalue text outputs remain available.
A four-process H4 example is `samples/hse_spatial/h4_scf.inp`; copy
`testsuites/pseudo/H_rps.dat` to its working directory.

H4 comparisons at density threshold 1e-10 gave a maximum energy difference
of 1.57e-9 eV and occupied eigenvalue difference of 2.09e-11 Ha against the
serial MLWF route. Iteration counts were 109, 49 and 476 on 1, 2 and 4 ranks,
respectively: this establishes final-state agreement, not parallel speedup or
identical SCF trajectories. Finite-temperature Gamma DC fragment SCF is supported as described below.

## Spatial HSE DC SCF with partial occupations

`yn_dc='y'`, `theory='dft'`, `xc='hse06'` can now use y/z spatial pencils
inside each fragment. The same unshifted Gamma, full-support and output
restrictions apply. FFT divisibility is checked on the **fragment** grid
`num_rgrid/num_fragment + 2*num_rgrid_buffer`, not the total grid.
`nproc_rgrid` describes each fragment; `nproc_rgrid_tot` describes the total
DC system. The total MPI count must accommodate both fragment and intra-fragment
parallelism. See `samples/hse_spatial/h4_dc_scf.inp` for two fragments with
two spatial ranks each.

Exchange factors are `Psi sqrt(f/2) U`, while temporal localization uses the
unweighted `Psi U` frame. This preserves the fractional density matrix;
zero occupations contribute no exchange, and an occupation-only update
invalidates cached exchange and ACE. DC core exchange includes spatial
reductions within each fragment. All positive-temperature HSE DC calculations,
including the legacy multi-k route, now use the same charge-converged Fermi
solver as PBEh. Insufficient weighted orbital capacity is explicitly rejected.
Gaussian centers are shared across hybrid spatial ranks, using the serial
initialization seed. Pure random initialization is unchanged.

Validation uses H4 at 10,000 K with six states per fragment on 1/2/4 spatial
ranks per fragment. Final printed total energies agree; core exchange differs
by at most 1.86e-11 Ha, and electron number errors are below 7e-14.
DC-LCFO output from every layout reconstructs mesh orbitals and runs fixed-ion
RT, with matching energy histories. The electronic SCF temperature is not
fitted or propagated during RT; the RT regression uses fixed occupations.

Convergence remains sensitive: the strict 1e-10 density test requires over
1,000 iterations in one layout. A 30,000 K stress case did not consistently
converge within 500 iterations, even on the legacy route; it is not certified
by these results. No speedup or general high-temperature convergence is claimed.
Legacy multi-k, finite-radius and projected implementations remain until their
spatial replacements are validated. HSE MD forces remain a separate task.

## Distribution of FFT and DC-LCFO projection work

Spatial exchange now keeps the forward FFT in Z pencils and applies the
reciprocal kernel there. The inverse returns directly to X pencils. Each
forward/inverse pair uses four pencil transposes instead of eight, and the
pair-density buffer is reused for the inverse result. This removes one local
complex FFT batch buffer without changing exchange normalization.

Complex DC-LCFO Hamiltonian projection no longer replicates two full-fragment
H-times-basis buffers on every rank. A single local-grid buffer replaces them;
each rank integrates its owned orbital columns and intersecting grid rows.
Only the small projected matrices are summed within the fragment. The
`DC_LCFO_HPSI local/global grid points` diagnostic exposes this layout.
For a fragment with N grid points, m states and s spins, these particular
buffers change from `32*N*m*s` bytes per rank to `16*Nlocal*m*s` bytes.
For the H4 test (N=1024, m=6, s=1), four spatial ranks reduce this workspace
from 192 KiB to 24 KiB per rank. This is not a claim about total process memory.

Remaining replication includes the representative rank's core basis, received
halo basis, dense LCFO diagonalization matrices, and all orbital columns in the spatial
MLWF/ACE representation. Halo data is still exchanged by fragment
representatives and broadcast within each fragment. Further orbital and
matrix distribution is needed for those parts; no fixed memory cap is imposed.

## Distributed core-basis construction

Complex DC-LCFO now restricts input orbitals to each rank's intersection with
the fragment core. Orbital and k-point peers combine only columns on the same
spatial rows. Overlap matrices, orthogonalization inner products and norms
are reduced across spatial peers; basis rotation acts on local rows.
Nonrepresentative ranks retain only their local core basis, including valid
zero-sized domains. `DC_LCFO_BASIS rank/local/stored grid points` reports this
layout once per fragment rank.

The fragment representative still holds the full core basis for the existing
binary writer and halo sender. Gathering now sends one actual spatial-domain orbital column at a time,
selecting one k/orbital representative per spatial domain. Senders use their
existing contiguous basis columns directly, with no extra send array. Root
receives into a buffer sized to the sending domain, then places it in the output
basis; its own domain is copied directly. Empty domains transfer no data.
Only small owner/extent metadata is collected over the fragment communicator.
`DC_LCFO_GATHER rank/send/receive grid points` reports the payload size and
maximum receive workspace (not total process memory).

The former full-core scratch columns on every process are gone. Root still
holds the complete output basis, so this is not a fully distributed I/O solution.
Received halo bases and global dense eigensolver matrices are unchanged.

Validation includes spatial HSE DC eigenvalues and reconstructed RT, a rotated
fragment geometry with ranks outside the core, and the existing four-k-point
Si LCFO references with k-point/orbital decomposition. No new memory cap or
input parameter is required.
