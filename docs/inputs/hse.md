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
| `hse_mlwf_interval` | `10` | Positive number of exchange refreshes between spread minimizations. Polar transport/reconstruction is performed on every changed source. |
| `hse_mlwf_maxiter` | `200` | Positive maximum number of unitary spread-descent iterations per minimization. |
| `hse_mlwf_tolerance` | `1d-6` | Positive finite spread-gradient norm tolerance, always in atomic units (bohr squared), independent of `unit_system`. |

The refresh count includes predictor and corrected states in RT, and is **not**
a count of physical time steps. Unchanged orbitals *and* occupations reuse the
cached exchange. Occupation changes invalidate the cache. `status=0` in
`HSE_WANNIER` diagnostics means the gradient criterion was met; `status=1` means
it was not. `status=2` denotes transport only (spread and gradient are -1, not
evaluated). A finite iteration limit does not guarantee an MLWF minimum.
Unconverged localization retains full-support exact exchange with an explicit
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

**Current performance boundary:** this is the exact full-support Wannier
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
