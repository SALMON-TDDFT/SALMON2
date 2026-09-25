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
