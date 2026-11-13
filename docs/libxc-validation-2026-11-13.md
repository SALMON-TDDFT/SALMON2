# Local Mac Libxc validation — 2026-11-13

The process-isolated TDDFT worker now permits `libxc_pz` and `libxc_pbe`
when built with `USE_LIBXC=ON`, in addition to the existing `pz`.
Only its input allowlist changed. GS, midpoint propagation, predictor-corrector,
nonlocal phases and current evaluation retain their existing implementation.
Other functionals remain rejected. The Python API is unchanged and A-driven.

## Build

Apple Silicon Mac, GNU Fortran, OpenBLAS, two OpenMP threads and one BLAS thread.
Libxc 5.2.3 was built in `/tmp/salmonpy-libxc-install` from its official archive:
https://gitlab.com/libxc/libxc/-/archive/5.2.3/libxc-5.2.3.tar.bz2
The installed Homebrew Libxc provides f03 but not the legacy f90 interface needed
by this SALMON implementation. `salmonpy/PACKAGING.md` records the compatible build.

Source and all application object files passed the forbidden-eigensolver audit.

| Executable | SHA256 |
|---|---|
| Libxc, before worker allowlist change | `de85bf8f275134a84a4764a4380b50d78b88a17f72b49d42bcc3d2dfa1bd0b89` |
| Libxc, after worker allowlist change | `3f6a01556838f4d7ac4456a1a1d849712f2a71131b1f31c2210ccd11cfe5aacf` |
| Libxc disabled, after change | `1b59b1000d56c7c3399ea96f8ee14cb5f5d1f2cd70dd2484a7739213f3a7061c` |

## Ground state and real-time tests

Eight-atom Si cell, lattice 10.26 a.u., 2×2×2 k points, 32 electrons/32 states.
For each case, a new GS was converged to a density residual below 1e-9 and its
restart used with the same functional and pseudopotential in ordinary TDDFT and
SALMONpy. RT: 500 steps, dt=0.05 a.u.; z-polarized A/c pulse amplitude 0.001,
frequency 0.2 and duration 20 a.u.

| XC | Grid | GS iterations | GS residual | Max ordinary/worker current difference (a.u.) |
|---|---|---:|---:|---:|
| PZ | 12³ | 37 | 6.0476828e-10 | 5.421010862427522e-20 |
| Libxc-PZ | 12³ | 37 | 6.0851142e-10 | 5.421010862427522e-20 |
| Libxc-PBE | 16³ | 60 | 8.2973377e-10 | 5.421010862427522e-20 |

Each case also passed 64 interleaved steps of two independent states driven with
opposite A, each compared with its separately computed trajectory.
Built-in PZ versus Libxc-PZ differed by 2.3032589524518698e-7 eV in converged total
energy and 2.4649698921559368e-14 a.u. in current. These are numerical comparisons,
not bitwise equivalence claims between the two XC implementations.

PZ and Libxc-PZ use the same repository `Si_rps.dat`, SHA256
`8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c`.
Its KY header does not establish the generating XC provenance; it is used here
for a controlled implementation comparison, not a certification of that provenance.

PBE uses the official PseudoDojo Si UPF:
https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard/Si.upf.gz
Uncompressed SHA256:
`39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`.
The inspected UPF 2.0.1 header specifies PBE, NC, scalar relativity, 4 valence
electrons, NLCC enabled, no SOC, no PAW and no ultrasoft augmentation.
The PBE run reported NLCC active. PBE and PZ use different potentials and grids;
their energies/currents must not be interpreted as a controlled XC-only comparison.

## Regression and rejection checks

Before changing the allowlist, both Libxc cases failed at the explicit PZ-only
worker guard. Afterwards `tools/test_tddft_entrypoint.py` passed for both XCs
with predictor-corrector off and on: ordinary pre/post-change current and energy
files were byte-identical; live current, independent states, invalid-A rejection
and killed-worker recovery passed. Those short synthetic-orbital tests check the
interface; the converged-GS runs above provide a separate physical initialization.
The regression helper now also rejects the unsupported `libxc_pw` worker input.
Libxc-disabled builds explicitly reject both Libxc choices before initialization;
the existing PZ integration still passes. All 34 Python tests passed.

## Limits and artifacts

This is local Mac validation, not new Windows/Linux CI evidence. No grid/k-point/
time convergence study, other materials, spin polarization, SOC, HSE, generic
Libxc functional support or PAW/USPP validation is claimed. The loader snapshots
and hashes input files but does not certify GS provenance or automatically detect
GS/RT functional and pseudopotential mismatches.

Inputs, GS/RT outputs, pseudopotential, summaries, scripts, logs and executable
copies are archived outside the code tree at:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/salmonpy-libxc-validation-20261113`.

Environment note: absolute paths in this report identify the historical validation host. They are not portable setup defaults or required user paths.
