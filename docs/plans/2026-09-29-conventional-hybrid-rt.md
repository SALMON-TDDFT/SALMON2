# Conventional hybrid GS to native RT implementation plan

**Goal:** Start PBE0/PBEh40/PBEh40+rVV10 fixed-ion native real-space RT from conventional GS files, without DC.

**Architecture:** Reuse checkpoint_gs/restart_rt and existing real-space wavefunction formats. Add versioned GS hybrid metadata and validate before wavefunction loading. Rebuild MLWF/ACE on initialization. Support Gamma spatial pencils and full uniform k meshes with k-only MPI. Preserve the DC route and existing HSE routes.

**Tech Stack:** Fortran, native SALMON I/O and communication wrappers, FFTW, Libxc, MPI/ScaLAPACK, numbered CTest.

User approved this design on 2026-09-29. RT continuation/restart, moving nuclei in the new route, reduced k meshes, and simultaneous multi-k spatial decomposition are out of scope.

## Tasks

1. Add conventional GS metadata: functional/fraction/screening/effective Coulomb radius/rVV10 settings, mesh/cell/k vectors/weights/occupations and ionic configuration. Reject missing/malformed/mismatched new-route metadata before read_bin. Validate stored occupation payload after loading. Do not change existing wavefunction format or legacy HSE compatibility.
2. Extend input gates for non-DC fixed-ion response/pulse and Gamma GS wfn export. Continue to reject RT checkpoint/resume, fractional/unoccupied states in this route, unsupported decomposition, fixed-radius RT and moving ions. Remove DC-source assumptions from native exchange refresh, field/energy updates and endpoint field sampling only for the newly admitted global hybrid route.
3. Add standard numbered GS→RT tests and small regression tests. Cover PBE0/PBEh40/rVV10, Gamma MPI layout, full multi-k mesh, impulse/Acos2, zero field stability, missing/wrong metadata and occupations. Check computation completes, norm conservation and finite current/energy; compare MPI and serial results.
4. Use separate hybrid-rules-build (MPI/HSE/Libxc/ScaLAPACK) and HSE serial build. Run new tests and existing DC 422–430 regression. Independently review interfaces and input guards. Document supported inputs, metadata/compatibility and untested cases. Commit verified result without touching running spectra binaries.

## Expected input

GS: yn_dc='n', theory='dft', write_gs_restart_data='wfn'.
RT: theory='tddft_response' or 'tddft_pulse', yn_conventional_from_dcdft='n', yn_restart='n', directory_read_data points to GS output.
Both use the same hybrid functional and physical system. A new GS must be written to obtain the metadata. Standard generic old GS files without hybrid metadata are not certified for this new route.

## Verification (2026-09-29)

- MPI/HSE/Libxc/ScaLAPACK build and HSE-enabled serial build succeeded.
- New numbered cases 431–437: all 21 preparation/run/verification stages passed.
- Existing HSE/DC cases 422–430: all 27 stages passed.
- Conventional integration suite: 8 tests passed, covering all three functionals,
  full-k impulse/Acos2, Gamma zero field/Acos2, adaptive 99.9% source ACE,
  one-rank/two-rank current and energy agreement, electron norm, and metadata/
  occupation rejection before wavefunction loading.
- Serial executable: all three functionals, Gamma/two-k meshes, zero field/Acos2
  (12 RT runs) completed with finite observables and conserved electron number.
- Independent static review found no required changes. The updated manual
  section parsed with docutils. These tests do not establish spectral or
  k-point convergence.
- HSE-disabled serial build also succeeded. Four focused HSE input/admission
  regression tests passed. An obsolete test expecting static fixed-radius SCF
  rejection was updated to check its already-supported admission instead.
- The optional full HSE spatial Python suite was stopped during its 16-rank
  fractional-occupation case to avoid competing with ongoing spectra runs;
  this is not counted as a completed suite. The numbered DC/HSE suite above
  completed, and the relevant input guards were rerun separately.
