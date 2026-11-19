# Computed ordinary GS on fast source — 2026-11-19

HSE-disabled native GS now supports SCF, convergence metadata, energy, eigenvalues,
forces, and the shared explicit export assembly path. HSE-enabled compute is not advertised.

Six host ranks form three two-rank solver groups. Orbital PZ/PBE and spatial
Libxc-PZ/PBE tests compare retained A+B with independent A/B; all five numerical
arrays agree exactly. Closing/reusing A preserves B; communicator resources release.
Unconverged SCF warns and retains its last iterate; the large-threshold B control
reports native convergence. Export arrays were inspected in RAM, without writing
GS archives, wavefunctions or restart datasets.

The prior fixed fast executable was reused after SHA256 verification. PZ and PBE
gs_info.data and gs_eigen.data are byte-identical to new native results (files were
compared locally; numerical GS datasets are not copied into this evidence store).
Old standard PBE differs by 3.37786711e-7 Hartree after four unconverged iterations;
that comparison failed and its log is retained. Its tolerance was not relaxed.
This is a regression/ownership test, not physical convergence certification.

Pending: fresh HSE-enabled build regression, computed-GS→RT lifecycle and installed
wheel verification for this slice, HSE computed GS and final source consolidation.
No Actions were run. Local OneDrive storage does not confirm cloud synchronization.

Source Python suite: 628 tests, 627 passed and 1 skipped. Manual JSON guide synchronization check passed. Full git diff whitespace inspection hit a OneDrive mmap timeout; no clean-diff claim is made.
