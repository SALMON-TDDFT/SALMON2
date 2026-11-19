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

## Computed GS→RT continuation

The existing retained-GS MPI tool now optionally starts RT from computed owners,
keeps A active while computing B and recycling GS slots, closes both source GS
owners, then continues their RT states. Six-rank orbital PZ/PBE and spatial
Libxc-PZ/PBE controls match separate execution exactly. Direct 8→16 restart
matches without replay. Temporary ordinary-RT checkpoint files are removed after
restore. No new GS restart/array dataset is collected.

The retained_gs example no longer requires import mode for --current-steps.
Its former rejection was reproduced before the sample change. Source and installed
MPI2 samples now run four SCF iterations then 16 current steps with three-component
A/c. This is ownership/numerical regression, not physical convergence.

A fresh HSE-enabled source/build passed the workspace source checker and all-object
audit. Learned full-k Taylor and reduced-k PT-CN MPI6 controls match the prior fixed
build exactly in currents, energy and history, including direct resume. HSE GS modes
remain import-only. Two earlier attempts reused generated execution bundles and
failed native initialization with both old and new binaries; those logs are retained.
The successful control uses the original fixed symmetry-history fixtures. The
specific generated-bundle failure mechanism was not established.

Source sdist→wheel→new environment installation succeeded. Installed GS→RT sample
and spatial MPI6 independent-state/direct-resume controls passed. Source and installed
Python suites each ran628 tests:627 passed,1 skipped. No new test methods were added.
Manual JSON synchronization and git diff whitespace checks passed.
No Actions or push were performed. OneDrive local storage does not certify sync.

Remaining: HSE-enabled computed GS, common ordinary/HSE registry, batch-report
completion and final removal of the old maintained source tree.
