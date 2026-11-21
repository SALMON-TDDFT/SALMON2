# R3 explicit MPI layout and mapper recovery

Reproduce prior validated runtime's STOP for unsupported Maxwell k/orbital factor; exit0 alone is not success (no completed validation record). Preserve old binaries/logs.

Return existing explicit layout/count errors1105 with the same check/broadcast. Add optional status to the actual init_communicator_dft error branches rather than whitelist names outside the mapper. Propagate unimplemented mapping1107 to cold-state/child cleanup. Default status-absent paths and implemented/Fujitsu/custom selection stay unchanged; no new global selector/collective/restriction. Extend the known cold-query protocol's accepted statuses only for guaranteed cleanup.

Extend the existing regression program, not separate per-case suites:27 combined checks per MPI1/2 layout including low-level k/orbital/count/mapping errors, Python failed constructors twice per category, prior local-grid checks, child resource release, caller MPI liveness and normal32-step direct resume. Compare against preceding audited runtime A/E/H, point files and checkpoint bytes. Explicitly exercise the alternate implemented orbital_sequential mapping. Source inventory/Git-only build/object audits and source→wheel verification required; no GS, heavy CI or arbitrary allocation/process/automatic-Fujitsu recovery claims.

Final verification: source and installed wheel679 tests, one optional skip. Installed SDK MPI2 repeats27 combined checks and direct resume; source/installed data and old control fields/output/checkpoint bytes match. Alternate implemented mapper passes. Guide/JSON sync and diff checks pass. Canonical records in calc2026/salmonpy-fdtd-mpi-layout-recovery-20261121; no manual CI/push/GS collection.
