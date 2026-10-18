# WSe2 HSE06 SOC NLCC preparation
Goal: prepare a separate orthogonal six-atom 1H WSe2 package, supported native HSE path, full uniform k grids, fixed nuclei, with official original fully relativistic NLCC VPS files. Build separately on Fugaku, preserve the old Si511 binary.
1. Verify official W/Se potential metadata and SHA; retain original files.
2. Record primary literature k meshes, method differences, cell folding and spacing.
3. Generate independent GS/RT cases for k, grid, vacuum, dt, field and time convergence. Reject unconverged GS before RT and pin a new binary hash.
4. Upload exact source snapshot, build separately on login node; only static binary inspection there. Prepare jobs, do not submit convergence production jobs in this preparation task.
5. Static preflight of all inputs/dependencies, geometry, electron/state counts, k weights/Gamma/K, scripts, source and binary provenance. Commit package and report precise remaining numerical checks.

User extended scope: implement SOC symmetry and validate GS/RT; polarized direction explicitly OFF. Bound implementation to signed-diagonal orthogonal operations with proper axial spin rotations; preserve legacy nonhybrid behavior. Full rectangular WSe2 group has eight operations including centering; x subgroup four, y subgroup eight. Separate gs-x avoids incompatible k-star restart.

Implementation commits a784add,3fbef68,7eca90d. Spec and independent quality reviews passed after gating density reconstruction to HSE06. A multi-k OpenMP ACE metric/application corruption was reproduced; sequential small ACE k loops retain BLAS teams and dominant Fock threading. Fifteen compiled/unit/MPI regressions pass. Final on/off GS+RT and new Fugaku build are recorded in the result ledger.

Preparation: 14 independent parameter cases, four stages each; official original W12/Se6 PBE19 VPS with finite NLCC and j-dependent channels, 48 occupied spinors. k spacing and literature evidence recorded; physical spectral convergence remains to be run, not asserted. Jobs are prepared only, not submitted. Existing Si511 executable SHA is protected before/after the independent build.
