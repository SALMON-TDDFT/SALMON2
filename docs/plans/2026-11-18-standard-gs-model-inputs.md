# Standard GS model-owned input entry

Goal: an internal sequential standard-GS entry that activates model-owned parsed
inputs before legacy main_dft, without changing default SALMON main or SCF arithmetic.
This is an incremental boundary, not retained GS wavefunction handles or concurrent
GS/TDDFT execution. Standard periodic unpolarized fixed-ion semilocal profile first.

1. Compile a red two-model load/run probe against the current audited runtime.
2. Add gs/gs_model.f90: owned parser/model values, load/run/clear, explicit status
   for empty/unsupported model, synchronous input view restoration, parked modulo
   and XC scratch. Preserve default entry and SCF warning policy.
3. Fresh source audit/build/object audit; MPI1/2 A→B→A, distinct grid/k/XC/cell/
   electron controls, compare individual old executable outputs. Do not collect GS
   bundles: no restart output, only selected diagnostic/eigen/energy evidence.
4. Regress TDDFT registry/default numerical comparisons, update manifests/checklist
   and guides/AI JSON with measured scope. No Actions or physical convergence claim.

## Completed measured scope

Internal GS model boundary, no retained orbitals. New model-only main_dft_model
retires request caches before local buffers; scoped children release afterwards.
Attached host MPI is required for this internal entry, while legacy main_dft
keeps its standalone path. Actual symmetry flags/all3 directions are swapped and
restored. Known unsupported load preserves its existing model.

Nine source and nine final-runtime/wheel cases pass: PZ/PBE, Libxc-PZ/TB-auto,
PBE/TB-fixed × MPI1/orbital2/spatial2. Different8³/12³ cells, atom/electron/k
controls; A→B→A status and independent old status/eigen/info/k bytes exact. Each
run leaves childlive0; six repeats keep warmed management slots stable. Warm
high-water is max(A,B), not only A's allocation. Preliminary failed count test
logs retained. A stale base-source manifest reference was synchronized; final
source and installed628 tests pass (one skip). Two source and one installed
TDDFT MPI6 cases keep currents/direct resume exact. Final legacy default GS
spatial2 also matches the old executable. Review reports no remaining blockers.

Four-iteration SCF diagnostics may be unconverged and continue with warnings.
No GS restart bundles written/collected, no material-convergence claim. No Actions.
Selected evidence at canonical calc2026/salmonpy-standard-gs-inputs-20261118;
local OneDrive writes verified, cloud synchronization unverified. GNU15 Apple
Silicon only; other compilers/platforms and retained/Python GS handles remain open.
