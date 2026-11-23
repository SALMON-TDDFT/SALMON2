# Ordinary heterogeneous HSE layer assembly

Preceding verified GS/registry/test-cleanup work committed as cd0841fe [skip ci].
This slice adds opt-in registry_model and ordinary ResponseSession group resources
for multiple HSE material IDs, plus assemble_native_layers model_groups admission.
Library/work root/symmetry permission and inner rank counts agree inside a group;
standard TDDFT/HSE native engines still use separate groups. Default plans and
same-model sharing are unchanged; custom backend contracts remain open.

The packaged mixed_hse_models example uses an offline material-ID input catalog
and layer order/thickness JSON, with field-local LD, SBE and two different HSE
models (full-k learned Taylor and identity-symmetry learned PT-CN), inner MPI2.
Grouped4 versus independent6 ranks produce bitwise-equal A/E/j for32 protected
steps on a40³ Maxwell grid and8³ electronic fixtures. Per-material native energy
output columns also agree exactly. Protected whole-field/LD/SBE/HSE save8 plus
24 further steps agrees bitwise with continuous32, without time replay.

Integration found fixed unit1001 in symmetry read_SymMat closed the first model's
live energy output during second-model initialization. A later model reused that
unit, so protected output capture correctly rejected the owners. NEWUNIT removes
the collision; registered initial output basenames are also unique. Existing
single-model basenames remain. Checks were not bypassed. The prior explicit
registry numerical tests did not establish this aggregate output-transaction path.
Use the newly built runtime for ordinary heterogeneous protected assembly.

An old-independent strict HSE MPI2 current/energy/history/direct-resume comparison
passes after the I/O fix. Source620 and installed620 tests pass (one skip each).
Installed-wheel mixed example matches source bitwise. JP/EN/AI guides, local input
and layout templates, active checklist and source inventories updated. Workspace
source/all linked application-object audits and wheel hashes pass. No CI, new SCF
or new GS calculation. This remains Apple Silicon/GNU/OpenMPI local interface
validation, not physical convergence or cross-platform certification.

Selected logs/JSON, development runtime/wheel and collision diagnostics are in
canonical calc2026/salmonpy-hse-layer-assembly-20261118. No GS restart/wavefunction
bundles were collected. Local OneDrive write verified; cloud sync unverified.
Remaining: retained-GS persistence/portable archives, broader native recovery and
remaining shared resources, larger-scale/platform/runtime distribution evidence.
