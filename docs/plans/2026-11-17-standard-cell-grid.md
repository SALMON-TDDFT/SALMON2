# Standard cell/grid initialization ownership

Continue the approved GS/TDDFT initialization extraction. Fast-PR HSE already
has explicit cell resolution and grid configuration. Port these exact value
interfaces/types to the standard solver and its optional RT-step structure file;
do not duplicate work or change valid grid arithmetic/parity/epsilon rules.
Add the same capture-only cell_config_from_input helper to both solver profiles.
Compatibility initialization preserves downstream al normalization; explicit cell
and grid functions read no globals and reject invalid input without target change.
No claim of complete heterogeneous GS/TDDFT or symmetry/communicator isolation.

Use existing cell/grid frozen-reference probes on standard and fast-PR, add cell
capture then poison globals. Compile against old standard first (expected missing
symbols), then new modules. Refresh complete source/overlay manifests including
RT-step structure before/after pins. New external sources/libraries; preserve old
binary/core provenance. Compare MPI GS/mixed/protected retry/direct restore and
source/wheel tests. Update JP/EN/AI/checklist, no Actions; budget warning-only.
