# Call-local HSE checkpoint sizing scratch

Remove the module SAVE array-limit variable. Keep the existing formula local
to rt_checkpoint and pass it explicitly through every generated serializer.
Preserve bounds/product-overflow admission, scalar IEEE behavior and stream
bytes. No new input/API or changed capacity restriction; native concurrency
remains unsupported.

Extend existing synthetic I/O fixture/runner, compiling the actual extracted
serializer procedures with the solver types. Small/large/small scalar and vector
roundtrips, lower bounds and oversize rejection; old writer byte comparison.
No real GS or electronic checkpoint collection. Real HSE current/energy/history
regressions use existing GS references and memory rollback only. Native audits,
source/wheel and JP/EN/AI/list sync, no CI.

## Acceptance

Actual extracted serializers pass synthetic source/wheel roundtrips and shape
rejection; all3 old/new byte streams exact. All7 serializer bodies normalize
exactly after threading the local bound. Source/installed638:637pass1skip.
MPI2/group full-Taylor/reduced-PT-CN heterogeneous physical arrays/history exact
vs prior independent controls; memory rollback and reused-model state pass.
No real electronic disk checkpoint collected or claimed; no GS copies/SCF/CI.
JP/EN/AI/list updated. calc2026/salmonpy-hse-checkpoint-scope-record-20261120.
