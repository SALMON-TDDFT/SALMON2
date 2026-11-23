# HSE model-owned initialization compatibility guard

Keep the existing fixed-settings predicate and numerical archive representation.
Move its pinned initialization copy/latch from the process selector into an
explicit s_rt_model_settings owner. Registered models own one each; the shared
same-model factory owns one; standalone/legacy calls retain the old fallback.
Select_model restores actual inputs but no longer overwrites a process-wide
guard. Only successful initialization commits the model guard. Failed opens
leave it unchanged. No C ABI, input restriction or checkpoint format change.

Verify alternating full Taylor/reduced PT-CN models, independent prior controls,
fresh state after model switching and rollback; existing GS references only.
Source/object audits, source/wheel and human JP/EN/AI/list sync; no CI.

Model generation reuse additionally clears only that slot's guard. Test reuse
of an A slot for material B with B still alive, then compare its first current,
energy and history to B's independent control. Failed evidence is preserved.

## Acceptance

MPI2/group heterogeneous full-Taylor/reduced-PT-CN current/energy/history exact
against prior independent controls; rollback, switched-model fresh creation and
slot reuse for a different material preserve live sibling counters and physics.
Intermediate missing-reset case fails; corrected runtime and final installed
heterogeneous/shared probes pass. Source/installed638:637pass1skip. Native
source/object audits and fixed-settings/archive-source parity verified. No GS
copies/SCF/CI. JP/EN/AI/list updated. Evidence in
calc2026/salmonpy-hse-model-guard-record-20261120; uncommitted.
