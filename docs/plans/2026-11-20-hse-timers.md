# Fast-PR HSE handle timer ownership

Reuse s_timer_context per native handle slot and a temporary host bank. Scope
first initialization, shared/registered model creation and selected operations;
restore host on returned failures and close. Reset on a new handle generation.
No numerical algorithm, direct-checkpoint bytes or symmetry/history change.
Expose optional protocol1 completed rank-local counters on shared/registered
Python states; singleton diagnostics remain unsupported. Old libraries retain
normal operations. Keep main-thread sequential calls.

Extend existing tests/probes, reuse existing GS read-only through references,
collect no GS copies, run no SCF or CI. Compare full/reduced learned Taylor and
PT-CN current/energy/history to audited prior independent controls. Exercise
protected memory rollback and parked sibling clocks. Native arithmetic stays
Fortran. Verify source/object audits and source/wheel/document parity.

## Acceptance

Full learned Taylor shared and reduced learned PT-CN registered handles use
2ranks/group,3groups;32steps. Current, energy and8 history counters exactly
match old independent controls.64 parked sibling timer checks per case;
memory rollback, close/reopen and output-collision failure pass. Existing GS
read by references; no copies/SCF/disk electronic checkpoints. Source638:
637pass1skip, audits pass. Test harness now sets group nproc_k and references
sym.dat; collision status uses the existing rank-count result contract.

Final wheel638:637pass1skip; installed reduced PT-CN registered probe passes.
Source/installed Python, adapter/header, guide/AI JSON and probe bytes match.
Evidence: calc2026/salmonpy-hse-timers-record-20261120. No CI; uncommitted.
