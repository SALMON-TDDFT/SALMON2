# Minimize outer/global communication

User instruction: keep whole-application communication minimal. No new world
collectives in timestep loops for diagnostics/receipts. Prefer owner/subgroup
communication and route A/j only to the point owners. Keep required physical
reductions, finite checks, consistent error decisions and retry guarantees.

First bounded changes: verification-only MPI example uses one setup broadcast
and one root-only final gather (no replicated allgather or health allreduce).
Per-worker receipt loads remain independent/no MPI calls.

NativeEHField.advance wrapper currently uses six subgroup collectives: input
admission, J equality, agreement on equality, native status/E finite, agreement
on status/E, agreement on integrated A. Fuse these into two ordered agreements:
(error,J) before native entry, then (code,E-finite,A-finite,integration-error)
after entry. All ranks see identical records and make identical decisions.
Numerical E/A arithmetic, native update order, recoverable status protocol and
commit point remain unchanged. Native internal physical/subgroup communication
is not removed. Preserve signed-zero equality and all ranks' failure behavior.

Tests: one consolidated synthetic wrapper test for two calls, mismatch before
entry, remote bad E, recoverable20 retry and overflow rejection; actual MPI2
native baseline/restore/overflow retries compare against existing retained
controls. Verification-only MPI4 checks root gather and error exit. Source/wheel,
paired guides/AI JSON, code ownership rules and TODO updates. No new GS or CI.

Remaining: targeted A/field-packet routing and reductions of redundant outer
agreements require separate design/validation with mixed-state rollback; do not
claim they are optimized by this bounded change.

Completed verification: source and fresh wheel668 tests each, one optional skip.
MPI4 communicator proxy measures exactly one setup bcast and one root gather
per rank in the verification example (measurement-only final gather excluded).
Native MPI2 field/LD400 steps: original control, new continuous run and direct
restored runs from source and wheel have bitwise identical A/E/H. Corrupt profile,
late checkpoint/model mismatch and overflowing trial reject/retry tests pass.
Known native rejection does not perform the Python A arithmetic before return.
Original numerical order and all finite checks retained; required internal native
collectives unchanged. No new GS, remote CI or multi-node/HSE/GPU claim.
Canonical record:calc2026/salmonpy-minimal-collectives-20261121/record.json.
