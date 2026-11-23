# Owner-only field input routing

User requested minimum global communication and authorized the next transfer
optimization. Existing J collection already gathers only each owner's local
currents and returns them in sorted global point order on field ranks. Keep it.

Replace leader Bcast(full A) by Scatterv in plan-group order, with one received
row per owned point; only material subgroup ranks receive their local rows.
Field-group members own no material points and need no input broadcast. Packet
mode scatters per-group dictionaries then broadcasts only within nonempty material
groups. Singleton subgroup broadcasts are skipped. Root may appear anywhere in
plan order; point IDs need not be contiguous/interleaved ownership is supported.

Cache topology indexes/counts once at session admission. Preallocate/send-pack
within the existing preflight error guard; keep existing group/leader failure
agreement and callback order. No new world collective or per-point diagnostic.

Keep response._A's global-indexed shape for checkpoint compatibility. Field
leader retains full accepted A; material ranks update owned rows only; nonowned
rows are private placeholders/not replicated global history. Field nonleaders
need no drive payload. Point adapters keep their own accepted A unchanged in
meaning. Old full-array checkpoint metadata remains accepted; owned rows must
match local material A. Protected rollback must preserve array aliases/values.

Tests: manual MPI6 owner-routing contract with field rank pair in middle and
interleaved points, typed/custom packet paths, local buffers/keys, finite rejection,
independent control currents and caller MPI health. Existing MPI4 LD and MPI3
packet failure tests, mixed SBE/LD checkpoint and protected rollback. No new GS.
Compare retained pre-change behavior via earlier wheel where useful; field
arithmetic/native libraries unchanged. Source+installed wheel suites and MPI
contracts, paired guides/AI JSON/checklist. No remote CI or multinode claim.

## Completed verification

Owner-only A/packet routing passed the six-rank unequal-ownership contract, including a middle-position field group and nontrivial return permutation. Source and installed-wheel suites each ran 668 tests with one optional skip and no failures. Native FDTD with field-local LD and two SBE points preserved A, E and J bitwise against the previous implementation, including new checkpoints, legacy checkpoints and installed-wheel resume. Protected restore preserved original state after injected failure and retried without electronic replay.

The initial Scatterv counts tuple was rejected by mpi4py; using a list resolved the buffer-specification ambiguity. Bounded local diagnostics are retained with the validation records. No remote CI or new GS calculation was run. Hierarchical error agreements and native field-subgroup spatial J distribution remain unchanged.

Validation record: `calc2026/salmonpy-owner-routing-20261121/record.json` in the canonical results store.
