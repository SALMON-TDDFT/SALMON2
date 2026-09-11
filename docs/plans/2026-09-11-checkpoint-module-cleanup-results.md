# Hybrid checkpoint responsibility cleanup

The ambiguous `src/rt/dg/rt_dg_hybrid_checkpoint.f90` is removed.

- `rt_dg_hybrid_checkpoint_v5.f90` now owns the v5 publication authorization
  type, collective validity/mapping checks, guarded publication endpoint, and
  existing v5 manifest/shard reader and writer.
- `rt_dg_hybrid_occupied_checkpoint.f90` owns only the auxiliary occupied-state
  stream and its I/O helpers. Its module name matches its filename.

The auxiliary stream is not consumed by production RT startup; that path reads
v5. It is still written by GS and read by production-evidence Python validators
and MPI compatibility tests, so it is retained rather than silently disabling
those checks. Eliminating this extra output would require migrating those
validators to v5 and is not part of this behavior-preserving rename/move.

All three moved v5 publication routine bodies were compared with their original
text and are unchanged. Occupied writer/reader/helper bodies are also unchanged.
The schema versions, magic strings, filenames of saved calculation data,
serialized layouts, and exact rank/fragment reuse policy are not changed.
Existing calculation files were not removed. Source is recoverable from the
preceding Git commit `7479ba8a`.

## Verification

The responsibility-separation test failed before the move because the ambiguous
file still existed, then passed. Updated imports, CMake source registration and
test compile lists refer to the new owners. The build completed with exit 0.

Passed on MPI 1/2/4/8 ranks:

- publication preconditions and occupied checkpoint round trips;
- v5 shards/atomic manifest, corruption and failure handling;
- v5 RT initializer and collective rejection of dense v3 input.

Static distributed-v5 architecture, continuation-alias publication,
self-consistent GS, terminal refinement, obsolete-route inventory, and
`git diff --check` also passed.

Logs are preserved under `/tmp/si8-route-cleanup-20260911/`:
`build-checkpoint-separation.log`, `occupied-separated-tests.log`,
`v5-separated-tests.log`, `initialization-separated-tests.log`.

No new physical GS calculation was necessary for this move. This checkpoint
does not certify the outstanding Si64 accuracy investigation.
