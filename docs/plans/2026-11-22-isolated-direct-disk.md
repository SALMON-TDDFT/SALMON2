# Isolated molecular direct disk restart

Continue the approved owned numerical savepoint into explicit disk persistence.
Fortran writes/reads all active numerical arrays through the same bounded codec
used by the periodic context. The helper source is shared, not copied. Retain the
live numerical reply E3, add SAVE_E/LOAD_E control operations using a fixed private
filename in each worker directory. LOAD_E is allowed only in an untouched worker;
read into an independent payload and install only after the complete trailer and
shape/timestep/cell checks. Output writers remain closed until Python restores
recorded stream bytes and explicitly releases the imported backup.

Python provides save_direct_checkpoint/restore_direct_checkpoint. Store model and
runtime fingerprint, accepted index/observations and optional explicit replay
history metadata, binary numerical state, original output streams, per-file hashes
and sizes in a new-directory manifest. No overwrite and no format fallback.
Restoration validates before allocating a candidate worker, loads without STEP_E,
then swaps ownership only after successful stream reconnection. The existing
worker remains unchanged when the candidate fails. Imported states never silently
recover by replay after a propagation failure. Explicit replay archives remain a
separate API. A checkpoint is tied to the same runtime/input/dt/capacity; no
cross-compiler portability or isolated inner-MPI claim.

Validation: a failing synthetic public process-exit test first; then synthetic
corruption/limit/foreign-setting rejection and candidate failure; real C2H2 PZ
steps0/7/8/20 with actual original process exit, zero new update counter after
import, exact remaining observations and output bytes. Reuse authorized GS only.
Verify regular worker reference/independence, shared periodic codec, source-only
tracked export build, prohibition audits, source and installed-wheel suites.
Update paired guides, AI task JSON, example, current checklist and evidence.
No added real-material CI; existing inventory/build/packaging checks remain.
