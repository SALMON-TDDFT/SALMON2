# Mixed checkpoint and output identity

Reuse existing point-run/1 and ResponseSession point records/maps. The ordinary
fixture declares a run ID (default mixed-example, user override --run-id), stable
point labels from layers or transport IDs, material definition fingerprints and
initial-state references. LD materials retain user IDs and constitutive-grid IDs.
Output paths/work directories and advisory memory budgets are not identities.

Identity-bearing native-direct/3 manifests embed the agreed point/ownership and
constitutive records. Compare every identity before numerical LOAD; also compare
identity digests across ranks before saving. Generic callers with legacy sessions
retain native-direct/2. Ordinary sample default is identified; explicit
--legacy-checkpoint-identity disables declarations for old unidentified archives,
without disabling protected restoration. It cannot bypass identities on version3.

Emit the same identity in normal results and identity-aware output archives. Keep
native payload encodings, update order and memory warning policy unchanged. Tests:
run/point/material/definition/initial/LD/ownership mismatches before any stage,
real MPI identity mismatch and preserved originals, retry equality, source/wheel
and fresh CLI resume with matching IDs; older archive mode explicit. Four cores,
advisory4GiB, no remote Actions. Broader archive migration and arbitrary backend
factories are not implied.
