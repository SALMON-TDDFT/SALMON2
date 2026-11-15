# Point-bound TDDFT disk checkpoint

Approved M6 continuation. Reuse existing TDDFT replay checkpoint, add strict outer
point identity envelope. PointState.save_checkpoint/restore_checkpoint support only
known TDDFTState (optionally wrapped by ConfiguredState). No custom deserialization,
LD/SBE new disk format or MPI whole-session resume.

Outer directory contains point.json and state/ existing native replay archive.
Envelope records schema/version, exact run/point/material/initial identity and
SHA256 of inner manifest. Same run ID for explicit resume; cross-run/point copies
rejected (branch API deferred). No model loaded from disk. User supplies validated
initial state/runtime; existing native loader checks data/settings and replays
atomically. Strict bounded JSON, duplicate-key rejection, symlink/extra-member
checks before backend restore. max_bytes retains inner loader semantics.

Tests first: mismatch rejected before restore, no overwrite, corrupt manifest,
malformed metadata, cleanup of failed save and unsupported backend. Native installed
Si restart/replay parity, independent point rejection, full tests/sdist->wheel.
Paired docs and AI JSON. Max two workers/one thread each, local only; no Actions,
commit/push. Memory peak not claimed.
