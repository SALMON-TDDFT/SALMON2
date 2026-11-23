# Point snapshot identity

Continue approved M6 identity propagation to process-local PointState snapshots.
Wrap the backend-owned opaque snapshot in an immutable, versioned point envelope;
restore must match run, point, material definition and optional initial-state
identity before invoking the backend. Preserve the backend's model/settings and
failure-atomicity checks. Do not convert native state, infer capabilities, add
cloning/branch semantics, or introduce disk formats. Raw backend tokens remain
usable on raw states, but are rejected by PointState.restore. No solver update
order changes. Verify custom callbacks plus real LD/SBE independent continuations,
closed-state rejection, malformed/cross-identity tokens and matching rewind.
