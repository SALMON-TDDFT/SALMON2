# Job-scoped runtime validation receipts

Goal: explicitly reuse one complete managed-bundle verification within a job.
Architecture: keep the existing runtime.py owner; add RuntimeValidation and
prepare_runtime_job, optional validation= on load_runtime/from_runtime and the
process-isolated TDDFT provider/backend. No implicit cache, environment activation,
network, build, numerical change or MPI initialization. Preserve default full path.

Contract: caller supplies an actual job ID and a fresh secret job_token (secrets).
One root validates everything; write a NEW external receipt bound to resolved
bundle root, complete manifest/hash and all file/directory stat fingerprints.
HMAC authenticates receipt fields using the separately distributed job token;
never store/log the token in receipt, engine report or calculation identities.
The token is ephemeral for one job, not a publisher signature or scientific proof.

Reuse verifies job/token, manifest SHA, file/directory fingerprints and ALWAYS
hashes executable bytes against recorded SHA256. No source/object startup audit.
Check metadata before/after use. This avoids hashing libraries/data and directory
walks, not all O(files) stat calls or the mandatory executable read. Treat bundle
as immutable/trusted during the job. OS metadata is not an adversarial attestation.
Cross-node stat identity/clock semantics require site validation. POSIX fast path;
other systems retain full validation rather than relying on weak change-time data.

Stale/changed/inapplicable receipt warns and falls back to full bundle validation;
valid calculations may continue. Corrupt actual runtime still fails its ordinary
integrity/capability checks. No hidden receipt rewriting or overwrite. Never trust
capabilities copied from an unauthenticated receipt; retain NLCC validation flow.
No arbitrary TTL: lifetime is explicit job ID and single-job secret. Do not reuse
that token across jobs. Missing/changed runtime files are never accepted from proof
alone. Input type errors remain normal API errors, not disabled scientific modules.

Tests: consolidated cases for fast I/O identity/secret redaction, job/signature/
metadata/manifest changes, same-size timestamp-restored tamper and extra/symlink
entries, safe fallback, no overwrite; NLCC scope parity. Measure hashes/read bytes
and metadata/walk counts rather than assume speedup. Actual caller-owned MPI4
GS-free synthetic receipt loading and a real prevalidated bundle where available.
Paired runtime guides, AI JSON, preview-first MPI sample, source and wheel suites.
No new GS, physics run, large-rank/site evidence or remote CI.

Implementation connected in existing runtime.py/fdtd.py/material_response.py;
source suite667 tests passed (one optional skip). Local MPI4 three loads/rank
retains caller MPI and uses receipt mode on all ranks. Auxiliary data/hash read
reduction measured on a DATA-ONLY18-file58MiB fixture; mandatory2MiB executable
stand-in hash remains. Full verification after stale receipt warns; actual byte
corruption still rejected. NLCC scope parity tested. Fresh wheel checks follow.
Receipt scope excludes embedded build/library loaders and automatic configured
catalog propagation; no new module required, no remote CI or physical GS data.

Final validation: source and installed wheel each667 tests (one optional skip,
no failures); installed MPI4 x3 loads all use receipt mode and host MPI survives.
Local warm-cache median per-load time about23ms full versus1.1ms receipt on a
data-only fixture; hash input60,817,396 bytes versus2,099,245 including manifest.
Instrumented Python stat/lstat/scandir counts:40/101/6 full,3/55/0 receipt. These
are local measurements, not FEFS/node scaling claims. No actual solver executed.
Canonical evidence and wheel/sdist retained in
calc2026/salmonpy-runtime-job-receipts-20261121/record.json.
