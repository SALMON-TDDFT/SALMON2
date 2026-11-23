# Public isolated live-savepoint validation — 2026-11-22

The public isolated PZ fixed-ion worker exposes one opaque live savepoint per
state. Numerical copying stays in Fortran. Python restores output positions and
resumes the existing process; no replay is used for MolecularFeedback cancellation.
Explicit history-replay snapshots and archives remain separate APIs. A live token
cannot survive process death and is not a durable checkpoint.

## Observed results

Existing authorized C2H2 GS reused, no SCF. Outer MPI6, two spatial field ranks,
LD environment, dt=0.003 au, 64 field samples. Source continuous and rejected-step30
retry, and final installed-wheel retry, match every saved history array and the
complete six-face energy record of box-grid-x2-source-continuous exactly.
Worker PID stays unchanged: 63 accepted updates, plus one rejected attempted update
in retry runs. Restoration itself adds zero updates. Independent public direct
continuation also matches exactly. No new physical-convergence claim is made.

Source and installed-wheel suites each ran 690 tests with one optional skip.
Packaged API, paired guides, AI JSON and sample bytes match the working source.
The worker receipt is revision4 and requires same-process/no-extra-evolution proof.
The final wheel SHA256 is c7b6a22e6f68db9c44e011169c3427efd26aafae181e476500dd9923b2cf6d93.

Synthetic tests cover stale/foreign tokens, temporary journal failure followed by
explicit completion, and refusal to step while restoration is incomplete. A real
MPI worker-loss injection reports unusable and original_preserved=False without
spawning a replacement. Field samples remain preserved; caller MPI survives.
The first loss-injection development harness aborted deliberately on exceptions;
the corrected public-example harness supplies the caller-survival evidence.

## Records and limits

Canonical local OneDrive campaign: calc2026/salmonpy-c2h2-feedback-20261122.
python-savepoint-summary.json records comparisons, artifact hash, package checks
and worker-loss details. Runs: python-savepoint-source-continuous,
python-savepoint-source-retry, python-savepoint-verified-wheel-retry and
python-savepoint-worker-loss-example. Native numerical-state coverage and audited
executable provenance are in isolated-savepoint-validation-2026-11-22.md.

One worker rank, one pending live token, PZ and fixed ions only in this campaign.
No internal isolated-MPI or durable direct disk-restoration claim. Local storage
does not establish completed cloud synchronization. No additional physical CI run.
