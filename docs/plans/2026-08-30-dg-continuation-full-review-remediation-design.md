# DG Continuation Full Review Remediation Design

## Scope

Resolve every Critical/Important finding from the review of
`0f59d9b2..23b73c63` before Task 10 stationarity work. Preserve legacy routes,
the current worktree, and unrelated dirty changes.

## Architecture

The checkpoint and RT path remains distributed end to end. Physical grid
points and retained basis rows have explicit exactly-once catalogs. Rank-count
changes redistribute serialized row/grid shards directly; no rank constructs a
global dense basis matrix, position cube, density grid, or projected local
potential. Startup invariants use distributed sparse actions and scalar
reductions.

Hybrid RT uses a dedicated environment initializer that allocates grids,
pseudopotentials, Hartree/XC state, and field samples but does not allocate the
ordinary RT orbital triplet. The density-dependent update reuses the existing
DC distributed density-to-potential path and projects only owned rows. The
periodic position operator comes from the authoritative wrapped length-gauge
convention and carries a convention fingerprint.

Checkpoint provenance records the physical grid extent, position convention,
pseudopotential identity, and the existing `s_dft_energy` decomposition needed
by the Task 10 evaluator. Publication rejects duplicate, missing, and
out-of-range row/grid IDs before hashing.

## Error handling

Every MPI result that controls allocation, displacement, or a later collective
is converted into a communicator-wide failure boundary before its output is
used. Invalid catalogs, scopes, provenance, or coupling envelopes are rejected
collectively with no rank-local early return into a later collective.

## Verification

Tests first reproduce duplicate/missing grid points, injected collective
failure, periodic-boundary position coupling, row-shard rank redistribution,
absence of ordinary RT orbital allocation, and real production Hartree/XC
startup. The obsolete DC contract is updated to the analyze/closure/freeze
route. All focused 1/2/4/8-rank runners and the full build run before review.
