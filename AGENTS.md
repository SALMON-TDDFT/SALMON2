# SALMONpy repository

## Start with the current instructions

Read README.md for the repository overview, then salmonpy/AGENTS.md and
salmonpy/AI_USAGE.md before using or changing the public interfaces. For an
installed-wheel first calculation, use `python -m salmonpy.manual` and
`python -m salmonpy.quickstart`; no historical plan or solver build is required.

Current implementation scope is documented in the topic guides under salmonpy/.
For remaining development work, read the current-status section and latest
records in docs/plans/2026-11-15-state-isolation-checklist.md; do not count
historical unchecked items as separate current tasks.

The dated docs/plans/2026-11-14-handoff.md is historical context, not the setup
procedure, current feature list or authorization checklist. Its pending source
import approval and fixed-upstream-plus-patches workflow were superseded by the
project lead's 2026-11-15 approval. Its section 7 has a dated status cross-reference.
Do not execute obsolete commands or restore the superseded source-import ban.
This does not relax the maintained eigensolver audits, scalar IEEE inquiries,
source inventories, licensing or evidence-preservation rules.
The handoff's six design principles remain the baseline; use current API and
execution-mode contracts for their supported guarantees.
Follow the user's current instructions and applicable AGENTS.md rules; confirm
supported conditions in the current topic guide and the implementation.

## Repository and source ownership

This repository is for independent SALMONpy distribution. Do not merge it into
SALMON upstream. The project lead authorized bundled solver sources on 2026-11-15.

The single maintained solver source lives in runtime/sources/solver; completed
changes are integrated there. Frozen prior sources live in runtime/source-history. Runtime preparation verifies manifests and copies to a NEW
external directory; do not apply compatibility patches or modify fixed binaries.
HSE support is the fast PR profile only; comparison baselines are not supported engines.
Preserve shared-kernel byte equality between native/ and the maintained solver.
Retain original commits, file inventories, licensing and provenance.

Build kernels explicitly; never compile on import. Use new output directories
outside the repository. Run Python tests and the source-distribution-to-wheel
check after changing packaging or layout. Historical docs/plans and validation
records may reference the old full-SALMON layout; current instructions are in
README.md and runtime/README.md.

After changing bundled source layout/inventories, verify every manifest asset is
Git-tracked and prepare/build a Git-only export before considering the import
complete. Local ignored assets must not mask a broken clone. Preserve nested
upstream .gitignore bytes; force-stage only hash-verified manifest paths, never
an entire work directory. Heavy CI avoidance does not replace this check.

Source/distribution integrity CI policy (2026-11-20 review instruction): commits
touching runtime/sources/**, runtime/source-trees/**, packaging/** or pyproject.toml
must not use [skip ci]. Keep required inventory verification, Git-only export
preparation/build and GS-free C/Fortran ABI query checks enabled. Do not dispatch
additional heavy CI manually without need. Put new per-instance native state in
owned contexts/registry slots, not additional process-wide selector flags.

CI scope clarification (2026-11-20): automated CI covers source inventories,
tracked-only builds, packaging, synthetic kernel checks and host/API contracts.
Do not add real-material physical-observable reproduction, convergence campaigns,
measured optical tables or GS/TDDFT physics calculations to required CI. Keep
those as separately requested validation tools; any retained heavy workflow is
manual-only. This does not waive source/object eigensolver audits.


MPI communication policy (project lead, 2026-11-21): minimize application-wide
collectives. Use them for necessary setup, consistent failure decisions and
explicit final aggregation; do not add per-point/per-load global diagnostics.
Prefer owning solver groups and point-owner routing. Coalesce compatible status
and data agreements rather than removing finite checks, required physical
reductions or native rollback/retry guarantees. Final reports belong on the
reporting rank unless the caller explicitly needs replicated results. Avoid
barriers or extra MPI-health reductions when an existing operation proves liveness.
