# DG Continuation Coding-Rules Compliance Design

## Scope

Bring only the Task 1--5 DG continuation additions into conformance with
`CODING_RULES.md`.  Do not change the DG mathematics, continuation state
machine, callback ordering, accepted-state provenance, or protected legacy
routes.

## Design

- Rename communicator and rank identifiers in the new public and private DG
  APIs to the SALMON `icomm_*`, `id_*`, and `nproc_*` conventions.
- Replace module-wide unrestricted MPI imports with explicit `only` imports.
  Retain direct MPI calls where the communication abstraction lacks the
  required integer-64 reduction or recoverable error contract.
- Enforce the 132-column free-form Fortran limit over all Task 1--5 source and
  fixture files.
- Preserve the approved focused `tests/dg` Python runners.  Moving these tests
  into `testsuites` would change the previously approved test workflow and is
  outside this behavior-preserving compliance pass.

## Verification

Use a static compliance runner as the failing-first test.  Then run all Task
1--5 focused MPI runners on their existing rank sets, compile a no-MPI probe
for every changed production module, and inspect the staged diff before the
single compliance commit.
