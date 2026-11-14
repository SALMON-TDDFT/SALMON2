# Proposal: field × material × execution

Status: approved design; first MPI Multiscale1D milestone in implementation. Existing numerical updates and
`salmonpy.time/1` declarations remain authoritative. No upstream source import or
merge is part of this work.

## First deliverable

Run `Multiscale1D` with distributed TDDFT: one field-role rank holds the Python
field state, each electronic solver group owns one material point and calls the
existing internally parallel native TDDFT backend. Preserve independent states,
point ordering, charge-current sign, initial feedback and time placement. Compare
against serial `maxwell.couple` using the same prescribed initial data.

## Boundaries

Separate three responsibilities instead of forcing MPI factories and a single
physical state into one class:

1. Material preparation creates reusable initial data; it does not propagate.
   Existing providers remain responsible for manifests and ownership.
2. Physical material state accepts A/c and produces charge-current density j.
   Public state contract: `step(A_new)`, `current`, `time`, `timing`, `close`.
   Old A is owned by the state/adapter, not supplied independently by the caller.
   Native J_matter is converted exactly once at its declared boundary.
3. Execution owns mapping and scheduling: serial or MPI. It creates independent
   point states, distributes A, reduces j as declared and returns point-ordered
   arrays. It never changes the field's discretization or sign convention.

Providers/factories still receive solver_comm and the point identity where
needed. All members of an MPI material group call collective native operations.
The factory is a lifecycle boundary, not a second physics protocol.

## Proposed adapters and migration

- Wrap existing crystal states directly; use their current, time and timing.
- Wrap legacy `advance(A_old, A_new)` states behind step, retaining accepted A in
  the wrapper. Commit that history only on successful advance. Do not negate j.
- Refactor ResponseSession into an MPI execution adapter without removing its
  existing entry point. New coupling code must call public adapter methods,
  not `_agree` or inspect material states directly.
- Provide a serial execution adapter with the same point-batch contract.
- Wrap uniform/Multiscale1D field states with a point map and field-execution
  adapter. For the first milestone permit exactly one field rank; reject other
  allocations instead of redundantly advancing a shared field.
- Wrap NativeEHField as a collective field implementation. Keep its spatial
  deposition, LD coefficients/Q ordering and ownership unchanged.
- Existing `maxwell.couple` and `CoupledStepper` become thin entry points to a
  shared schedule only after parity tests pass. ResponseSession remains usable
  independently for prescribed-A workloads.

No generic multiplication of execution choices is promised: Python 1D field
stepping is initially serial, native E/H is collective, and electronic native
states retain one lifetime per process. Reject unsupported combinations before
state advancement.

## Time and initialization

The schedule is j(previous) -> field step -> A(new) -> material step -> j(new).
Check every field/material port with TimeSpec before running. Native E/H keeps
its explicit +1/2 electronic clock adapter. Multiscale1D uses endpoint clock 0;
never apply the E/H half step universally.

Keep the material's physical initial current as an observable, but do not feed it
back where the field declares uses_t0=False. In particular both SALMON-based
field paths begin with zero feedback. No calculation of initial current may
advance a state. Declare whether current is physical or zero explicitly.

Distribute field A(0) before creating material states. The current native MPI
backend accepts zero initial A only; first milestone rejects nonzero material-site
A(0) before native setup. Report the point ID and value and suggest increasing
the left vacuum cells or delaying the pulse. Never silently zero A(0). Pulses initially supported lie in vacuum so material A(0)
is zero. Supporting initially illuminated matter requires an explicit extension
and comparison with the native initialization, not silently replacing A(0).

## Failure guarantees

Local input validation must precede mutation and collective native calls.
Successful point steps do not make an MPI group transaction atomic: a later
point/group can fail after another has evolved. The existing native field and
TDDFT do not provide rollback. A wrapper cannot repair that by retaining old A.

For the first migration, preflight errors leave states unchanged; any failure
following evolution invalidates the entire coupled session and requires restart
from initial data. Report this restriction explicitly. The project goal of
atomic failed steps requires checkpoint/restore or prepare/commit operations in
*both* field and material components, plus a collective commit decision. Do not
advertise that guarantee before implementing and testing those capabilities.

## Acceptance tests and order

1. Same manufactured material, same field: serial execution and MPI execution
   give equal point-ordered histories. Two heterogeneous points, nonconsecutive
   point IDs, multiple material ranks, wrong shape/time/sign declarations.
2. Multiscale1D + distributed TDDFT: zero feedback, two independent points from
   the same GS, comparison to serial workers and saved-A replay. Test rank-count
   parity, initialization rejection and failure invalidation. Extend
   `tools/test_longitudinal_salmon.py` ms1d to at least two points and compare
   each point's j plus the entire field A grid against SALMON itself. Align
   native saved feedback j(n-1) with endpoint histories explicitly.
3. E/H migration: preserve all existing 13 coupling regressions, both clocks,
   LD coexistence and all-axis MPI-boundary deposition.
4. Keep physical convergence tests separate from interface parity. First add
   source-level/small fixtures to normal CI; converged material runs have pinned
   input/potential provenance and a separately declared execution budget.

Only then retire duplicate scheduling code. Persistent Fortran handles, GS v2
and copy reduction are separate milestones; do not enlarge this first migration
into a simultaneous runtime rewrite.

## Communication measurements

Measure A distribution and j collection separately from material execution,
recording per-rank elapsed times, call counts and logical payload bytes. MPI
call elapsed times may include waiting; do not claim isolated network latency.
Collection/reduction scope must be named explicitly in the report. Measurement
errors mark metrics unavailable without turning a successful calculation into
a failed calculation. First acceptance records one-host results for two material
points and multiple inner rank counts; multi-node scaling remains unverified.
