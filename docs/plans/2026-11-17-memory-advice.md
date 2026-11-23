# Advisory memory budgets

The user replaced budget admission stops with warnings. Budget excess is advisory:
never cancel computation/restoration merely because an estimate exceeds a budget.

Expose immutable MemoryEstimate/MemoryReport and memory_advice for user-owned
backends, plus a lazy standard-library process-memory probe. No MPI initialization,
allocation probes, downloads or compilation on import. Aggregate supplied rank
estimates over a host communicator; budget means total across those ranks, not a
per-rank or scheduler-enforced limit. Explicitly distinguish unavailable estimates
from a numerical zero. Reporting failures do not invalidate scientific states.

Ordinary mixed sample: optional --memory-budget-gib. Before initialization, report
heuristic field/model/input-size estimates and initialization unknowns. After
initialization, report observed process peak RSS plus reserve. Before protected
restore, include current rank peak RSS, candidate coexistence/read/copy allowance
(four times required local payload sizes), decoded Python metadata and per-rank
reserve. Native payload sizes are an estimation basis, not a proven heap bound.
Keep stage/commit arithmetic and update order unchanged. Save serializable advice
for inspection; warn on root only, allow computation/retry even with tiny budget.

Verify report arithmetic, partial unknowns, unavailable probes, warning delivery
errors, collective budget disagreement, over-budget warning without stopping,
real four-rank mixed restore including injected late failure and exact retry, and
fresh CLI resume. Use at most four ranks/OMP1; verify the 4GiB report without
claiming an OS-enforced RSS cap. Sync JP/EN/AI and source-to-wheel checks.
