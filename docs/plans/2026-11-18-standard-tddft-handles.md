# Same-model standard TDDFT native handles

Continue approved native state isolation. First standard profile: one immutable
non-symmetric periodic scalar semilocal model, common dt/capacity/decomposition,
independent wavefunctions/clock/three-component A/current, transactions, output
journals and pending direct restore. Legacy singleton ABI and numerical order stay.
Native token ABI first; Python shared-model factory and different input/dt scopes
remain subsequent work, not implicit claims.

Move accepted A, step_saved/XC/cval/output metadata, pending/cval and serializer
bound into rt_context. Use a pointer selector only during synchronous calls and
an opaque generation registry owning pointers (never copy owning contexts).
Read input once; fresh states call existing initialization from original GS with
unique construction-only output labels and per-state MPI ownership. Model outer
timers close only with the last state, permit idle recreation. All subgroup ranks
call same sequence on main thread, attachment lives until states close. Native
checkpoint encoding unchanged, public handle boundary returns charge current.

TDD: missing-ABI red; shared states vs old singleton controls, x/y/circular drive,
MPI orbital/rgrid2 and serial1, all supported XC. Fresh direct restore equality,
simultaneous transaction/pending and close independence, stale token and MPI
argument mismatch, bounded MPI children. Test source/installed/package; paired
JP/EN/AI and active queue. Explicit new external source/build only, no Actions.
