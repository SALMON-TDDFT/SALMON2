# Common protected step, first implementation

Continue the approved failure-atomicity objective with explicit backend-owned
capture/rollback/release callbacks. Tokens remain local/opaque; no serialization,
cloning, native-type guessing or hidden checkpoint replay. A public StepTarget
binds callbacks to one live owner. guarded_step executes one quiescent operation,
agrees errors across solver groups and leaders, rolls back every captured target
in reverse local order, then releases backups. Capture and release may be
collective inside the supplied group; callbacks must never use outer leaders.
Any incomplete rollback is reported explicitly and cannot revive a failed run.

CoupledStepper.step_protected requires exactly one explicit target per local field
or material owner, plus internal clock/current/adapter metadata backup. No native
backend is automatically certified: callers must cover its complete numerical and
external side effects. Normal step remains unchanged. First acceptance: real
Multiscale1D + independent LD/SBE handles plus custom collective model, uneven
point counts, late one-rank mutation/failure, strict original state and retry
agreement; capture/rollback/release failure negatives. Native aggregate output/
transaction support is a following extension, not claimed by this first entry.

Implement fail-first unit tests, common serial/MPI runner, coupled coverage and
metadata integration, real MPI example/tests, paired/AI docs, source/wheel tests
and persistent evidence. Warning-only memory budgets, no Actions.
