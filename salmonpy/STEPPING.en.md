# Protected steps with explicit backend callbacks

[日本語](STEPPING.md) · Topic: `protected_step` · Paired revision: `2026-11-17.1`

`salmonpy.stepping` lets a host protect one operation across explicit owners.
`StepTarget(key, owner, capture, rollback, release=None)` declares a complete
backend-owned backup contract. Tokens remain local and opaque: SALMONpy does not
copy, pickle, transmit or interpret them. `capture()` must be read-only.
`rollback(token)` must restore all numerical state, clocks, caches, failure flags
and side effects necessary for retry. `release(token)` only retires backup
resources. After partial capture failure, every rank that entered that target
still calls release; the failed rank supplies None. Therefore release must accept
None and clean any partial preparation idempotently. None is also a valid token. If release is omitted, capture must clean its own
failed preparation and tokens must need no explicit retirement.

```python
from salmonpy.stepping import StepTarget, StepRecoveryError, guarded_step

# Your backend explicitly promises complete rollback through these methods.
target = StepTarget("material", state, state.capture, state.rollback, state.release)
try:
    result = guarded_step(lambda: state.step(A), [target])
except StepRecoveryError as error:
    if not error.original_preserved or error.cleanup_errors:
        raise  # Do not retry an unsafe lifetime.
    # Accepted state is restored. Fix the cause before explicitly trying again.
```

For MPI pass `context=parallel_context` on every participating rank. All ranks in
a solver group must declare the same target key order and release availability.
Different groups may have different target counts. Callbacks may communicate
inside their solver group, never through the host leaders communicator. Both
the operation and callbacks must return quiescently on all ranks: this is not
recovery from deadlock, MPI/process termination, allocation exhaustion or a
native STOP. A declaration cannot certify hidden mutable state or external I/O.

Admission completes before capture; every capture slot is agreed before starting
the next. If the operation fails, all captured targets roll back in reverse local
order, including targets whose operations succeeded. Tokens are then released.
Errors are agreed through solver groups and their leaders without sending tokens
or results. A rollback error makes original_preserved false; all remaining
rollback/cleanup callbacks still run. A release error after successful execution
reports phase `release` and original_preserved false: the step may already have
been accepted. Keep these cases separate from a safely rejected operation.

## Mixed field/material stepping

`CoupledStepper.step_protected(targets)` requires exactly one explicit target for
each live local field or material object. Each target.owner must be that actual
object. Missing, duplicate or extra owners are rejected collectively before any
capture or field update. The coordinator additionally protects its own current,
clocks, response A/failure flag, endpoint adapters and measurement metadata.

This protects physical component state only to the extent guaranteed by the
registered callbacks. On incomplete rollback or cleanup the response is marked
failed, and the field is also invalidated where present. The ordinary `step()`
keeps its previous failure behavior. There is no automatic fallback, timestep
change or retry.

Run the installed example with six ranks and OMP/BLAS1:

```sh
mpiexec -n 6 python -m salmonpy.examples.mpi_protected_step NEW_OUTPUT
```

It joins one Multiscale1D rank, one LD rank, one SBE rank and a three-rank custom
group owning two independent points. Explicit callbacks protect the full wrapper
and its raw LD/SBE handle, not just A/current. First/warm/final late errors, a
failure after field update and an LD failure preserve every checked state;
80 accepted steps match ordinary execution exactly, without propagation replay.
Admission, partial capture, rollback and cleanup failure cases are included.
These are synthetic interface checks, not a physical-material prediction.

**Native E/H/TDDFT/fast-PR HSE whole-step protection is not provided by these
example callbacks.** Their individual numerical recovery/checkpoint protocols do
not cover accepted output streams or a whole mixed transaction. No native targets
are generated automatically; without explicit complete targets the call rejects
before stepping. Native aggregate transaction/output support remains separate.
Backup memory/time depends on the backend; budgets remain advisory. Importing
this module does not initialize MPI, load kernels, compile or download anything.
