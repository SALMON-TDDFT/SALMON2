# Standard singleton native timer ownership

Keep init/step/close in a legacy-state timer bank, preserving host counters and
running latch across returned paths. Shared handles keep their existing bank
and bypass this wrapper; numerical initialization/propagation/close bodies are
unchanged. Add optional legacy timer query; old runtime calculations remain
valid and only the new query requires its capability. No numerical checkpoint
fields, no concurrency/reopen guarantee. HSE singleton remains separate.

Use existing timer/native IO fixtures to seed host counters and compare before
initialization/after every step/after close. MPI independent controls, returned
input rejection/retry and host output survival; GS references only, no copies
or disk electronic checkpoints. Native audits, source/wheel, JP/EN/AI/list;CI0.

## Acceptance

Source PBE/rgrid2 and installed TB-mBJ/orbital2,3groups of2ranks: singleton
currents exact vs prior independent controls,36 host counter checks per case,
stale-A rejection/retry and close/detach preserve host. Shared TB-mBJ controls
and rollback still pass. Corrected old-runtime harness reproduces counter
mutation; earlier macOS helper loader failures were unrelated harness defects.
Source/installed638:637pass1skip. Audits and JP/EN/AI/list/installed parity pass.
No GS copies/SCF/CI. calc2026/salmonpy-standard-legacy-timer-record-20261120.
