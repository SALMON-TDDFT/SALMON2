# HSE continuation validation

Version-controlled copy of the validated Nr16 continuation checker and its tests.
This is an offline validation template, not a job submitter. It does not modify
existing production inputs, executables, submission records, or locks.

The old post-run checker required at least `nt // 8` SVD records even when the
input used exact interval K=16. A valid checkpoint could therefore be rejected.
`check_endpoints` now checks the exact teacher schedule using the configured K
and the absolute saved start step. The expected SVD count is
`end // K - start // K`; the initial step-0 teacher is included only for a fresh
run. Missing/extra teachers and incorrect counts remain errors.

The accompanying Nr16 post-run template retains all-rank completion, placement,
memory, finite values, Ne, predictor flags, pulse phase, GS projection, and
checkpoint integrity checks. Its checkpoint schema assumes Nr16, 80 k owners,
16 occupied orbitals and K16/warmup160. Pulse parameters are the existing Si pump
parameters. It must not be copied blindly to other physical configurations.
The Nr20 continuation uses the same schedule helper with its own grid/schema
and its authorized projection-off checker. No SALMON propagation code changes
or new executable builds are part of this commit.

Run from this directory:

```sh
python3 test_continuation.py
python3 test_gate.py
```

The first suite checks schedule and pulse acceptance/rejection; the second uses
synthetic all-rank evidence and checks a valid restart plus seven rejection cases.
Actual production restart compatibility is not established by these fixtures.
