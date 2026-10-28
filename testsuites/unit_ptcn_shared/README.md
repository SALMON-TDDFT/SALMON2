# Shared native PT-CN regression tests

`run.py` accepts an existing MPI + ScaLAPACK build and explicit platform link flags.
It audits the source, all SALMON application objects, and new probe objects before execution.

- Fresh exchange: stale-inner convergence is rejected; acceptance uses a newly built physical operator.
- Rectangular FFT: nonorthogonal finite-difference kinetic inverse, plus scalar NaN/Infinity rejection.
- Distributed spinor: unequal occupations, analytic CN reference, below-tolerance plateau acceptance and above-tolerance rejection, on one and two MPI ranks.

These tests do not certify large time steps or HSE SOC/thermal input configurations.

- Learned inner predictor fixture: correct prediction reduces three fresh actions
  to two; deliberately wrong prediction requires exact correction and agrees with
  the ordinary result. See `learned_exchange_probe.f90`.

- Lowdin: serial multi-k orthogonality, scalar NaN/Infinity/invalid-metric rejection,
  and orbital-distributed1/2/4MPI comparison with a sequential eigensolver reference.
- `check_orthogonalization.py LOG` asserts the HSE hook, pre-root gate and restored
  post-Gram accuracy. The post-retraction CN residual is independently reported,
  and this integration test does not silently substitute it for the root residual.
