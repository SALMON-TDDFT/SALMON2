# Public isolated E-worker entry

User direction: stop extra convergence studies; implement remaining component APIs.
Implement load_isolated_runtime(binary, validation) with a binary-bound bounded
local validation receipt. No download/build/import side effects or silent A-worker
fallback. Existing releases without E1 evidence are rejected. This is explicit
local runtime loading, not a newly published catalog asset.

IsolatedTDDFTInitialState snapshots native GS bundle using existing ownership/
lease machinery. new_state(dt,position) creates independent E1 worker. Initial E
and time0 are zero in this first supported contract. PZ/fixed-ion/scalar/isolated/
midpoint restrictions remain in native. step_fields validates E time/unit/location;
step(E) explicitly uses canonical atomic units. Named DipoleObservables exposes
charge dipole and electronic change, never an invented current density.

Reuse TDDFT process transport and candidate-worker replay through private protocol
hooks; defaults retain existing periodic behavior. E reply has six components;
input history always three. Public wrapper exposes in-memory snapshot/restore,
not inherited periodic disk/checkpoint/output APIs. Restore validates in a second
worker before swap. Failed step replay restores accepted history or closes state.
No whole-MPI restore, native memory checkpoints, or new convergence tests.

Tests first: receipt rejection/binary drift, independent state replay, wrong
snapshot/dt rejection, lease lifetime, E packet rejection, failed-worker recovery.
Then actual C2H2 20-step parity and snapshot continuation with at most two children.
Run full source tests, public docs/sample and installed-wheel tests. No Actions.

Implemented and verified in source/installed wheel. Evidence:
`../isolated-public-api-validation-2026-11-15.md`.
