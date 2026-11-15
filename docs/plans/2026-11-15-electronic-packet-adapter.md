# Explicit electronic binding to field packets

Status: implemented and locally validated; see ../electronic-packet-adapter-validation-2026-11-15.md.

## Scope and choice

Wrap an existing `MaterialBinding`, preserving its provider, execution mode,
reduction and material identity. Do not modify TDDFT/SBE numerical kernels or
teach `ResponseSession` to guess a model's driving field. This permits both
`step(A)` and legacy `advance(A_old, A_new)` through one explicit boundary.

Alternatives are changing every backend (duplicates clock handling), or implicit
adaptation in the session (hides physical phase selection). An explicit binding
wrapper is preferred and follows the existing explicit-adapter decision.

Proposed entry: `packetize_A_binding(binding, *, field_phase)` in
`salmonpy.material_response`. Require `field_phase` to be exactly `endpoint` or
`yee_half_step`; do not provide a default. Keep position and per-point contracts
in the caller's `field_contracts`, not in a shared provider.

## Clocks, units and signs

The wrapped state retains its local clock n*dt. The wrapper exposes point A/c and
charge-current j ports, both offset 0 (`endpoint`) or +1/2 (`yee_half_step`).
The latter is the existing `CoupledStepper` electronic-clock mapping, not
interpolation, an extra integration step, or a change to the worker's time origin.
Record original and mapped TimeSpec plus the phase choice in result metadata.

Admit only explicit original A/j endpoint point-sampling declarations with
t0=clock_offset=0, lag=0 and matching dt. Reject J_matter, E/H requirements,
interval averages and nonzero initial A. Existing bindings already normalize
charge current; do not negate it again. Initial field feedback remains zero.
Document that a model depending on absolute time must support this mapping;
the adapter cannot certify arbitrary callbacks' physical validity.

## Lifetime and failure

Provider preparation remains cached by the original key. Every create_state
constructs an independent wrapper and accepted-A history. Preserve reduction. Restrict this increment to leader_only execution and one
rank per material group. Collective backends are rejected before preparation;
rank-local validation cannot safely precede their collective callbacks.
If validation fails after creation, close the newly created raw state before
raising. Closing the wrapper closes its owned state once, never its provider.

Use FieldContract/FieldAdapter validation for packet units, position and time.
Before each callback check any available raw clock and its unchanged declaration.
Only update accepted A and local step count after success. Failure after a
callback invalidates the wrapper/session; do not promise collective rollback.
Expose initial_current when supplied by the raw state without treating it as
initial feedback. Preserve separate per-point histories even for shared bindings.

## Verification and implementation sequence

1. Add failing tests for both callback forms and both phases, invalid clocks,
   nonzero A, sign preservation, independent histories, cleanup and failure.
2. Implement the wrapper, publish its API inventory entry, and test real SBE
   states against unwrapped states for identical accepted A histories.
3. Extend the three-rank 1D/native packet sample to select the explicit wrapper.
   Compare A/E/j against legacy coupling, including multiple independent points.
4. Run the managed TDDFT provider with two independent material groups, at most
   four total cores and two worker processes. Compare the same accepted A/j
   sequence with legacy coupling. A synthetic test alone does not close this item.
5. Update paired FIELD_INPUT guides, AI JSON and runnable examples; build from
   sdist and repeat tests using the installed wheel. Do not run GitHub Actions.

No dipole-to-current conversion, local/self-field policy, general subcycling,
magnetic response or new physical convergence sweep belongs to this adapter.
