# Common A/E/H field input

User approved A/E/H common input with None distinct from zero; keep crystal A-only,
LD E-only, and explicit conversions. Implement FieldValue (finite real vector,
absolute atomic time, position and declared units/sampling), FieldPacket optional
A/E/H, FieldContract using TimeSpec input ports and units/position/gauge, and
FieldAdapter.step_fields. New module field_input.py, no native numerical edits.

A means A/c. Unit labels are strict caller declarations with no conversion; H is
never interpreted as B. Each component may have a distinct time, using TimeSpec
clock/offset/lag. One model sample position; Yee-grid interpolation is explicit
upstream, not hidden here. Required-port validation before backend invocation.
Additional components allowed but unused unless declared. Missing is not zero.

Adapter either explicitly delegates one A/E vector to legacy step, preserving its
return value, or calls a custom state.step_fields(packet) for declared multiport
models. It borrows state. Caller supplies current accepted step index (default0),
keeps clocks aligned and does not independently advance/restore while attached.
Validation failures leave adapter unchanged; backend exceptions invalidate adapter,
without pretending that custom backend state rolled back. No MPI initialization,
E/H communication, isolated native TDDFT claim or implicit A->E/H->B conversion.

TDD: immutable ownership, missing/complex/nonfinite values, times/units/positions,
Yee-like H offset, gauge requirements, custom callback, failure behavior, exact
native LD/SBE legacy parity. Paired guide, AI JSON, runnable example; full tests,
sdist->wheel and installed example. Local only, no Actions/commit/push.

## F3 source audit after b6ae34a0

At the initial audit this was investigation only. See the F3a implementation
record below for subsequent experimental native validation. The installed public
packet/adapter alone does not imply new backend physics.

Observed in the bundled standard runtime:

- `src/rt/tddft_worker.f90` accepts periodic, transverse, midpoint A-driven
  calculations. Its protocol returns `rt%curr` via the main TDDFT loop.
- `src/rt/time_evolution_step.f90` calls `calcVbox` for `iperiodic=0`, before
  propagation. The periodic branch separately sets midpoint A and computes
  endpoint current. The isolated branch does not populate current in this path.
- `src/rt/em_field.f90:calcVbox` reads `rt%E_ext(:,itt_t)` and adds E dot r
  to the local scalar potential. Preserve this sample and update order; do not
  substitute the periodic midpoint-A implementation.
- `src/rt/dip.f90:subdip` computes electronic/ionic dipoles and their electronic
  change. Its isolated output is dipole change, not periodic current density.
  Verify the actual sign against density integration; comments alone are not a
  sufficient sign contract.

Implementation sequence:

1. Establish a small isolated GS/reference calculation and record its sampled E,
   electronic dipole, charge, and energy. Confirm the E sample time through the
   complete propagation path, including predictor/corrector evaluations.
2. Add a distinct, versioned E-worker capability and protocol. Reject an A-only
   executable before state creation. Keep the existing periodic worker protocol
   and numerical order unchanged. Retain explicit runtime manifests/provenance.
3. Add an isolated initial-state loader and state with `step_fields` requiring E.
   A/H availability must not imply that the native solver consumes those ports.
   Preserve independent workspaces and include E history and gauge in replay.
4. Return named dipole observables initially. Do not label a dipole difference
   as instantaneous charge current density. Any field coupling must explicitly
   define sign, differentiation time, effective volume and physical normalization
   and be tested separately before exposing a current-density adapter.
5. Compare the worker and unmodified reference path on identical E samples;
   interleave two independent states; restore and continue a third execution.
   Re-run periodic A-worker regression, then installed-wheel tests. Use local
   runs within 4 GiB / 4 cores; do not trigger GitHub Actions.

Alternatives considered: relaxing only the periodic guard would return the wrong
observable; reconstructing E through artificial A changes the input contract and
can change temporal sampling. A dedicated E path is therefore the intended next
implementation, with observable normalization resolved before Maxwell coupling.

This audit does not claim isolated runtime support, current normalization,
physical convergence, or a completed F3 acceptance test.

### Reference prerequisite executed

`tools/test_isolated_reference.py` now runs a small converged C2H2 GS and two
independent midpoint batch trajectories. See
`../isolated-reference-validation-2026-11-15.md` for measured results and the
isolated E_ext output-column caveat. F3 itself remains open. AETRS requires
look-ahead E and is excluded from the first worker implementation.

### F3a implemented and locally checked

The standard runtime E1 protocol and validation driver now pass direct reference,
interleaved independence, history replay and periodic regression checks. See
`../isolated-worker-validation-2026-11-15.md`. F3b public loader/capability/recovery
work remains open. The RT-step-config overlay explicitly rejects E1 because its
isolated path has not received this numerical validation.
