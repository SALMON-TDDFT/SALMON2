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
