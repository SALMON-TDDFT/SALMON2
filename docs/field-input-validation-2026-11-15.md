# Common A/E/H field input validation (2026-11-15)

User-approved common local input is implemented in salmonpy.field_input:
FieldValue (immutable finite real vector, atomic time, position, unit/sampling),
FieldPacket (optional A/E/H; None distinct from zero), FieldContract (TimeSpec,
units, gauge, position) and FieldAdapter.step_fields. A means A/c; H is not B.
No numerical kernel or MPI payload was changed.

Legacy A/E adapters explicitly select one endpoint vector in canonical atomic
units and preserve original return values. Custom multiport callbacks receive
only declared components. Per-port times include clock offset, port offset and
lag; spatial interpolation and all unit/field conversions remain explicit upstream.
Gauge declarations constrain required input but do not certify backend physics.

Ten focused tests cover immutable ownership, missing versus zero inputs, invalid
numeric types, units/sampling/positions/time, staggered H, gauge/port validation,
legacy unit/endpoint restrictions, undeclared-field removal, failure behavior,
backend declaration/clock drift and native LD/SBE direct-call parity.

Review found that initial code checked backend timing only on the legacy path.
The final implementation checks exposed timing/dt/time for both paths before and
after stepping, including clock_offset; explicit contradictions are rejected.
A post-call mismatch invalidates the adapter, without promising backend rollback.
Models exposing none of this metadata still rely on caller declarations.

Evidence:

- Source suite: 399 tests, one filesystem-condition skip.
- Fresh installed wheel suite: 399 tests, two skips (filesystem and optional plot).
- sdist -> wheel succeeded (max two build jobs); isolated installed imports.
- Native LD and SBE wrapped calls match direct calls exactly for supplied four-step
  fixtures. Installed sample checks 20 native LD steps and custom staggered E/H
  delivery. H is recorded independently, not treated as a magnetic material law.
- Paired FIELD_INPUT guides and AI JSON synchronized; diff whitespace check clean.

Local logs: `/tmp/salmonpy-field-input-*.log`.
Installed sample: `python -I -m salmonpy.examples.field_inputs`.
Environment: Apple Silicon/macOS, Python 3.14; one process/thread for examples.
No peak-memory or physical-convergence claims. No native solver rebuild, GitHub
Actions, commit or push.

Remaining: native isolated-system TDDFT length-gauge E connection, MPI A/E/H
transport and capability checks, field solver collocation/interpolation adapters,
magnetic-interaction backend support. Existing crystal A-only convention remains.
This local interface does not imply that any current native solver consumes E/H
or that MPI ResponseSession now transfers them. Input values, clocks and unit
labels must reflect the backend's actual implemented equations.
