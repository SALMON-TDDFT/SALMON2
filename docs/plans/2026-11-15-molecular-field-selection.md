# Explicit molecular drive selection: next implementation boundary

Status: design, not a public implemented API. Builds on
2026-11-15-molecular-self-field.md; does not enable two-way TDDFT by itself.

## First increment

Keep E->dipole, dipole->current and field selection separate. Introduce a pure
selection operation over already sampled electric FieldValue objects. It must
not start solvers, interpolate, infer volumes or alter the FDTD update order.

Require an explicit policy, with no default:
- `provided_drive`: use the caller's already selected drive, with a required
  provenance record describing its physical meaning.
- `environment_scattered_only`: total environment field minus this molecule's
  own-source vacuum-reference field. External illumination is present only in
  the environment run; the reference has no external illumination.

A raw total-field feedback policy is not automatically selected. A finite
physical footprint/cavity prescription can be supplied through provided_drive;
this operation cannot certify such a model's physical validity.

## Data contract

Every selected drive records target molecule ID, field sample time and position,
policy, input values, and source provenance. The subtraction input must identify
exactly one source molecule, matching the target. Multiple source IDs and a
summed vacuum response are rejected. This preserves other molecules' direct
vacuum interactions.

Use existing FieldValue E:au, bohr and point-sampling contracts. Subtraction
requires matching time, position, units and sampling. No zero-filled missing
references, silent time alignment or fallback to total E. Return an immutable
FieldValue plus immutable metadata; inputs remain unchanged. Reject nonfinite
arithmetic results.

Reference provenance must declare identical deposition footprint and sample map
identifiers and vacuum/no-incident conditions, with matching source history
identity. These are caller declarations, not a cryptographic certification that
separate Maxwell calculations used identical inputs. Future solver helpers can
construct and verify these records from their own manifests.

## Regression sequence

1. Analytic additive fields: own-only vacuum feedback cancels, prescribed drive
   survives, and a second molecule's contribution is not removed.
2. Reject target/source mismatch, multiple-source reference, different positions,
   time/units/footprints/sample maps and incident illumination in the reference.
3. Verify immutable results, metadata independence and finite-value guards.
4. Supply a small runnable sample, paired guides and searchable AI task entry.
5. Before closing R4b, separately implement and test a causal native/reference
   schedule with oscillator source-work accounting; only then attach isolated
   TDDFT. A pure subtraction test does not complete two-way molecular coupling.

The environment_scattered_only policy removes the vacuum radiation-reaction
contribution as well as the grid self contribution. It is not full radiative
molecular dynamics or spontaneous emission. Production physical convergence is
user-owned; developers verify signs, clocks, causality and bookkeeping.
