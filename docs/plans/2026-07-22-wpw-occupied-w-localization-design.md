# Occupied-W Localization Design

## Purpose

Measure the localization of the 128 projected occupied Wannier functions with
the same dimensional convention as Wannier90, then localize them within the
occupied subspace if their widths are materially larger than the earlier
`sqrt(Omega) ~= 1.2 Angstrom` result.  Do not infer localization from the norm
on the edge of a fragment buffer: once the unwrapped buffer crosses a periodic
cell, that edge can approach another periodic image of the Wannier center.

## Diagnostic contract

For each occupied W, construct the three first-harmonic periodic moments

`z_a = integral |W(r)|^2 exp(i 2 pi r_a/L_a) dr`.

The phase of `z_a` gives a branch-safe periodic center.  Relative to that
center, use the minimum-image displacement on the physical 24-point periodic
cell to accumulate the real-space second moment

`Omega = integral |W(r)|^2 |Delta r_min|^2 dr / integral |W(r)|^2 dr`.

Report `Omega` in Angstrom squared and `sqrt(Omega)` in Angstrom.  Publish the
minimum, mean, maximum, number above 1.2 Angstrom, and the global W identity,
fragment, center, norm, and width of the maximum.  Also report per-W radial
cumulative norm at fixed physical radii so a genuinely broad W can be
distinguished from a periodic-image shell.

The diagnostic is evaluated on one canonical physical cell, not on P.  Every
canonical grid point contributes exactly once.  MPI ranks reduce moments and
norms collectively.  Non-finite values, non-positive norms, ambiguous source
ownership, or a center with a vanishing first harmonic are collective
failures.

## Localization

Projected Wannier construction remains the deterministic initial gauge.  A
localization rotation, if needed, is restricted to the 128-dimensional
occupied subspace and is unitary.  Use the projected periodic-position
matrices `U_a = <W|exp(i 2 pi r_a/L_a)|W>` and minimize the total periodic
spread by unitary Jacobi rotations.  This is the Gamma-point specialization of
the Marzari-Vanderbilt localization problem and does not require an external
Wannier90 run.

The accepted rotation must reduce (or preserve within roundoff) total spread
and must preserve:

- W orthonormality;
- the occupied density/projector;
- the number and stable identity of W rows;
- the analytic-gradient transformation;
- core values and periodic P construction under the same rotation.

Failure to converge does not silently fall back to a partially rotated gauge.
The last accepted unitary state is restored and the bootstrap stops
collectively.

## Staging

1. Add and validate the standalone periodic-center/spread diagnostic.
2. Run B=6 and compare `sqrt(Omega)` directly with the earlier 1.2 Angstrom
   reference.
3. Only if the diagnostic confirms excessive spread, add the occupied-subspace
   localization rotation.
4. Replace the current outer-shell pass/fail gate only after radial-tail and
   spread behavior supply a physically defined buffer convergence criterion.

No tolerance is relaxed, and fixed-H, checkpoint publication, and RT remain
blocked while localization or buffer sufficiency is unresolved.

## Verification

Unit tests use analytic localized functions placed across a periodic seam and
verify their center and spread, translation invariance, normalization
invariance, and rejection of invalid moments.  A two-state rotated fixture
verifies that localization lowers total spread while preserving its projector.
MPI tests split the physical cell and require the same results as the serial
fixture.  The Si64 B=6 diagnostic records all 128 widths before any production
gate is changed.
