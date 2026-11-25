# Radial-complete molecular field input validation (2026-11-25)

The E/H input uses real r^(l+2k)Y_lm terms through total degree8. Fortran owns
fitting, evaluation, electric/magnetic coupling and numerical checkpoint data.
An explicit center is an arbitrary real coordinate, not a required grid point.
FieldV2 and direct binaryV5 identify the coefficient representation without fallback.

- Independent Fortran oracle:400 fitting points and53 new off-grid positions;
  degrees2/4/6/8. Maximum E/H error1.30e-15au; scalar potential1.18e-15,
  vector potential A/c6.85e-19. A separate Cartesian-polynomial SVD implementation
  cross-checked the radial-complete basis on saved FDTD samples.
- C2H2, existing PZ GS, fixed ions, CPU1rank/1thread, dt=.003au,20updates;
  off-grid center(.2,-.3,.1)bohr, nonzero quadratic radial E/H terms.
  Prescribed/live input agrees exactly; independent control state is unchanged.
- Direct checkpoint at7, process termination, restore and updates8..20:
  all9observables bit-exact with continuous execution;0extra evolution on restore.
  Savepoint/output-stream checks also pass. Electron loss fraction1.12e-14.
- Normal uniform-field response still agrees with the batch reference within
  3.86e-12e*bohr. Git-tracked-only source preparation/build and application-object
  eigensolver audit pass. Source tests31pass; installed-wheel tests41pass.
- All measurements are local Mac validation, not cross-platform or physical
  convergence evidence. No new SCF. Tip/Bow Tie required order, FDTD grid convergence,
  actual molecular-response convergence and bidirectional feedback remain separate.

Detailed local records: `calc2026/salmonpy-c2h2-nearfield-20261125/` under the configured
results store: `radial-final-validation.json`, `radial-analytic-validation.log`,
`radial-final-wheel-molecule-validation/`, `radial-git-wheel-uniform-validation/`,
`radial-complete-basis-crosscheck.json`. Local OneDrive storage is not cloud-sync confirmation.

The final receipt explicitly records nonzero radial-term testing (`radial_terms=true`).
Normal uniform-field results before/after the basis change were bit-exact for all9observables
(Fortran comparison maximum absolute difference0). No CI or remote physics campaign was launched.
