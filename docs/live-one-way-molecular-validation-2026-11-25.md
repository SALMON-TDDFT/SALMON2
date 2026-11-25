# Single-launch spatial FDTD -> molecule (2026-11-25)

Public entry: `python -m salmonpy.examples.molecular_fdtd_response --help`.
The shared runner in examples/_grid_coupling.py is also used by the existing
manual validation tool; its numerical implementation is not duplicated.

Scope: one MPI launch, two FDTD solver ranks plus one coordination rank,
one separate single-thread molecular worker. Fixed40³ FDTD grid,2bohr spacing,
synthetic LD box,24³ E/H patch, off-grid C2H2 at(-2.54,1.62,1.94)bohr,
dt0.003au,96molecular updates /97field updates. Molecular current never enters
FDTD, there is no own-vacuum FDTD and no charge-density or magnetization query.
No new GS was computed. Existing authorized reference data were reused.

Unchanged validated binaries:
- FDTD library SHA256:5366b7005297fa9b5ed9fc1141f313a423b9a7ed9c2812f8fa28e3e007ccb600
- Molecular executable SHA256:14557305ee88f46a99d3e15352f11ffb88c7c1c80626b93b21b2659a1b035897

Checks:
- Live field histories match the prior illumination-only FDTD control bit-exactly.
- Live molecular dipole change, charge dipole and electron number match the saved-field
  route bit-exactly when both use the same actual field-step numbers.
- The old illumination history began at field step1, not0. A first comparison
  incorrectly treated its first row as step0 and disagreed. That comparison is
  not evidence of numerical disagreement between methods. The new sample explicitly
  records field step0; aligned saved histories prepend this measured frame.
  Failed trial logs remain available. Do not silently relabel field-step1 as0.
- Save at molecular48/field49; separate MPI launch imports with0extra evolution
  and completes48new molecular updates. All histories match continuous execution.
- Installed wheel --stop-after48 exits with those actual indices; another launch
  resumes. All history arrays and four ordinary field/molecular .data streams match
  the source continuous run bit-for-bit. Resumed history includes the full prefix.
- The previous8-step electric feedback validation still matches its prior record
  bit-for-bit after moving the common runner. This is a code regression, not a
  recommendation to use feedback for the current molecular workflow.
- Related tests:45source and45installed. The existing test was extended; test
  method count was not increased. Help works without MPI initialization.
- Source distribution -> wheel passed, including retained source/object
  eigensolver prohibition checks. Make was used because Ninja was unavailable.
- Day-to-day guides and manual.json synchronization and whitespace checks pass.
  No additional CI run was manually dispatched.

Canonical local records under calc2026/salmonpy-c2h2-nearfield-20261125:
- one-way-live-ld-final/comparison.json (continuous and aligned saved-field check)
- one-way-live-ld-resumed/summary.json (source process-exit restore)
- one-way-live-installed-first/summary.json (installed stop/save)
- one-way-live-installed-resumed/comparison.json (installed restart/stream check)
- one-way-refactor-feedback-regression/summary.json (existing-route regression)

No physical-convergence, strong-back-action, general geometry, cross-layout restart
or arbitrary MPI-failure recovery claim is made. Radiation and molecular-to-field
feedback remain outside the current one-way scope.
