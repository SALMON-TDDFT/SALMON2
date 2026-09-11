# Active DG route cleanup

## Approved scope

Keep conventional DC/LCFO and the current fragment-local Hybrid GS, terminal
LCFO refinement, v5 handoff and Exp RT. Remove obsolete alternatives rather
than silently translating their input into a different calculation.
Preserve the current integration checkout, saved calculations and caches.

## Checkpoint 1: unused internal routines (implemented)

Removed seven internal routines from main_dft.f90 with no source or test
callers: form_dg_hybrid_coefficient_actions,
measure_dg_hybrid_operator_covariance, measure_dg_hybrid_projector_covariance,
evaluate_dg_hybrid_distributed_low_energy_symmetry,
measure_dg_hybrid_core_density_covariance, ow_hybrid_density_to_dc and
ow_hybrid_density_from_dc. Their bodies total 294 lines. This does not remove
the active symmetry validation or change the production numerical path.

The terminal-route test first failed on an obsolete helper, then passed after
removal. The production CMake build exited 0. Five source contracts passed:
obsolete routes, terminal refinement, continuation publication, divided DC
controls, and overlapping Wannier. git diff --check passed. No fresh MPI
calculation was performed for this checkpoint. No commit or push performed.

## Checkpoint 2

Checkpoint 2: removed the unreachable localization-first arm (264 lines
including its conditional wrapper) from the legacy GS entry and its unused
DG_W90_UNCONSTRAINED import. Retained the constrained legacy arm pending full
legacy-entry retirement; the separated fragment-Wannier production routine
was not changed. The OW route test failed before deletion and passed after
it. Production build, OW route, terminal refinement, divided DC controls and
obsolete-route contracts passed. An additional fragment-Wannier source gate
fails because it forbids every dc_lcfo call, including the existing diagnostic
call in the active routine. That call also exists in HEAD before cleanup;
the gate needs a scoped diagnostic-versus-production assertion. No runtime
MPI regression has been run for this checkpoint.

## Checkpoint 3

Removed the unreachable divided/continuation preparation arm (186 lines)
and its now-unreferenced external prepare_dg_hybrid_divided_production_basis
routine (100 lines). The new absence assertions failed before removal.
Migrated divided DC controls from the dead arm to the separated active entry
and its local phase. The OW occupation test now checks the retained LCFO
electron-count validation that remains reachable instead of an assignment
that existed only in the removed arm. The fragment-Wannier gate now permits
exactly one ordinary LCFO call inside the optional density-diagnostic block,
requires its core-density output and disabled file output, and forbids the
call elsewhere in that entry.

Production build exited 0. Source contracts passed: terminal refinement,
fragment Wannier, divided DC controls, overlapping Wannier, continuation
publication, obsolete routes. git diff --check passed. No MPI calculation,
commit or push in this checkpoint. Deleted source is recoverable from HEAD
f94985e8; existing logs, caches and unrelated files are preserved.

## Remaining work (not implemented)

Checkpoint 4: removed the selectable legacy global Hybrid SCF arm and the
V3 overlapping-Wannier RT main driver/dispatch, including their direct
imports. The three legacy selectors remain parseable solely to report
explicit retirement errors; no automatic v5 conversion is performed.
The input check rejects old RT even when the new RT selector is also set.
The retired-entry source contract failed before changes and passes now.
Updated old production assertions to test retirement while retaining shared
library tests. Build exited 0; terminal refinement, fragment Wannier, divided
controls, OW route, old SCF library/retirement, retired entries, continuation
publication and obsolete-route source contracts passed. No runtime MPI
input-rejection test or physical regression was run at this checkpoint.
The legacy helper libraries/callbacks and redundant old input validation
below the retirement errors still require pruning; full cleanup is not done.

Checkpoint 5: removed all six old ow_hybrid SCF callbacks (92 lines) and
their exclusive workspace declarations/density import. Removed the unused
dg_hybrid_scf and rt_dg_overlapping_wannier libraries from production CMake
source lists; their sources remain available to standalone tests pending
test-support cleanup. Removed unreachable legacy selector requirements and
V3 RT-specific validation after the retirement errors. Length gauge now
explicitly requires the current Hybrid RT selector. The retired-entry
absence test failed before callback removal and passed after it.
The final production build exited 0 and eight source contracts passed
(retired entries, OW route, terminal refinement, fragment Wannier, divided
controls, continuation publication, old SCF test support, obsolete routes).
git diff --check passed. No MPI regression, commit or push was performed.

Checkpoint 6: extracted the candidate receipt, v5 rank-policy validator and
its collective/fingerprint helpers into dg_hybrid_publication_policy. Main
GS imports that module directly. The old continuation controller reuses the
same implementation for standalone tests and is no longer compiled into
the production executable. Updated the three standalone test compile lists.
The retired-entry assertion failed before the split and passed afterward.
Production build and active source contracts passed. v5 checkpoint tests
and continuation-controller compatibility tests passed on 1/2/4/8 ranks
after permitting local MPI communication (the sandbox-only attempt failed
at socket startup, before numerical execution).

Fresh cached eight-rank production run:
/tmp/dg-active-route-regression.eNtfPW/run.log, exit 0, eight end SALMON
receipts. DC seed and fragment WF were read-only hits. All four printed
terminal density-change/energy records exactly match
/tmp/si8-route-cleanup-20260911/run.log. Final eigensystem residual is
8.69e-13 versus 9.33e-13 before cleanup; electron defect is zero. Both runs
exhaust the three additional refinement steps without density convergence;
this is regression agreement, not a claim of physical convergence or Si64
accuracy. No existing calculation was overwritten. No commit or push.

Checkpoint 7: removed the unreachable continuation scope-receipt arm in the
legacy GS entry and its local type/import. Removed continuation_state,
divided_scf and real_space_residual from production CMake sources after
checking source consumers. Standalone diagnostic/test sources are retained.
Removed unused DG ONLY imports in main_dft. A continuation-line formatting
error during import cleanup was caught by compilation and corrected; the
final production build exited 0. Retired entries, OW route, terminal
refinement, divided controls and fragment-Wannier source contracts passed,
as did git diff --check. No fresh physical run in checkpoint 7; checkpoint 6
records the most recent cached regression. No commit or push.

Checkpoint 8: removed whole unused declaration statements from main_dft,
including old divided LCFO matrices, face/interior workspaces, selection and
closure receipts, SCF fingerprints and certified-localizer arrays. Each
removed variable had only its declaration occurrence in the full source;
mixed used/unused declaration statements were intentionally left intact.
Representative absence assertions failed before removal and passed after.
The final production build exited 0 and eight source contracts passed
(retired entries, OW route, terminal refinement, fragment Wannier, divided
controls, continuation publication, old SCF test support, obsolete routes).
git diff --check passed. No numerical statements changed in checkpoint 8;
the most recent runtime regression remains the checkpoint 6 cached run.
No commit or push; existing logs and caches preserved.

Checkpoint 9: pruned unused variables from mixed declaration statements in
main_dft. Retained variables preserve their type and array shape. The
selection required exactly one identifier occurrence in the entire source,
namely the declaration; declarations with initializers were excluded.
Representative absence assertions failed before changes and passed after.
The production build exited 0 and the same eight source contracts passed.
git diff --check passed. No numerical executable statements changed and no
new MPI calculation was run at this checkpoint. No commit or push.

Checkpoint 10: moved six non-production reference modules to
tests/dg/legacy_support, with a README explicitly excluding production use.
Updated all explicit Python test compile/read paths. The retired-entry
contract first failed on old source locations and passed after migration.
No numerical implementation changed in these moved modules (only a test-only
header was added). Production build and eight route/source contracts passed.
All six standalone MPI runners passed: SCF/divided SCF/controller/V3 RT on
1/2/4/8 ranks; real-space residual and continuation state on 1/2/4 ranks.
git diff --check passed. Sources were moved, not discarded; saved calculations
and caches were untouched. No commit or push.

Checkpoint 11 (audit, no numerical changes): production no longer imports
legacy_support and the retired SCF/RT entry contract passes. However the
GS dispatcher still calls run_dg_overlapping_wannier_ground_state_for_main
when yn_dg_dc_overlapping_wannier=y and both active Hybrid selectors are n.
This is the retained bare one-shot constrained-Wannier/EigenExa route, not
ordinary DC+LCFO. It still writes the old OW checkpoint. Consequently the
entire cleanup must NOT be described as current-route-only or complete.
The OW source gate currently verifies this retained arm, so removing it
requires separating its standalone numerical checks from entry wiring.

Checkpoint 12: removed the approximately 2,000-line bare one-shot OW GS
driver and replaced its dispatch with a defensive retirement error. Input
rejects bare OW selection regardless of EigenExa availability. Ordinary
DC/LCFO, the divided entry and its compatibility alias are retained.
Replaced the retired driver's extensive static source gate with current
dispatch/ordinary-LCFO preservation checks and active route gates; standalone
numerical module tests remain. Migrated rank, PP provenance, seed ordering
and v2.3 integration checks to the active entry/publication path. Build,
active dispatch, DC controls, DC seed, canonical PP, continuation, obsolete
route and v2.3 integration checks passed. git diff --check passed.

Fresh cached 8-rank run /tmp/dg-bare-gs-retirement.hNay9h/run.log exited 0
with eight end SALMON records. All four terminal density/energy records
exactly match the pre-cleanup /tmp/si8-route-cleanup-20260911/run.log.
The pre-existing refinement exhaustion remains; this is a regression check,
not physical convergence. Existing caches/logs were not overwritten.
No commit or push.

Checkpoint 13: removed 40 unreferenced internal subroutines and five
unreferenced internal functions, 2,881 lines, in successive dependency
passes. No reference to each routine remained outside its own definition
when removed. Build and source checks passed (current dispatch, terminal,
fragment WF, DC controls, seed, canonical PP, integration, continuation,
old SCF test support and obsolete routes).

RUNTIME REGRESSION NOT PASSED: /tmp/dg-helper-retirement.v45xBr/run.log
exited 1 before Hybrid GS, reporting status=invalid cause=ownership_mapping
detail=rank_owned_grid_mismatch. Its inputfile is byte-identical to the
successful checkpoint 12 input. Existing cache files were not rewritten.
The failure is the aggregate ownership fingerprint comparison in
dg_dc_seed_checkpoint.f90, not the fragment-id comparison. It fingerprints
total/fragment grids, orbital ownership, fragment maps and projector grids;
the reported diagnostic alone does not identify which component differs.
Root cause remains unresolved. Do not relax reuse guards or claim numerical
regression success for checkpoint 13. Changes remain uncommitted/unpushed.

Checkpoint 14: ownership regression root cause identified and fixed.
init_ps assigned ppg%jxyz_min=ppg%nps before calc_nps initialized nps.
The fresh fragment grid carried 62992832 in every jxyz_min entry, while the
total-system grid carried zero. These otherwise-unused bounds are included
in the DC seed ownership fingerprint. Removing uncalled routines exposed
this pre-existing uninitialized read through changed memory layout.

Changed only the fresh jxyz_min initialization to zero; kept both min/max
arrays in the fingerprint and did not change any reuse guard, MPI-size or
rank-fragment requirement. A source regression test failed before the fix
and passed afterward. Diagnostic comparison showed only ownership_map(7),
the fragment projector-grid component, changed on all eight ranks.
The stored mapping fingerprint 398352827599223578 then matched again.

/tmp/dg-ownership-trace.r0mSNh/bounds-trace.log records failure and raw bounds;
fixed-trace.log records successful read-only DC/WF cache hits, exit 0 and
eight end SALMON records. All four terminal density/energy records match
the pre-cleanup baseline exactly. Existing refinement exhaustion remains;
this is regression agreement, not physical convergence. Temporary trace
writes were removed after diagnosis; the final non-trace build exited 0
and the initializer/current-route/DC-seed/PP/integration source gates passed.
No cache data changed and no commit or push was performed. A cache produced
with nonzero uninitialized metadata may still correctly fail the unchanged
compatibility check; no unsafe automatic migration was added.

1. Continue final unused imports/state audit after the resolved regression.
2. Remove the retired drivers' exclusive callbacks and obsolete input checks
   after checking shared consumers. Entry retirement itself is implemented.
3. Separate the live v5 publication policy from obsolete continuation control;
   remove production compilation of test-only drivers where appropriate.
4. Run build, active GS/RT contracts and cached MPI regression calculations.
   Do not claim completion of the full route cleanup at checkpoint 1.
