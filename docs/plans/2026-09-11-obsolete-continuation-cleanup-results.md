# Obsolete continuation cleanup checkpoint

## Deleted scope

The top-level GS dispatcher already sends both divided and continuation
selectors to `run_dg_hybrid_divided_ground_state_for_main` (the latter via an
input-compatible alias). Therefore the continuation-only branch inside the
legacy overlapping-Wannier driver cannot be entered through these selectors.

Removed that branch and `run_dg_hybrid_concrete_continuation`, which performed
complete diagonalization within its local stage/SCF loop. Also removed its
unreferenced internal helpers `checkpoint_grid_real_fingerprint`,
`selection_added_members`, and `set_dg_hybrid_trial_state`, associated output
arrays, and now-unused imports. `main_dft.f90` has a net reduction of 542 lines.

The conventional DC/LCFO paths, local CG/terminal LCFO route, Exp RT, and
checkpoint formats are unchanged. The shared controller module is retained:
its candidate-acceptance/publication-rank policy is still used by the v5
publisher. Its standalone schedule tests remain, without asserting reachability
of the removed production driver. This is not a complete dead-code audit of
every remaining module or legacy selector.

## Test migration

Old driver-dependent assertions now follow the compatibility alias to the
terminal refinement route and verify publication after that loop. The new
absence assertion failed before deletion and passed afterward.

Two pre-existing audit mismatches were corrected: the LCFO signature assertion
now recognizes the appended optional core-density output, and three explicitly
identified Hybrid tests are exempt from a blanket obsolete-WPW filename glob.
No active PW projection test was deleted. Historical monolithic WPW references
in the August 20 plan are labeled accordingly.

## Verification

Passed: terminal route, continuation alias/publication, divided DC controls,
overlapping-Wannier route, obsolete-route inventory, density snapshot route,
density decomposition unit tests, Fourier unit test, and `git diff --check`.
The production build exited 0. Checkpoint/publication MPI tests passed on
1/2/4/8 ranks; the retained standalone continuation schedule runner also passed.

Fresh production run: `/tmp/si8-route-cleanup-20260911`, eight ranks, exit 0,
eight `end SALMON` receipts. DC seed and WF caches were read-only hits. All four
terminal density-change/total-energy log records match the earlier diagnostic
run exactly at printed precision. The final algebraic residual is 9.33e-13,
orthogonality defect 4.33e-15, electron defect zero. This confirms numerical
regression stability for this fixture, not physical accuracy or Si64 validation.

Run log SHA256:
`72acda01789701411a61f710cf5aca5fe857f1e13efe026acb2ce91678783665`.
Build log SHA256:
`f2c9c8f4e2e4c4545141f0892317bd530e24fb2313cd7b1eaf3d874ded77f599`.

No original-worktree files, saved calculations, or reusable caches were deleted.
Deleted source can be recovered from Git history (preceding commit `46dccf8a`).
