# Native field spatial current routing

The field leader keeps ordered global J for public output and checkpoints. Each other field rank receives only point currents whose transpose sample-map deposits intersect its owned Yee edges, including cross-boundary deposits. A native query computes this map from the initialized mg bounds, never a Python guess of spatial decomposition. A separate packed-current step ABI preserves the legacy replicated-current ABI. Missing capabilities retain legacy routing.

Response routing is selected once during coupling initialization. The existing native update order, edge weights, A/E sampling, rollback, and checkpoint formats remain unchanged. Rank-local histories retain global shape with unspecified non-owned rows. Legacy checkpoints are consumed by selecting only the required rows. No additional application-wide collectives per step are introduced.

Verify first with the missing-ABI failure, then a freshly audited Git-only runtime build. Compare old/full and new/local execution bitwise with points on and away from decomposition boundaries and empty owners. Verify saved/restarted execution and protected-restore failure/retry. Run source and installed-wheel regressions; sync Japanese/English docs and AI JSON. No new GS calculations or manually dispatched heavy CI.

## Verification completed

- Prior changes committed/pushed as `a9ba9212`. New routing changes remain reviewable in the worktree.
- Git-only source export: all 740 manifest assets match; fresh embedded CPU/MPI runtime built with retained source/object eigensolver audits. Scalar IEEE checks retained.
- Source and source-distribution→wheel installed suites: 668 tests, one optional skip, zero failures.
- Four-rank field decomposition 2×1×1: point ownership [0,1]/[1]. Six-rank total with field decomposition 4×1×1: []/[0,1]/[1]/[]. Each owned point carries three current components. The field reporting leader retains full ordered J even if it owns no deposit edges.
- A/E/J bitwise equality against the preceding binary/wheel and against the new runtime's legacy replicated-current path. A nonzero z-current crossing a z decomposition boundary was tested separately. Native packet coupling also matches the replicated-current control.
- Owned and legacy full-row checkpoints resume bitwise identically, with compatibility checks preserved. Protected restoration injected into field rank 1 preserves all original bytes, replays zero steps and matches three future retry steps.
- Direct native NaN, Infinity and count mismatches are rejected collectively before advance; subsequent valid execution matches the continuous control.
- Final installed wheel: six-rank owned-current/invalid-input test and four-rank direct resume passed. Japanese/English guides, example notices and AI JSON synchronized.

Manual invocation (ordinary sample arguments; use the SBE-only layer fixture for GS-free tests):

```sh
mpiexec -n 4 python tools/test_spatial_current_mpi.py BUILD RUNTIME BUNDLE NEW_OUTPUT --layers LAYOUT --field-grid 2 1 1 --tddft-ranks 1 --steps 96
mpiexec -n 4 python tools/test_mixed_protected_restore_mpi.py --failure-rank 1 BUILD RUNTIME BUNDLE NEW_OUTPUT --layers LAYOUT --field-grid 2 1 1 --tddft-ranks 1 --steps 96 --direct-save-at 48
```

Canonical validation record: `calc2026/salmonpy-spatial-current-20261121/record.json`; scripts, initial failed fixture log and corrected runs retained there. No new GS collection or manual CI dispatch. Local Mac evidence does not establish multi-node scaling or Fujitsu/GPU support. Remaining communication task: coalesce outer failure agreements only after consistent failure/retry ordering is proven.
