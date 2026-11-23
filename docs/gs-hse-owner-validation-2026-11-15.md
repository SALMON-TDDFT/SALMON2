# GS-owned HSE lifecycle validation

Local Apple Silicon/Mac verification on 2026-11-15, following the SCF flag routing stage. No GitHub Actions were used.

## Change and limits

`main_dft` creates a local `s_exx_rt_state` and calls `exx_bind_gs_system` after system/grid initialization and before electronic initialization. Eligibility reuses the existing native-resource predicate and additionally requires ordinary DFT without geometry optimization. Successful binding emits `EXX_GS_OWNER explicit`. The owner stays bound through all result/checkpoint output and is unbound/cleared at normal completion.

The scope is periodic scalar unpolarized HSE06, full-support CPU exchange, no Wannier/symmetry reduction, integer-occupation temperature setting and one rank in each electronic decomposition. Unsupported modes return status 2 without binding and keep the legacy path. Status 1 indicates a binding conflict. DC, distributed/thermal/localized exchange and geometry optimization are not newly supported independent-owner configurations.

This is not a new Python GS context, complete independent GS state, concurrency support, or failure-atomic cleanup after native `error stop`. Remaining shared solver settings and DC resources are unchanged.

## Verification

- `test_pre_scf_runtime.py --require-gs-owner` failed against the preceding build because the explicit-owner branch was absent: `/tmp/salmonpy-gs-owner-red-20261115`.
- Serial and MPI builds completed; source/pre-link/all-object forbidden-eigensolver audits passed. Audited records: `/tmp/salmonpy-gs-owner-{serial,mpi}-build-20261115.json`.
- Native owner probe passed, including eligible GS binding, rejection without mutation for spatial decomposition, independent resource/flag behavior and localization routing: `/tmp/salmonpy-gs-owner-probe-20261115`.
- Four tiny serial GS cases with required explicit ownership match original immutable fast-HSE phase/energy diagnostics and information/eigenvalue/k files: `/tmp/salmonpy-gs-owner-gs-20261115/result.json`.
- Two-/four-rank DC, four pre-SCF configurations each, remain on the legacy owner and match original diagnostics/scalar files: `/tmp/salmonpy-gs-owner-dc2-20261115/result.json`, `/tmp/salmonpy-gs-owner-dc4-20261115/result.json`.
- Eight TDDFT independent-context comparisons remain identical: `/tmp/salmonpy-gs-owner-contexts-20261115/result.json`.
- Python source suite: 305 passed. sdist → wheel → fresh install: 305 tests, one optional plotting skip, no failures. The installed tool independently repeated all four GS cases with explicit ownership required.
- Installed prepared Fortran bytes and runtime/AI guides match the tested checkout. Source line-length and manual synchronization checks passed.

Local artifacts: `/tmp/salmonpy-gs-owner-dist-20261115`, `/tmp/salmonpy-gs-owner-installed-20261115`, `/tmp/salmonpy-gs-owner-wheel-gs-20261115/result.json`, `/tmp/salmonpy-gs-owner-python-tests.log`, `/tmp/salmonpy-gs-owner-installed-tests.log`.

These are coarse/synthetic same-platform implementation regressions, not converged materials validation. Full checkpoint payload equivalence and repeated complete GS solves in one process are not claimed. The original fast-HSE comparison commit remains `c9e8ce431f43eefd4a88998134adb4796f8ec837`; its audited sources/binaries were not modified. Aggregate build/execution parallelism was limited to four.
