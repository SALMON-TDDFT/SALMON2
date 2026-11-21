# Coalesce outer coupling admission

Combine the completed field operation's status with ResponseSession input admission, before any A/packet scatter or material callback. Keep the pre-field clock/usability vote and the post-material vote. Do not change numerical kernels, native physical communication or transactional snapshot/rollback ordering.

Packet coupling without H can combine field/sample errors, packet construction and response preflight into one decision. With H, retain the explicit post-sampling vote before entering the H collective, then combine H errors, packet construction and response admission. Recipient packet validation remains separate.

Use an explicit private upstream-error argument, not mutable process/session flags. Standalone response APIs keep their signatures. Custom public response-step overrides must remain honored via the prior separate agreement path.
Negotiate the dispatch once across every group before stepping; if any rank has a
custom entry, use the legacy path everywhere. Stable original callable identities
cover subclass, instance and class-level overrides. Freeze response-entry dispatch
when constructing the coupling; do not decide vote counts independently per rank.

Expected hierarchical agreement counts per successful ordinary step: A 4→3, packet without H 7→4, packet with H 7→5. These counts exclude native subgroup votes, point reductions, setup, checkpoints and protected-step transaction votes. Each hierarchical agreement still uses group/leader/group collectives.

Test a failing count contract first, then successful counts, single-rank field/sample/H/backend failures, protected original-state/alias/side-effect preservation and retry, legacy/custom response overrides, unchanged native currents/fields and direct resume. Source and installed-wheel suites, Japanese/English manuals, AI JSON and checklist must agree. No new GS or manually dispatched heavy CI.

## Verification complete

- Prior spatial-current work committed/pushed as `ea42d24c`; this follow-on remains in the worktree. No Fortran/runtime source change or full-runtime rebuild was needed; each run verifies the unchanged audited library SHA256.
- Red contract measured4 A votes against requested3. Final MPI4 normal A counts3; packet counts7→4 without H and7→5 with H, measured against the preceding installed wheel. Counts exclude setup and native/transaction votes.
- A response subclass, one-rank-only subclass and one-rank class-method replacement preserve public dispatch and globally select the legacy4-vote path. The class-method test initially exposed bypassed dispatch; stable original callable identities fixed it.
- MPI4 native field+LD+two SBE points: late field/SBE failures and deep field overflow preserve all original bytes, aliases and output files; protected retry matches continuous execution with zero replay. Direct protected restore and resume also match.
- Native packet faults: nonleader A/E sampling failure with/without H, quiescent H failure, and leader packet-construction error reject before any material callback. Pre-H sample failure executes no H collective. Ordinary packet step remains nontransactional.
- MPI6 protected LD/SBE/custom state preservation, MPI3 MS1D admission/backend failures and standalone packet regressions pass.
- Final source and sdist-built installed-wheel suites each run668 tests with one optional skip and no failures. Installed wheel additionally passes MPI4 normal/packet counts, class override, protected numerical failure/retry, direct resume and pre/post-H failures. Fourteen A/E/J comparisons are bitwise identical to preceding controls.
- Japanese/English guides, AI JSON and checklist synchronized. Canonical evidence: `calc2026/salmonpy-outer-admission-20261121/record.json` with distributions, scripts and initial failure records. No new GS or manual CI dispatch. No multi-node, GPU or new TDDFT/HSE material-physics claim.

Manual count invocation using the ordinary sample arguments:

```sh
mpiexec -n 4 python tools/test_outer_agreement_mpi.py --expected-votes 3 --report NEW_REPORT.json BUILD RUNTIME BUNDLE NEW_OUTPUT --layers LAYOUT --field-grid 2 1 1 --tddft-ranks 1 --steps 96
mpiexec -n 4 python tools/test_native_protected_step_mpi.py BUILD RUNTIME BUNDLE NEW_OUTPUT --layers LAYOUT --field-grid 2 1 1 --tddft-ranks 1 --steps 96
```

Use an explicit SBE-only layer layout for GS-free validation. The optional count driver accepts `--driver SCRIPT` for packet/custom contract drivers. Existing protected transaction votes and native physical communication are deliberately retained.
