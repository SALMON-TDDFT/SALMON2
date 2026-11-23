# SCF nonconvergence review pause (2026-11-21)

User decision: pause after terminal unconverged GS, warn and ask to inspect residual evolution.
This supersedes warning-only automatic continuation. Preserve native update order,
thresholds and normal GS saving; no additional GS collection for validation.

- DFT.run retains complete manifest/output before raising SCFReviewRequired.
- Exception.result retains the batch result/live native owner, report_directory
  points to CSV/optional PNG/metadata. Report failure cannot bypass review.
- Batch diagnostics parse actual native printed rho_dne history; older output
  without history is explicitly labeled final-status-only, never invented.
- Native run/continue use existing native residual history and subgroup agreements.
  No extra per-step application-wide collective. Intentional paused budgets are valid.
- Saved result inspection and warm starts remain available; batch GS→band/TDDFT
  checks known convergence before copying/initializing downstream data.
- Explicit allow_unconverged=True only after human review; never auto-retry it.
- ASE inherits pause without publishing unconverged/old force results.

Validation: synthetic regression includes converged success, cap pause, full
residual CSV, optional plotting failure, report disk failure, native owner survival,
partial budgets, reviewed override and prelaunch boolean rejection. Source suite
672 tests, 1 optional skip. Actual audited-executable Si8 8-iteration run collects
no GS restart: six input/numerical output files match the preceding wheel bytewise;
CSV has all 8 iterations and PNG is saved. This is local control-flow regression,
not a physical convergence or multi-node MPI claim. Record is in the canonical
calc2026/salmonpy-scf-review-20261121/record.json; cloud sync not asserted.

Final sdist→wheel→fresh environment verification: 672 tests, 1 optional skip;
installed wheel real Si8 pause/8-point CSV/PNG and six-file byte comparison passed.
Logs and distribution-validation.json are retained alongside the source record.
