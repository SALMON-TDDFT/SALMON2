# Synthetic LSDA versus corrected noncollinear MGGA trajectory comparison

This is an operator/trajectory workflow test using deterministic non-ground-state
spinors. It contains no material pseudopotential, no ionic SOC, no physical ground
state and no collected GS data. Both functionals use the **same initial state**,
so this explicitly compares functional quenches, not response from each
functional's self-consistent ground state. No physical demagnetization or material
accuracy should be inferred from these trajectories.

The program uses the actual SALMON `init_xc`, `exchange_correlation`, `hpsi`,
Taylor propagator and `write_rt_spin`. LSDA follows the existing noncollinear
density rotation and XC matrix path; MGGA uses its full local, tau and current
payload. Both use the same midpoint predictor/corrector procedure, 6th-order
Taylor expansion, time step, grid, occupations and spin-independent scalar field.
The stream includes time zero and 20 propagated samples. Spin and norm
conservation are checked before accepting outputs. This tests the components
called by RT propagation; it is not a complete input-driven `main_tddft` or
restart test.

Example (adjust library flags for the existing build):

```sh
MPLCONFIGDIR=/private/tmp/nc-mgga23-matplotlib \
python3 testsuites/unit_nc_mgga23_comparison/run.py \
  --build-dir /private/tmp/salmon-ncmgga23-serial \
  --checker /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/tools/check_forbidden_eigensolvers.py \
  --libs /opt/homebrew/opt/openblas/lib/libopenblas.dylib \
  --output /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/nc-mgga23-validation-20261024/synthetic-comparison \
  --plot
```

The destination must not already exist. The runner audits the full source and
all linked SALMON application objects plus its new test object before running.
It writes native spin/norm files, validation bounds, executable SHA256, method
metadata and a run log, then invokes `tools/compare_spin.py` to produce JSON/CSV
and optionally PNG plots. Files at the canonical OneDrive path demonstrate local
storage only; cloud synchronization must be confirmed separately.
