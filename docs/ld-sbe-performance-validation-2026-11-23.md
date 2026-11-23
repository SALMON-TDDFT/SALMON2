# Standalone LD/SBE performance mode — 2026-11-23

`execution_mode="performance"` is fixed per state; `checked` remains default.
Both single-step and batch entry points support it. Physical Fortran update
arithmetic is unchanged. LD omits finite scans and pole backup transfers. SBE
allocates only the evolved density, omitting trial density and per-step scans.
Initial/library/shape and explicit snapshot/restore/output checks remain.
There is no automatic failed-update rollback guarantee in performance mode.

Synthetic scalar/tensor LD P/J and SBE current/density were exactly equal to
checked mode. SBE covered correction orders 0–3 and nonlocal correction on/off.
Interleaved independent states and explicit restore/continuation also matched.
NaN propagation in performance mode was exercised; explicit valid restoration
still succeeded. This is numerical/ownership validation, not a physical
convergence study or a performance benchmark.

Source and separately installed sdist-built wheel suites each ran 699 tests:
698 passed, one optional test skipped. Installed scalar/tensor LD and SBE
examples produced identical numerical stdout in checked/performance modes.
The C/Fortran application objects passed the forbidden-eigensolver audit.
No GS calculation, remote CI or heavy runtime rebuild was needed.

日英LD/SBEガイド、AI向けJSON、同梱サンプルと作業リストを更新済み。
Records and wheel: `calc2026/salmonpy-ld-sbe-performance-20261123`.
Local OneDrive storage does not confirm cloud synchronization.

Wheel SHA256: `fe930eeac6251843908e7bd6f65b40f11350bad2d6cb47ce0d592691ba1598ed`.
