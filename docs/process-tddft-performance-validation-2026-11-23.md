# Periodic subprocess TDDFT performance mode — 2026-11-23

TDDFTInitialState.new_state, crystal.respond and TDDFTResponseBackend now accept
execution_mode="performance". Checked remains default. A dedicated INIT_PERFORMANCE
and SALMONPY_TDDFT_PERFORMANCE_V1 reply prove child support; old workers fail
instead of silently accepting the request. The child-only environment leaves
host/embedded and molecular paths unchanged. Setup/hash/protocol/capacity and
explicit output validation remain. Python/Fortran input/current finite checks
and automatic failed-worker respawn/replay are omitted. Physical stepping,
requested electron-number output and solver failure handling remain. Failed
performance states close without a retry/previous-state guarantee.

Explicit snapshots still retain history; restore spawns a candidate and replays
that history, checking saved currents. It is NOT direct wavefunction restoration.
Cross-mode restore preserves the receiving mode. No new numerical Python solver.

PZ/TB-mBJ/Libxc-PZ/PBE from source and separately installed sdist-built wheel:
16-step current and measured electron number exactly match checked and previous
validated executable. Both checked/performance states restored the other's step8
snapshot and reproduced currents/electron number9–16. Forced worker termination
closed the performance state with no automatic respawn. An old executable was
rejected during initialization. Installed/source arrays also matched.

To honor the no-GS-collection instruction, only the validation storage copier was
patched to make read links to pre-existing synthetic GS. Production copier is
unchanged. Actual initialization/protocol/stepping/restore use production code.
GS bytes/names remained unchanged; no new GS or electronic disk save. These
regressions do not establish real-material convergence or speed gains.

Git-only source preparation/new executable build and source/application-object
prohibition audits passed. Source and installed-wheel suites: 705 tests each,
704 passed, one optional skip. Wheel code, samples, paired guides and AI JSON
match source. No remote CI dispatched. Molecular/HSE subprocess and outer MPI
performance mode remain separate work. Records/artifacts reside at
calc2026/salmonpy-process-tddft-performance-20261123. Local OneDrive storage does
not confirm cloud synchronization.
Worker SHA256: `0fcc04e75b79ebc094a7610e2825a22c30006b7de2be62bb163148dd5ffc35d1`.
Wheel SHA256: `52f6d506a9d1c00a90284eb94255f60c17aa027eb342d6ede180b025e91eee47`.
