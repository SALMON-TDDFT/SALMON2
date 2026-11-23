# Outer MPI performance validation — 2026-11-23

{
  "schema": "salmonpy.outer-performance-validation/1",
  "ranks": [
    3,
    5
  ],
  "steps": 24,
  "source_and_wheel": true,
  "checked_performance": "array_equal",
  "explicit_protected_failure_rollback_retry": true,
  "tests": {
    "total": 707,
    "passed": 706,
    "skipped": 1
  },
  "wheel_sha256": "bcce4a924ae0f75defeeee4a2bae7b9fedd90a1654b418d705e2cb2eef0785bf",
  "GS_collected": false,
  "scope": "Multiscale1D plus two collective custom A-to-j groups; synthetic test, not material convergence, native EH/electronic coupling or speed benchmark"
}

Records: calc2026/salmonpy-outer-performance-20261123. Local storage verified; cloud synchronization not confirmed. Necessary exception coordination and physical reductions remain. Performance mode does not promise automatic rollback; explicit protected steps retain their existing contract.
