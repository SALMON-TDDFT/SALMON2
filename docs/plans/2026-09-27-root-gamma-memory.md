# Root Gamma localization memory

Goal: reduce simultaneous root link storage without changing the MV objective, tolerance, or Jacobi pair update.

Design: factor the MV functional evaluation on already-rotated links from the existing gauge functional. Add a consuming Gamma entry point which rotates the six input links in place, evaluates that working representation, and applies the existing common Jacobi sweeps. Keep the non-consuming API as the reference, using original-link reconstruction for gradient checks. The consuming API has no raw-link copy or separate six-link evaluator array. Snapshot is written before input links are consumed.

Also reduce root lifetimes: form the links only after pivoted seed construction, write the coefficient prefix of the initial snapshot and free gathered coefficients before MV sweeps; append the remaining same-format records once localization succeeds. Release temporary U/Gram as previously. This avoids QR and links overlapping and avoids keeping gathered coefficients through MV. QR gather itself remains a production limitation.

Tests: RED missing consuming API; compare original/in-place on commuting and noncommuting links, nonidentity initial U, independent original-link gradient and unitarity. Probe root process peak RSS at fixed N with maxiter0 (construction/evaluation only, not convergence benchmark). Native MPI2/4 and C128 before/after compare currents/density/energy plus initial snapshots. Static review. Do not claim 8³/10³ production readiness without all-phase memory validation.
