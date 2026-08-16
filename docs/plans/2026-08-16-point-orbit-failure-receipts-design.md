# Point-Orbit Failure Receipts Design

## Goal

Distinguish DC/fragment information loss from an over-strict point-orbit
construction assumption without changing any acceptance threshold.

## Design

`build_point_orbit_blocks` will attach a stage-specific message before each
failure return:

- incomplete candidate/cover construction: candidate count, covered rank, and
  required dimension;
- rank-deficient orbit Gram matrix: minimum eigenvalue and rank threshold;
- point-operation cluster leakage: maximum leakage and center tolerance.

The caller will preserve this detailed message and use the existing generic
message only if no internal stage supplied one.  No arrays, communication,
thresholds, or success-path arithmetic change.

