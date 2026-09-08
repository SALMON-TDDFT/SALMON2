# Defer Periodic Centers Until Character Inversion Design

## Context

The Si64 reference translation-character sector has dimension 12.  Projecting
the supercell periodic-position phase into that one sector loses the
off-diagonal couplings to the other 31 translation characters.  The resulting
point-orbit cluster leakage is 0.702, versus a 0.00316 threshold, despite DC
boundary and symmetry defects near 1e-13.

## Decision

Remove the production-only single-character periodic-position tuple and joint
point-orbit canonicalization.  Keep the W90/LCFO reference-sector anchor,
materialize that anchored reference once, transport its gauge to every other
character with the prepared translation intertwiner, apply Gamma sewing, and
perform the inverse character transform.

Center measurement and affine center-orbit validation remain after inversion,
where the full translation-character space is available.  The standalone
periodic-position primitives and their focused tests remain available; only
their mathematically invalid production placement is removed.

## Memory and Provenance

This removes the point-overlap rows, all 48 small point representations,
position tuple, weighted reference copy, canonical weighted output, and joint
center arrays.  The post-gauge provenance remains bound to the operator anchor,
reference materialization, per-character alignment, Gamma sewing, inverse
orbit payload, and final point-cogroup proof.

