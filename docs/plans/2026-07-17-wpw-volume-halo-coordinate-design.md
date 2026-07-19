# WPW Piecewise-Wannier Volume Correction Design

## Problem

The production SCF density currently uses local and remote Wannier values on
each uniquely owned fragment core.  The assembled WW metric is nevertheless
`orthonormal_ww`.  In the two-fragment
smoke test, the integrated WW density from local plus halo values is much larger
than the coefficient-space WW norm, while the owned-only diagnostic agrees with
that norm.  The fragment basis is orthonormalized on its own core and represents
a piecewise DG basis function there.  Therefore extending the same W row into a
neighbor core changes its norm and contradicts `orthonormal_ww` even when the
halo coordinates are correct.

## Accepted direction

Keep the frozen `orthonormal_ww` metric and enforce piecewise-Wannier volume
semantics.  On a spatial fragment core, WW, WP, and density volume terms use
only W rows owned by that fragment.  Neighbor W rows couple through canonical
SIPG face traces, not through a second volume representation.  P rows retain
their bounded multi-fragment window support.

## Data-flow contract

At each uniquely owned core point, the W vector has exactly the fragment-owned
rows and the P vector has all bounded support rows.  WW potential updates are
owned-by-owned, WP volume candidates are owned-W by support-P, and density uses
the same vectors as the published overlap.  Remote W data remains available
only to the canonical face scanner.  No density renormalization is permitted.

## Implementation boundary

Add a focused source/linked fixture requiring the production volume accumulator,
WP publication, and SCF density path to use owned W rows only, while the face
scanner continues to use support W rows.  Remove volume-W halo preparation from
the production consumer chain without changing the independently tested halo
provider API.  Do not change the WW metric or P support contract.

## Verification

Run the focused volume-layout fixture, face and halo fixtures, rank-local
quadrature tests, matrix-free SCF test, and a full build.  Then
rerun the two-fragment GS smoke and require the WW real-space norm to agree with
the coefficient-space WW norm within quadrature tolerance, finite charge and
energy through all configured SCF iterations, and no algebra collapse.  Remove
temporary success-path diagnostics after the invariant is covered by tests;
retain concise failure diagnostics that identify the failed stage.
