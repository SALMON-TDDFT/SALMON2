# Upstream correctness backport design

## Goal

Keep the SALMON2_RTDG history and DG/Wannier implementation intact while
backporting applicable correctness fixes from SALMON-TDDFT/SALMON2
`develop-2.0.0`.

## Scope

The first batch includes bug fixes for DC input/restart handling, RT restart,
non-orthogonal cells, NLCC, TBmBJ, and automatic MPI decomposition. New
features and dependencies, including the SLEPc DC-LCFO solver, are excluded.

## Integration strategy

Upstream and the local repository have unrelated Git histories, so commits are
not merged wholesale. Each upstream patch is compared with the corresponding
local implementation, adapted only where the local API differs, and recorded
with its upstream commit ID. Existing Task 9 changes are stashed and restored
around the backport commit.

## Verification

Every behavioral change receives a focused RED test where the local test
infrastructure permits it. The resulting tree is built from a clean parent
prerequisite overlay, relevant focused tests are run, and the final diff is
reviewed for DG/Wannier and normal-DC regressions.
