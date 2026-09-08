# LCM occupation plumbing cleanup design

## Scope

Remove occupation-weight data flow made obsolete by the sharp all-electron LCM
contract, without changing projector construction, multiplicities, parallel
ownership, or output.

## Design

Keep `occ_w` because the sharp-occupation validator consumes the actual SALMON
occupations.  Remove `local_occ_w` and the occupation arrays from
`build_occ_distribution_cache`.  Delete the unused `local_occ_weight` helper
and the `local_occ_index` and `local_occ_global_io` helpers that have no other
callers.  Retain `local_occ_position` as the fallback owner lookup.  Apply the
same cleanup to scalar and SOI modules.

## Verification

Add a source-contract regression check that rejects the dead symbols and
weight arguments while requiring the sharp validator to retain `occ_w`.  Run
the complete focused LCM test set and rebuild the MPI/EigenExa target before
committing and pushing to both configured remotes.

