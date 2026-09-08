## DG RT next optimization steps

### Goal

Reduce the two remaining dominant costs after block-pruning of `H0`:

1. overlap solve
2. density projection

### Step 1: overlap solve pruning

Current `solve_overlap_operator_batch` solves all rows in `rhs(:, :)` with column-by-column PCG.
In DG RT, only owner rows are ultimately consumed, so the first reduction is:

- build local overlap row list from `coef_owner`
- build local overlap block ids from `S_mat_prop_blocks` / `S_mat_blocks`
- solve only on owner rows

This should reduce both vector length and each `apply_overlap_operator` cost before any batched-PCG work.

### Step 2: blocked density projection

Current density reconstruction loops over:

- fragment grid point
- occupied state
- basis function

and explicitly forms `psi_val` point by point.

Replace the fragment-basis contribution with blocked projection:

- flatten a fragment grid block into `ngrid_blk`
- build `A(ngrid_blk, nbasis_blk)` from `phi_frag`
- use `B(nbasis_blk, nocc_blk)` from `coef`
- compute `Psi_blk = A * B`
- immediately accumulate row-wise `sum(abs(Psi_blk)**2)` into density

Do not materialize the full `Psi`.
Keep plane-wave contribution separate at first.

### Order

1. overlap row-owner/pruned solve
2. blocked density projection for fragment-basis part
3. later: batched overlap iterations and PW-side projection blocking
