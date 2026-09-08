# DG Nonlocal H Block Design

## Goal

Replace the dense fragment-fragment nonlocal pseudopotential cache `H_nl_cache` in the RT hot path with a direct block representation so the nonlocal `H0` contribution can be applied through row-owner-local block lists without keeping an always-live dense fragment matrix.

## Current Problem

The DG RT `H0` path now passes the local Hamiltonian block apply, but still adds the nonlocal PP contribution through a dense fragment matrix multiply:

- `matmul(dg_frag%H_nl_cache(1:n_frag,1:n_frag,ispin), coef_all(1:n_frag,:))`

This has two costs:

- compute cost: every rank performs a dense fragment-fragment matmul after the block path
- memory cost: `H_nl_cache` is stored as a dense `n_mat_max x n_mat_max x nspin` array

That defeats the row-owner pruning already added for `H_mat_blocks`.

## Chosen Approach

Use a direct block build for the nonlocal PP matrix.

The new design is:

- add `H_nl_blocks(:)` and `H_nl_block_map(:,:)` to the DG fragment state
- use the same fragment-pair topology as `H_mat_blocks`
- build nonlocal PP contributions directly into `H_nl_blocks`
- apply the nonlocal contribution through `apply_matrix_blocks_batch(..., block_ids=...)`
- reuse the existing row-owner-local block list concept with a dedicated local list for nonlocal blocks if needed

Dense `H_nl_cache` should no longer be required for the normal DG fragment-only RT path once the block path is complete. A temporary dense fallback may remain only where mixed/dense paths still need it.

## Data Model

Add to `s_dg_fragment_rt`:

- `type(matrix_block_info), allocatable :: H_nl_blocks(:)`
- `integer, allocatable :: H_nl_block_map(:,:)`
- `integer :: n_H_nl_blocks = 0`
- optionally `integer, allocatable :: H_nl_local_block_ids(:)` if we keep local lists per operator instead of reusing `H_local_block_ids`

Preferred topology rule:

- `H_nl_blocks` should use the same fragment-pair block topology as `H_mat_blocks`
- this keeps ownership and local block pruning consistent
- the simplest safe rule is to initialize with the same `init_matrix_blocks(...)`

## Build Path

Today `build_nonlocal_pp_matrix_A` writes into a dense output array.

The direct-block design changes that flow:

1. initialize zeroed `H_nl_blocks`
2. loop over local fragment contributions exactly as today
3. for each `(ifrag_row, ifrag_col, io, jo, ispin)` contribution:
   - find the corresponding block id with `H_nl_block_map`
   - accumulate directly into `H_nl_blocks(iblk)%val(io, jo, ispin)`
4. reduce block data across MPI with `reduce_matrix_blocks`

This keeps the communication model aligned with the existing Hamiltonian block path.

## Apply Path

In `calculate_time_derivative`, the fragment-only `H0` path becomes:

1. apply `H_mat_blocks` through the owner-local block list
2. apply `H_nl_blocks` through the same owner-local block list shape
3. add the `A^2/2` diagonal term
4. multiply by `-i`

This removes the dense nonlocal matmul from the RT block path.

## Observable Path

The same nonlocal operator also appears in observables.

That path should be updated to use `H_nl_blocks` so we do not keep one dense nonlocal path in RT propagation and another in observables. The goal is consistent operator storage and application for both:

- RT derivative
- observables using fragment coefficients

## Memory Policy

The main reason for choosing direct block build is memory.

So the intended policy is:

- fragment-only DG RT: use `H_nl_blocks` as the primary storage
- dense `H_nl_cache`: avoid allocating unless a dense-only path explicitly requires it

If mixed basis or another fallback still needs dense data, allocate it lazily and treat it as a fallback representation, not the default runtime state.

## Risks

### Block topology mismatch

If nonlocal PP creates couplings outside the current fragment-pair block topology, direct block build would silently drop terms.

Mitigation:

- verify that nonlocal PP couplings are limited to the same fragment-neighbor pattern already used for `H_mat_blocks`
- add a debug assertion that every `(ifrag_row, ifrag_col)` nonlocal contribution finds a valid block id

### Rebuild timing

`H_nl` depends on `A(t)`, so it must be rebuilt whenever the current nonlocal cache logic would rebuild.

Mitigation:

- keep rebuild policy inside `ensure_nonlocal_pp_matrix_A`
- treat block rebuild as the canonical nonlocal refresh path

### Mixed and dense fallbacks

Some paths still depend on dense matrices.

Mitigation:

- keep dense fallback only where strictly needed
- gate it behind existing mixed/dense conditions

## Verification Strategy

The first checkpoint is functional, not performance-perfect:

- build succeeds
- RT logs advance past the current stop point after `after-apply-block-h`
- no new ownership or block-map errors appear

Then compare:

- `before-h0 -> after-h0` progress
- memory footprint trends
- whether the next bottleneck moves into overlap solve or another stage
