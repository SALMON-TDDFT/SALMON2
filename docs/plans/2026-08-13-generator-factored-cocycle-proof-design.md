# Generator-Factored Point-Cogroup Proof Design

## Goal

Reduce the pre-RT Si64 point-cogroup proof from all 48 squared dense retained-space products to a generator-complete proof, without weakening the affine cocycle contract.

## Mathematical contract

The affine catalog defines

```text
r_left r_right = t_cocycle(left,right) r_product(left,right).
```

Spatial maps are pullbacks. Their inexpensive integer action is therefore checked for every ordered point pair. This establishes the supplied product table, cocycle orientation, and map payload exactly on the full spatial grid.

For the retained-space matrices, choose a deterministic generating set of the point cogroup. Verify the cocycle matrix relation for every generator and every point element, in both generator-left and generator-right directions. Verify by closure traversal that these generators reach all point elements. These relations recursively determine every representative matrix from the identity, so all-pair dense products are redundant once the full integer action has passed.

## Implementation

- Reuse `select_dg_group_generators` on `point_product`.
- Keep the current all-pair integer map composition gate.
- Assemble and multiply retained-space matrices only when either operand is a selected generator.
- Count checked ordered matrix pairs and publish generator count and checked-pair count in the diagnostic receipt.
- Keep the same tolerance, cocycle orientation, row ownership, MPI consensus, and conservative workspace accounting.
- Do not retain all 48 dense matrices. Continue streaming one left/right/expected operation at a time.

For a 48-element point cogroup with 3 generators, the dense checks fall from 2304 to at most 282 unique ordered pairs, while persistent memory remains unchanged.

## Failure handling

Reject collectively if generator selection fails, does not generate every point element, metadata disagree, or any checked generator relation exceeds tolerance. Allocation and MPI failures continue through the existing collective cleanup path.

## Tests

- Add a point group with elements not all selected as generators.
- Corrupt a generator relation and require rejection.
- Corrupt a non-generator spatial map/cocycle relation and require rejection by the retained all-pair integer gate.
- Assert the dense checked-pair count is below `npoint*npoint` and is rank independent.
- Run construction MPI tests on 1/2/4/8 ranks, W90 MPI tests, route checks, and a clean overlay build.
- Only then rerun ideal Si64 MPI 8 / OMP 1 from a persistent result directory.

