# Hybrid LCFO Low-Energy Symmetry Design

## Goal

Allow the Hybrid LCFO generalized eigensystem to use the complete retained
WF+PW variational basis while requiring symmetry only of the physically used
low-energy eigenspace: the occupied 128 states plus the existing LCFO target
empty-state allowance, currently 384 states in total.  If the upper boundary
cuts a degenerate cluster, include the complete cluster.

## Motivation

The current pre-solve gate constructs a projected symmetry action on all 1560
retained WF+PW basis functions and requires that action to be unitary in the
retained metric.  That is a valid test for invariance of the entire retained
basis span, but it is stronger than the production requirement.  The PW
component is intended to absorb symmetry content missing from the localized
WF component, and high-energy directions outside the LCFO target window do
not need to form a closed representation.

The present gate therefore can reject a basis that is sufficient to produce a
symmetry-preserving ground state and the required low-energy excitations.  It
also does not by itself prove that the final occupied projector or density is
symmetric.

## Chosen Approach

Keep the all-basis closure calculation as a diagnostic receipt, but never use
its failure alone to stop the continuation.  Run the LCFO generalized
eigensystem, form the final low-energy target projector, and validate that
projector under every retained physical symmetry operation.  Validate the
occupied projector and final real-space density independently.

Do not explicitly symmetrize the Hamiltonian, coefficients, or density.  The
first production result must show whether the WF+PW variational space is
sufficient naturally.  Explicit symmetry projection remains a later remedy
only if the post-solve physical validation fails.

## Target Window

Use the existing LCFO target rank, currently 384 states, as the minimum target
window.  This comprises 128 occupied states and 256 low empty states for the
Si64 acceptance case.  Inspect the final generalized eigenvalues at the upper
boundary and extend the window through every state in the same numerical
degeneracy cluster.  A fixed count must never split a degenerate multiplet.

The occupied window is selected from the final occupations rather than from
individual eigenvector labels.  Symmetry is assessed on projectors because
an eigensolver may return any unitary rotation within a degenerate subspace.

## Symmetry Validation

For each physical symmetry operation, apply the existing real-space pencil
map to the reconstructed target states and measure the component outside the
target eigenspace in the retained metric.  Equivalently, for target
coefficients `C_t` normalized by `C_t^H S C_t = I`, construct the induced
target action and require its projector leakage and metric-unitarity defect to
remain below the configured symmetry tolerance.

Report at least:

- requested and degeneracy-extended target ranks;
- the operation with the largest target-projector leakage;
- the maximum target metric-unitarity defect;
- the corresponding occupied-projector defects;
- the final real-space density symmetry defect.

The production acceptance gate fails only when the occupied projector, the
degeneracy-complete target projector, or the final density exceeds tolerance.
Failure of the full 1560-dimensional basis closure is retained as a warning
and diagnostic value.

## Control Flow

1. Build the retained WF+PW basis, metric, and operator rows as today.
2. Attempt the full-basis representation analysis and emit its defect without
   stopping when the span is not closed.
3. Continue through the LCFO generalized eigensystem and final SCF refresh.
4. Select the occupied and low-energy target windows, extending the latter at
   a degenerate upper boundary.
5. Measure occupied-projector, target-projector, and density symmetry.
6. Publish the measurements in the final Hybrid GS acceptance receipt.
7. Stop before checkpoint publication only if a physical post-solve symmetry
   condition fails.

## Error Handling

Structural failures remain fatal: invalid dimensions, incomplete distributed
rows, singular retained metric, failed collectives, or non-finite values.
Only a finite, well-formed full-basis closure defect above tolerance is
downgraded to a diagnostic.  Post-solve occupied, target, and density defects
remain fatal acceptance failures.

## Testing

Add a focused MPI regression where the complete retained basis is not symmetry
closed but the selected low-energy eigenspace is closed.  It must continue and
pass.  Add failures for target-window leakage and for splitting a degenerate
upper cluster.  Protect the existing retained-basis structural-error tests,
variational payload tests, continuation route contract, and final Si64 8-rank
GS-to-RT acceptance calculation.
