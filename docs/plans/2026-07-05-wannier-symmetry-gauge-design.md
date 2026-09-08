# Wannier Symmetry Gauge Design

## Goal

Make the DG mixed-Z Wannier basis symmetry-consistent inside each DG fragment without hard-coding C, Si, sp3 bonds, or diamond-specific inversion pairs.

## Motivation

The C64 laser tests showed that the mixed-Z `WW` position block can generate nonphysical even-order HHG. Disabling the `WW` block removes the asymmetry, which means the issue is not the external field itself but the gauge/branch used for the Wannier position block. The right fix is not to delete `WW`, but to make the Wannier gauge respect the symmetries that close inside each DG fragment. Global crystal symmetry is only relevant when its restriction maps the fragment-owned atoms and Wannier functions back into the same fragment-local block.

## Design

Add a general fragment-local symmetry-gauge post-processing layer after Wannier90 output is imported and before SALMON writes `wannier90_global_basis.bin` and `wannier_flux_eigen_seed.bin`.

The layer has four parts:

1. Detect spatial symmetries from each fragment's local atomic environment.
   Candidate operations are affine operations `{R | tau}` in the fragment-local coordinate system. `R` is an orthogonal operation compatible with the local metric, and `tau` is accepted only if it maps atoms owned by that fragment to atoms of the same species inside the same fragment. The global crystal may have more symmetry than a fragment, but DG only requires the symmetry that closes inside the fragment-local block.

2. Build the action of each accepted operation on Wannier centers owned by the fragment.
   For each owner fragment, and for each center `r_i` in that fragment, find `g r_i = r_j` within the same owner block. This gives a fragment-local permutation and residual error. If an operation maps a center outside the owner block, the operation is not used for that fragment.

3. Build the unitary representation of the symmetry in each owner-fragment Wannier subspace.
   This must not be inferred from centers alone. For each operation, SALMON reconstructs the real-space Wannier functions from the DC wavefunction seed and the Wannier transform, applies the spatial operation on the real-space grid, and computes the overlap matrix with the original Wannier functions. A polar decomposition gives the nearest unitary representation `D_g`.

4. Rotate the Wannier gauge consistently.
   The rotation is block diagonal in owner-fragment space. The same block-diagonal unitary is applied to:
   - `v_matrix` from Wannier90 checkpoint data
   - Wannier centers
   - `AA_R` / position matrix
   - `wannier_flux_eigen_seed.bin`
   - owner/cluster metadata after centers are updated

The initial production target is fragment-local inversion symmetry because it directly controls even-order HHG in the DG local response of centrosymmetric fragment blocks. The infrastructure should represent local symmetry operations generally so rotations and mirrors that close inside a fragment can be added without changing the data model.

## Non-Goals

- Do not hard-code diamond, silicon, carbon, sp3 hybrids, atom-index pair tables, or global space-group operations.
- Do not rely on `WW=0` as the final physics route. That switch remains a diagnostic/control path.
- Do not modify Wannier90 internals or require SALMON to call Wannier90 during RT.
- Do not force symmetries that map one fragment to a different fragment. Inter-fragment symmetry equivalence can be diagnosed later, but the gauge correction is local to each owner block.
- Do not require the full supercell to be centrosymmetric. The only required test is fragment-local closure.

## Verification

Use C64 as the first regression because the current failure is clear:

- Before symmetry gauge: at least one owner-fragment WF block fails local inversion-pair closure and `include_ww=y` gives strong even-order HHG.
- After symmetry gauge: center and representation residuals should fall inside each fragment-local block, `include_ww=y` should keep H2/H1 and H4/H1 near the `include_ww=n` control, and H1 should remain comparable.

Then repeat on Si64 impulse and laser cases. A successful route must not depend on material-specific code paths.
