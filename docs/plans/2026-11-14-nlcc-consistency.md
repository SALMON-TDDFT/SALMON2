# NLCC XC consistency and force implementation plan

**Goal:** Correct the approved unpolarized scalar NC NLCC path in the external SALMONpy runtime patch, then verify discrete derivatives and force/energy consistency. Existing builds and jobs remain unchanged. The Calculator force guard stays until validated runtime capability can be established.

**Architecture:** Supply total density gradients to Libxc and weight its energy density by total density. Use the derivative of the existing piecewise-linear sampled core density for the explicit NLCC force, with outer-atom OpenMP. Reduce its grid-only contribution over the spatial communicator, without multiplying by k/orbital replicas. Add the contribution before force symmetry.

1. Prepare fresh external source from the pinned upstream plus current compatibility patch. Preserve all source/configure/incremental/object audits. Construct native XC variational tests and verify the existing code fails the core-density energy derivative test.
2. Correct XC weighting/gradient dataflow. Keep no-NLCC arithmetic unchanged. Compile in a new directory, audit linked application objects, record hashes, test energy variations and builtin-PZ/Libxc-PZ equivalence with NLCC.
3. Add explicit core force using the same images, coordinates, radial interpolation, and cutoff as calc_nlcc. Test frozen-density core displacements, skew cells, outer OMP and MPI replica counting. Preserve unsupported spin-Libxc guard and do not claim meta-GGA force support.
4. Regenerate runtime compatibility patch against pinned upstream, update its provenance hash, and verify fresh preparation reproduces the modified source. No full solver source is imported into the repository.
5. Compare original/XC-fixed/final scalar-only Si diagnostics and no-NLCC/KY regressions. Save requested scalar results only to canonical calc2026. Ephemeral native eigen/k/info outputs are deleted with scratch. Run Python/diagnostic tests and independent review.
