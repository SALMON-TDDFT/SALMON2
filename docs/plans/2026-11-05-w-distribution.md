# W preparation efficiency

Approved scope: preserve RPA, microscopic cutoff, occupation factors, eta/q convergence guards and HXR1/HXT1/HXP1 layouts. No changes to production jobs or GS collection.

1. Exact separable finite-grid box factors, bounded phase/snapshot reuse. Compare with explicit real-grid sums and uncached solvers.
2. Keep target-k Sternheimer distribution; reduce each q to an owner, assigning q/-q pairs together. Invert dielectric matrices and project boxes on q owners. Preserve the original root-return API for reference tests.
3. Gather only metallic-head stencils and small diagnostics. Transpose q-owned box matrices to pair owners; slabwise interpolation and normalized FFT; bounded scratch/memmap, reciprocal/finite validation and streaming output with unchanged layout.
4. Serial/MPI1/2/3 comparisons against preserved baseline, including independent q meshes, integer/fractional snapshots, empty owners and rejection. Benchmark isolated box factors and complete synthetic W pipeline, preserving failed artifacts. Register tests and document limits/measurements.
