# Native three-dimensional domain EXX

Approved direction: keep SALMON wavefunctions, MLWF support and exchange action in the original n1 x n2 x n3 spatial partition. Use locality to reduce pairs and communications; do not require a y/z pencil partition or replicate the global box.

Rejected: forcing the simulation mesh to the FFT input layout; globally replicating a complete pair-density box.

Implementation: local ownership is the Cartesian block (global grid, local extent, origin). For extended pairs, transform each axis through its existing axis communicator in bounded line batches: send segments to distributed line owners, FFT owned lines, return segments immediately to the original block. Neither wavefunctions nor persistent spectral arrays change ownership. Temporary storage is local-sized plus line batching padding, not a replicated global box. For compact pairs, construct the subgroup touching source support or required discrete-kernel displacements; compact arrays exist only in that subgroup. Preserve the existing discrete HSE/Coulomb kernel including G=0 and compact padding.

Verification: analytic Fourier modes and round trips on 8x1x1/2x2x2/1x2x4, nontrivial Cartesian coordinates, MPI and serial builds, HSE and Coulomb actions against the existing numerical reference, localized supports crossing periodic/rank boundaries. Run tests before resuming long laser measurements. Keep old 2D test interfaces compatible while migrating production to 3D.
