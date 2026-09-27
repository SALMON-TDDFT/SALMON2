# Streamed Gamma seed coefficients

Goal: eliminate simultaneous root full_coeff and its QR transpose without changing global pivoted QR or initial snapshot format.

1. Extend the existing tiled root gather with optional adjoint layout and coefficient stream output. Default API unchanged. Stream each full coefficient column in original order using O(nbasis) scratch, while retaining only the conjugate-transposed global matrix.
2. Add gauge seed pivot and polar-finish stages, preserving LAPACK calls/work sizes. Release the QR matrix before collecting selected original rows.
3. New lcfo_seed module orchestrates gather -> root QR -> status/pivot broadcast -> tiled selected-row reduce -> root SVD. Nonroot stores no global coefficients or selected N² matrix; scratch limited to N*64. Root SVD status is broadcast. Caller opens/writes header and closes snapshot prefix; no full_coeff remains.
4. MPI2/4 regression with empty root and unequal row counts, tile tails, exact streamed data, reference U and invalid/singular fallback. Isolated RSS benchmark in separate processes, same matrices. Freeze previous executable and run native RT and sequential C128 regression.
5. Review, document measured scope and residual root QR/U/ACE limitations, publish notes and portable patch. No 3D production job.
