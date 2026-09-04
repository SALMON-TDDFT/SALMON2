# Saved Si64 DC core-density handoff audit

## Scope

Read-only audit of the existing conventional DC checkpoint; no DC, Wannier90,
SCF or material calculation was rerun. This is a necessary-condition check
for C3/C5 admission, not an observed failure of a newly connected main route.

Input directory (relative to the user-fixed worktree):
`verification-si64-localization-first-lcfo-20260902-final/dc-seed`.
Publication: `8CBCD7D8424BFE03`, eight shards, rank 0..7 maps to fragment 1..8.
Manifest SHA-256:
`565f54dfeff1972021682ecc29fc81b3a972685f2a7253fcbc09b2901a0d7549`.
This audit is not a replacement for the production checkpoint integrity and
exact MPI/rank--fragment compatibility validator.

## Method

Decode the existing stream layout in `write_shard_raw` in
`src/gs/dc/dg_dc_seed_checkpoint.f90`: 32-byte magic, four int32, four int64,
38 int32 array bounds, two float64 and one int32 (252-byte header).
Read the five float64 arrays in Fortran order using those bounds; require
the calculated final byte offset to equal the shard file size. Observed rwf
shape is `(28,29,28,1,400,1,1)`; do not infer this shape from state-count defaults.

The saved input has grid 32^3, cell side 20.52 bohr, uniform 2x2x2 fragments,
and fragment optimization disabled. The density assembly in `dcdft.f90` combines local grid
indices 1..16 in each direction. The integration weight `(20.52/32)^3` equals
the manifest density weight, 0.263683001953125. Use this same core and weight
to evaluate `G = Psi_core^T W Psi_core`, `rho = sum_i f_i Psi_i^2`, and
`N_core = sum_i f_i G_ii`. Read f from the stored rocc array without modifying it.
The occupied range below uses f > 1e-8; electron totals include all f unless
explicitly labeled selected. No new state-count cap is applied.

## Results

All eight fragments agree to the displayed precision:

| Quantity | Per fragment |
| --- | ---: |
| Orbitals with f > 1e-8 | 88 |
| Stored-grid occupied orbital norms | 1 within 1.1e-15 |
| Occupied core norm range | 0.0147253 .. 0.4806794 |
| Core electrons from orbitals | 32.00000000002 .. 32.00000000007 |
| Core electrons from saved density | 31.99999999879 |
| Sum of occupations for f > 1e-8 | 168.08551131164 .. 168.08551131172 |
| Maximum off-diagonal core overlap (all 400 orbitals) | 0.3144073 |

Total core electrons from orbitals: 256.0000000003457.
Total from saved density: 255.99999999031618.
Manifest expected electrons: 256.
Sum of occupations selected at f > 1e-8: 1344.6840904934284 over eight fragments.

## Consequence and checkpoint

The current initializer enforces `C^dagger S C = I` and reconstructs density
using the unchanged original occupations for the selected seed IDs. If the
occupied seed set above is retained and orthonormalization succeeds, its
integrated core density is therefore sum(f), about 168.09 per fragment, not
32. Exact WF+PW reconstruction before initialization cannot remove this
discrepancy. Additional unoccupied guards do not solve it. Failure earlier
due to insufficient span or dependence is also possible; this audit does not
claim to have run that combined path on the checkpoint.

The appreciable off-diagonal core overlaps additionally rule out assuming
that independent orbital rescaling generally preserves the density under a
subsequent orthogonalizing rotation. A density-preserving occupation-matrix
transformation, or a separately approved change to the starting-density
policy, needs explicit design. Such a change must also revisit post-initializer
support comparisons, which currently compare individual raw seed columns.

The approved C3 plan explicitly requires a pause when representative DC data
needs this treatment. C5/C6/main promotion is paused pending that decision.
Independent read-only review reproduced the eight-shard norms, electron counts,
header extents and physical core mapping and confirmed this stop condition.
Do not raise cutoff, rescale density, drop occupied seeds or relax the gate to
conceal the discrepancy. Existing source changes and verification files remain
untouched by this audit.
