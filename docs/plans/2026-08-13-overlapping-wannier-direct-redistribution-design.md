# Overlapping-Wannier Direct Redistribution Design

## Goal

Remove the post-Wannier replicated core value/gradient copies and the intermediate orbital-owned full-grid value array without changing the final center-fragment data.

## Design

The final buffer arrays remain the canonical spatial storage. Routines that need core data receive the buffer arrays plus the core-to-buffer position map and stream bounded core tiles instead of requiring `ow_core_values` and `ow_core_gradients` copies.

The existing two-stage spatial-to-orbital transpose and orbital-to-center redistribution is replaced in production by one fused routine. It streams bounded orbital batches from the local buffer/core view, routes each orbital directly to its final center-owner rank with checked `MPI_Alltoallv`, and returns only the center-local result. No `owned_values(local_orbitals,global_core)` intermediate is retained.

The old primitives remain available for focused compatibility tests. Production route and memory-lifetime checks prohibit the intermediate arrays.

## Safety and verification

- Preserve canonical orbital numbering and global physical-ID ordering.
- Check replicated dimensions, row ownership, MPI counts/displacements, allocation status, and finite output collectively.
- Compare the fused result with the established two-stage result on MPI 1/2/4/8.
- Compile the full MPI/EigenExa/Wannier90 configuration.
- Do not restart or alter the running Si64 observation binary; it uses an immutable source snapshot.
