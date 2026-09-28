# LCFO halo tile distribution

Continue approved memory/compute distribution. Keep inter-fragment representative exchange, then replace whole-halo broadcast with point-to-point delivery of each receiving rank's halo/domain intersection. Skip ranks without the current k point or ket columns and empty intersections. Copy noncontiguous outgoing slices into explicit temporary tiles, wait before release, and retain existing small projected-matrix reductions. Root still holds full inter-fragment halos; do not claim that cost removed.

Verify local tile sizes, HSE spatial/empty-core/reconstructed RT cases, Si combined k/orbital references, PBEh water and pulse, ON/OFF builds, independent review.

## Scope clarification and validation

This change concerns assembling the global LCFO Hamiltonian, not distributing its eigensolver. It also does not distribute MLWF/ACE orbital columns in native SCF/RT.

- HSE regression covers tile-size scaling, inactive k-point ranks, empty-core layouts, spectra and reconstructed RT.
- Si reference checks pass for four k points/192 eigenvalues with k distribution and combined k/orbital distribution.
- PBEh water mesh propagation and spatial pulse tests pass. ON/OFF builds pass.
- Independent review found no blocking issue; explicit empty orbital-range guard added.
- Root inter-fragment send/receive halo buffers, global eigensolver matrices and EXX orbital-column distribution remain unchanged.
- Final HSE suite: all eight tests pass after the ownership guard. Diff whitespace check passes.
