# Distribute complex LCFO basis construction

Continue the approved distribution work. Construct restricted fragment orbitals only on each rank's intersection with the core. Gather columns only across k/orbital peers sharing the same grid; reduce overlap matrices and Gram-Schmidt inner products across spatial peers. Store local basis rows on ordinary ranks; retain a full core basis only on the fragment representative for existing binary output/halo sends. Gather one basis column at a time, avoiding all-band global scratch. Retain existing received halo broadcast for now.

Validate basis orthogonality, HSE spatial DC eigenvalues/reconstructed RT, legacy multi-k and multi-orbital Si references, workspace sizes, ON/OFF build and PBEh regression. Root basis and one-column gather scratch remain explicit limitations.

## Results

- HSE suite: eight tests pass, including local/stored basis-size assertions and a rotated two-fragment case with zero-core-domain ranks. LCFO spectra and reconstructed RT retain parity.
- Si: four k points and 192 eigenvalues match existing references on four ranks (k distribution) and eight ranks (combined k/orbital distribution).
- ON/OFF builds pass; independent reviewer found no blocking issue in communicator scopes, per-fragment loop counts, root-only gather or zero-domain handling.
- Root full-core basis, one full-core send column per rank, root receive column, halo replication and dense matrices remain explicitly documented.
- Final PBEh Ehrenfest regression: all 16 tests pass. Combined HSE/PBEh suite covers 24 tests. Diff whitespace check passes.
