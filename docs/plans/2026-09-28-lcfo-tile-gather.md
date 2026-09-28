# LCFO basis tile gather

Continue the approved basis-memory distribution. Replace the full-core one-column MPI_Reduce with point-to-point transfers from one k/orbital representative per spatial domain. Exchange only small extent/owner metadata collectively. Root receives one actual domain column at a time into a domain-sized buffer and places it in the existing root-only basis. Empty domains send nothing; local root data is copied directly. Preserve output and halo formats.

Verify transfer-workspace diagnostics, zero-domain and combined k/orbital ownership, HSE/PBEh reconstruction and Si reference spectra; ON/OFF builds and independent review. Root full basis and halo/dense-matrix replication remain future work.

## Results

- HSE regression: eight tests pass, including zero-core ranks, domain-sized receive-workspace assertions, LCFO spectra and reconstructed mesh RT.
- Si reference: four physical k points and 192 eigenvalues pass on both four ranks (k split) and eight ranks (k/orbital split).
- Independent review found no blocking issues in request lifetime, contiguous column transfer, explicit source/tag ordering or duplicate-owner selection.
- ON/OFF builds and diff whitespace check pass. Root output basis, received halos and dense eigensolver matrices remain unchanged; this phase removes only the full-core gather scratch arrays.
- Targeted PBEh regressions: water mesh propagation and spatial mesh pulse both pass; water serial/spatial energy difference 7.93e-12 eV. No total-process memory or speedup claim.
