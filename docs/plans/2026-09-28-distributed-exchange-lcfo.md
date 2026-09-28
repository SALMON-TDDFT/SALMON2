# Exchange and LCFO distribution plan

Continue the approved computation/memory distribution. First remove avoidable FFT transposes by retaining Z spectral pencils, evaluate kernels in that layout and reuse the density buffer for inverse output. Then remove replicated full-fragment Hpsi arrays from complex DC-LCFO: retain local grid rows, compute diagonal/halo projections over local intersections, reduce only band matrices. Halo basis exchange remains between fragment representatives followed by fragment broadcast. Basis and dense diagonalization replication remain documented future work; do not claim full orbital distribution.

Validate screened/fractional exchange oracle across ranks, HSE DC-LCFO spectra and reconstructed RT parity, PBEh regressions and ON/OFF builds. Add local-workspace diagnostics/assertions demonstrating Hpsi grid storage scales with spatial decomposition. Independent review before commit.

## Validation

- Screened/Coulomb, fractional and zero-occupation exchange compared with serial Wannier on 1/2/4 ranks; counter asserts exactly four pencil transposes per source/batch pair.
- HSE regression includes assertions that local Hpsi grid points times spatial rank count equal fragment grid points, LCFO eigenvalue parity, reconstructed RT, and legacy multi-k thermal charge.
- Existing Si complex DC-LCFO reference test: four physical k points, 192 eigenvalues and total energy -1702.2269 eV verified on 4 ranks (k split) and 8 ranks (k plus orbital split).
- Independent review found no blockers in Z-pencil ordering/normalization or distributed projection/root ownership; requested additional orbital ownership coverage was completed with Si.
- Remaining full-basis, halo and dense-matrix replication is explicitly documented; this phase does not implement full orbital-column distribution.
- Final results: HSE 7 tests and PBEh Ehrenfest 16 tests pass. All 36 exchange rank/kernel/occupation comparisons pass (maximum action difference 5.18e-15). ON/OFF builds pass. Whitespace diff check passes.
- H4 Hpsi-only workspace: 192 KiB old per-rank full-grid buffers -> 24 KiB local buffer at four spatial ranks; basis/halo/other solver memory excluded.
