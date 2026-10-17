# Global hybrid symmetry FFT reuse

Approved goal: extend the existing dense-k space-group FFT reuse from HSE to PBE0/PBEh, preserving SR as HSE-only. Reuse the existing isotropic finite-radius Coulomb kernel and its zero-frequency value; no new long-range approximation.

1. Red test: run PBE0 symmetry/reuse input against current binary and confirm the HSE-only input guard rejects it.
2. Broaden the reuse selector to is_hybrid and bypass automatic Wannier route selection only for the supported global-hybrid symmetric dense-k native path. Preserve CPU, all-q, integer occupation, k-only distribution, cubic/orthogonal restrictions and field validation.
3. Test omega=0 exchange oracle (shifted and unshifted, MPI), PBE0/PBEh converged GS and short RT on/off, and standard HSE regressions.
4. Document finite Coulomb radius conventions and remaining restrictions. Review final changes and leave uncommitted pending user instruction.
