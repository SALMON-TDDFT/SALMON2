# Dense-k Hermitian real-space pair reuse implementation plan

Goal: automatically use Hermitian real-space pair exchange for native dense k without a crystal symmetry file. Keep a namelist opt-out and the reference path.

Architecture: on the full uniform k mesh, pack the upper triangle (r<=r') of the phase-adjusted occupied density matrix in bounded pair tiles. Fourier-transform each retained pair over the k mesh with the existing Coulomb/HSE multiplier and q=0 prescription. Return k-local entries, accumulate V(r,r') target(r') and its conjugate partner V(r',r) target(r), then remove Bloch phases. Targets may differ from sources. Preserve all q points and CPU-only initial coverage; q downsampling, crystal symmetry and GPU retain their existing paths. Tail FFT plans avoid transforming padded pair slots and are owned/destroyed by the existing kernel object.

1. Write a failing independent full-vs-pair action oracle using complex sources and distinct targets, shifted/unshifted k meshes, HSE/global kernels, uneven/empty k ownership and invalid layouts. Confirm missing entry point before implementation.
2. Add src/xc/exx_k_pairs.f90 and bounded kernel-owned tail plans/counters in exx_k_exchange. Validate collective conditions before variable-size communications, use int64 pair indices, no full N_r^2 matrix.
3. Add yn_exx_k_pair_symmetry='y' to functional declarations/defaults/broadcast/lowercase/log/validation. Select only CPU all-q full dense-k with no crystal symmetry; 'n' uses reference. Leave crystal-symmetry selector unchanged.
4. Validate pair mapping diagonal/off-diagonal updates, arbitrary target counts, exact pair FFT counts, Hermiticity and bounds with 1/2/4+ MPI. Add a small numbered CPU/MPI regression, keeping restart probes outside testsuites.
5. Test GS/RT on/off for HSE/PBE0, existing crystal symmetry and q-reduction fallbacks, and serial/MPI builds. Update both repository docs and SALMON-DOCS; review and leave uncommitted until requested.
