# Dense-k Hermitian real-space pairs

Native CPU exchange on a full uniform k mesh automatically uses Hermitian
real-space pairs when crystal symmetry is not enabled. Set
`yn_exx_k_pair_symmetry='n'` in `&functional` to use the full-pair reference.
The default is `'y'`. HSE06 and the supported global hybrids share this route.
Space-group reuse (`yn_exx_fft_orbits`) takes precedence when crystal symmetry
is enabled. Gamma, q-downsampled, GPU, spatially distributed and Wannier routes
retain their existing implementation. No symmetry file is required.

After the existing Bloch-phase adjustment, form the occupied projector

\[
P_k(r,r')=\sum_n s_{nk}(r)s_{nk}^*(r').
\]

The existing k-mesh forward FFT, real periodic kernel multiplication and backward
FFT produce an exchange matrix `Q_k(r,r')`. The real even kernel preserves
`Q_k(r',r)=conjg(Q_k(r,r'))`. This is Hermiticity, not a point-group assumption
about atoms or the external field. A target need not belong to the source set.
For each retained `r <= r'`, accumulate

\[
(A_k t)(r)\mathrel{+}=-Q_k(r,r')t_k(r')/N_k,
\qquad
(A_k t)(r')\mathrel{+}=-Q_k^*(r,r')t_k(r)/N_k.
\]

The second update is omitted on the diagonal. The final Bloch phase is removed
as in the reference route. The Coulomb cutoff, HSE screening, q=0 value,
occupations, hybrid fraction and integration weights are unchanged.

Pair indices use 64-bit arithmetic. Each rank owns bounded tiles containing at
most `hse_block_rows * N_r` pairs. There is no stored `N_r * N_r` matrix or full
pair-index table. Both MPI transposes retain all k points; only the pair axis is
reduced. Exact tail-batch FFT plans avoid transforming padded final pairs and are
cached and destroyed with the kernel. Each direction transforms
`N_r * (N_r + 1) / 2` pairs rather than `N_r * N_r`. This count is not a claim
that total wall time halves. Packing, both action updates and communication still
contribute and require machine-specific measurement.

Collective validation precedes variable-length communication. Uneven and empty
k owners are supported. Kernel and ownership layouts must agree across ranks.
The application prints `EXX_K_PAIR_SYMMETRY unique/full pairs` when selected;
`variables.log` records the switch. Per-stage profiling remains available.

Validation includes complex independent source and target orbitals, shifted and
permuted meshes, HSE and Coulomb kernels, both FFT layouts, partial tiles, empty
k owners, collective rejection of inconsistent ownership, and GS/short-RT
comparisons. Runtime scratch and tail plans are not checkpoint state, so no
restart file format changes are introduced.
