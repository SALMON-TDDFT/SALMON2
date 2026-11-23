# TDDFT workspace ownership implementation plan

Continue the approved extraction using existing isolated worktree and local tests.
XC 3D and PZ work arrays belong to s_dft_system; calc_xc accepts an optional
explicit PZ workspace and otherwise uses call-local storage. This preserves the
immutable functional/Libxc handle contract and avoids aliasing a parent mutable
functional with its output-buffer components. Resize only that owner's arrays
when shape/spin changes. Stencil modulo tables and the optimized-callable flag
belong to s_stencil and are passed to numerical consumers explicitly.

1. Add a failing native probe for two grid shapes, PZ evaluation, workspace resize,
   independent modulo tables, and deallocation of one owner.
2. Move workspace definitions into structures, route exchange_correlation and
   builtin PZ, then route RT/DFT stencil initialization and Hamiltonian consumers.
3. Fresh audited native build; lower-level probe, full-context interleaving at
   independent dt, and original fast-HSE driver comparisons.
4. Synchronize Japanese/English/AI docs, run Python and sdist-to-wheel tests.

Do not remove common-grid/XC guards: global input/parallel flags still need work.
No heterogeneous full-TDDFT, MPI, GPU, concurrent host calls or physical-convergence
claim. No GitHub Actions. Preserve earlier immutable validated builds and fixtures.
