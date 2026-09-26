# Fragment x orbital MPI implementation plan

User approved Gamma-point two-level MPI and Si64 8x1x1 MPI8 versus MPI16.
Reuse native icomm_r (same orbitals across fragments) and icomm_o (same core
across orbital groups). Native orbitals, Hamiltonian actions, density/current
remain distributed by SALMON's io_s:io_e ranges. Configure one core per r rank.
Distribute local LCFO basis to orbital replicas during initial reconstruction.
Gather coefficient columns to orbital group zero for refresh only; construct
MLWF sources/Hx/ACE there once, then broadcast coefficient operators along o.
ACE application reduces overlaps along r independently for each orbital block.
Full Hx and ACE coefficient factors remain replicated; refresh masters still
hold all occupied coefficients/source WFs. This is not fully distributed ACE
factorization. Preserve native Hartree/XC updates and impulse-first rebuilding.

1. Add MPI2x2 regression reusing a single GS; confirm old layout rejection.
2. Adapt LCFO setup, local target packing, refresh collection and state broadcast.
3. Compare current, energy, density and rebuild schedule; test unequal columns.
4. Use converged Si64 8x1x1 LCFO input for MPI8 and MPI16 sequential benchmarks.
   OMP1/BLAS1, same GS/initial gauge, radius9, ACE4. Record elapsed and RSS.
5. Review, document limitations, commit and push the existing GitHub branch.
