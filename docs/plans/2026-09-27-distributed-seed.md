# Distributed Gamma seed QR

Goal: remove the remaining root Noccupied × Nbasis QR allocation.
Architecture: opt-in SALMON_LCFO_SEED_DISTRIBUTED=1; one-row BLACS grid,
32-column cyclic distribution. Transfer conjugated coefficient row blocks directly
from their original owners to QR owners. Write snapshot columns using bounded
Gatherv buffers. PZGEQPF returns distributed global pivots; recover only selected
original rows through the existing seed SVD path. Default root ZGEQP3 unchanged.
Tech stack: Fortran, MPI, ScaLAPACK, existing GNU/OpenBLAS test environment.

Alternatives: root QR preserves exact pivot ordering but cannot remove root storage;
2D QR distribution balances tall matrices better but this matrix is wide, and a
one-row grid distributes basis columns without introducing another row redistribution.
Do not change localization tolerances or introduce coefficient truncation.

1. Extend MPI seed tests with explicit backend detection, snapshots, empty and
   unequal owners, block/tile tails, singular/invalid dimensions. Demonstrate failure.
2. Implement bounded redistribution, snapshot streaming, distributed pivot QR,
   collective status handling. Unsupported builds must reject explicit requests.
3. Compare serial/root/distributed paths, measure isolated RSS in separate processes.
4. Build and compare C128/MPI16 RT from identical GS; inspect localization and current
   if pivots differ. Keep opt-in unless physical equivalence has adequate evidence.
5. Review, record limitations and measured results, provide archive-friendly patch,
   push existing development branch. Fugaku compiler/numerical validation remains external.
