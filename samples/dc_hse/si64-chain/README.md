# Si64 8x1x1 chain experiment

The user confirmed rearrangement of the same 64 atoms into eight conventional
Si cells along x, not replication of the original Si64 cell. Lattice constant
10.26 bohr, total box 82.08 x 10.26 x 10.26 bohr, grid 128 x 16 x 16,
eight DC fragments, each core 16^3 / 8 Si atoms. Periodic in all directions;
this is a bulk supercell, not an isolated nanowire.

The committed input is a 20-iteration pilot, not a converged reference.
x buffer 8 gives a 32 x 16 x 16 periodic fragment (2a x a x a), with 64
retained states. Proposed buffers 16 and 24 give periods 3a and 4a. Keep
zero y/z buffers for unsplit periodic directions. Retained-state convergence
must be checked when increasing the fragment.

Run only one simulation at a time. Current independent-pair comparison jobs
were cancelled after CPU contention; their timings are not scaling evidence.
The 20-step optimized and no-loop-vectorization pilots both had large SCF
residuals. Saved source-frame orthogonality defects were about 1e-14.

The recorded Netlib ZLARF1L issue is already guarded in the fallback build.
Current binaries use OpenBLAS; the saved LAPACK reproducer passes at 2.05e-15
relative residual. Do not describe SALMON's no-loop build as fixing that
separate dependency issue.

WF sizes will be measured using the minimum-image x coordinate around each
periodic center, integrating over the periodic y/z cross-section. Axial tail
norms, not spherical radii, are the intended metric. Compare Phi (orthonormal
frame) separately from Q (fractionally occupied density factors). A small norm
tail alone does not certify energy, force, or time-propagation accuracy.
No acceptable numerical cutoff radius has been established yet.

The user clarified that "full calculation" means the eight-fragment DC system
with no support truncation and no pair pruning. An undivided MPI1 Si64 trial was
a mistaken interpretation and must not be presented as the requested baseline.
The current requested baseline uses MPI8, OpenMP2 and BLAS1, with a single job.
Its convergence must be established before reducing support.
