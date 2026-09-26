# Active WF reconstruction implementation plan

**Goal:** Reduce repeated real-space WF reconstruction without adding a physical approximation beyond the existing spherical source mask.

**Architecture:** Precompute the fixed periodic point masks after the initial centers are known. Reconstruct only WF columns whose support intersects each fragment/core (always include protected WFs). Return compact fragment sources. For the diagnostic, evaluate total core norm with the cached basis Gram matrix and subtract the directly reconstructed retained norm. Keep dense coefficient-frame rotation and polar U transport unchanged in this first step.

**Tech Stack:** Fortran, BLAS, existing MPI/OMP and frozen native reference trajectories.

The user's existing authorization covers WF reconstruction/communication efficiency and local-mask development; the latest steering prioritizes the WF bottleneck. Keep the dedicated branch and freeze the ongoing33dbea5a weak-scaling binary. No numerical regression or compilation during its performance runs.

1. Add standalone algebra test against dense real-space reconstruction followed by the original point mask. Cover complex nonorthogonal basis, distant/protected WFs, periodic edges, empty support, full range, and discarded-norm identity. Confirm missing-module compile failure before implementation.
2. Add `lcfo_wf_support.f90` with a reusable support plan, compact reconstruction, and Gram-based total norm. Register CMake and standalone fixtures.
3. Integrate two cached plans (fragment/core) in `lcfo_rt_wannier.f90`. Preserve initial gauge, U transport, impulse/laser schedule, periodic distance, protected WFs, and global diagnostic definition. Report reconstructed/possible columns.
4. Once frozen weak scaling completes, run algebra, transport/sphere/native MPI2/MPI4 suites and build. Request static review.
5. Re-run only the C64/C128 R6 conditions with the same GS/OMP1/BLAS1/ACE1/16 steps and compare current, density, energy, initial dump, loss diagnostic, phase times and total RT. One numerical job at a time. Record remaining dense U costs; commit/push only after verification.
