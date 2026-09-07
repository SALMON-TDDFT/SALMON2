# Wannier90 Native Line-Search Convergence Design

## Purpose

Test the simplest history-free adaptive-step route before introducing a custom Pulay/DIIS controller. The Si64 reference must use the corrected Bohr geometry and the production global/spectral seed, not the random projection used by the failed 10,000-iteration Hybrid-SCF fixture.

## Design

SALMON writes `trial_step = 2.0d0` and `num_cg_steps = 0` to the Gamma-point Wannier90 input. It does not write `fixed_step`. Wannier90 therefore evaluates its native parabolic line search from the standard large trial step on every iteration, while rebuilding the steepest-descent direction from the current gradient without conjugate-gradient history.

No custom step schedule, Pulay history, gauge restart, or post-hoc acceptance of an unconverged transform is added. The existing convergence validator remains fail-closed.

The Si64 run uses `dg_ow_w90_initial_projection='spectral'`, `wannier_num_iter=10000`, MPI 8, OMP 1, and no wall-time cutoff. The exported A/DMN/MMN/EIG bundle is retained. Evidence records spread, gradient, accepted line-search step, symmetry, memory, and whether downstream Hybrid-SCF is reached.

## Acceptance

- Generated `.win` contains `trial_step = 2.0d0` and `num_cg_steps = 0` and contains no `fixed_step`.
- Correct Bohr geometry and spectral initialization are present.
- No custom adaptive-step patch or Pulay/CG history is active.
- Wannier90 itself reaches its convergence condition before downstream TDDFT validation begins.
- Failure to converge is reported without weakening the convergence gate.

