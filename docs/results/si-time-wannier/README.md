# Si: propagated Bloch orbitals and time-dependent Wannier localization

These are offline diagnostics of the existing alpha0=1, K0=0 polarization-closure trajectories. No Wannier-to-alpha formula, production alpha-update cadence, or screening feedback has been added. Existing alpha is sampled every ten steps in the dense analysis. In particular, these trajectories do **not** use fixed alpha=0.2.

## Setup and meaning

Si8, 12³ real-space grid, 4³ shifted k mesh, 16 doubly occupied orbitals; dt=0.08 au, 1600 steps = 3.09617 fs. Pulse frequency0.2 au, envelope60 au (1.45133 fs), intensities0, 1e8, 1e13 W/cm². No k convergence scan. The original instantaneous bare-K alpha estimator and its known field-gate dependence remain unchanged.

The user requested reusing the preceding localization U(k) as the next initial guess. `analyze_dense.py` minimizes the discrete Marzari–Vanderbilt spread every10 steps (0.0193511 fs), starting from the preceding converged U transported into the current occupied Bloch basis. Specifically, Q(k)=polar[u(k,t)† u(k,t_previous) U(k,t_previous)] is the initial guess; the overlap includes the real-space integration weight. A pure occupied-basis rotation therefore leaves the starting Wannier orbitals unchanged. Subsequent localization still uses the current occupied subspace and does not remove its physical evolution. The initial gauge is a projection onto16 Si bond Gaussians followed by polar orthonormalization and spread minimization. Localization uses all6 nearest k neighbors including reciprocal-cell boundary phases. This is a locally converged MLWF gauge; global optimality is not guaranteed. Orbital overlaps and periodic center displacements monitor continuity. Reusing U does not time-average any observable.

`mlwf.py` reports the finite-mesh spread and its gauge-invariant part. Its unitary descent uses analytic derivatives, with an Armijo search on the gauge-dependent part alone to avoid loss of significance against a large invariant constant. Direct U reuse without temporal transport is retained as a no-pump control in `none_direct_u_mlwf.json`: it converges but typically needs 219.09 iterations on average to undo the band phase evolution. Polar transport reduces this to zero additional localization iterations at every no-pump snapshot, with at most8.1e-12 Å² difference from the direct-U minimized spread. The stopping criterion is Frobenius gradient norm <1e-6 in atomic units. Spread is a one-particle localization diagnostic, not an exciton relative-coordinate radius, dielectric constant, or screened Coulomb interaction.

The occupied time-evolved subspace remains fully occupied. Carrier excitation is diagnosed separately by projection onto the ground-state reference bands with SALMON's `projection_option='gs'` (Houston field shift included). The finite32-band reference window slightly undercounts excited electrons relative to holes.

## Every-ten-step MLWF results

All161 snapshots of each of the3 cases converged with gradient norm below1e-6. No orbital label changes occurred. The final time is3.09617 fs.

| Case | Mean MV spread (Å²/orbital) | Invariant part (Å²/orbital) | Mean / max iterations after initial frame | Minimum adjacent orbital overlap |
|---|---:|---:|---:|---:|
| No pump |2.057679|1.885606|0 / 0|0.99999999994|
| Weak |2.057711|1.885635|103.3 / 130|0.9999999894|
| Strong |4.936765|4.350975|240.4 / 286|0.9989425|

The no-pump spread varies by less than4e-9 Å². Strong-pump maximum center displacement between adjacent snapshots is0.0400 Å, and the minimum temporal occupied-subspace singular value is0.99565. A separate initial-ground-gauge start at the final strong snapshot converges in489 iterations to the same spread as the sequential warm start within2.4e-13 Å². This supports the warm-start result for this trajectory, without proving global optimality for arbitrary states.

Strong-pump broadening persists after localization, including an increase of the gauge-invariant contribution from1.88561 to4.35097 Å². It is therefore not solely a consequence of an unoptimized occupied gauge. A transient occupied-projector localization change does not by itself establish metallic screening or an exciton radius.

The original estimator's final alpha remains1 (no pump),0.0129912 (weak), and0.000206291 (strong); these are logged values of the old model, **not** inferred from MLWF. Weak MLWF spread barely changes despite the old estimator's large alpha reduction. No screening coefficient is calibrated from this observation.

See `mlwf_metrics.json`, `{none,weak,strong}_mlwf.json`, and `mlwf_dynamics.png` / `.pdf`. The `.npz` files retain centers at every snapshot and the final U; intermediate U values are reproducible by sequential analysis of the raw checkpoints.

## Fixed-gauge controls

`analyze.py` keeps the initial projected Wannier U fixed, and additionally removes the same no-pump occupied rotation from all cases as a control. This control is **not** MLWF minimization. Its moments are real-space minimum-image moments in the 4³-cell supercell, length21.7174 Å; they are not numerically identical to the finite-mesh MV functional.

At3.09617 fs:

| Case | Fixed initial U: spread (Å²/orbital) | Field-free-aligned: spread (Å²/orbital) | Excited electrons / Si8 |
|---|---:|---:|---:|
| No pump |93.9881|2.25856|~0|
| Weak |93.9881|2.25864|0.00002014|
| Strong |95.2155|9.69603|1.78874|

The strong-pump projection gives1.83997 holes/cell and1.78874 electrons/cell (about1.12e22 cm⁻³). The missing0.05124 electron reflects the finite reference-band window. This is substantial excitation; an earlier assumption that the pulse produces only a small excited population is not supported by this diagnostic. This population is not automatically the free-carrier density.

No-pump fixed-U broadening is dominated by orbital phase evolution; removing that rotation restores the initial localization. Fixed-U late-time boundary weights reach0.27, so its very large widths cannot be treated as unbounded real-space radii on this finite mesh. Strong aligned boundary weight reaches0.028. See `metrics.json`, `dynamics.png`, and `density.png`.

## Verification and reproduction

Independent review confirmed normalization, boundary phases, the analytic gradient, the gauge-dependent line-search objective and the temporal transport orientation. Nine tests cover temporal-basis transport, projected-gauge covariance, singular projections, supercell norm and density, periodic moments, the MLWF gradient, an exactly localized mesh-boundary case and occupied-gauge invariance. Coarse reconstructed density errors are below1.3e-16; orbital norm errors below1.1e-9. Diagnostic current differences from the preceding bare-K runs are below1.1e-17.

From the repository root, use Python with NumPy and Matplotlib, `OPENBLAS_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`:

```sh
python3 -m unittest discover -s docs/results/si-time-wannier -p 'test_*.py'
python3 docs/results/si-time-wannier/run.py none
python3 docs/results/si-time-wannier/run.py weak
python3 docs/results/si-time-wannier/run.py strong
python3 docs/results/si-time-wannier/analyze.py
python3 docs/results/si-time-wannier/initialize_mlwf.py
python3 docs/results/si-time-wannier/run_dense.py none
python3 docs/results/si-time-wannier/run_dense.py weak
python3 docs/results/si-time-wannier/run_dense.py strong
python3 docs/results/si-time-wannier/analyze_dense.py
```

Run scripts use the existing MPI executable `/private/tmp/salmon-tdcdft-mpi/salmon` and ground-state restart under `calculations/si_tdcdft_k4/gs`. Raw wavefunction checkpoints are in `/private/tmp/salmon-si-time-wannier{,-dense}` and are not committed. Dense checkpoints occupy about13 GB; regenerate them if removed. Inputs, completion statuses and reduced results are retained here.

Reference: [Marzari and Vanderbilt, PRB56, 12847 (1997)](https://doi.org/10.1103/PhysRevB.56.12847). The spread minimization is over occupied-band unitary rotations; dielectric screening requires additional physical information beyond those rotations.

## Subsequent alpha coupling (not implemented here)

The requested integration point is t_j=10*j*dt: obtain the current localized basis, evaluate an explicitly chosen screened interaction / response in that basis, then update alpha. Alpha can be held over the intervening RT steps while P and the polarization-consistent XC field continue evolving each step. Previous U is only a numerical initial guess; no time averaging of alpha, current or field is required by this scheme. Restart would need both the held alpha and the most recent U / update step.

What remains undetermined is the physical functional mapping the nonequilibrium state to the screened interaction and then alpha. Neither an occupied-band unitary transform nor its minimized spread alone provides it. Using a width ratio as alpha would introduce a new empirical model and is intentionally not done in these diagnostics.
