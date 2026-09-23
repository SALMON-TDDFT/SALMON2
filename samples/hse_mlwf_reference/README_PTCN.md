# Parallel-transport Crank–Nicolson + ACE

`ptcn.py` implements Eq.12 of [Jia and Lin](https://arxiv.org/abs/1809.09609), separately from `ace_rt.py` (original-gauge midpoint). With column notation and grid-weighted products, define F(U)=H[U]U-U(U†H[U]U). The accepted endpoint solves

    V + i dt F(V)/2 = U - i dt F(U)/2.

Arrays in code store orbitals as rows. `pt_residual` performs the equivalent row operations independently for each k. Initial orbitals are orthonormal; nonlinear iterates are not forced onto the orthonormal manifold. Neither a Gram-inverse replacement of the projection nor post-solve normalization is applied, since either would change the discrete equation. Trapezoidal PT is not asserted to preserve endpoint norms exactly. Endpoint Gram and electron-count errors are monitored.

An inner iteration uses kinetic preconditioning, current Hartree/semilocal terms and a fixed ACE. Each outer iteration constructs full MLWF exchange from the trial endpoint. Acceptance requires full nonlinear relative residual <1e-10; inner threshold is 1e-12. If a sufficiently small inner residual stagnates, the solver can hand control to a fresh full-exchange check without changing the full tolerance. It records these inexact-inner exits; failed outer convergence or an unsolved large inner residual raises an error. No pair truncation or ACE rank truncation is used.

The driver caches the accepted endpoint's full exchange and compressed operator. At the next step this exact initial-endpoint action supplies the right-hand side; the compressed operator is only an initial seed for the new endpoint. This is correct for the implemented constant vector potential after an impulse. A laser implementation must separately set A(t_n) for the right-hand side and A(t_n+1) for endpoint iterations; the current driver does not implement a time-dependent pulse.

The propagated PT orbital gauge and the MLWF representation gauge serve different purposes. The existing MLWF gauge is transported/relocalized every10 steps and used only inside full exchange evaluation. Physical observables are calculated directly from the propagated orbitals and the native current operator. No extra gauge connection term is added to the already transformed PT equation.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 samples/hse_mlwf_reference/benchmark_ptcn.py calculations/si_hse_reference/export calculations/si_hse_reference/scf/state.npz calculations/si_hse_reference/ptcn_dt032 1 0.32 0.0001
```

The driver requires the previous converged all-pair Si HSE ground state and its adjacent result metadata. It reports bootstrap cost, per-step full builds and inner iterations, energy/norm/current diagnostics, and actual BLAS environment/single-thread FFTW. Large checkpoints remain ignored by Git. These short runs establish neither a production time step nor an optical spectrum.
