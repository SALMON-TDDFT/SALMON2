# ACE reference extension

`ace.py` constructs a separate negative Hermitian low-rank exchange for each k point. The input action contains the full screened exchange, without the HSE 0.25 coefficient. Inner products include the real-space quadrature volume. The coefficient is applied once by the Hamiltonian caller. ACE interpolates its construction subspace; this is not a guarantee of accuracy on arbitrary changed orbitals.

`ace_rt.py` implements a **self-consistent implicit midpoint method in the original Bloch gauge**, not parallel-transport Crank–Nicolson. For midpoint orbitals X and initial orbitals U, it solves

    X - U + i dt/2 H[X] X = 0,    U_new = 2 X - U.

Local density-dependent terms change on every inner iteration. Exchange is compressed and reused in the inner loop. The outer loop rebuilds the full, untruncated exchange from X, checks the full nonlinear residual, and refreshes ACE if needed. Acceptance requires the full residual, not merely convergence under a stale ACE. Midpoint orbitals are not orthonormalized: this would change the implicit equation. At convergence the Hermitian midpoint Hamiltonian preserves endpoint orbital overlaps up to solver tolerance. Energy conservation and time-discretization error must still be measured. The code rejects inner/outer nonconvergence.

`benchmark_ace.py` measures full exchange, ACE construction, seven applications after warm-up, and actual short RT steps including full rebuilds. Additional full energy evaluations and standalone setup timings are outside per-step wall time and included in total runtime. FFTW/BLAS run in one thread. A constant vector potential represents the impulse; no pump waveform or optical spectrum is implemented here. The existing SCF driver remains unchanged in algorithm; only its exchange-action method was separated for reuse.

From the repository root, after producing the earlier converged Si reference:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 samples/hse_mlwf_reference/benchmark_ace.py calculations/si_hse_reference/export calculations/si_hse_reference/scf/state.npz calculations/si_hse_reference/ace_plus 1 0.08 0.0001
```

MLWF gauge minimization has its existing 10-step cadence. ACE rebuilds are determined separately by the nonlinear solve. The complete all-pair RK4 implementation remains available for comparison. PT gauge, larger-step convergence, long-time stability, and pump–probe spectra are subsequent work.

For steps after the first, the previous accepted midpoint's ACE is used as an initial operator. This skips an initial full rebuild only; every accepted step still performs fresh full-exchange residual checks, and a stale seed cannot bypass the convergence gate. No fixed multi-step exchange freezing is used. This seeded variant is recorded separately from the first single-step/time-step pilot timings.
