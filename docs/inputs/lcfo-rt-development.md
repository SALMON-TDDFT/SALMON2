# LCFO RT development status

The DC input guard remains unchanged: `yn_dc=y` currently accepts `theory=dft`
only. An experimental fixed-LCFO-subspace RT adapter now reuses native
`initialization_rt`, density, Hartree, semilocal XC, pseudopotential, Taylor4
predictor/corrector and current routines. No additional SCF is performed.
Si128/16-fragment self-consistent RT has passed a four-step integration check;
there is still no converged dielectric spectrum.

## Experimental native path

Set environment `SALMON_LCFO_RT=1` and use `theory='tddft_response'`,
`yn_dc='n'`, `yn_conventional_from_dcdft='y'` in `&calculation`.
Use Gamma, unpolarized HSE06, `yn_hse_wannier='n'`, default `hse_taylor4`,
one real-space MPI rank per fragment, and no orbital/k distribution.
Each rank's grid must coincide exactly with its fragment core. The existing
`./data_dcdft/fragments` LCFO records provide both the initial orbitals and fixed
orthonormal core bases. Native `hpsi` is projected after all Hamiltonian terms;
Hartree and density are evaluated on the whole physical system by existing code.

The current implementation requires occupied-only fixed occupations (for Si128,
256 states, 512 electrons, omit temperature). It loads the first256 saved LCFO
states without reconverging the accepted DC-to-LCFO density difference.
Restart input/output, checkpoints and time_shutdown are rejected until LCFO
basis and exchange state metadata are supported.

`hse_lcfo_rt.f90` reconstructs periodic fragment+buffer bases and density factors,
uses the existing screened FFT exchange kernel, core-weights and Hermitianizes
the projected exchange, and sums contributions once across spatial ranks.
The global trace energy is passed to the existing native energy bookkeeping.
Existing ACE is built in coefficient space; a rejected indefinite/singular
metric falls back explicitly to the full projected operator without clipping.
The predictor/corrector averages endpoint operators. Full support is used:
density eigenfactors are NOT MLWFs. A relative density eigenvalue threshold of
1e-14 only removes numerical null modes, with discarded trace logged.

The environment switch is a development opt-in, not a new production input.
Without it the ordinary reconstruction and native RT paths are unchanged.


`src/rt/lcfo_rt_core.f90` supplies a fixed orthonormal complex basis kernel:

- `lcfo_cayley_step(H,C,dt,next,status)` solves
  `(I+i dt H/2) next=(I-i dt H/2) C`. It does not diagonalize H or normalize C.
- `lcfo_density(C,f,P,status)` forms `P=C diag(f) C†`, retaining complex coherence
  and fractional occupations. f is the physical occupation; no implicit factor2.
- `lcfo_midpoint_step` iterates a Hamiltonian callback on `(P_start+P_end)/2`
  at `time+dt/2` to the specified relative Frobenius density residual. Each trial
  propagates from the unchanged starting C. It is a small reference solver;
  fixed-point convergence is not guaranteed for arbitrary dt or functionals.
- `lcfo_grid_density` reconstructs density on a core's grid from its basis and
  the corresponding coefficient rows of ALL global occupied states.
- `lcfo_project_potential` projects one real grid potential using grid quadrature.

Outputs use separate storage from inputs. Failure retains the input C in `next`
when its shape permits; status1 invalid input,2 non-Hermitian H,3 linear solve
failure,4 fixed point not converged,5 callback failure. A caller must reject the
step on any nonzero status. No silent symmetrization of caller H is performed.

The independent Python reader `samples/dc_hse/lcfo_rt_reference.py` accepts only
Gamma, unpolarized V1 complex LCFO data. It verifies completion/footer lengths,
run identity, geometry, periodic halo mapping, orthonormal core bases and saved
eigenvectors. It reconstructs the native symmetrized matrix, retaining a separate
raw-directed anti-Hermitian diagnostic. Basis arrays use Fortran grid ordering.
The `lcfo_frozen_probe.f90` helper tests real saved LCFO eigenstates against exact
Cayley phases and time reversal; this is deliberately NOT a TDHSE driver.

Still required: MLWF U transport and integration-support controls in this native
adapter, longer-time/dt/LCFO-basis convergence, and dielectric comparisons.
The accepted initial density mismatch is retained; stationarity is not a gate
for starting RT.
A Hermitian matrix alone establishes neither energy-functional consistency nor
optical-response accuracy. In particular, the static source-cutoff scripts are
not silently promoted into a variational time-dependent functional.

Tests: CTest `lcfo_rt_core`; `OPENBLAS_NUM_THREADS=1 python3
testsuites/unit_lcfo_rt/test_reference.py`. The standalone check.py compiles
against the local Homebrew BLAS for the current development machine; CTest uses
the build's selected BLAS/LAPACK and is the portable verification path.

Native integration regression: `python3 testsuites/unit_lcfo_rt/test_native.py
--binary /absolute/path/to/salmon --pseudo /absolute/path/to/H_rps.dat`.
It runs MPI2 jobs sequentially in a fresh temporary directory: Gamma DC-SCF,
LCFO RT, half-dt RT, and rejection of unsupported restart. It checks normal
completion, finite output, endpoint-current agreement and post-impulse energy
width. This small integration test is not a production convergence criterion.
