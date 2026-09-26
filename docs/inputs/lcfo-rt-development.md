# LCFO RT development status

The DC input guard remains unchanged: `yn_dc=y` currently accepts `theory=dft`
only. There is no Si128 self-consistent TDHSE dielectric spectrum from this work
and no new user input switch enabling one.

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

Still required: variationally consistent global exchange assembly, global
Hartree/semilocal updates, U transport and support controls, ACE integration,
self-consistent initial state, electromagnetic coupling/current, and spectra.
A Hermitian matrix alone establishes neither energy-functional consistency nor
optical-response accuracy. In particular, the static source-cutoff scripts are
not silently promoted into a variational time-dependent functional.

Tests: CTest `lcfo_rt_core`; `OPENBLAS_NUM_THREADS=1 python3
testsuites/unit_lcfo_rt/test_reference.py`. The standalone check.py compiles
against the local Homebrew BLAS for the current development machine; CTest uses
the build's selected BLAS/LAPACK and is the portable verification path.
