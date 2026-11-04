# Occupied-only real-space periodic PZ–KLI smoke test

Based on the committed native Wannier development line at b3388e45e34fec8c65ee354ad9c338ff9cfa27cd. Si uses two atoms, a diamond primitive lattice, the KY Si potential and 6³ real-space grid from regression 131, and a full equal-weight 2³ k mesh. Eight electrons occupy exactly four bands per k. LSDA and static KLI use four additional empty bands in the GS subspace diagonalization. Wannier localization and real-time propagation use the four occupied bands. The occupied grid wavefunctions are propagated by SALMON's existing midpoint Taylor predictor/corrector. Wannier and projected Hamiltonian rotations evaluate a potential; they never replace the propagation variables.

Build with native Wannier/FFTW support; at runtime select built-in PZ, not hybrid exchange. The implementation is disabled unless selected in the functional namelist. For the static run, start from the converged LSDA WFN checkpoint and diagonalize from the first iteration. A fresh RT run reads the converged KLI GS restart directory, including kli-v1.bin; RT-to-RT resume is not supported. The initial occupied projector, crystal, grid and occupations are checked. Use shared restart files, fixed ions, integer occupations, no spatial decomposition, and no orbital decomposition for RT. The RT adapter rejects rectangular spin allocations which would include empty orbitals. k-point decomposition gathers occupied states before cross-k reconstruction; the small solver is replicated, not distributed for large systems.

The seed and transport are SALMON's native wannier_localize/gauge_transport and the BvK synthesis uses its forward/backward transforms. The forward phase is evaluated in fractional coordinates with the common k-origin phase removed (it cancels in densities and overlaps); physical quadrature uses the primitive-cell volume per point. The home Wannier has 1/Nk normalization; its translates span the full occupied projector. The discrete gradient is the anti-Hermitian part of W(k)†[uw](k); all translation coefficients are checked. The optimizer descends the SIC energy by unitary Cayley steps with backtracking. This verifies stationarity, not global optimality. Self-Hartree uses a tracked image cut fixed within each minimization and an enlarged FFT box with the full nonorthogonal metric and analytic spherical cutoff kernel (including G=0). Neighbor images cannot enter the cutoff. The kernel and FFT plans are cached. Refining padding does not improve the finite-k orbital. The auxiliary check compares a second box with 25% larger linear grid dimensions; the density-weighted RMS change in the potential must pass kli_hartree_tol (Ha).

For KLI, p_i is primitive-normalized and pq_i=Re(conjg(psi_i)q_i). Solve in the zero-mean D gauge, reject extra null modes, and set the Slater and response terms to zero on numerically empty points (rho <= tiny), rejecting nonzero numerators there. GS uses SALMON's Ritz basis, including its choice within exact degeneracies. RT uses the projected endpoint Hamiltonian with endpoint Hartree, LSDA and vector potential and the last committed local SIC. Eigenvalue blocks within 1e-7 Ha are aligned to the previous frame by polar overlap. The predictor cannot commit a trial gauge. Instantaneous KLI is an adiabatic model, not time-dependent OEP. Frozen mode holds the GS correction fixed; the constant GS SIC energy recorded there is not a recomputed instantaneous SIC energy functional.

**Run** from the new source root. Set CHECKER to the parent workspace's tools/check_forbidden_eigensolvers.py, BUILD to a fresh external build directory, and SCRATCH to a nonexistent temporary directory. Example macOS GNU/OpenMPI configuration:

```sh
python3 "$CHECKER" --source-root "$PWD/src"
cmake -S . -B "$BUILD" -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc \
  -DUSE_MPI=ON -DUSE_HSE=ON -DUSE_LIBXC=OFF -DCMAKE_PREFIX_PATH=/opt/homebrew \
  -DCMAKE_Fortran_FLAGS="-O2 -fallow-argument-mismatch -ffree-line-length-none"
cmake --build "$BUILD" --target salmon -j8
# After the post-link object audit passes, record the validated executable digest.
shasum -a 256 "$BUILD/salmon" > "$BUILD/validated.sha256"
python3 testsuites/kli_static/formulas.py --checker "$CHECKER" \
  --fftw-includes=-I/opt/homebrew/include --fftw-libs='-L/opt/homebrew/lib -lfftw3' \
  --lapack-flags='-L/opt/homebrew/opt/openblas/lib -lopenblas'
python3 testsuites/kli_static/smoke.py "$BUILD/salmon" "$(cut -d' ' -f1 "$BUILD/validated.sha256")" "$SCRATCH" --ranks 1
```

Repeat with a fresh SCRATCH and --ranks 2 for k decomposition. smoke.py requires NumPy. It records run times, invokes each calculation with one OpenMP thread, and verifies completion and the physical state count. The independent diagonal-KLI probe checks D=(+1/8,-1/8), the potential, rejection of additional null modes and NaN/Infinity. The normalized Gaussian probe checks the central potential and energy with the center both inside the box and on its boundary, plus NaN/Infinity rejection.

**Pass criteria.** LSDA energy within the 0.001 eV regression 131 tolerance and the occupied splitting within 4e-5 Ha. The static SCF density threshold is 1e-9; the pair stationarity threshold is 1e-5 Ha and the algebra/normalization threshold 1e-8. In RT, allow at most 1e-8 in the maximum absolute zero-field density/current/potential and occupied orthogonality defects; this is a smoke bound above the GS residual and finite-step error. Require 32 zero-field steps, no extra PZ rotations, a nonzero weak-pulse response, and maximum matched-time current change under 1% of the finer-step peak at 0.04 versus 0.02 a.u. The pulse has peak intensity 1e8 W/cm² and duration 1.28 a.u. These measurements do not establish convergence in k, real-space resolution or optical spectrum. Fujitsu, spin-asymmetric cases and large grids require separate validation.

The driven instantaneous case uses the same GS restart and a 1e14 W/cm² pulse to exercise conditional PZ rotations. Each predictor and each corrector reports the number of rotations, the PZ residual after gauge transport, the residual after minimization, and the minimum transport singular value. The test checks the decision for every evaluation against the 1e-5 Ha threshold, checks that both branches occur in each stage, and checks the final residual and occupied orthogonality.

## Optional representative k points

For the Si primitive cell and Cartesian z polarization used here, copy
`sym.dat.si-xy` to `sym.dat` in **each** run directory, and set
`yn_symmetry='yyy'` in `&system` in the LSDA start, static KLI and RT
inputs. This file explicitly lists identity and the x/y mirror (which
interchanges primitive axes 2 and 3). The default `yn_symmetry='n'` uses
the complete mesh. Both routes use the same occupied-only Wannier solver.

The mirror reduces the standard shifted 4×4×4 mesh from 64 points to 40.
Density and current use the star weights and SALMON's symmetry operations.
KLI expands the occupied orbitals and, during RT, the Hamiltonian actions
with the grid permutation and Bloch phase before gauge transport and
the full-mesh Wannier transform. The existing map accepts grid-compatible
signed permutations of primitive axes and commensurate translations that
preserve the crystal. A single applied field must be invariant under every
retained operation. The potential is averaged over the group only when its
maximum deviation from that average is within `kli_local_tol` (Ha).
Reduced checkpoints store the operations and full k ordering. Use the
same reduction throughout GS and RT.

In the KLI evaluation, each Wannier orbital's forward transform,
self-Hartree solve, and reverse transform run on its owner rank. The
reverse-transformed action is broadcast to form the KLI system, and
the density and Slater checks sum the orbital owners' contributions.

Static KLI always checks the enlarged self-Hartree box. For instantaneous RT the diagnostic is off by default. Set `yn_kli_rt_box_check='y'` in `&functional` to repeat it at each predictor and corrector. The final box-audit column is -1 when skipped.

Input energy tolerances (`kli_local_tol`, `kli_hartree_tol`) follow
`unit_system`; their defaults are equivalent to 1e-5 and 1e-3 Ha,
respectively. When migrating an explicit Hartree value in an existing
`A_eV_fs` input, multiply it by 27.21138505 to specify the same value in eV.
The acceptance threshold `kli_solve_tol` retains its numerical convention.
The 463_kli_input_units CTest covers both a supplied eV value and a physical
default. The kli_static_formulas CTest runs the two numerical probes when
installed FFTW and LAPACK are available in a native HSE build.
