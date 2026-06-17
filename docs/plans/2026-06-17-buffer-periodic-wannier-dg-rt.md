# Buffer-Periodic Wannier DG-RT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a DG-RT path where local Wannier bases are generated from buffered DC fragment wavefunctions with periodic boundary conditions, then RT starts from eigenstates of the Flux-including DG Hamiltonian without keeping full-system density matrices.

**Architecture:** Each fragment constructs a local periodic cell from `core + buffer`, builds an overcomplete local Wannier candidate space including occupied and excited states, removes linear dependence, projects the Flux-DG Hamiltonian into the retained basis, diagonalizes that local Hamiltonian, and uses only center-owned states for RT. Density is reconstructed from owned Wannier states on the buffered real-space support and sent to destination fragments when the Wannier tail lies outside the owner core.

**Tech Stack:** Fortran 90/2003, SALMON DC-LCFO/DG-RT modules, MPI neighbor communication, BLAS/LAPACK/EigenExa, optional Wannier90 MPI executable.

---

### Task 1: Define the New Seed Format

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`

**Plan:**
Add a new file beside `local_wannier_basis.bin`, for example `buffer_periodic_wannier_basis.bin`, rather than changing the current file format in place.

The seed must store:
- fragment id
- core grid shape
- buffer grid shape
- periodic buffer-cell grid shape
- retained basis count
- candidate basis count before pruning
- overlap eigenvalues
- retained candidate indices
- Wannier coefficients in the buffered DG basis
- Wannier centers
- spread/tail diagnostics
- projected `H_flux`
- projected `r`
- projected velocity/current matrix if available

Acceptance:
- Old `local_wannier_basis.bin` path still reads.
- New RT path refuses to run if the seed lacks `H_flux`.

### Task 2: Build Buffered Periodic Candidate Orbitals

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Plan:**
Use the already exported buffered basis support from the DC fragment. Treat `core + buffer` as a small periodic cell when evaluating projection overlaps and position moments. The periodicity is only for local Wannier construction; DG Flux remains the physical boundary coupling.

Candidate basis should include:
- occupied-like local subspace
- enough excited local states from `nstate_frag`
- optional projection families beyond `C:sp3` later

Acceptance:
- Log candidate count, retained count, minimum overlap eigenvalue, maximum spread, and tail outside core.
- For the 2x2x2 diamond case, candidate count must exceed occupied-local count.

### Task 3: Remove Overcompleteness

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Plan:**
Construct the candidate overlap matrix `S_W`. Diagonalize it and keep only modes with eigenvalue above a configurable threshold. Then orthonormalize:

```text
W_orth = W U_keep lambda_keep^{-1/2}
```

Use the same retained basis for `H_flux`, density, and RT coefficients.

Acceptance:
- Abort if retained count is zero.
- Warn if `lambda_min_keep / lambda_max` is close to the threshold.
- Log discarded count and condition number.

### Task 4: Project the Flux-Including DG Hamiltonian

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Possibly modify: `src/rt/dg/rt_dg_fragment_hamiltonian.f90`

**Plan:**
Project the already built DG Hamiltonian, including surface Flux blocks, onto the retained buffered-periodic Wannier basis:

```text
H_W = W_orth^T H_DG_flux W_orth
```

For local RT ownership, diagonalize the fragment-local `H_W`. The initial RT states are occupied eigenvectors of this `H_W`, not raw projection/Wannier columns.

Acceptance:
- `H_W` Hermiticity error is logged and below tolerance.
- Eigenvalues are finite and sorted.
- RT can initialize from identity in the `H_W` eigenbasis.

### Task 5: Assign Owner by Wannier Center

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`

**Plan:**
Map each retained local Wannier/eigenstate center to a core fragment. Only the owner fragment stores time-dependent coefficients for that state. Non-owner images may keep enough real-space shape information for density/tail evaluation if needed, but not a propagated coefficient.

Acceptance:
- Logs show `center_owned` and `center_total`.
- No full-system Wannier density matrix is allocated.
- Coefficient arrays scale with locally owned Wannier count, not total fragment count.

### Task 6: Density From Owned States Plus Tail Exchange

**Files:**
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Possibly modify: `src/rt/dg/rt_dg_fragment_types.f90`

**Plan:**
For each owner fragment:

```text
psi_n(r) = sum_i C_i,n W_i(r)
rho_owner(r) = sum_n occ_n |psi_n(r)|^2
```

Accumulate core-local density directly. For grid points outside the owner core but inside another fragment core, send the density contribution to the destination fragment via point-to-point or neighbor collective communication.

Acceptance:
- No per-Wannier diagonal-only density approximation.
- Cross terms between owned Wannier functions are included.
- Integrated electron count is conserved after tail exchange.

### Task 7: RT Integrator Uses Fixed H During Taylor

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_taylor.f90`
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Plan:**
The Taylor expansion must apply a fixed Hamiltonian for the whole step:

```text
exp(-i H(t) dt)
```

Rebuild density/H/Flux only between steps, never inside Taylor terms.

Acceptance:
- No call path updates Flux or density inside Taylor sub-terms.
- Existing derivative cache remains valid for one step.

### Task 8: Verification Cases

**Files:**
- Test inputs under: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/`

**Plan:**
Run:
- 2x2x2 diamond GS with buffered periodic Wannier seed
- RT no-field: polarization/current should remain near stationary baseline
- impulse RT: `J_tot,z` should drop from initial kick response and oscillate around zero
- basis convergence: vary local excited count and overlap cutoff
- buffer convergence: vary buffer width

Acceptance:
- GS finishes with Flux-DG eigenstate seed.
- No-field drift is small.
- Impulse response is qualitatively insulating.
- Increasing excited local basis improves `J_para` response without overcomplete instability.

---

## Critical Decisions

1. The RT initial state must be the occupied eigenstate set of the Flux-including projected DG Hamiltonian.
2. The local Wannier construction may use periodic boundary conditions on `core + buffer`, but the physical inter-fragment coupling remains the DG Flux term.
3. Overcompleteness is removed before Hamiltonian projection.
4. Density construction includes cross terms in the retained local basis.
5. Full-system density matrices are not allowed in the production path.
