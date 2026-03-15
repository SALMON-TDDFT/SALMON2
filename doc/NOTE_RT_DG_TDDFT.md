# RT-DG-TDDFT Implementation Note

## Scope

This note is a developer-oriented overview of the real-time DG-Fragment TDDFT implementation in SALMON.

The focus is not user input syntax but internal structure:

- what state is propagated,
- how the real-time loop is organized,
- how fragment basis data and MPI distribution are handled,
- where density, Hamiltonian, and observables are updated,
- how adaptive basis updates and mixed fragment+PW space fit into the flow.

Relevant entry points:

- [src/rt/main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90)
- [src/rt/rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90)
- [src/rt/rt_dg_fragment_soi.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_soi.f90)
- [src/rt/rt_dg_fragment_types.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_types.f90)

## Method Summary

The DG-Fragment RT path replaces direct real-space time propagation of Kohn-Sham orbitals with propagation in a fragment basis coefficient space.

Conceptually, the time-dependent orbitals are expanded as

```text
|psi_n(t)> = sum_i c_in(t) |phi_i>
```

where:

- `|phi_i>` are fragment-local basis functions, usually imported from a prior DC-LCFO calculation,
- `c_in(t)` are the propagated coefficients,
- the Hamiltonian is represented in this basis, optionally augmented by plane-wave components.

The implementation is currently in velocity gauge form. The module header in [rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90) summarizes it as

```text
H(t) = H_0 - i A(t) · nabla + A(t)^2 / 2
```

with additional nonlocal pseudopotential contributions handled separately.

## High-Level Execution Flow

The DG-Fragment RT branch is selected in [main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90#L113).

The top-level flow is:

1. `main_tddft` checks `yn_dg_fragment_rt`.
2. `time_evolution_dg_fragment(...)` is called instead of the conventional orbital RT loop.
3. `init_dg_fragment_rt(...)` initializes fragment data structures and imports basis data.
4. `calculate_hamiltonian_matrix(...)` builds the initial coefficient-space Hamiltonian.
5. The RT loop calls `tddft_dg_fragment_iteration(...)` every step.
6. Observables are mapped back into the usual RT output path where possible.
7. `finalize_dg_fragment_rt(...)` releases fragment data.

This keeps DG-Fragment RT inside the same global RT driver while replacing the propagated state representation and the internal update logic.

## State Representation

The core state container is `type(s_dg_fragment_rt)` in [rt_dg_fragment_types.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_types.f90).

Important groups of members are:

- Basis coefficients
  - `coef(:,:,:)`
  - `coef_new(:,:,:)`
  - `coef_work(:,:,:)`
- Basis-space matrices
  - `H_mat`, `H_mat_c`
  - `S_mat`, `S_mat_prop`
  - `S_mat_c`, `S_mat_prop_c`
- Basis metadata
  - `n_basis`
  - `index_basis`
  - `n_mat`, `n_mat_max`
- Real-space fragment basis data
  - `phi_frag`
  - `phi_frag_c`
- Fragment geometry / MPI ownership
  - `n_frag`
  - `num_fragment`
  - `nxyz_domain`
  - `ixyz_frag`
  - `ifrag_start`, `ifrag_end`
  - `id_array`
  - `halo(:)`, `n_halo`
- Self-consistent RT fields on the global grid
  - `rho_frag`
  - `rho_s_frag`
  - `Vh_frag`
  - `Vxc_frag`
- Adaptive basis update state
  - `H_mat_old`
  - `H_mat_kinetic`
  - `hamiltonian_change_norm`
  - `basis_update_threshold`
  - `nbasis_update_count`
  - `last_basis_update_step`
- Mixed fragment + plane-wave basis
  - `use_plane_wave_basis`
  - `n_plane_waves`
  - `coef_pw`
  - `H_mat_mixed`
  - `S_mat_frag_pw`

In practice, the minimum mental model is:

- `phi_frag` is the static or adaptively updated fragment basis in real space,
- `coef` is the propagated electronic state in that basis,
- `H_mat` is the coefficient-space Hamiltonian built from current fields and potentials.

## Initialization Path

Initialization is performed in [rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90#L114) by `init_dg_fragment_rt`.

The main tasks are:

1. Set global sizes and communicator
   - `n_frag = product(num_fragment)`
   - `nstate_frag`
   - `nstate_tot = system%no`
   - `dg_frag%icomm = info%icomm_rko`
2. Enforce the current parallel restriction
   - `np <= n_frag`
3. Decide whether fragment basis data are loaded from DC-LCFO
   - `yn_dg_fragment_from_dcdft`
4. Distribute fragments over MPI ranks
   - `distribute_fragments(...)`
5. Read fragment geometry and basis data
   - `read_fragment_basis_data(...)`
6. Initialize halo communication if real-space basis functions are available
   - `init_halo_communication(...)`
   - `exchange_phi_frag_halo(...)`
7. Initialize integrator coefficients and working arrays

The expected data source for production runs is the DC-LCFO fragment data directory:

```text
./data_dcdft/fragments/
```

Without this imported basis, the DG RT path is not the intended working mode.

## Fragment Distribution And MPI

The fragment ownership logic is in [rt_dg_fragment_halo.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_halo.f90).

The current distribution policy is block distribution over fragment index:

```text
fragment indices -> contiguous chunks -> MPI ranks
```

Key routines:

- `distribute_fragments`
- `get_fragment_range_for_rank`
- `get_fragment_owner_rank`

Each rank owns

```text
ifrag = ifrag_start : ifrag_end
```

and `id_array(ifrag)` stores the owner rank for each fragment.

### Halo exchange

When `phi_frag` is represented in fragment-local real space with ghost cells, neighboring fragment boundaries are synchronized by halo exchange.

Key points:

- halo width is currently tied to the real-space stencil
  - `nxyz_buffer = 4` for the 4th-order stencil
- periodic wrapping is assumed in halo neighbor construction
- each halo entry stores
  - source/destination ranks,
  - source/destination fragment indices,
  - direction vector,
  - extents and displacements,
  - send/receive buffers

This machinery is required for fragment-boundary consistency when basis functions are used in local real-space operators.

## Real-Time Loop

The DG RT loop is implemented in [main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90#L289).

At each time step:

1. `tddft_dg_fragment_iteration(...)` advances coefficients.
2. DG-side observables are stored in `dg_frag`.
3. `energy%E_tot` is updated from `dg_frag%total_energy`.
4. If single-scale Maxwell coupling is enabled, the DG current is fed into the Maxwell update path.
5. Standard RT output writers are reused as much as possible.

The design intent is:

- coefficient-space propagation stays inside DG modules,
- global SALMON RT infrastructure still handles external-field coupling, Maxwell coupling, and common outputs.

## One-Step DG Iteration

The exact integrator path depends on `time_integrator_dg_fragment` and on whether SOI is active.

The per-step logic conceptually contains:

1. determine the current external field / `Ac_tot`,
2. propagate `coef`,
3. reconstruct density from fragment coefficients,
4. solve Hartree and XC using the updated density,
5. rebuild the coefficient-space Hamiltonian,
6. compute observables,
7. optionally trigger basis update.

The self-consistent density/Hamiltonian stage update is factored in [rt_dg_integrator_stage_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_integrator_stage_update.f90) via `update_density_hamiltonian_stage`.

That routine performs:

```text
calculate_density_from_fragments
-> hartree
-> exchange_correlation
-> (+U if active)
-> reconstruct_hamiltonian_matrix
-> mixed-basis update if PW mixing is enabled
```

This is the key self-consistent bridge from coefficient space back to grid-based KS ingredients.

## Density Reconstruction

The DG method does not directly propagate full real-space orbitals on the global grid. To update KS potentials, the electron density must be reconstructed from fragment coefficients and basis functions.

The reconstruction path is entered from `calculate_density_from_fragments(...)`.

Practical meaning:

- `coef` carries time dependence,
- `phi_frag` carries spatial structure,
- `rho(r,t)` is rebuilt from the combination of both.

Once `rho` and `rho_s` are available on the SALMON grids, the existing Hartree and XC infrastructure can be reused.

## Hamiltonian Construction

The initial Hamiltonian is built right after DG initialization in [main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90#L285).

During RT, the Hamiltonian is refreshed through `reconstruct_hamiltonian_matrix(...)`.

The Hamiltonian contains several contributions:

- kinetic term
- local potentials
  - `Vh`
  - `Vxc`
  - `Vpsl`
- velocity-gauge coupling through `Ac_tot`
- nonlocal pseudopotential contributions
- optional HSE / RI pieces
- optional fragment-PW coupling when the mixed basis is active

The code distinguishes:

- `H_mat_kinetic`
  - cached kinetic contribution
- `H_mat_old`
  - previous step Hamiltonian for adaptive-basis trigger logic
- `H_mat_mixed`
  - enlarged Hamiltonian for fragment + PW space

## Time Integrators

`time_integrator_dg_fragment` is parsed into an internal integer code:

- `ssprk3`
- `aetrs`
- `rk4`

and stored in `dg_frag%time_integrator`.

The implementation goal is to keep the propagation in coefficient space while preserving a conventional RT-TDDFT-like self-consistent update pattern around it.

In practice, when tracing integrator-specific behavior, start from:

- [rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90)
- [rt_dg_fragment_soi.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_soi.f90)
- [rt_dg_integrator_stage_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_integrator_stage_update.f90)

## Adaptive Basis Update

Adaptive basis update is one of the main extensions in this branch.

The purpose is to detect when the initial fragment basis has become insufficient under strong driving and then refresh the basis.

The logic is documented extensively in [rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90) and implemented further in [rt_dg_fragment_basis_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_basis_update.f90).

### Trigger

The trigger criterion is based on the Frobenius norm of the Hamiltonian change:

```text
||H_new - H_old||_F
```

see `check_hamiltonian_change_fragments(...)`.

This is accumulated rank-locally, globally reduced, and compared against

```text
basis_update_threshold
```

so that the decision is MPI-size independent.

### Update modes

Two update modes exist:

1. DC-CG / DC-LCFO-based update
   - preferred
   - can expand the basis space
2. Simple diagonalization fallback
   - rotates within the existing space
   - does not enlarge the basis

The selection is controlled by `yn_dc_cg_basis_update`.

### Projection

After a basis refresh, the old propagated state must be mapped to the new basis.

That projection machinery is split out into:

- [rt_dg_basis_projection.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_basis_projection.f90)
- [rt_dg_fragment_basis_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_basis_update.f90)

The core idea is overlap-based projection of old coefficients onto the new basis.

## Mixed Fragment + Plane-Wave Basis

The branch also supports adding plane waves to the fragment basis for metallic or delocalized components.

Important state:

- `use_plane_wave_basis`
- `n_plane_waves`
- `coef_pw`
- `S_mat_frag_pw`
- `H_mat_mixed`

Mixed-basis support is integrated into the self-consistent stage update:

```text
compute_fragment_pw_overlap
compute_fragment_pw_hamiltonian
build_mixed_hamiltonian
```

This path is especially important when a purely fragment-local basis is too restrictive for delocalized or conduction-like states.

## Observables

DG-Fragment RT computes observables in coefficient space and then feeds them back into the usual SALMON output path as much as possible.

Important quantities stored in `dg_frag`:

- `dipole_moment(3)`
- `current(3)`
- `total_energy`

In [main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90#L312), these are written through existing RT output routines.

Current limitations:

- response-spectrum postprocessing hooks are still marked as not implemented for DG RT
- some newer observables added for the conventional RT path are not yet connected to DG RT

## Single-Scale Maxwell Coupling

DG-Fragment RT is wired into the single-scale Maxwell path.

When

```text
theory == 'single_scale_maxwell_tddft'
```

the DG current is forwarded into:

- `fdtd_singlescale(...)`, or
- `fourier_singlescale(...)`

depending on the single-scale mode.

This means:

- DG RT provides the electronic current,
- the Maxwell solver still lives in the existing single-scale infrastructure,
- field feedback remains compatible with the rest of SALMON's RT machinery.

## File Map

The minimum file map to keep in mind is:

- Driver / entry
  - [src/rt/main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90)
- Main DG RT implementation
  - [src/rt/rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90)
- SOI variant
  - [src/rt/rt_dg_fragment_soi.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_soi.f90)
- Core DG data structures
  - [src/rt/rt_dg_fragment_types.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_types.f90)
- Halo and fragment ownership
  - [src/rt/rt_dg_fragment_halo.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_halo.f90)
- Density / Hamiltonian stage update
  - [src/rt/rt_dg_integrator_stage_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_integrator_stage_update.f90)
- Basis update
  - [src/rt/rt_dg_fragment_basis_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_basis_update.f90)
- Basis projection between old/new fragment sets
  - [src/rt/rt_dg_basis_projection.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_basis_projection.f90)
- Fragment-space operators and current
  - [src/rt/rt_dg_fragment_ops.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_ops.f90)
- Plane-wave mixing
  - [src/rt/rt_dg_plane_wave.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_plane_wave.f90)

## Current Assumptions And Restrictions

At the current stage, developers should assume the following unless they verify otherwise in code:

- DG-Fragment RT expects fragment basis data from DC-LCFO in normal use.
- The MPI model is fragment-parallel and requires `np <= n_frag`.
- Halo handling assumes the fragment geometry metadata are correct and periodic wrapping is available.
- Some observables remain less complete than in the conventional RT path.
- DFT+U support is present as framework, but the comment in initialization explicitly warns that full RT coupling is not yet fully wired.
- Mixed fragment+PW support exists but should still be treated as an advanced path needing careful validation.
- Adaptive basis update is implemented, but practical stability depends on threshold choice and on the quality of the DC-CG refresh path.

## Practical Reading Order

For someone modifying the DG RT implementation, the most efficient reading order is:

1. [src/rt/main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90)
2. [src/rt/rt_dg_fragment_types.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_types.f90)
3. [src/rt/rt_dg_fragment.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90)
4. [src/rt/rt_dg_integrator_stage_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_integrator_stage_update.f90)
5. [src/rt/rt_dg_fragment_halo.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_halo.f90)
6. [src/rt/rt_dg_fragment_basis_update.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment_basis_update.f90)

That order is usually enough to understand where to patch:

- initialization,
- per-step propagation,
- density/Hamiltonian refresh,
- MPI fragment ownership,
- basis refresh and projection.
