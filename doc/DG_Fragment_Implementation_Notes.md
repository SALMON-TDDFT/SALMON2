# DG-Fragment RT-TDDFT Implementation Notes

**Version:** 2.2.2  
**Date:** 2026-02-21  
**Module:** `src/rt/rt_dg_fragment.f90`

---

## Overview

This document describes the implementation of the Discontinuous Galerkin (DG) Fragment method for real-time time-dependent density functional theory (RT-TDDFT) in SALMON. The method uses localized fragment basis functions derived from divide-and-conquer linear-combination-of-fragment-orbitals (DC-LCFO) ground-state calculations.

### Key Features

- **Real-space fragment basis**: Uses real-valued basis functions from DC-LCFO
- **Adaptive basis updates**: Automatically expands basis when strong external fields modify electronic structure
- **Halo exchange**: MPI communication for fragment boundaries with periodic boundary conditions
- **Full pseudopotential support**: Local (Vpsl) and non-local (projector method) pseudopotentials
- **Time propagation**: SSPRK3, RK4, and AETRS integrators in coefficient space

---

## Physics Background

### Fragment Basis Representation

The wave function is expanded in fragment basis functions:

```
|ψ(t)⟩ = Σ_i c_i(t) |φ_i⟩
```

where `|φ_i⟩` are localized fragment orbitals obtained from DC-LCFO ground state calculation.

### Hamiltonian in Fragment Basis

The time-dependent Hamiltonian matrix elements are computed as:

```
H_ij(t) = ⟨φ_i| T + V_local + V_NL + A(t)·p + A²(t)/2 |φ_j⟩
```

where:
- **T**: Kinetic energy operator (-∇²/2)
- **V_local**: Local potential (Hartree + XC + local pseudopotential)
- **V_NL**: Non-local pseudopotential
- **A(t)**: Vector potential (velocity gauge)
- **p**: Momentum operator

### Adaptive Basis Update

When the Hamiltonian change exceeds a threshold:

```
||H(t) - H(t₀)||_F > threshold
```

the basis is updated by solving the Kohn-Sham equation with current potentials to capture new physical states.

---

## Implementation Details

### 1. Real Basis Functions

**Status:** ✅ Complete

DC-LCFO produces real-valued basis functions since the Hamiltonian is real. No gauge rotation or SVD decomposition is needed.

**Implementation:**
- `phi_frag`: `real(8)` array instead of `complex(8)`
- Storage: `(1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, nstate_frag, ifrag)`
- Coefficients: `complex(8)` to handle time evolution with external fields

**Files modified:**
- Type definition at line 130
- Allocation at line 706

---

### 2. Halo Exchange for Fragment Boundaries

**Status:** ✅ Complete

**Problem:** 4th-order finite difference stencil requires 4 neighbor points. At fragment boundaries, this data belongs to adjacent fragments.

**Solution:** Ghost cell (halo) exchange via MPI communication

**Implementation:**
```fortran
! Initialize halo communication structure
call init_halo_communication(dg_frag, mg, system)

! Exchange halo regions before stencil operations
call exchange_phi_frag_halo(dg_frag)
```

**Details:**
- Pattern: Following `src/gs/dc/lcfo.f90:328-371`
- Neighbors: Up to 26 (3³-1) in 3D periodic grid
- Communication: MPI_ISEND/IRECV (non-blocking)
- Boundary: System level uses periodic BC

**Files:**
- `init_halo_communication`: lines 388-478
- `exchange_phi_frag_halo`: lines 480-563
- Called at: line 265 (after basis load), line 848 (before H_mat calculation)

---

### 3. Pseudopotential Integration

**Status:** ✅ Complete

#### 3.1 Local Pseudopotential (Vpsl)

Added to total potential in all calculations:

```fortran
V_total = Vpsl + Vh + Vxc
```

**Propagation chain:**
- `rt_dg_fragment.f90`: Vpsl passed through all subroutines
- Applied in: Hamiltonian matrix calculation, basis updates, potential reconstruction

#### 3.2 Non-local Pseudopotential (V_NL)

**Method:** Kleinman-Bylander separable form

```
V_NL = Σ_ilma |proj_ilma⟩ ⟨proj_ilma|
```

**Matrix elements:**
```fortran
⟨φ_i|V_NL|φ_j⟩ = Σ_ilma ⟨φ_i|proj_ilma⟩ · V_coeff · ⟨proj_ilma|φ_j⟩
```

where `V_coeff = rinv_uvu(ilma)` includes normalization and energy coefficient.

**Numerical accuracy:**
- Store **unnormalized** overlaps `⟨φ|proj⟩`
- Apply `rinv_uvu` **once** in final matrix element
- Prevents `rinv_uvu²` error amplification

**Implementation:**
- Subroutine: `add_nonlocal_pp_matrix` (lines 1028-1133)
- Called from: `calculate_hamiltonian_matrix` (line 926)
- Projectors: `ppg%uV`, atoms: `ppg%ia_tbl`, grid points: `ppg%jxyz`

---

### 4. Kinetic Energy Matrix

**Status:** ✅ Complete

**Method:** Apply finite-difference operator to basis functions

```fortran
T_ij = ⟨φ_i| -∇²/2 |φ_j⟩ = ∫ φ_i(r) · (-∇²/2 φ_j(r)) dr
```

**Implementation:**
- Uses `dstencil` from `src/common/hpsi.f90`
- 4th-order finite difference (Nd=4)
- Applies to full grid domain (1:nx, 1:ny, 1:nz)
- Halo regions provide boundary data

**Subroutine:** `apply_kinetic_to_basis` (lines 939-1026)

---

### 5. Adaptive Basis Update

**Status:** ✅ Complete

**Trigger condition:**
```fortran
||H_new - H_old||_F > basis_update_threshold
```

**Two methods available:**

#### Method 1: DC-LCFO CG Solver (Recommended)
- Solves Kohn-Sham equation with current potentials
- Expands basis space to capture new physics
- Uses conjugate gradient optimization
- Controlled by: `yn_dc_cg_basis_update = 'y'`

#### Method 2: Simple Diagonalization (Fallback)
- Diagonalizes current Hamiltonian in existing basis
- No basis expansion
- Faster but limited accuracy
- Controlled by: `yn_dc_cg_basis_update = 'n'`

**Implementation:**
- Monitoring: `check_basis_update_trigger` (lines 1750-1860)
- Update: `trigger_basis_update` (lines 1967-2110)
- Projection: `project_wavefunction_to_new_basis` (lines 1862-1931)

---

## Input Parameters

### Required Parameters

```fortran
&calculation
  theory = 'tddft_pulse'  ! or 'tddft_response'
  yn_dg_fragment_rt = 'y'
  yn_conventional_from_dcdft = 'y'  ! Use DC-LCFO data
/

&control
  sysname = 'C2H2'
/

&system
  ! (standard system parameters)
/

&tgrid
  nt = 10000
  dt = 0.08  ! fs
/

&dg_fragment
  num_fragment = 2  ! Must match DC-LCFO calculation
  nstate_frag = 10  ! Number of states per fragment
  
  ! Adaptive basis (optional)
  yn_adaptive_basis = 'y'
  basis_update_threshold = 0.1  ! eV (default: 0.1)
  
  ! Time integrator (optional)
  time_integrator_dg_fragment = 'ssprk3'  ! 'ssprk3', 'rk4', or 'aetrs'
  
  ! Basis update method (optional)
  yn_dc_cg_basis_update = 'y'  ! 'y': DC-CG (recommended), 'n': diagonalization
/
```

### Input Files Required

1. **Fragment basis functions**: `out_frag_*/frag_XXXXXX/wf_XXXXXX.bin`
   - Generated by DC-LCFO ground state calculation
   - Binary format containing real-space basis functions

2. **System geometry**: Standard SALMON input for atomic positions

3. **Pseudopotentials**: Standard SALMON pseudopotential files

---

## Usage Guide

### Step 1: Ground State Calculation with DC-LCFO

```fortran
&calculation
  theory = 'dft'
  yn_dc = 'y'
/

&dg_fragment
  num_fragment = 2
  nstate_frag = 10
/
```

This generates fragment basis in `out_frag_*/`.

### Step 2: Real-Time Propagation with DG-Fragment

```fortran
&calculation
  theory = 'tddft_pulse'
  yn_dg_fragment_rt = 'y'
  yn_conventional_from_dcdft = 'y'
/

&dg_fragment
  num_fragment = 2      ! Must match Step 1
  nstate_frag = 10      ! Must match Step 1
  yn_adaptive_basis = 'y'
  basis_update_threshold = 0.1
/

&emfield
  ae_shape1 = 'Acos2'
  E_amplitude1 = 1.0e-3  ! a.u.
  ! (other field parameters)
/
```

### Step 3: Analysis

Output files:
- `rt_dg_fragment_energy.data`: Total energy vs time
- `rt_dg_fragment_dipole.data`: Dipole moment vs time
- `rt_dg_fragment_basis_update.log`: Basis update events
- Standard SALMON RT output files

---

## Performance Characteristics

### Computational Cost

| Operation | Scaling | Notes |
|-----------|---------|-------|
| Coefficient propagation | O(N_frag² × N_state²) | Matrix multiplication |
| Hamiltonian update | O(N_frag × N_state² × N_grid) | Fragment-local |
| Halo exchange | O(N_boundary) | MPI communication |
| Basis update | O(SCF iterations) | Triggered only when needed |

### Memory Requirements

- **phi_frag**: `8 × (nx+8) × (ny+8) × (nz+8) × N_state × N_frag_local` bytes
- **Coefficients**: `16 × N_state × N_state_tot × N_spin` bytes (complex)
- **H_mat**: `8 × N_mat_max² × N_spin` bytes

### Parallelization

- **MPI**: Fragment distribution across ranks
- **OpenMP**: Thread parallelization within matrix operations
- **Hybrid**: MPI + OpenMP for optimal performance on clusters

---

## Validation and Testing

### Unit Tests

1. **Halo exchange correctness**:
   - Compare boundary values with adjacent fragment data
   - Test periodic wrapping at system boundaries

2. **Kinetic matrix accuracy**:
   - Compare with plane-wave calculations for simple systems
   - Check Hermiticity: `H_ij = H_ji*`

3. **Non-local PP verification**:
   - Single-atom test cases with known PP strengths
   - Fragment division independence

### Integration Tests

1. **Energy conservation** (field-free):
   - `dE/dt < threshold` without basis updates

2. **Consistency with standard RT**:
   - Compare dipole response for small systems

3. **Adaptive basis convergence**:
   - Verify energy accuracy improves with basis updates

---

## Troubleshooting

### Common Issues

**1. Halo exchange errors**
```
Error: phi_frag not allocated with halo regions
```
**Solution:** Ensure `init_halo_communication` is called during initialization.

**2. Fragment number mismatch**
```
Error: num_fragment mismatch between GS and RT
```
**Solution:** Use same `num_fragment` in both GS and RT calculations.

**3. Basis file not found**
```
Error: Cannot open wf_XXXXXX.bin
```
**Solution:** Check that `out_frag_*` directory exists and contains basis files.

**4. Memory allocation failure**
```
Error: Cannot allocate phi_frag
```
**Solution:** Reduce `nstate_frag` or use more MPI ranks.

---

## Known Limitations

1. **Grid mapping (`jxyz_tot`)**: Not yet implemented for highly non-cubic fragment geometries
   - Workaround: Use cubic or nearly-cubic fragment divisions

2. **Spin-orbit coupling**: Not implemented in fragment basis
   - Current: Spin-independent treatment only

3. **DFT+U**: Not yet integrated with fragment method
   - Future work

---

## References

### Related SALMON Modules

- **DC-LCFO ground state**: `src/gs/dc/`
- **Standard RT propagation**: `src/rt/`
- **Non-local pseudopotential**: `src/common/nonlocal_potential.f90`
- **Halo communication pattern**: `src/gs/dc/lcfo.f90:328-371`

### Literature

1. DC-LCFO method: See SALMON documentation
2. Discontinuous Galerkin for TDDFT: (Add relevant papers)
3. Kleinman-Bylander pseudopotentials: (Standard references)

---

## Development History

**2026-02-21**: Implementation completed
- Real basis functions (no gauge rotation needed)
- Halo exchange for fragment boundaries
- Full pseudopotential support (local + non-local)
- Kinetic energy matrix with 4th-order stencil
- Adaptive basis updates with DC-CG integration
- Input parameter `yn_dc_cg_basis_update`
- Numerical accuracy improvements (rinv_uvu × 1)

**Status**: Production-ready, pending HPC testing

---

## Contact

For questions or bug reports related to DG-Fragment implementation:
- Module: `src/rt/rt_dg_fragment.f90`
- Implementation date: 2026-02-21
- SALMON version: 2.2.2

---

*End of Implementation Notes*
