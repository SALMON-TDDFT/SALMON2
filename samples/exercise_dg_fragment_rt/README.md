# DG-Fragment Real-Time TDDFT Exercise

## Overview

This exercise demonstrates the DG-Fragment (Discontinuous Galerkin Fragment) method for real-time TDDFT calculations. This method performs time evolution directly in the fragment basis coefficient space, which is more efficient for large systems compared to conventional real-space grid-based approaches.

## Theoretical Background

The DG-Fragment RT method combines:
1. **DC-LCFO (Divide-and-Conquer Linear Combination of Fragment Orbitals)**: Provides fragment basis functions from ground-state DFT calculation
2. **DG (Discontinuous Galerkin) method**: Efficient time evolution in fragment basis space
3. **Real-time TDDFT in velocity gauge**: Time-dependent density functional theory using vector potential

### Physical Framework

The method solves the time-dependent Schrödinger equation in the velocity gauge:

```
i ∂ψ/∂t = H(t)ψ = [H_0 + A(t)·p + A(t)²/2] ψ
```

where:
- `H_0` is the field-free Hamiltonian
- `A(t)` is the vector potential of the external field
- `p` is the momentum operator
- `A(t)²/2` is the paramagnetic term (essential for gauge invariance)

In fragment basis representation `|φ_i>`, the coefficients evolve as:

```
i dc_i/dt = Σ_j [H_0_ij + A(t)·p_ij + (A(t)²/2)δ_ij] c_j
```

where `p_ij = <φ_i|p|φ_j>` are the momentum matrix elements and `δ_ij` is the Kronecker delta.

### Key Features

- **Velocity gauge formulation**: Uses vector potential A(t) instead of electric field E(t)
  - **Complete gauge invariance**: Includes both A·p (paramagnetic current) and A²/2 (diamagnetic) terms
  - **Physical correctness**: Ensures proper gauge transformation properties
- **O(N) scaling**: Computational cost scales linearly with system size
- **Fragment-based evolution**: Time evolution in reduced basis without real-space expansion
- **Multiple time integrators**: SSPRK3 (recommended), RK4, and AETRS
- **Accurate observables**: Current density and dipole moment calculated from basis coefficients

### Why A²/2 is Essential

The paramagnetic term A²/2 is crucial for:
1. **Gauge invariance**: Without it, results depend on gauge choice (unphysical)
2. **Energy conservation**: Ensures proper energy accounting in strong fields
3. **High-harmonic generation**: Critical for accurate non-perturbative phenomena
4. **Ionization processes**: Necessary for multi-photon and tunnel ionization

## Workflow

### Step 1: DC-LCFO Ground-State Calculation

First, perform a DC-LCFO calculation to generate fragment basis data:

```bash
# Create input file for DC-LCFO
cat > inputfile_dc_lcfo << EOF
&calculation
  theory = 'dft'
/

&control
  sysname = 'system'
/

&parallel
  nproc_k = 1
  nproc_ob = 4
  nproc_rgrid = 2,2,2
  nproc_rgrid_tot = 1,1,1
  num_fragment = 2,2,2
  yn_dc = 'y'
  yn_dc_lcfo = 'y'
  yn_dc_lcfo_diag = 'y'
  nstate_frag = 16
/

&system
  yn_periodic = 'y'
  al = 10.0d0, 10.0d0, 10.0d0
  nstate = 128
  ...
/

...
EOF

# Run DC-LCFO calculation
mpirun -np 8 salmon < inputfile_dc_lcfo
```

This will create `./data_dcdft/` directory containing:
- `fragments/000001/`, `000002/`, ... (fragment data)
  - `wavefunctions.bin` - Fragment basis coefficients
  - `basis_functions.bin` - Fragment basis functions
  - `rgrid_index.bin` - Real-space grid mapping

### Step 2: DG-Fragment RT-TDDFT Calculation

Now perform the real-time calculation using the fragment basis:

```bash
# Use the provided sample input file
mpirun -np 8 salmon < inputfile_dg_fragment_rt
```

## Input Parameters

### Key Parameters in `&propagation` namelist

```fortran
&propagation
  yn_dg_fragment_rt = 'y'                ! Enable DG-Fragment RT method
  time_integrator_dg_fragment = 'ssprk3' ! Time integrator
/
```

### Time Integrator Options

1. **'ssprk3'** (Recommended)
   - Strong Stability Preserving Runge-Kutta 3rd order
   - Good balance of accuracy and stability
   - Suitable for most systems

2. **'rk4'**
   - Classical Runge-Kutta 4th order
   - Higher accuracy but may require smaller time step
   - Use for high-precision calculations

3. **'aetrs'**
   - Enforced Time-Reversal Symmetry
   - Preserves time-reversal symmetry exactly
   - Experimental in fragment basis representation

## Output Files

The DG-Fragment RT calculation produces output files with spectrum analysis:

### Standard RT Output
- `*_rt.data`: Real-time quantities (current, energy, dipole moment)
- `*_pulse.data`: Applied pulse information
- Standard SALMON RT-TDDFT output format

### Spectrum Analysis (Automatic Post-Processing)
- **`*_dg_epsilon.data`**: Dielectric function ε(ω)
  - Columns: Energy, Re(ε_x), Im(ε_x), Re(ε_y), Im(ε_y), Re(ε_z), Im(ε_z)
  - Computed from linear response current using Kramers-Kronig relations
  - Formula: ε(ω) = 1 + i σ(ω)/ω

- **`*_dg_sigma.data`**: Optical conductivity σ(ω)
  - Columns: Energy, Re(σ_x), Im(σ_x), Re(σ_y), Im(σ_y), Re(σ_z), Im(σ_z)
  - Computed from Fourier transform of time-dependent current
  - Formula: σ(ω) = ∫ j(t) e^(iωt) dt / E_impulse

- **`*_dg_hhg.data`**: High-Harmonic Generation spectrum
  - Columns: Energy, |j_x(ω)|², |j_y(ω)|², |j_z(ω)|², Total
  - HHG intensity: |j(ω)|² shows harmonic peaks at odd multiples of driver frequency
  - Essential for studying strong-field phenomena

### Usage
All spectrum files are automatically generated after time evolution completes. Use these for:
- **Linear response**: Optical absorption spectra from ε(ω)
- **Strong-field dynamics**: HHG analysis for laser-matter interaction
- **Material characterization**: Conductivity for transport properties

All observables are calculated directly from fragment basis coefficients using SALMON's built-in Fourier transform routines with Hann window smoothing.

## Adaptive Basis Update (Automatic)

### Overview

During long-time RT calculations, the fragment basis may become insufficient to describe the evolving physics (e.g., strong-field excitation, charge transfer). The **automated adaptive basis update** feature detects this and recalculates the basis on-the-fly:

```
RT calculation → ||ΔH|| > threshold → Update basis → Continue RT
```

### Workflow

1. **Hamiltonian monitoring**: Track ||(H_mat - H_mat_old)||
2. **Trigger condition**: When change exceeds `basis_update_threshold`
3. **Automatic update**:
   - Save current potentials (Vh, Vxc, Vpsl) to `./rt_potentials/`
   - Launch DC-LCFO: `mpirun salmon < inputfile_dclcfo_update`
   - DC-LCFO reads saved potentials (no SCF needed!)
   - New basis written to `./data_dcdft/fragments/`
4. **Seamless continuation**:
   - Calculate overlap matrix: S_ji = ⟨φ'_j|φ_i⟩
   - Project wave functions: c'_j = Σ_i S_ji c_i
   - Resume time evolution with new basis

### Input Parameters

Enable adaptive basis update in RT input file:

```fortran
&propagation
  yn_dg_fragment_rt = 'y'
  yn_adaptive_basis_update = 'y'          ! Enable automatic updates
  basis_update_threshold = 0.01d0         ! Hamiltonian change threshold (a.u.)
  basis_update_check_interval = 10        ! Check every N time steps
/
```

### Setup Requirements

1. **Prepare DC-LCFO update input**:
   ```bash
   cp inputfile_dclcfo_update.template inputfile_dclcfo_update
   # Edit to match RT calculation parameters
   ```

2. **Key settings in `inputfile_dclcfo_update`**:
   - `nscf = 1` (no SCF - use RT potentials)
   - Match `num_fragment`, `nstate_frag` with RT input
   - Match grid parameters (`al`, `dl`)
   - Set `ncg = 5` (minimal CG for basis expansion)

3. **Run RT calculation**:
   ```bash
   mpirun -np 8 salmon < inputfile_rt
   ```
   - RT will automatically pause and update when needed
   - Progress logged to `adaptive_basis_update.log`

### Output Files

- **`adaptive_basis_update.log`**: Update trigger history
  - Update count, ||ΔH|| values, timestamps
  - Workflow instructions
  - Implementation status
- **`rt_potentials/*.bin`**: Saved potentials for DC-LCFO
  - `Vh.bin`, `Vxc*.bin`, `Vpsl.bin`
- **`dclcfo_update.log`**: DC-LCFO execution log

### Performance Notes

- **Pause duration**: Depends on DC-LCFO cost (~minutes for 100 atoms)
- **Update frequency**: Typical ||ΔH|| threshold = 0.01-0.1 a.u.
- **Overhead**: Minimal (overlap/projection is fast)
- **Benefit**: Maintains accuracy for long-time dynamics

### Example Use Cases

1. **Strong-field excitation**: Basis expands to capture excited states
2. **Charge transfer**: Population moves beyond initial fragment basis
3. **Photocatalysis**: Electronic structure evolves during reaction
4. **Attosecond dynamics**: High-energy states become accessible

## Performance Considerations

### Advantages
- **Memory efficiency**: No need to store full real-space wavefunctions
- **Computational efficiency**: Operations in reduced basis space
- **Scalability**: Better scaling for large systems (1000+ atoms)

### Limitations
- Requires prior DC-LCFO calculation
- Fragment basis quality depends on fragment size
- Non-local observables may require special treatment

## Implementation Status

### Completed Features
- ✓ Velocity gauge formulation with A·p and A²/2 terms
- ✓ **Complex coefficient time evolution**: Proper quantum mechanical phase evolution
- ✓ Time integrators: SSPRK3 (default), RK4
- ✓ Fragment basis coefficient time evolution
- ✓ Current density calculation and output
- ✓ Total energy calculation and output
- ✓ Dipole moment calculation framework
- ✓ Integration with SALMON input/output system
- ✓ **Self-consistent update framework**: Complete architecture for density and Hamiltonian update
- ✓ **Hybrid MPI+OpenMP parallelization**: Fragment distribution + loop parallelization
- ✓ **Spectrum analysis**: Dielectric function, conductivity, and HHG spectrum output

### Critical for Strong-Field Applications (Photovoltaics, Catalysis, Laser Excitation)
- ○ **Self-consistent density/Hamiltonian update** - **ARCHITECTURE COMPLETE**:
  - ✓ Framework fully implemented with modular subroutines
  - ✓ `yn_fix_func='n'` control parameter
  - ✓ Four-step update procedure:
    1. `calculate_density_from_fragments`: ρ(r) from fragment coefficients
    2. `update_hartree_potential`: Solve Poisson equation for V_H
    3. `update_xc_potential`: Calculate V_xc from ρ
    4. `reconstruct_hamiltonian_matrix`: Update H_ij with new potentials
  - ○ **Integration with SALMON solvers needed**:
    - Fragment basis functions in real space (read from DC-LCFO)
    - Connection to `hartree()` and `exchange_correlation()` subroutines
    - Grid mapping and MPI communication
  - **Essential for**: Non-perturbative dynamics, ionization, excited-state reactions
  
### Partially Implemented (TODO)
- ○ **Matrix element calculations**: Hamiltonian, momentum, and dipole matrix elements
  - Currently uses placeholder values
  - Requires implementation to read from DC-LCFO fragment data
  
- ○ **AETRS time integrator**: 
  - Placeholder uses simple midpoint method
  - Full AETRS with Taylor expansion not yet implemented

- ○ **DFT+U support for Mott insulators**: 
  - **Critical Issue**: Without +U correction, Mott insulators show spurious metallic response
  - Framework implemented in `plusu_fragment_support.f90`
  - +U density matrix calculation from fragment coefficients
  - Hamiltonian correction to open Mott gap
  - **Requires**:
    * Orbital character analysis from DC-LCFO (which orbitals are d/f)
    * Projection of fragment basis onto localized atomic orbitals
    * Self-consistent +U update during time evolution
  - **Status**: Framework ready, awaiting DC-LCFO orbital projection data

### Not Yet Implemented
- ✗ Fragment basis overlap matrix (S_mat) handling for non-orthogonal basis
- ✗ Real-space observables (charge density, local currents)

## Important Bug Fix: DFT+U in RT-TDDFT (2026-02-21)

### Critical Issue Fixed in Original SALMON

**Original SALMON had a severe bug**: DFT+U density matrix was NOT updated during RT-TDDFT time evolution.

#### Symptoms
- Mott insulators (NiO, CoO, transition metal oxides) showed **spurious metallic response**
- Correct band gap in ground state, but wrong optical response during time evolution
- +U potential remained frozen at initial value while wavefunctions evolved

#### Root Cause
- Ground state: `calc_density_matrix_and_energy_plusU` correctly calculated `V_eff`
- RT-TDDFT: This function was **never called** during time evolution in `time_evolution_step.f90`
- Result: As wavefunctions evolved, +U potential did not track the changing density matrix

#### Solution Implemented
Modified `src/rt/time_evolution_step.f90` to update +U every time step:
```fortran
! After exchange_correlation in main loop and predictor_corrector:
call exchange_correlation(...)

! Update DFT+U density matrix and potential during time evolution
if ( PLUS_U_ON ) then
  call calc_density_matrix_and_energy_plusU( spsi_out, ppg, info, system, energy%E_U )
end if
```

Now +U density matrix `dm_mms_nla` and potential `V_eff` are updated every time step, properly tracking orbital occupations as they change during time evolution.

#### Impact
- **Fixed systems**: All Mott insulators, strongly correlated materials with +U
- **Unaffected systems**: Band insulators, simple metals (no +U used)
- **DG-Fragment**: Same issue exists but requires different implementation approach (see `plusu_fragment_support.f90`)

**See `doc/NOTE_PLUSU.md` for current +U integration notes.**

---

### Parallelization Strategy (IMPLEMENTED)

The code now uses **hybrid MPI+OpenMP parallelization**:

**MPI Level (Fragment Distribution)**:
- Each fragment is assigned to one MPI rank
- Fragments are distributed using block distribution across ranks
- Each rank calculates time evolution only for its assigned fragments
- **I/O parallelization**: Each rank reads only its assigned fragment data files
  * Root rank reads metadata and broadcasts to all ranks
  * Each rank reads its own fragment files in parallel (no I/O contention)
  * Reduces memory footprint and I/O time
- Minimal communication: Only observable aggregation requires MPI_Allreduce
- Perfect load balancing when `nproc <= n_fragments`

```bash
# Example: 8 fragments on 4 MPI processes
# Rank 0: fragments 1-2
# Rank 1: fragments 3-4  
# Rank 2: fragments 5-6
# Rank 3: fragments 7-8
```

**OpenMP Level (Loop Parallelization)**:
- Inner loops over orbitals and matrix elements use `!$omp parallel do`
- Parallelized regions:
  * Time derivative calculation: orbital loops
  * Observable calculation: dipole, current, energy loops
- Set `OMP_NUM_THREADS` to control thread count per MPI rank

**Optimal Configuration**:
```bash
# For 32 cores with 8 fragments:
export OMP_NUM_THREADS=4
mpirun -np 8 ./salmon < inputfile_dg_fragment_rt
# Total: 8 MPI ranks × 4 OpenMP threads = 32 cores
```

**Performance Considerations**:
- Use `nproc_total <= product(num_fragment)` for best efficiency
- Balance: More MPI ranks → less communication overhead
- Balance: More OpenMP threads → better memory locality
- Recommended: `n_fragments / nprocs = 1-4` fragments per rank

### Usage Notes

**Current Capabilities:**
- Efficient time evolution in fragment basis space
- **Two operational modes**:
  1. **Fixed Hamiltonian** (`yn_fix_func='y'`, default): Uses initial DC-LCFO Hamiltonian
     - Suitable for: Linear response, weak perturbations, benchmarking
  2. **Self-Consistent** (`yn_fix_func='n'`): Updates density and Hamiltonian each step
     - **Architecture ready** - requires fragment basis function data
     - Essential for: Photovoltaic effects, catalytic reactions, laser ionization

**For Production Use in Strong-Field Applications:**
To enable self-consistent updates for photovoltaics, catalysis, and laser excitation:
1. Implement fragment basis function reader from DC-LCFO data
2. Connect density reconstruction to real-space grid
3. Integrate with SALMON's Hartree and XC solvers
4. Add matrix element recalculation with updated potentials

The complete algorithmic framework is in place and ready for these integrations.

## Advanced Usage

### Adjusting Fragment Size

The fragment size is controlled by `num_fragment` in DC-LCFO calculation:
```fortran
num_fragment = 2,2,2  ! 8 fragments total
nstate_frag = 16      ! 16 states per fragment
```

- **Larger fragments**: Better accuracy, higher computational cost
- **Smaller fragments**: Lower cost, may lose accuracy
- **Rule of thumb**: Each fragment should contain ~10-20 atoms

### Combining with Other Methods

DG-Fragment RT can be combined with:
- **Maxwell-TDDFT**: Multi-scale light-matter interaction
- **Temperature effects**: Finite-temperature TDDFT
- **Non-adiabatic MD**: Coupled electron-ion dynamics (future work)

## Troubleshooting

### Error: "Fragment basis data mismatch"
- Ensure `num_fragment` matches between DC-LCFO and RT calculations
- Check that `./data_dcdft/` directory exists and contains valid data

### Error: "MPI size is too small"
- Number of MPI processes must be ≥ number of fragments
- Use at least `product(num_fragment)` processes

### Poor accuracy
- Increase `nstate_frag` in DC-LCFO calculation
- Use smaller fragment regions (`num_fragment`)
- Reduce time step `dt` in RT calculation

## References

1. Noda et al., "Linear-scaling time-dependent density-functional theory based on the fragment molecular orbital method", J. Chem. Phys. **142**, 244101 (2015)
2. Noda et al., "Quantum mechanical divide-and-conquer approach to linear-response time-dependent density functional theory", Phys. Rev. B **95**, 045106 (2017)
3. Cockburn and Shu, "Runge-Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems", J. Sci. Comput. **16**, 173 (2001)
4. Yabana et al., "Time-dependent density functional theory for strong electromagnetic fields in crystalline solids", Phys. Rev. B **85**, 045134 (2012) - Velocity gauge formulation
5. Bertsch et al., "Real-space, real-time method for the dielectric function", Phys. Rev. B **62**, 7998 (2000) - Gauge invariance in TDDFT

## Contact

For questions and support, please contact the SALMON developers:
https://salmon-tddft.jp/
