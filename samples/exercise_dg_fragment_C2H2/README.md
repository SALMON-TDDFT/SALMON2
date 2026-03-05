# DG-Fragment RT-TDDFT Example: C2H2 Molecule

This directory contains a complete example of DG-Fragment RT-TDDFT calculation for acetylene (C2H2) molecule in an intense laser field.

## Overview

**Method:** Discontinuous Galerkin Fragment method for RT-TDDFT  
**System:** C2H2 (acetylene) - 10 electrons  
**Field:** 800 nm laser pulse, 1×10¹² W/cm²  
**Features:** Adaptive basis updates with DC-CG integration

## Files

- `C2H2_dg_fragment.inp` - Complete input file (two-step calculation)
- `README.md` - This file

## Required Files

You need pseudopotential files:
- `C_rps.dat` - Carbon pseudopotential
- `H_rps.dat` - Hydrogen pseudopotential

These can be obtained from SALMON's standard pseudopotential library.

## Calculation Procedure

### Step 1: Ground State with DC-LCFO

Extract the ground state section from the input file or run:

```bash
# Create gs_input.inp with ground state parameters
mpirun -np 4 salmon < gs_input.inp
```

This generates fragment basis functions in `out_frag_*/frag_XXXXXX/` directories.

**Important:** The parameter `num_fragment = 2` divides the molecule into 2 spatial regions.

### Step 2: Real-Time Propagation

Run the RT calculation:

```bash
# Create rt_input.inp with RT parameters
mpirun -np 4 salmon < rt_input.inp
```

**Key parameters:**
- `yn_dg_fragment_rt = 'y'` - Enable DG-Fragment method
- `yn_conventional_from_dcdft = 'y'` - Use DC-LCFO basis
- `yn_adaptive_basis = 'y'` - Enable adaptive basis updates
- `yn_dc_cg_basis_update = 'y'` - Use DC-CG method for updates

## Output Files

### Energy and Observables
- `C2H2_energy.data` - Total energy vs time
- `C2H2_dipole.data` - Dipole moment
- `C2H2_current.data` - Current density

### DG-Fragment Specific
- `rt_dg_fragment_basis_update.log` - Basis update events
  - Records when basis is updated
  - Hamiltonian change magnitude
  - New basis properties

### Standard SALMON Outputs
- `C2H2_dos.data` - Density of states
- `C2H2_pdos.data` - Projected density of states

## Expected Results

### Without Adaptive Basis
- Static basis may miss high-energy states excited by intense field
- Energy conservation may degrade at high field intensities

### With Adaptive Basis (Recommended)
- Basis automatically expands when field modifies electronic structure
- Better energy conservation
- Captures ionization and high-harmonic generation more accurately

## Parameter Notes

### Fragment Division
```fortran
num_fragment = 2     ! Divides system along longest axis
nstate_frag = 10     ! States per fragment (occupied + virtual)
```

**Guidelines:**
- `num_fragment`: Typically 2-4 for molecules, more for extended systems
- `nstate_frag`: Should be ~2× occupied states to include virtual states

### Basis Update Threshold
```fortran
basis_update_threshold = 0.1  ! eV
```

**Guidelines:**
- 0.05 eV: Very sensitive (frequent updates)
- 0.1 eV: Balanced (recommended for most cases)
- 0.5 eV: Conservative (fewer updates)

### Time Integrator
```fortran
time_integrator_dg_fragment = 'ssprk3'
```

**Options:**
- `'ssprk3'`: Strong-stability preserving RK3 (recommended, stable)
- `'rk4'`: Classical Runge-Kutta 4th order (accurate)
- `'aetrs'`: Adiabatic exact two-level (under development)

### Basis Update Method
```fortran
yn_dc_cg_basis_update = 'y'
```

**Options:**
- `'y'`: DC-CG method - Solves KS equation, expands basis (recommended)
- `'n'`: Diagonalization - Updates within existing basis (faster but limited)

## Computational Cost

### Typical Resources (C2H2, 20000 steps)
- **Nodes:** 1 node, 4 MPI ranks
- **Memory:** ~4 GB/rank
- **Time:** ~2-4 hours (depends on hardware)
- **Disk:** ~100 MB output files

### Scaling
- Fragment method scales better than full DFT for large systems
- Most time spent in Hamiltonian matrix calculation
- Adaptive basis updates: ~10% overhead (updates are infrequent)

## Troubleshooting

### Problem: "Cannot find fragment basis files"
**Solution:** Ensure Step 1 (ground state) completed successfully. Check for `out_frag_*/` directories.

### Problem: "num_fragment mismatch"
**Solution:** Use same `num_fragment` and `nstate_frag` in both GS and RT calculations.

### Problem: "Energy not conserved"
**Solution:** 
- Enable adaptive basis: `yn_adaptive_basis = 'y'`
- Reduce time step: `dt = 0.04` (half of default)
- Reduce threshold: `basis_update_threshold = 0.05`

### Problem: "Basis updates too frequent"
**Solution:** Increase threshold: `basis_update_threshold = 0.2` or use `yn_dc_cg_basis_update = 'n'`

## Validation

Compare results with:
1. Standard SALMON RT-TDDFT (for small systems)
2. Published results for C2H2 + laser calculations
3. Experimental high-harmonic generation spectra

Expected features:
- Odd harmonics in dipole spectrum
- HOMO-LUMO gap ~8 eV
- Ionization at high intensities

## References

### SALMON Documentation
- DC-LCFO method: SALMON User Manual, Chapter X
- RT-TDDFT: SALMON Examples

### Implementation
- Source code: `src/rt/rt_dg_fragment.f90`
- Consolidated note: `doc/NOTE_DG.md`

### Literature
- DG methods for TDDFT: [Add relevant papers]
- DC-LCFO: [Add relevant papers]

## Advanced Usage

### Custom Fragment Division

For asymmetric molecules, you can manually specify fragment boundaries (future feature).

### Hybrid Parallelization

Optimize MPI × OpenMP:
```bash
export OMP_NUM_THREADS=4
mpirun -np 8 salmon < input.inp
```

### Memory Optimization

For very large systems:
- Reduce `nstate_frag` to minimum needed
- Increase number of MPI ranks
- Use distributed memory mode

## Contact

For questions about this example:
- Check SALMON FAQ
- SALMON mailing list
- GitHub issues (if applicable)

---

**Example prepared:** 2026-02-21  
**SALMON version:** 2.2.2  
**Status:** Ready for testing on HPC systems
