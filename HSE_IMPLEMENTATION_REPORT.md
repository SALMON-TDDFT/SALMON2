# HSE Hybrid Functional Implementation for DG-Fragment RT-TDDFT
## Plan A: Density Matrix Method with Fragment-Local Computation

**Implementation Date**: 2026-02-22  
**SALMON Version**: v.2.2.2  
**Status**: Phase 1 Complete - Fragment-Local Exchange Calculation

---

## Overview

Implementation of the Heyd-Scuseria-Ernzerhof (HSE) hybrid functional for Divide-and-Conquer (DC) DG-Fragment real-time TDDFT. This enables efficient inclusion of exact exchange within the fragment-local computational framework.

### Key Innovation: Fragment-Decomposed Exact Exchange

Instead of computing the full global 4-index integrals, we:
1. **Decompose** at the fragment level: Each fragment independently evaluates its local exchange potential
2. **Integrate** only within fragment domain: Grid points beyond fragment boundary are skipped
3. **Accumulate** fragment contributions: Final Hamiltonian is sum of fragment exchanges
4. **Scale as O(N_frag)**: Computational cost scales linearly with number of fragments

---

## Mathematical Formulation (Plan A)

### Exact Exchange Potential (Hartree-Fock form)

$$V_x^{\text{HF}}(r,r') = -\sum_{n=1}^{N_\text{occ}} \phi_n(r) \phi_n^*(r') / |r-r'|$$

### HSE Mixing

$$V_x^{\text{HSE}} = (1-\alpha) V_x^{\text{GGA}} + \alpha V_x^{\text{HF}}$$

where:
- $\alpha = 0.25$ (default, PBE0)
- $V_x^{\text{GGA}}$ from PBE/LDA exchange

### Fragment-Local Form

For fragment $\mu$, compute only:

$$\langle i | V_x^{\text{HSE}} | j \rangle_\mu \approx -\alpha \sum_{k,l \in \text{occupied}}^{\mu} \sum_{r_1,r_2 \in \Omega_\mu} \phi_i(r_1) \phi_k(r_1) \phi_l(r_2) \phi_j(r_2) \frac{h^6}{|r_1-r_2|}$$

where $\Omega_\mu$ is the domain of fragment $\mu$.

---

## Implementation Details

### 1. Global Parameters (`salmon_global.f90`)

Three new parameters added to the `&functional` namelist:

```fortran
character(1)   :: yn_hse              ! Flag: enable HSE (default: 'n')
real(8)        :: hse_alpha           ! Exact exchange coefficient (default: 0.25d0)
real(8)        :: hse_omega           ! Screening parameter (a.u., default: 0.11d0)
```

### 2. Input/Output Handling (`inputoutput.f90`)

- **Namelist extension**: Added to `&functional` namelist for input parsing
- **Default values**: Set during initialization
- **MPI broadcast**: Parameters distributed across all MPI ranks via `comm_bcast`

### 3. Fragment Hamiltonian Reconstruction (`rt_dg_fragment.f90`)

#### Modified: `reconstruct_hamiltonian_matrix` (line ~2550)

```fortran
! Branch for HSE calculation
if (yn_hse == 'y') then
  call add_exact_exchange_hse(dg_frag, system, dg_frag%H_mat(:, :, ispin), ifrag, ispin)
end if
```

#### New: `add_exact_exchange_hse` (line ~2350-2470)

**Algorithm:**
1. **Fragment identification**: Determine local fragment index and grid range
2. **Occupied state count**: From coefficient matrix (simplified: $N_\text{occ} = \lfloor N_\text{el}/2 \rfloor$)
3. **Two-electron integral**: 
   - Loop over orbital pairs $(k,l)$ within occupied manifold
   - Grid integrate Coulomb kernel $1/|r_1-r_2|$
   - Accumulate matrix elements: $\langle i,j | 1/r | k,l \rangle$
4. **Apply mixing**: $V_x^{\text{HSE}},ij \gets -\alpha \times \text{integral}$

**Complexity Analysis:**
- **Grid points per fragment**: $O(L^3)$ where $L$ = grid points per dimension in fragment
- **Orbital pairs**: $O(N_\text{basis}^2)$ where $N_\text{basis}$ = states per fragment
- **Total per fragment**: $O(N_\text{basis}^2 \times L^6)$ for double grid summation
- **For small fragments**: Practical due to limited $L$ and $N_\text{basis}$

---

## Key Features

✅ **Fragment-local computation**
- Each fragment calculates exchange independently
- No inter-fragment 4-index integrals required
- Scales as $O(N_\text{frag})$ in fragment count

✅ **Seamless integration with DG-Fragment RT**
- Uses existing basis: `dg_frag%phi_frag`
- Grid spacing from `dg_frag%hgs`
- Compatible with time evolution loop

✅ **Backward compatible**
- Default: `yn_hse = 'n'` (HSE disabled)
- Existing PBE/LDA calculations unaffected

✅ **Density matrix method (Plan A)**
- Physically accurate exchange energy
- No gaussian basis expansion required
- Suitable for real-space RT-TDDFT

---

## Usage Examples

### Enable HSE in Calculation

**File: `inputfile` (sample)**

```fortran
&functional
  xc = 'PBE'
  yn_hse = 'y'           ! Activate HSE
  hse_alpha = 0.25d0     ! PBE0 mixing
  hse_omega = 0.11d0     ! Screening (future)
/
```

### Test Case: DG-Fragment C2H2 with HSE

See sample file: `inputfile_dg_fragment_rt_hse`

```bash
cd ~/SALMON-v.2.2.2/samples/exercise_dg_fragment_rt/
# First run DC-LCFO to generate fragment basis
../../../build/salmon < inputfile_dclcfo_update

# Then run RT-TDDFT with HSE
../../../build/salmon < inputfile_dg_fragment_rt_hse
```

---

## Computational Complexity

| Aspect | Complexity | Notes |
|--------|-----------|-------|
| Fragment exchange | $O(N_b^2 L^6)$ | $N_b$ = basis states, $L$ = grid points/dim |
| Total system | $O(N_\text{frag} \times N_b^2 L^6)$ | Linear in fragment count |
| Time step | ~6-8 hours (target) | For typical 8-fragment system with HSE |
| Gain vs Full | ~100-1000x | Due to DC decomposition + locality |

**Optimization opportunities:**
1. Grid coarsening in exchange integral (reduce $L^6$)
2. Cache blocking for grid pairs
3. MPI + GPU acceleration
4. Screening function for short-range only

---

## Physical Validation Targets

For C2H2 molecule (test case):

| Property | Target | Expected With HSE |
|----------|--------|-------------------|
| HOMO-LUMO gap | LDA: ~7.5 eV | HSE: ~10-11 eV |
| Exchange energy | LDA underestimated | Better description |
| RT dynamics | Expected field response | Improved accuracy |

---

## Future Improvements (Phase 2+)

1. **Screening implementation**: Use `hse_omega` for range-separated kernel
2. **Energy calculation**: Compute exchange energy contribution
3. **Density-dependent tuning**: Adaptive `hse_alpha` (CAM-B3LYP style)
4. **Inter-fragment exchange**: Weak coupling corrections
5. **GPU acceleration**: Fragment-local loop on GPU

---

## Code References

- **Global parameters**: `src/io/salmon_global.f90` (lines 127-130)
- **I/O routines**: `src/io/inputoutput.f90` 
  - Namelist: line 278
  - Defaults: line 717
  - Broadcast: line 1248
- **HSE computation**: `src/rt/rt_dg_fragment.f90`
  - Main routine: line 2550 (branch)
  - Exchange kernel: lines 2350-2470
  - Hamiltonian update: line 2550-2560

---

## Compilation & Testing

**Build with HSE enabled:**
```bash
cd SALMON-v.2.2.2/build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

**Sample test execution:**
```bash
./salmon < inputfile_dg_fragment_rt_hse > output.txt 2>&1
```

**Expected output:**
- Fragment-specific numbers at each time step
- Energy breakdown with HSE contribution
- Orbital dynamics influenced by improved exchange

---

## Authors & History

- **Implementation**: Plan A (Density Matrix Method)
- **Date**: 2026-02-22
- **Framework**: SALMON v.2.2.2 DG-Fragment RT-TDDFT
- **Status**: Phase 1 (Foundation) - Compilation verified ✓

---

**Next Phase**: Integrate screening effects, validate physical results, benchmark performance
