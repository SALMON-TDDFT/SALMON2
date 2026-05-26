# SALMON2-SBE — SALMON fork with SBE dephasing

This repository is a fork of the original SALMON project, an open-source
software package for ab-initio quantum-mechanical calculations of light-matter
interactions.

For the upstream SALMON project, documentation, and publications, see:

http://salmon-tddft.jp/

## Fork features

This fork implements the **gauge-covariant Velocity Gauge formulation** of the Semiconductor Bloch Equations (SBE) based on arXiv:2012.00994v1. Key improvements include:

1. **Gauge-covariant decoherence** (Eq. 7): Replaces the standard `-ρ_ij/T₂` term with a double-commutator form `D = -1/(T₂·E_g²) [H_eff, [H_eff, ρ]]` that preserves gauge invariance and trace conservation.

2. **Exact current operator** (Eq. 15): Computes the total current as `J = Tr[(p - A·I) ρ]` without perturbative expansions, ensuring proper inter/intra-band compensation in the Velocity Gauge.

3. **Hermiticity stabilization**: Enforces `ρ = ρ†` at each time step to maintain real-valued currents and FFT stability.

### Configuration parameters

The `&sbe` namelist now accepts:

```fortran
&sbe
  ! ... existing SBE parameters ...
  
  ! Dephasing time (femtoseconds)
  t2_sbe_fs = 50.0d0
  
  ! Minimum bandgap for gauge-covariant decoherence (eV)
  ! If > 0, uses this value; if <= 0, auto-calculates from band structure
  eg_ev = -1.0d0  ! Default: auto-calculate
/
```

Parameter details:

- **`t2_sbe_fs`**
  - **Units:** femtoseconds (`fs`)
  - **Default:** `1.0d10` (disables dephasing)
  - **Active range:** `0 < t2_sbe_fs < 1.0d9` enables gauge-covariant decoherence
  - **Physical effect:** Off-diagonal density-matrix elements decay via the gauge-invariant double-commutator form (Eq. 7 of arXiv:2012.00994v1)

- **`eg_ev`**
  - **Units:** electron-volts (`eV`)
  - **Default:** `1.5d0` (typical semiconductor bandgap)
  - **Special value:** `-1.0d0` triggers automatic calculation from band structure
  - **Usage:** Set to a positive value to use a fixed bandgap; set to `-1.0d0` for auto-calculation
  - **Physical effect:** Sets the energy scale `E_g` in the decoherence prefactor `1/(T₂·E_g²)`
  - **Conversion:** Automatically converted to atomic units internally using `au_ev` (= 27.211386245988 eV/Hartree)

### Physical background

Standard decoherence terms break gauge invariance, leading to artifacts in high-harmonic generation (HHG) spectra and incorrect Light-Matter interaction physics. The implementation follows Sec. II.A and Sec. III.A of arXiv:2012.00994v1:

- **Eq. 4:** Velocity Gauge Hamiltonian `H_eff = H₀ + A·p`
- **Eq. 7:** Gauge-covariant decoherence superoperator
- **Eq. 15:** Exact current operator in Velocity Gauge

The Taylor expansion (4th order) combined with gauge-covariant decoherence preserves:
- Trace conservation: `Tr[ρ(t)] = const`
- Hermiticity: `ρ(t) = ρ†(t)`
- Gauge covariance: Results transform correctly between Length and Velocity gauges

### Input propagation

The new parameters are:

- Declared globally as `t2_sbe_fs` and `eg_ev`;
- Read from the `&sbe` namelist;
- Initialized to defaults (`t2_sbe_fs = 1.0d10`, `eg_ev = 1.5d0`);
- Broadcast to all MPI ranks before SBE propagation;
- `eg_ev` is converted from eV to atomic units internally (`eg_au = eg_ev / au_ev`).
- If `eg_ev = -1.0d0`, the minimum bandgap is automatically calculated from the band structure.

## Minimal SBE input example

```fortran
&calculation
  theory = 'sbe'
/

&sbe
  ! Other required SBE system parameters must be supplied as usual.
  
  ! Enable gauge-covariant dephasing with T2 = 50 fs
  t2_sbe_fs = 50.0d0
  
  ! Option 1: Use default bandgap (1.5 eV)
  ! (eg_ev is optional, defaults to 1.5 eV)
  
  ! Option 2: Auto-calculate minimum bandgap from band structure or Override with a fixed bandgap value (e.g., experimental value)
  eg_ev = -1.0d0
  
  ! Option 3: Active bands close to E_F (±15 eV)
  frozen_core_threshold_ev = -15.0d0,
  frozen_free_threshold_ev =  15.0d0,
/
```

Set `t2_sbe_fs` back to the default value, or any value greater than or equal to
`1.0d9`, to recover the original no-dephasing behavior.

Set `eg_ev = -1.0d0` to use the automatically calculated minimum bandgap from the
band structure. Set `eg_ev > 0` to override with a specific value in eV (default: 1.5 eV).

### ETDRK4 Solver for Semiconductor Bloch Equations

The `dt_evolve_bloch_etdrk4` subroutine implements a **4th-order Exponential Time Differencing Runge-Kutta (ETDRK4)** scheme for solving the Semiconductor Bloch Equations (SBE) in the velocity gauge, following the methodology of Kassam & Trefethen (SIAM J. Sci. Comput., 2005).

#### Operator Splitting

The method exploits the structure of the SBE by splitting the Liouville-von Neumann equation into linear and nonlinear parts:

```
dρ/dt = L·ρ + N(ρ, t)
```

**Linear operator L (integrated exactly):**
- Diagonal in the band basis: `L_{nm} = -i(ε_n - ε_m) - Γ_{nm}`
- Contains the unperturbed Hamiltonian evolution `-i[H₀, ρ]`
- Includes static decoherence: `Γ_{nm} = (ε_n - ε_m)² / (T₂·E_g²)` for active bands
- Solved analytically via matrix exponential: `ρ → exp(L·Δt)·ρ`

**Nonlinear operator N (integrated via RK4):**
- Contains the light-matter interaction: `N = -i[V(t), ρ]` where `V = A(t)·p`
- Includes dynamic decoherence: `D_dynamic = -1/(T₂·E_g²)·[H_eff, [H_eff, ρ]] - L_decoh`
- Evaluated at 4 intermediate stages using Runge-Kutta quadrature

#### Why ETDRK4?

1. **Exact treatment of fast oscillations**: The linear part `exp(-i·ΔE·Δt)` is computed exactly, eliminating phase errors for large energy gaps (deep valence bands, high conduction bands).

2. **Larger time steps**: Stable for `Δt ~ 1-5 a.u.` (compared to `Δt < 0.5 a.u.` for explicit Taylor-4), providing **3-5× speedup** at the same accuracy.

3. **Contour integral stabilization**: Coefficients `Q, f₁, f₂, f₃` are computed via complex contour integration to avoid catastrophic cancellation when `|L·Δt| ≈ 0` (critical for near-degenerate bands).

#### Frozen Core Approximation

For systems with many deep bands (e.g., 80 bands with 60 below -20 eV), the **selective nonlinear term** optimization can be enabled:

```fortran
! In input file:
frozen_core_threshold_ev = -15.0  ! Freeze bands below E_F - 15 eV
frozen_free_threshold_ev = +20.0  ! Freeze bands above E_F + 20 eV
```

**How it works:**
- Frozen bands evolve only under the linear operator (exact phase oscillation).
- Nonlinear terms `N(ρ)` are computed only for the active subspace (e.g., 20×20 instead of 80×80 matrices).
- After each step, frozen bands are reset to ground state: `ρ_{nn} = 1` (occupied) or `0` (empty).
  
## License

SALMON is available under Apache License version 2.0.

    Copyright 2017-2025 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
