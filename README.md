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
  
  ! Option 2: Auto-calculate minimum bandgap from band structure
  eg_ev = -1.0d0
  
  ! Option 3: Override with a fixed bandgap value (e.g., experimental value)
  ! eg_ev = 1.519d0  ! GaAs bandgap at room temperature
  ! Активные зоны около E_F (±15 eV)
    frozen_core_threshold_ev = -15.0d0,
    frozen_free_threshold_ev =  15.0d0,
/
```

Set `t2_sbe_fs` back to the default value, or any value greater than or equal to
`1.0d9`, to recover the original no-dephasing behavior.

Set `eg_ev = -1.0d0` to use the automatically calculated minimum bandgap from the
band structure. Set `eg_ev > 0` to override with a specific value in eV (default: 1.5 eV).

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
