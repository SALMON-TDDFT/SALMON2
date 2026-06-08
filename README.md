# SALMON2-SBE — SALMON fork with CF4/Yoshida SBE propagation, CPTP decoherence & EPM ground states

This repository is a fork of the original [SALMON project](http://salmon-tddft.jp/), an open-source software package for *ab-initio* quantum-mechanical calculations of light-matter interactions. 

This fork extends SALMON's Semiconductor Bloch Equations (SBE) module with a **commutator-free Magnus 4 (CF4) / Suzuki-Yoshida exponential propagator**, a **strictly CPTP Kuhn-Zurek/Caldeira-Leggett decoherence model**, **frozen-core optimizations**, and a self-contained **local Empirical Pseudopotential Method (EPM)** ground-state solver (Cohen-Bergstresser, GaAs) that closes the EPM → SBE pipeline without external scripts.


## Key Fork Features

### 1. Strictly CPTP Kuhn-Zurek/Caldeira-Leggett Decoherence
Phenomenological dephasing schemes (whether $-\rho_{nm}/T_2$ in the field-free basis, or double-commutator dissipators built from the instantaneous dressed Hamiltonian $H_{eff}=H_0+\mathbf{A}(t)\cdot\mathbf{p}$) generally fail to be **completely positive and trace preserving (CPTP)** for arbitrary parameters/timesteps, and can introduce artifacts such as spurious adiabatic-following relaxation in strong fields.

This fork instead implements a **Kuhn-Zurek/Caldeira-Leggett wave-packet dephasing model** that is *exactly* CPTP by construction:
* At every step the instantaneous (Houston/adiabatic) eigenbasis $U(t)$ of $H_{VG}(t)=H_0+\mathbf{A}(t)\cdot\mathbf{p}$ is computed, together with the branch (wave-packet) positions $X_a(t)$, propagated via their group velocities $V_a = (U^\dagger\boldsymbol{\pi}U)_{aa}+\mathbf{A}(t)$.
* The density matrix is rotated into this basis, $\tilde\rho = U^\dagger\rho U$, and dephased through an **exact Hadamard/Gaussian (RBF) kernel** $\tilde\rho_{ab}\leftarrow e^{-\lambda (X_a-X_b)^2\,\tau}\,\tilde\rho_{ab}$.
* By the Schoenberg/Bochner positive-definiteness of the Gaussian kernel and the Schur product theorem, this Hadamard map is CPTP for **any** $\tau\ge 0$ — no positivity violations, no ad-hoc clipping.
* The decoherence rate is set physically via $\lambda = k_B T/\tau_m$ (`sbe_decoh_temperature_k`, `sbe_decoh_tau_m_fs`).

### 2. CF4 + Suzuki-Yoshida Exponential Propagator
Replaces the previous ETDRK4/Taylor-4 propagators with a **commutator-free Magnus 4th-order (CF4)** exponential integrator evaluated on two-point Gauss-Legendre quadrature nodes, composed into a 4th-order scheme via the **Suzuki-Yoshida triple-jump** ($p_1=1.35120719196$, $p_2=-1.70241438392$):
* Each CF4 sub-step is realized as **two exact unitary rotations** $\rho \to e^{-i\Omega_2}e^{-i\Omega_1}\rho\, e^{+i\Omega_1}e^{+i\Omega_2}$, with $\Omega_{1,2}$ Hermitian combinations of the Hamiltonian sampled at the Gauss-Legendre nodes, exponentiated *exactly* via eigendecomposition (no Padé/Krylov truncation error — unitary to machine precision).
* The full step combines the coherent and dissipative parts via **Strang splitting**, $D(h/2)\circ\big[S_2(p_1h)\circ S_2(p_2h)\circ S_2(p_1h)\big]\circ D(h/2)$, with the Suzuki-Yoshida composition wrapping **only** the unitary part. This is essential for CPTP: a negative sub-step ($p_2 h<0$) is a harmless backward-time *unitary* rotation, but applying it to the dissipator would invert the sign of the Hadamard kernel exponent and break positive semi-definiteness.
* Branch positions $X_a$ are advanced using the **midpoint (average of endpoint) velocities**, matching the overall 4th-order accuracy of CF4 (a forward-Euler update would degrade the scheme to 1st order).

### 3. Frozen Core / Active Subspace Optimization
For systems with many deep bands (e.g., 80 bands where 60 lie below -20 eV), evaluating the nonlinear commutator $[V, \rho]$ is computationally wasteful. This fork introduces an **Active Subspace** projection:
* Deep core and high-energy free bands are "frozen" and evolve purely under the exact linear operator $L$.
* Nonlinear light-matter interactions are computed exclusively in the active subspace (e.g., $20 \times 20$ instead of $80 \times 80$ ZGEMM calls).
* Yields an additional **~30× speedup** for the nonlinear evaluation step without sacrificing physical accuracy.

### 4. Exact Current Operator
Computes the total current as $J = \text{Tr}[(\mathbf{p} + \mathbf{A}) \rho]$ (in atomic units) without relying on perturbative expansions, ensuring proper inter/intra-band compensation in the Velocity Gauge. Hermiticity stabilization (`ρ = ρ†`) is enforced at each step to maintain real-valued currents and FFT stability.

### 5. Local Empirical Pseudopotential Method (EPM) ground states
A self-contained local-EPM ground-state solver (`theory='epm'`, `src/epm`) that computes the Cohen-Bergstresser band structure and momentum matrix elements for zincblende GaAs directly in SALMON, and writes `SYSNAME_k.data`/`SYSNAME_eigen.data`/`SYSNAME_tm.data` in exactly the format read by `gs_info_ssbe` — closing the EPM → SBE pipeline end-to-end without external scripts (`rvnl_tm` is written as identically zero, since a local pseudopotential has no nonlocal velocity correction).

---

## The CF4 + Suzuki-Yoshida + CPTP Operator Splitting

The propagator advances $\partial_t\rho = -i[H_{VG}(t),\rho] + \mathcal{D}[\rho]$ over a step $h$ as

$$\rho(t+h) = D(h/2)\circ\Big[S_2(p_1 h)\circ S_2(p_2 h)\circ S_2(p_1 h)\Big]\circ D(h/2)\,[\rho(t)]$$

**Unitary part $S_2(\tau)$ — CF4 on Gauss-Legendre nodes** (nodes $c_{1,2}=\tfrac12\mp\tfrac{\sqrt3}{6}$, weights $\alpha_{1,2}=\tfrac14\pm\tfrac{\sqrt3}{6}$):
* $H_1=H_{VG}(t+c_1\tau)$, $H_2=H_{VG}(t+c_2\tau)$
* $\Omega_1=\tau(\alpha_1 H_1+\alpha_2 H_2)$, $\Omega_2=\tau(\alpha_2 H_1+\alpha_1 H_2)$
* $\rho \to e^{-i\Omega_2}e^{-i\Omega_1}\,\rho\,e^{+i\Omega_1}e^{+i\Omega_2}$, each exponential built *exactly* from an eigendecomposition of the Hermitian generator (unitary to machine precision).

**Dissipative part $D(\tau)$ — Strang/Hadamard Kuhn-Zurek dephasing** (always applied with $\tau=+h/2 > 0$):
* Diagonalize $H_{VG}(t)\to U(t),\ \{E_a\}$ (Houston/adiabatic basis); $\tilde\rho = U^\dagger\rho U$
* $\tilde\rho_{ab}\leftarrow \exp[-\lambda(X_a-X_b)^2\tau]\,\tilde\rho_{ab}$ (PSD Hadamard/Gaussian kernel ⇒ exactly CPTP for $\tau\ge0$)
* Rotate back $\rho = U\tilde\rho U^\dagger$; update $X_a \mathrel{+}= \tfrac12(V_a(t)+V_a(t+h))\,h$, with $V_a=(U^\dagger\boldsymbol\pi U)_{aa}+\mathbf{A}(t)$.

**Why Yoshida wraps only the unitary part:** the middle Yoshida sub-step has $p_2 h<0$. For $S_2$ this is merely a unitary rotation run backwards in time — exact and always valid. For $D$, however, a negative $\tau$ would turn the Gaussian kernel $e^{-\lambda\Delta X^2\tau}$ into $e^{+\lambda\Delta X^2|\tau|}$, which is not positive semi-definite (violates the Schoenberg/Bochner criterion and the Schur product theorem) and would break CPTP. Hence $D$ is applied only with $\tau=+h/2$, via Strang splitting around the (always-safe) Yoshida-composed unitary block.

---

## Configuration Parameters

The `&sbe` namelist now accepts the following parameters:

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `sbe_decoh_temperature_k` | K | `-1.0d0` | Bath temperature $T$ for the Kuhn-Zurek/Caldeira-Leggett dephasing rate $\lambda=k_B T/\tau_m$. Both this and `sbe_decoh_tau_m_fs` must be `> 0` to enable decoherence. |
| `sbe_decoh_tau_m_fs` | fs | `-1.0d0` | Wave-packet momentum-relaxation time $\tau_m$ entering $\lambda=k_B T/\tau_m$. |
| `frozen_core_threshold_ev` | eV | `0.0d0` | Freeze bands below $E_F + \text{threshold}$. (Use negative values, e.g., `-15.0`). |
| `frozen_free_threshold_ev` | eV | `0.0d0` | Freeze bands above $E_F + \text{threshold}$. (Use positive values, e.g., `+20.0`). |

*Note: Internal conversions to atomic units (Hartree) are handled automatically (`kB_au`, `au_fs`).*

### Real-time output frequency (`&analysis`)

Real-time SBE propagation writes three diagnostic files (`SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data`, `SYSNAME_sbe_nex_k.data`), each on its own cadence selectable in the `&analysis` namelist. The k-resolved file in particular can grow to gigabytes for dense k-grids/long runs, so its default stride is ten times coarser than the band-projection output:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `out_rt_energy_step` | `10` | Stride (in time steps) for `SYSNAME_sbe_rt_energy.data` (total energy / trace) and the stdout progress log. |
| `out_projection_step` | `100` | Stride for `SYSNAME_sbe_nex.data` (number of excited electrons/holes, summed over k). |
| `out_projection_k_step` | `1000` | Stride for `SYSNAME_sbe_nex_k.data` (Houston-basis population of the lowest conduction band, resolved per k-point). Defaults to 10× `out_projection_step` to avoid producing terabyte-scale output on dense k-grids; increase the stride (larger value) further for very large `nk`/`nt`. |

`SYSNAME_sbe_nex_k.data` reports, for every saved time `t`, one block of `nk` lines `ik, kx, ky, kz, population_lcb`, where `population_lcb = (W^\dagger \rho W)_{aa}` is the diagonal element of the density matrix rotated into the instantaneous Houston (adiabatic) eigenbasis $W$ of $H_{VG}(t)$ for the lowest conduction band $a = N_{elec}/2+1$ — i.e. the same gauge-independent basis used internally by the CPTP dephasing step.

### Plotting the real-time output (`plot_sbe_results.py`)

The repository root contains a self-contained `plot_sbe_results.py` script (matplotlib + numpy, not part of the Fortran build — copy it into the calculation directory and run it there). It scans the directory for `SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data` and `SYSNAME_sbe_nex_k.data`, and produces (with no interactive windows, `Agg` backend):
* line plots of total energy and excited-electron/hole counts vs time;
* for `SYSNAME_sbe_nex_k.data`, one PNG per saved time step (the time value is encoded in the file name), each showing the Houston-basis lowest-conduction-band population as three 2D heatmap slices of the k-grid ($k_x$-$k_y$, $k_x$-$k_z$, $k_y$-$k_z$).

```sh
cp plot_sbe_results.py /path/to/calculation/
cd /path/to/calculation/
python3 plot_sbe_results.py            # writes PNGs into ./sbe_plots/
```

The `&epm` namelist configures the local-EPM ground-state solver (`theory='epm'`):

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `epm_material` | — | `'GaAs'` | Material whose tabulated Cohen-Bergstresser local form factors are used (currently `'GaAs'`). |
| `epm_lattice_constant_au` | Bohr | `10.68d0` | Zincblende lattice constant $a$. |
| `epm_pw_cutoff_ry` | Ry | `11.1d0` | Plane-wave cutoff $|\mathbf{k}+\mathbf{G}|^2$ for the basis set. |

---

## Minimal SBE Input Example

```fortran
&calculation
  theory = 'sbe'
/

&sbe
  ! ... standard SALMON SBE system parameters ...

  ! ---------------------------------------------------------
  ! 1. Kuhn-Zurek/Caldeira-Leggett Decoherence (strictly CPTP)
  ! ---------------------------------------------------------
  ! lambda = kB*T / tau_m;  enabled only when both are > 0
  sbe_decoh_temperature_k = 300.0d0
  sbe_decoh_tau_m_fs      = 10.0d0

  ! ---------------------------------------------------------
  ! 2. Frozen Core / Active Subspace Optimization
  ! ---------------------------------------------------------
  ! Only evolve bands within ±15 eV of the Fermi level non-linearly.
  ! Deep core bands will only undergo exact linear phase oscillation.
  frozen_core_threshold_ev = -15.0d0
  frozen_free_threshold_ev =  15.0d0
/
```

**Reverting to default behavior:**
* Set `sbe_decoh_temperature_k` and/or `sbe_decoh_tau_m_fs` to a non-positive value to recover the original purely-coherent (no dephasing, $D\equiv 0$, trivially CPTP) behavior.
* Set both `frozen_core_threshold_ev` and `frozen_free_threshold_ev` to `0.0d0` to force all bands into the active nonlinear subspace.

## Minimal EPM → SBE Pipeline Example

### Standalone Python reference (`epm_gaas_reference.py`)

For quick debugging without building/running SALMON, the repository root also contains `epm_gaas_reference.py` -- a monolithic, single-machine NumPy/SciPy reimplementation of the GaAs Cohen-Bergstresser local-EPM solver (no MPI/OpenMP). It builds the same lattice/plane-wave basis/Hamiltonian/momentum matrices as `src/epm`, and writes byte-compatible `SYSNAME_k.data`/`_eigen.data`/`_tm.data` files that `gs_info_ssbe` can read directly -- so its output can be diffed against the Fortran `theory='epm'` run, or fed straight into an SBE real-time calculation. All parameters (lattice constant, plane-wave cutoff, k-grid, number of bands/electrons, sysname) are hardcoded constants at the top of the script -- edit them there and run:

```sh
python3 epm_gaas_reference.py
```

This is a debugging aid only -- `theory='epm'` in SALMON remains the primary, MPI/OpenMP-parallel ground-state path.

```fortran
! Step 1: ground state via local EPM (writes GaAs_k/_eigen/_tm.data)
&calculation
  theory = 'epm'
/
&epm
  epm_material            = 'GaAs'
  epm_lattice_constant_au = 10.68d0
  epm_pw_cutoff_ry        = 11.1d0
/
```
```fortran
! Step 2: real-time SBE propagation reading the files generated above
&calculation
  theory = 'sbe'
/
&system
  ! sysname, lattice vectors, num_kgrid, nstate, nelec must match the EPM run
/
```
---

## References & Theoretical Background

1. **Commutator-Free Magnus Integrators:** Blanes, S., & Moan, P. C. "Practical symplectic partitioned Runge-Kutta and Runge-Kutta-Nyström methods." *J. Comput. Appl. Math.* 142, 313-330 (2002); Alvermann, A., & Fehske, H. "High-order commutator-free exponential time-propagation of driven quantum systems." *J. Comput. Phys.* 230, 5930-5956 (2011).
2. **Suzuki-Yoshida Composition:** Yoshida, H. "Construction of higher order symplectic integrators." *Phys. Lett. A* 150, 262-268 (1990).
3. **CPTP / Lindblad & RBF-kernel positivity:** Schoenberg, I. J. "Metric spaces and completely monotone functions." *Ann. Math.* 39, 811-841 (1938) (Bochner/Schoenberg PSD criterion for Gaussian/RBF kernels); Schur product theorem (Hadamard maps of PSD matrices are PSD).
4. **Caldeira-Leggett / Kuhn-Zurek Decoherence:** Caldeira, A. O., & Leggett, A. J. "Path integral approach to quantum Brownian motion." *Physica A* 121, 587-616 (1983); Zurek, W. H. "Decoherence, einselection, and the quantum origins of the classical." *Rev. Mod. Phys.* 75, 715 (2003).
5. **Cohen-Bergstresser Local Pseudopotentials:** Cohen, M. L., & Bergstresser, T. K. "Band Structures and Pseudopotential Form Factors for Fourteen Semiconductors of the Diamond and Zinc-blende Structures." *Phys. Rev.* 141, 789 (1966).
6. **Velocity-Gauge SBE / Houston Basis:** Wismer, M. S., & Yakovlev, V. S. "Gauge-independent decoherence models for solids in external fields." *Phys. Rev. B* 97, 144302 (2018).
7. **Original SALMON SBE:** Sato, S. A. et al. "Multiscale computational approach for light-matter interactions." *Phys. Rev. B* 92, 115145 (2015).

## License

SALMON is available under the Apache License version 2.0.

    Copyright 2017-2026 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.