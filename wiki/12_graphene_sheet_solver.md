# Graphene sheet solver: gapless-cone dissipation with collisional memory, a two-temperature Coulomb sector, and a self-consistent sheet field

**Status: implemented and unit-tested (2026-09-04); calc-validation in exercise x14. Maintained together with `wiki/10` §8.11.**
*Written in the form of a methods paper so that it can be lifted into a manuscript; every equation is the one the code evaluates, every constant is cited, every claim carries its test.*

---

## Abstract

We describe the graphene branch of the SALMON2-TROUT semiconductor-Bloch-equation (SBE) solver, built to study field-induced transparency / absorption of a monolayer graphene sheet under 1–100 kV/cm single-cycle THz transients (the maintainer's DAST source) and near-infrared pulses. Three elements were added to the general velocity-gauge, completely-positive (GKLS) SBE machinery: (i) a **collisional-memory (non-Markovian) treatment of the ring dissipators adapted to the gapless Dirac cone** — the electron–phonon sectors keep their phonon-line kernels while the Coulomb (Auger / carrier-multiplication) sector, which on the cone is a *global* rate model, receives the 2D Dirac-plasmon line of the instantaneous electron–hole plasma as its memory kernel; (ii) a **two-temperature description of the Coulomb sector**, in which a carrier temperature $T_e$ and the quasi-Fermi levels are read from the first two moments of the gathered Dirac-cone populations while the lattice stays at the phonon-bath temperature and cools the carriers through the phonon channel; (iii) a **self-consistent sheet field** (radiation reaction) in the single-cell driver, so that the total field written by the solver *is* the transmitted field and the transmission coefficient follows from the field before and after the sheet; and (iv) a **doped / finite-temperature initial occupation**, which turns the intrinsic semimetal into the metal a real sample is and gives it a Drude sector, with the velocity-gauge pure-gauge reference kept undoped so that the physical intraband current survives the correction of point (v). A further element turned out to be indispensable at THz (v): a **parameter-free pure-gauge restoration of the velocity-gauge current** — a truncated basis cancels only part of the diamagnetic current of the filled π band (70 % for two bands), the remainder, proportional to $A=E/\omega$, turns the sheet into a plasma mirror, and the fitted linear correction first tried over-corrects at high field and makes the self-consistent sheet unstable; subtracting the adiabatic ground-state current of the same truncated Hamiltonian removes the artifact exactly at every field with no adjustable quantity (§6a). We give the equations, the numerical realization, the cost scaling ($O(N_k^2)$ for the graphene ring), the k-mesh rules (resonance-shell resolution; Dirac point on the half-shifted Monkhorst–Pack mesh only for odd multiples of 3; Zener excursion $A_0$ versus mesh spacing), the unit tests that pin each piece, and the calculation-level validations. A level check performed on the way exposed and removed a spurious 0.21 eV gap at the Dirac point of the previously used 7-plane-wave empirical pseudopotential basis.

---

## 1. Scope and notation

Monolayer graphene is represented by the π/π* pair of the Ramanujam local empirical pseudopotential (EPM) in the 2-atom hexagonal primitive cell embedded in a 20 Å vacuum slab (strictly two-dimensional plane-wave basis, no $G_z$ components and no dispersion along $z$: an isolated sheet, not graphite); the SBE runs on `nstate = 2` bands (the pure-gauge restoration of §6a makes the result independent of `nstate`). Hartree atomic units are used throughout ($e=\hbar=m_e=1$, $c = 137.036$); the sheet lies in the $xy$ plane, the driving field is in-plane. $N_k$ denotes the number of k-points of the $N\times N\times1$ Monkhorst–Pack (MP) mesh; $A_{2D}=(\sqrt3/2)a^2$ is the primitive-cell area, $L_z$ the slab height along the vacuum axis, so the cell volume is $V=A_{2D}L_z$.

The ground-state basis is the shell-complete set $|\mathbf G|^2\le 29.4$ a.u. (43 plane waves). The earlier 7-plane-wave set ($|\mathbf G|^2\le 2.94$) is **not closed under the little group $C_{3v}$ of K** (the rotation about K maps a first-shell vector onto the second-shell vector $\mathbf b_2-\mathbf b_1$ absent from the set); the symmetry protection of the Dirac degeneracy is then lost and a spurious gap of 0.2125 eV opens at K (identical in the Python reference and in the Fortran band path). With 43 plane waves the gap is $<10^{-5}$ eV, $v_F = 0.960\times10^6$ m/s, and the Γ-bottom / M-dip values fall inside the thesis acceptance windows (`tests/test_graphene_dirac_levels.py`).

## 2. Master equation and channels

The density matrix $\rho_{\mathbf k}(t)$ obeys the velocity-gauge SBE in GKLS form,
$$
\dot\rho_{\mathbf k} = -i\,[H_{\mathbf k}(t),\rho_{\mathbf k}] + \sum_c \mathcal D_c[\rho],\qquad
H_{\mathbf k}(t) = \varepsilon_{\mathbf k} + \mathbf A(t)\cdot\boldsymbol\pi_{\mathbf k} + \tfrac12 A^2,
$$
propagated with the fourth-order commutator-free Magnus (CF4) scheme in Strang splitting with the dissipators $\mathcal D_c$, which act in the instantaneous Houston (field-dressed) basis of $H_{\mathbf k}(t)$ (`wiki/03`, `wiki/08`). On graphene the "ring" (inter-k) channels are:

* **electron–phonon**, the two Kohn-anomaly optical modes $E_{2g}$ (Γ, 196 meV) and $A_1'$ (K, 160 meV) with $\langle g^2_\Gamma\rangle=0.0405$, $\langle g^2_K\rangle=0.0994$ eV² and the GW factor 2 on the K mode [1,2], plus the quasi-elastic acoustic deformation-potential mode ($D=16$ eV, $v_{ph}=2\times10^6$ cm/s [3], grid-resolved $q$, Thomas–Fermi screened);
* **Coulomb (Auger recombination / carrier multiplication)**, the Rana rate model [4]: CCCV/CVVV recombination $R$ and their CVCC/VCCC time-reverses $G$ evaluated on the instantaneous quasi-Fermi levels of the gathered sheet densities $n$, $p$, applied as a uniform-fractional, trace-exact, bounded population transfer $R-G$. Impact ionization with a gap threshold, the carrier–carrier Fermi–Dirac fit and Kuhn–Zurek dephasing are not defined on a gapless cone and are refused by the code.

The Coulomb balance $R=G$ holds iff the electron and hole quasi-Fermi levels coincide; for a symmetric pair population this is $\mu=0$, i.e. the intrinsic density of the cone at the temperature $T$ the rates are evaluated at,
$$
n_0(T)=n_i(T)=\frac{\pi}{6}\left(\frac{k_BT}{\hbar v_F}\right)^2 \;=\; 8.08\times10^{10}\ \mathrm{cm^{-2}}\ \text{at 300 K}. \tag{1}
$$
Below $n_i$ the plasma net-multiplies, above it net-recombines: the pair population saturates at $n_i(T)$ (`tests/test_rana_saturation.f90`: root of $R-G$ at $n_i$ to $10^{-4}$, two-sided monotone CPTP relaxation, $T^2$ law).

## 3. Collisional memory on the gapless cone (the "2D colmem analog")

### 3.1 Motivation
The Markovian ring dissipators scatter the reversible field-induced admixture of conduction character ("dressing") as if it were real population (`wiki/10` §1). For gapped materials this was cured in three sectors (§8.7–8.10 of `wiki/10`): a memory kernel in the coherence damping, a memory filter on the populations that feed the collision kernels, and a dressed reference for the carrier measure. Graphene had been excluded by a guard. Two features of the cone required a dedicated version: its population channel is the **Coulomb rate model on the global densities $n,p$**, so the virtual share inflates the very quantities the R07 rates are evaluated on; and the memory line of a Coulomb collision is not a phonon energy but the **plasma response** — screening builds up on the inverse plasma frequency [5].

### 3.2 Kernels
All memory filters have the form used in `wiki/10` §8.6–8.9: a set of Lorentzian lines $\{c_j,\mu_j\}$, $\mu_j = 1/\tau_c \pm i\omega_j$ with the thermal split $(N_j+1):N_j$, the common width $1/\tau_c=\sigma_E$ (the ring's energy-matching width), and the discrete Markov anchor $R(0)=1$, so that a constant quantity is a machine-exact fixed point (calibrated rates untouched) while a modulation at frequency $\omega$ is transmitted with $|R(\omega)|=|\sum_j c_j/(\mu_j+i\omega)|$.

| sector | line set |
|---|---|
| e-ph coherence damping (`yn_sbe_colmem`) | graphene phonon table: $E_{2g}$, $A_1'$, acoustic |
| e-ph population source (`yn_sbe_colmem_pop`) | same |
| **Coulomb source densities $n,p$** | **2D Dirac plasmon** $\omega_{pl}(n,p)$ at $q=Q_{TF}$ |

The long-wavelength plasmon of the Dirac gas [6], $\omega_{pl}^2(q)=2e^2E_Fq/(\kappa\hbar^2)$ (Gaussian units), is generalized to the two-component electron–hole plasma by adding the Drude weights, each branch's $E_F$ becoming the finite-temperature intraband Drude weight [7]
$$
W(\mu)=2k_BT\,\ln\!\big[2\cosh(\mu/2k_BT)\big]\;\to\;|\mu|\ (\text{degenerate}),\quad 2k_BT\ln2\ (\text{intrinsic}),
$$
and it is evaluated at the collision's own screening momentum, the Thomas–Fermi vector of R07 Eq. (13),
$$
\omega_{pl}^2=\frac{2\,[W(\mu_c)+W(\mu_h)]\,Q_{TF}}{\varepsilon_r},\qquad
Q_{TF}=\frac{4k_BT}{\varepsilon_r v_F^2}\ln\!\big[(e^{\mu_c/k_BT}+1)(e^{\mu_h/k_BT}+1)\big]. \tag{2}
$$
At 300 K and $\varepsilon_r=10$: $\hbar\omega_{pl}=31$ meV (intrinsic), 35 meV at $n=p=10^{11}$, 133 meV at $10^{12}$ cm⁻² — the phonon scale, all $\ll 2\hbar\omega_{\rm laser}$ for near-IR drives (`tests/test_colmem_2d.f90`: degenerate limit $\omega_{pl}^2=2E_FQ_{TF}/\varepsilon_r$ to $10^{-6}$; a $2\omega$ breathing at 0.8 eV transmitted at $0.17=|R(2\omega)|$). No new input variables and no new free parameters enter: $T$, $\varepsilon_r$ and the $\mu$'s are the Rana channel's own.

### 3.3 Composition
Per ring step: dressed-reference subtraction (basis level) → ring gate → the raw Coulomb source moments $(n,p,\varepsilon)$ are captured → the phonon-line filter replaces the gathered populations for the e-ph kernel → the Coulomb channel filters $n,p$ with the plasmon line and evaluates $R-G$ on them (`rana_auger_dpop(n2d_in,p2d_in)`), while its transfer stencils and the CPTP limiter still use the instantaneous populations. Trace is exact and positivity unchanged. Kuhn–Zurek dephasing remains forbidden on graphene (many-body coherence loss).

## 4. Two-temperature Coulomb sector

### 4.1 Model
The R07 rates are quasi-equilibrium expressions in $(\mu_c,\mu_h,T)$. Evaluating them at the lattice temperature under-estimates the balance density of a hot plasma by $(T_e/T_L)^2$ (Eq. 1). We therefore read a **common carrier temperature** from the gathered populations. With the Dirac-point energy $E_D$ (midpoint of the instantaneous band edges, exact by electron–hole symmetry), the moments per unit area are
$$
n=\frac1{N_kA_{2D}}\sum_{\mathbf k} f_{c\mathbf k},\quad
p=\frac1{N_kA_{2D}}\sum_{\mathbf k}(2-f_{v\mathbf k}),\quad
\varepsilon=\frac1{N_kA_{2D}}\Big[\sum_{\mathbf k} f_{c\mathbf k}(\epsilon_{c\mathbf k}-E_D)+\sum_{\mathbf k}(2-f_{v\mathbf k})(E_D-\epsilon_{v\mathbf k})\Big],
$$
and $(k_BT_e,\mu_c,\mu_h)$ solve
$$
n=n_D(\mu_c,T_e),\quad p=n_D(\mu_h,T_e),\quad \varepsilon=\varepsilon_D(\mu_c,T_e)+\varepsilon_D(\mu_h,T_e), \tag{3}
$$
$$
n_D(\mu,T)=\frac{g}{2\pi}\!\int_0^\infty\!k\,f\!\Big(\frac{v_Fk-\mu}{k_BT}\Big)dk,\qquad
\varepsilon_D(\mu,T)=\frac{g}{2\pi}\!\int_0^\infty\!k\,(v_Fk)\,f\!\Big(\frac{v_Fk-\mu}{k_BT}\Big)dk,\quad g=4,
$$
with the closed forms $n_D(0,T)=\frac{\pi}{6}(k_BT/v_F)^2$ and $\varepsilon_D(0,T)=\frac{2}{\pi}\cdot\frac{3\zeta(3)}{2}\,(k_BT)^3/v_F^2$. At fixed $(n,p)$ the energy is monotone in $T$, so Eq. (3) is solved by a bisection in $\ln T$ with the density inversions nested inside (`dirac_fit_te`). $T_e$ is clamped from below at the lattice temperature (carriers are not colder than the bath) and the fit falls back to the bath when there are no carriers.

### 4.2 Use
$k_BT_e$ replaces the bath temperature in every carrier-side quantity of the Coulomb sector: the R07 integrals, $Q_{TF}$, the plasmon line and its thermal split. The phonon Bose factors keep the lattice temperature. **No separate rate equation for $T_e$ is integrated**: heating by the field and cooling by optical/acoustic phonon emission are already contained in the SBE kinetics; the two-temperature model is realized by *reading* $T_e$ from the distribution at every ring step. The time series $(T_e,\mu_c,\mu_h,n,p,n_i(T_e))$ is written to `*_sbe_te.data`. Interpretation: for a hot, non-thermal photo-excited shell the fit returns an effective temperature at which the R07 channel *generates* pairs toward $n_i(T_e)$ — the carrier-multiplication regime of graphene [8,9] — and switches to recombination as phonon cooling lowers $T_e$.

Unit test `tests/test_dirac_te_fit.f90`: the closed form of $\varepsilon_D(0,T)$; explicit 2D-mesh moments of thermal populations reproduce $n_D,\varepsilon_D$ to $0.5\%$ and the fit recovers $T$ to $10^{-4}$ and $\mu$ to $10^{-4}k_BT$; the degenerate case $E_F=0.3$ eV at 300 K; fallbacks; monotonicity.

## 4a. The doped sheet: initial occupation, the Drude sector, and the field scales of transparency

### 4a.0 Choosing the doping: it is read off the measurement, not fitted
$E_F$ is a property of the sample, and a THz transmission measurement determines
it together with the momentum-relaxation time $\tau$. The chain is closed and has
no free parameter left over.

**Step 1 -- from the measured transmission to a sheet conductance.** A measurement
that has not had the substrate divided out contains the substrate's own Fresnel
loss. For a slab of index $n_s$ with two incoherent faces,
$$
T_{\rm bare}=\Big[\frac{4n_s}{(1+n_s)^2}\Big]^{2},\qquad
\frac{T_{\rm meas}}{T_{\rm bare}}=\Big[\frac{1+n_s}{1+n_s+Z_0\sigma}\Big]^{2}
\ \Longrightarrow\
Z_0\sigma=(1+n_s)\Big[\sqrt{T_{\rm bare}/T_{\rm meas}}-1\Big]. \tag{4a.1}
$$
PET, $n_s=1.65$: $T_{\rm bare}=0.883$, so $T_{\rm meas}=0.60\to\sigma=24.7\,\sigma_{\rm univ}$
and $0.70\to14.3\,\sigma_{\rm univ}$ (`drude_check.py --t-meas`, pinned in
`tests/test_doped_drude.py`). Use the *low-field* value: that is the linear
conductance the doping alone produces.

**Step 2 -- split $\sigma_{dc}=D\tau/\pi$ into $E_F$ and $\tau$.** One equation, two
unknowns; any one of three independent handles closes it.

* *An independent measurement of the doping* (gate voltage, Hall, the Raman
  2D/G ratio) gives $E_F$, and then $\tau=\pi\sigma_{dc}/E_F$.
* *A resolved THz spectrum*: the Drude roll-off frequency is $1/\tau$, and $D$
  follows from the low-frequency plateau.
* *The field at which the transmission starts to rise* -- available from the very
  scan being modelled. Saturation begins where the vector-potential excursion
  reaches the Fermi radius, $A_0(E_{\rm sat})=k_F$; with the scaled DAST transient
  $A_0=6.213\times10^{-4}\,\text{a.u.}\times E_0[\mathrm{kV/cm}]$, so
  $$
  k_F = 6.213\times10^{-4}\,E_{\rm sat}[\mathrm{kV/cm}],\qquad E_F=\hbar v_Fk_F. \tag{4a.2}
  $$
  An onset near 30 kV/cm therefore means $E_F\simeq0.22$ eV, and $\tau$ follows from
  step 1. This is the self-contained route when nothing but the transmission scan
  is available.

**Is the answer a plausible sample?** Convert the same conductance to the units a
transport measurement quotes: $\sigma=24.7\,\sigma_{\rm univ}=1.50$ mS/sq, i.e.
$R_s=666\ \Omega/\square$ (the high-field value $14.3\,\sigma_{\rm univ}$ is
$1153\ \Omega/\square$). That is an ordinary as-transferred CVD monolayer. Splitting
it across dopings, each with the $\tau$ that reproduces the same $\sigma$:

| $E_F$ [eV] | $n_{2D}$ [cm⁻²] | $\tau$ [fs] | mobility [cm²/V s] | mean free path [nm] | verdict |
|---|---|---|---|---|---|
| 0.1 | 8.0×10¹¹ | 127 | 11700 | 122 | too clean for CVD on polymer |
| **0.2** | **3.2×10¹²** | **64** | **2900** | **61** | **typical as-transferred CVD** |
| 0.3 | 7.2×10¹² | 43 | 1300 | 41 | typical |
| 0.4 | 1.3×10¹³ | 32 | 730 | 31 | typical chemically doped (HNO₃, AuCl₃) |
| 0.6 | 2.9×10¹³ | 21 | 330 | 20 | heavily doped; low mobility, at the edge |

$E_F=0.2$–$0.4$ eV is where an ordinary sample sits, and Eq. (4a.2) discriminates
inside that window through the onset field (27 / 54 / 81 kV/cm for 0.2 / 0.4 /
0.6 eV).

This is corroborated independently of the optics. Hassanpour Amiri *et al.* [18]
identify the unintentional doping of transferred CVD graphene as ionic residue from
the copper etch, with a **typical density of $4\times10^{12}$ cm⁻²**, and show that an
aqueous-ammonia wash removes most of it (the Dirac voltage returns near zero and the
geometrically normalised mobility exceeds $2.4\times10^{4}$ cm²/V s). On the Dirac
cone that density is
$$
k_F=\sqrt{\pi n}=0.0188\ \text{a.u.},\qquad E_F=\hbar v_Fk_F=0.224\ \text{eV},
\qquad E_{\rm sat}=30\ \text{kV/cm}, \tag{4a.2b}
$$
i.e. within 12 % of the $E_F=0.2$ eV the transmission measurement gives on its own,
and it places the current-saturation onset squarely inside the 1–100 kV/cm DAST
range. A sample that has *not* had such a wash is therefore expected to sit at
$E_F\simeq0.22$ eV and to start brightening near 30 kV/cm; an ammonia-washed
(doping-free) sample should behave like the intrinsic curve of §4a.3 instead --
darkening with field, not brightening. That is a sharp, cheap experimental test of
the mechanism. The $E_F=0.6$ eV used in the $48^2$ runs of §4a.3 is therefore a
*mesh-affordable proxy*, higher than a typical sample, and its curve maps onto the
sample by the rescaling of step 3.

**Temperatures.** Three distinct ones enter and must not be conflated.
`sbe_temp_init_k` sets the occupation the run starts from, `sbe_eph_temperature_k`
the phonon bath the dissipators relax into -- both 300 K for a room-temperature
measurement -- while the *carrier* temperature is an output. It follows from the
change of the electronic energy $\Delta E_{\rm all}$, which for a Dirac gas at fixed
density fixes $T_e$ through $\varepsilon(\mu,T_e)$:

| $E_0$ [kV/cm] | $\Delta E_{\rm all}$ [eV/cell] | $T_e$ (peak) | $T_e$ (384 fs) | $\tau$ [fs] | mean free path [nm] |
|---|---|---|---|---|---|
| 10 | $2.5\times10^{-5}$ | 361 K | 361 K | 141 | 136 |
| 100 | $2.28\times10^{-3}$ | **2050 K** | 2038 K | **60** | 58 |

(48², $E_F=0.6$ eV, `diss`, lattice at 300 K throughout; the $\tau$ values carry the
mesh caveat of §4a.5.5.) At 2050 K the Drude weight of §4a.3 has moved by about 10 %,
while $\tau$ falls by a factor 2.3 within this mesh -- the statement that the
bleaching is carried by the scattering time and the drift rather than by the Drude
weight rests on that *ratio*, which is far larger than the ~10 % the Drude weight can
supply, not on the absolute $\tau$. With
`yn_sbe_rana_te = 'y'` (the `mem` variant) the same $T_e$ is fitted from the
distribution each ring step and written to `*_sbe_te.data`.

*Caveat on the channel ledger.* For a doped initial state the per-channel column
`dE_eph` of `*_sbe_channels.data` is a **gross** exchange counter, not a net loss: it
grows at $\approx7\times10^{-6}$ eV/cell/fs from $t=0$ in both runs above, i.e.
identically with and without an appreciable field. That this is bookkeeping and not
physics is settled by the total electronic energy, which at 10 kV/cm moves by only
$+2.5\times10^{-5}$ eV/cell over 384 fs: the doped Fermi sea is stationary under the
dissipators to one part in $10^{5}$ of its energy, and its carrier number is conserved
to 0.14 % (1.3 % at 100 kV/cm). Use $\Delta E_{\rm all}$, not `dE_eph`, for the energy
balance of a doped run.

**Step 3 -- check the doping against the mesh you can afford, and rescale if not.**
§4a.2 requires $k_F\gtrsim3\,|\mathbf b|/N$, i.e.
$$
N \gtrsim 3\,|\mathbf b|/k_F = 4.68/k_F \quad\Longrightarrow\quad
N\gtrsim280\ (E_F=0.2\ \text{eV}),\quad 140\ (0.4),\quad 93\ (0.6). \tag{4a.3}
$$
If the sample's doping is out of reach, run a *larger* $E_F$ on the mesh you have
and map back. In the collisionless limit the displaced Fermi disc gives
$J = (k_F^2v_F/\pi)\,g(A/k_F)$ with a single universal $g$: the **shape** of
$\sigma(E_0)/\sigma(0)$ is a function of $A_0/k_F$ alone, so a run at $E_F'$
reproduces the sample's curve with the field axis scaled by $k_F/k_F'=E_F/E_F'$,
while the conductivity scale itself goes as $E_F$. The absolute transmission does
*not* transfer (it depends on $Z_0\sigma$ through Eq. 4), so compare
$\sigma(E_0)/\sigma_{\max}$, not $T$, when rescaling. The 48², $E_F=0.6$ eV scan of
§4a.3 is exactly such a proxy: its onset at 81 kV/cm maps to 27 kV/cm at the
sample's 0.2 eV.

| $E_F$ [eV] | $n_{2D}$ [cm⁻²] | $E_{\rm sat}$ [kV/cm] | $N$ for $k_F\ge3\Delta k$ |
|---|---|---|---|
| 0.1 | 8.0×10¹¹ | 13 | 560 |
| 0.2 | 3.2×10¹² | 27 | 280 |
| 0.4 | 1.3×10¹³ | 54 | 140 |
| 0.6 | 2.9×10¹³ | 81 | 93 |

### 4a.1 Initial occupation
Everything above describes an *intrinsic* sheet: integer filling, $f_v=2$, $f_c=0$.
A real CVD sample is a metal. `sbe_ef_ev` (the Fermi level measured from the
undoped one -- the Dirac point here, mid-gap in a semiconductor) and
`sbe_temp_init_k` replace the integer filling by
$$
f_n(\mathbf k)=\text{occ}_{\max}\,f_{\rm FD}\!\big(\varepsilon_n(\mathbf k);\ \mu,\ T_{\rm init}\big),
\qquad \mu=E_F^{\rm undoped}+\texttt{sbe\_ef\_ev}, \tag{4a.4}
$$
with the added charge left uncompensated (a gated or adsorbate-doped sheet) and
reported per cell and as a sheet density. Three couplings in the solver had to
follow.

**(i) The pure-gauge reference must stay undoped.** The restoration of §6a
subtracts the adiabatic ground-state current of the truncated $H_{\mathbf k}(\mathbf A)$.
For a *filled* band that current is the truncation artifact and nothing else. For
a *partially filled* band the same object is the physical intraband response --
the shifted Fermi sea *is* the Drude current -- so subtracting it would delete the
quantity a doped run exists to compute. The solver therefore keeps the undoped
$T\to0$ filling in `gs%occup_ref` and uses it, and only it, in Eq. (10). What
remains uncorrected is the doped carriers' own truncation error: each carries an
uncancelled diamagnetic unit against a Drude weight per carrier
$\langle\partial^2\varepsilon/\partial k_a^2\rangle=v_F/k_F$, i.e. a relative
error $k_F/v_F\approx3.8\,\%$ at $E_F=0.2$ eV, falling as $1/E_F$.

**(ii) The dressed reference follows the doped baseline.** `dressed_ref_delta`
now accepts the initial diagonal $f^{(0)}$ and returns
$\delta_a=\sum_b f^{(0)}_b|W_{ba}|^2-f^{(0)}_a$ -- the adiabatically rotated
initial state minus the initial state. Trace-neutral by unitarity, zero at
$\mathbf A\to0$, and identical to the old form when $f^{(0)}$ is
$\{\text{occ},\dots,\text{occ},0,\dots\}$.

**(iii) The basis-edge monitor measures an excess.** $P_{\rm top}$ is now
$\max_{\mathbf k}[f_{\rm top}(\mathbf k)-f^{(0)}_{\rm top}(\mathbf k)]$: a metal
legitimately fills its top band at every $\mathbf k$ inside the Fermi surface,
which is not a velocity-gauge basis failure.

### 4a.2 The mesh a Fermi surface needs
A uniform mesh represents a Fermi sea only if $k_F=E_F/\hbar v_F$ exceeds several
spacings $|\mathbf b|/N$. The number of mesh points inside the Fermi disc per
valley is $\pi k_F^2/(\Delta k^2\sqrt3/2)$:

| $E_F$ [eV] | $k_F$ [a.u.] | $n_{2D}$ [cm⁻²] | $N=48$ | $N=147$ | $N=300$ | $N=600$ |
|---|---|---|---|---|---|---|
| 0.2 | 0.0167 | 3.2×10¹² | 0.9 | 9.0 | 37.6 | 150 |
| 0.4 | 0.0335 | 1.3×10¹³ | 3.7 | 36.1 | 150 | 601 |
| 0.6 | 0.0502 | 2.9×10¹³ | 8.4 | 81.2 | 338 | 1353 |

Measured: at $N=147$, $E_F=0.2$ eV, $T=300$ K the solver reports
$n=3.27\times10^{12}$ cm⁻² against the analytic $3.36\times10^{12}$ (2.6 %); at
$N=24$ it reports $1.8\times10^{10}$ -- the Fermi circle contains no mesh point at
all. The start-up banner counts the partially occupied points and warns below 20.

![initial level occupation, undoped and doped, on the 147x147 mesh](figures/graphene_doped_levels.png)

*The initial density matrix the solver starts from, on the same $147^2$ mesh
(`plot_occupation.py`, which applies Eq. (4a.4) to the ground-state files).* **Left** --
the undoped filling: the valence cone full, the conduction cone empty, nothing
partially occupied, $n_{2D}=0$. **Middle** -- $E_F=0.2$ eV, $T_{\rm init}=300$ K:
the conduction cone is filled up to $E_F$, which adds $1.71\times10^{-3}$
electrons per cell, i.e. $n_{2D}=3.27\times10^{12}$ cm⁻² (analytic
$3.36\times10^{12}$), and leaves 36 partially occupied k-points -- these carry the
whole intraband response. **Right** -- the radial profile of the conduction
occupation: at this doping the Fermi disc holds one full mesh shell plus a
partially filled second one, which is why the density is good to 3 % while the
Drude weight, weighted by $\partial^2\varepsilon/\partial k^2\propto1/k$, is still
$\approx30\,\%$ low. Run this picture before any doped production run.
The Drude *weight* converges more slowly than the density, because
$\partial^2\varepsilon/\partial k_a^2=v_F\sin^2\theta/k$ weights the innermost
shells. Measured with the trajectory fit of §4a.3 at 1 kV/cm (linear regime):

| mesh | $E_F$ [eV] | $k_F/\Delta k$ | partially occupied points | $n_{2D}$ vs analytic | $D_{\rm fit}/E_F$ | $D_{\rm spec}/E_F$ |
|---|---|---|---|---|---|---|
| $147^2$ | 0.2 | 1.58 | 36 | $3.27$ vs $3.36\times10^{12}$ (−2.6 %) | 0.659 | — |
| $300^2$ | 0.2 | 3.22 | 116 | $3.345$ vs $3.36\times10^{12}$ (−0.5 %) | **0.930** | — |
| $147^2$ | 0.4 | 3.15 | 72 | $1.284\times10^{13}$ | 0.894 | — |
| $72^2$ | 0.4 | 1.54 | — | — | 0.846 | 0.858 |
| $48^2$ | 0.6 | 1.54 | 8 | $3.04$ vs $2.89\times10^{13}$ (+5 %) | 0.888 | 0.900 |
| $72^2$ | 0.6 | 2.32 | 22 | — | 0.961 | 0.976 |
| $111^2$ | 0.6 | 3.57 | — | — | 0.954 | 0.970 |
| $147^2$ | 0.6 | 4.73 | — | — | **0.998** | **1.013** |

$D_{\rm spec}$ is the same weight read off the reactive response with **nothing
fitted**: a collisionless metal has $\sigma(\omega)=iD/\pi\omega$, so each bin of the
driven band gives $D(\omega)=-\pi\omega\,\mathrm{Im}\,\sigma(\omega)$ directly
(`drude_check.py`, column `D_spec`). It matters because for a coherent run the
trajectory fit of §4a.3 has to determine an essentially infinite $\tau$ at the same
time, and $D_{\rm fit}$ rides on that; the two nevertheless agree to $\sim1\,\%$
throughout, which is a check on both.

The density is good to a few per cent as soon as the Fermi circle contains a shell;
the Drude weight needs three or four. Read against $k_F/\Delta k$ rather than against
the mesh, the deficit collapses onto one curve regardless of doping — $0.85$–$0.90$ at
$1.5$ spacings, $0.97$ at $2.3$–$3.6$, and $1.00$ at $4.7$. Quantitative Drude work
therefore wants $k_F\gtrsim3\,\Delta k$ (Eq. 4a.3); below that the deficit is a
field-independent scale error, so the *shape* of $\sigma(E_0)$ survives while its
absolute value does not.

**How much of that last row is a real convergence, though?** Every number in the table
is `nstate = 2`. Repeating the $147^2$, $E_F=0.6$ eV run at `nstate = 4` (1 kV/cm,
$\Delta t=0.1$ fs, energy ledger $5\times10^{-9}$ eV per cell, no leakage) gives
$D_{\rm spec}=0.548$ eV instead of $0.608$, $\mathrm{Re}\,\sigma=23.8$ instead of
$27.0\,\sigma_{\rm univ}$, and $T=0.735$ instead of $0.700$ — an $11\,\%$ change in the
sheet response from the basis alone, at a field so small that the vector potential
cannot mix anything. Two things follow.

* The $D/E_F\to1$ agreement at `nstate = 2` is a **consistency check, not an
  independent confirmation**. Near K the two-band EPM *is* the Dirac model
  (Eq. 4a.9), so the adiabatic band derivative $\partial E/\partial A$ that the
  restoration leaves behind for the doped carriers is the exact Dirac velocity, and
  $D=E_F$ follows almost by construction once the disc is sampled well enough.
* A **doped** sheet is therefore not `nstate`-converged in the way the intrinsic one
  is. The restoration of §6a makes the *undoped* filling basis-independent to $10^{-6}$
  in $T$, because it subtracts the adiabatic ground-state current of exactly the same
  truncated Hamiltonian. The carriers *added* on top of the reference have no such
  subtraction, and the natural reading of the $11\,\%$ is the second-order repulsion of
  the $\pi^*$ level by the bands above it, $\Delta E\propto A^2$, which is precisely a
  shift of the Drude weight. Absolute doped conductances quoted from this solver carry
  that systematic; ratios and shapes at fixed `nstate` do not.

### 4a.3 What can bleach a doped sheet
The sheet conductance of a Drude metal is $\sigma_{dc}=D\tau/\pi$ with the
Dirac-cone Drude weight
$$
D(\mu,T)=2k_BT\,\ln\!\big[2\cosh(\mu/2k_BT)\big]\ \longrightarrow\ |\mu|\quad(\mu\gg k_BT). \tag{4a.7}
$$
A field-induced *rise* of transmission must therefore reduce $D$, reduce $\tau$, or
break the linear relation between current and field. All three are separable:

1. **Drude weight, heating at fixed density.** As $T_e$ rises, $\mu$ falls, but
   thermally generated pairs add weight; $D$ passes a shallow minimum and returns.
   Over 300–3000 K the deepest excursion is $D/D_0=0.88$ at $10^{13}$ cm⁻² and
   $0.90$ at $3\times10^{12}$ cm⁻² (`tests/test_doped_drude.py`). **Heating alone
   is worth ~10 %.**
2. **Momentum relaxation.** $v_FA_0=0.74$ eV at 100 kV/cm, far above the
   optical-phonon thresholds (E$_{2g}$ 196 meV, A$_1'$ 160 meV): twice per cycle
   the whole distribution is pushed over the emission threshold, and every
   emission randomises momentum. This is the channel the solver already has
   (§2); it needs a dissipative run on a Fermi-surface-resolving mesh.
3. **Current saturation on the cone.** The Fermi sea is displaced by $\mathbf A(t)$.
   When the excursion $A_0$ exceeds $k_F$, the displaced sea is no longer a small
   perturbation: the drift velocity saturates at $v_F$ and the differential
   conductivity falls as $k_F/A_0$. For the DAST transient
   $A_0=6.21\times10^{-4}\,$a.u. per kV/cm, so
   $$
   A_0=k_F \iff E_0 = 27\ \text{kV/cm}\ (E_F=0.2\ \text{eV}),\quad 81\ \text{kV/cm}\ (E_F=0.6\ \text{eV}). \tag{7a}
   $$
   The collisional excursion $eE\tau/\hbar$ crosses $k_F$ at the same place for
   $\tau\simeq60$ fs, so at THz the two scales coincide. Saturation is a
   *coherent* effect and needs no dissipator at all. The calculated curve (doped
   coherent runs, $N=48$, $E_F=0.6$ eV, $T_{\rm init}=300$ K, self-consistent
   sheet; $D_{\rm fit}$ from the run's own current):

   | $E_0$ [kV/cm] | 1 | 10 | 30 | 100 | 200 | 300 | 500 | 1000 |
   |---|---|---|---|---|---|---|---|---|
   | $T$ | 0.7286 | 0.7257 | 0.6990 | 0.6622 | 0.7202 | 0.7716 | 0.8122 | 0.8263 |
   | $D_{\rm fit}$ [eV] | 0.533 | 0.537 | 0.580 | 0.664 | 0.546 | 0.424 | 0.325 | 0.290 |
   | $\mathrm{Re}\,\sigma/\sigma_{\rm univ}$ | 23.6 | 23.8 | 25.8 | 30.3 | 25.4 | 20.5 | 15.8 | 13.4 |

   $T$ is flat through the linear regime, dips near 60 kV/cm, then rises monotonically
   once $A_0>k_F$ (81 kV/cm at this doping). **Most of the dip is not physics**:
   §4a.5.5 refines this run on three finer meshes and shows that three quarters of the
   darkening is the discretization bump of a Fermi disc holding only 8 mesh points —
   the dip falls from $11.8\,\%$ here to $3.3\,\%$ once the disc holds $3.6$ spacings,
   after which it stops falling. A residual few-per-cent darkening near $A_0\simeq k_F$
   therefore survives mesh refinement and is not yet explained. The rise itself is
   robust: $\sigma$ falls from $30.3$ to
   $13.4\,\sigma_{\rm univ}$, $-56\,\%$.

   ![strong-doping proxy on the 48x48 mesh](figures/graphene_doped_proxy_ef06.png)

   *Strong-doping proxy, $48^2$, $E_F=0.6$ eV. The same mesh, the same pulse, the
   same solver settings; only the initial occupation differs.* **Left** --
   transmission against peak field: the intrinsic sheet (dark squares) darkens
   monotonically as Landau-Zener pairs are created (§7), while the doped sheet (red
   circles) is flat in the linear regime, dips into the mesh artifact of §4a.5.5,
   and then **brightens** past the shaded region, which begins at
   $A_0=k_F$. Dotted green: the transmission of the measured sample with the PET
   Fresnel loss divided out, $0.68\to0.79$. **Middle** -- the extinction the doping
   carriers alone contribute, $1-T_{\rm doped}/T_{\rm intrinsic}$: dividing by the
   intrinsic curve removes the Landau-Zener darkening and leaves the Drude response,
   which peaks at the saturation field and collapses. **Right** -- the sheet
   conductivity against the two conductances Eq. (4a.1) extracts from the measured
   transmissions, $24.7$ and $14.3\,\sigma_{\rm univ}$: calculation and measurement
   compared as sheet conductances in the same units, nothing fitted in between, the
   calculated fall across saturation ($-56\,\%$) against the measured $-42\,\%$.

**The same scan on four meshes.** The $48^2$ figure above is the coarsest member of a
series. Repeating it at $72^2$, $111^2$ and $147^2$ — same doping, same transient, same
solver settings, only $\Delta k$ changing — separates what the mesh does from what the
sheet does:

![T(E0) of a doped Dirac sheet on four k-meshes, against the continuum drift-saturation curve](figures/graphene_T_of_field_mesh.png)

*$E_F=0.6$ eV, coherent, self-consistent sheet, 100 fs ring-down; the legend gives
$k_F/\Delta k$, the number of mesh spacings inside the Fermi disc.* **Left** —
transmission against peak field. The four curves agree above $A_0=k_F$ (shaded) to
better than 0.02 and disagree below it, which is the signature of a discretization
error rather than a field scale: the dip at 60 kV/cm deepens from $3.3\,\%$ to
$11.8\,\%$ as the disc empties from 4.7 spacings to 1.5. Dotted: the continuum
drift-saturation prediction $T(\sigma_{\rm lin}G(u)/u)$ of §4a.5.3, anchored on each
series' own lowest-field point and drawn only where the cone picture behind it holds
($u\le2$); the gap to it at $u\lesssim1$ is the artifact, and it closes with the mesh.
**Right** — the sheet conductivity, with the two conductances Eq. (4a.1) extracts from
the measured transmissions. The spurious pre-saturation peak falls from $+32\,\%$ to
$+7.6\,\%$ over the series while the collapse above saturation, $-56\,\%$, does not move.
Reproduce with `field_scan_plot.py --continuum --series ...` (x14 README §7.12).

**How the substrate and the layer count enter these figures — read before comparing
with a measurement.**

* *The calculation has no substrate.* The sheet is **free-standing**: vacuum on both
  sides, $n_s=1$ in Eq. (4), and $T$, $R$, $A$ are the fluence ratios of the bare
  sheet. Nothing about PET, glass or silicon enters the propagation.
* *The measured numbers plotted beside them have the substrate divided OUT, not added.*
  An unreferenced measurement still contains the substrate's own Fresnel loss, and
  Eq. (4a.1) removes it: $T_{\rm sheet}=T_{\rm meas}/T_{\rm bare}$ with
  $T_{\rm bare}=[4n_s/(1+n_s)^2]^2=0.883$ for PET ($n_s=1.65$, two incoherent air/PET
  faces, no etalon). The dotted green lines are therefore $0.60/0.883=0.68$ and
  $0.70/0.883=0.79$, **not** 0.60 and 0.70. The opposite convention — putting the
  substrate into the calculation instead — is available through Eq. (4) with
  $n_s\neq1$ (`transmission.py --n-sub 1.65`, `drude_check.py --n-sub`), and then the
  raw 0.60/0.70 are the right comparison. The two must not be mixed; every number in
  this page uses the first.
* *One layer.* Unless stated otherwise the calculation is a **single monolayer**,
  `sbe_sheet_nlayers = 1`. For $N$ electronically decoupled layers within
  $d\ll\lambda$ the driver adds their currents in the same local field,
  $J_s\to NL_zJ_m$ (§5) — the incoherent or large-angle-twisted stack. Bernal and
  small-angle moiré bilayers change the electronic structure itself and are outside
  the model.

**One layer against two, at two dopings.** Same $72^2$ mesh, same transient with its
100 fs ring-down, coherent, self-consistent sheet; only `sbe_sheet_nlayers` and
`sbe_ef_ev` differ. $u=A_0/k_F$ is the drift parameter of §4a.5.1.

![one layer against two at two Fermi levels](figures/graphene_layers_1_vs_2.png)

*Left* — the two occupations, each with one and two layers, and the naive $T_1^2$
dotted. The dotted curve lies **above** the propagated two-layer curve everywhere, at
both dopings. *Right* — the size of that error, $100(T_1^2/T_2-1)$: $+15\,\%$ at
$E_F=0.6$ eV and $+11\,\%$ at $0.4$ eV in the linear regime, never negative. Dashed
verticals mark each doping's own saturation field $A_0=k_F$ (81 and 54 kV/cm), where
the transmission of both stacks turns up. Reproduce with `layers_plot.py`.

| $E_0$ [kV/cm] | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|
| **$E_F=0.6$ eV** ($E_{\rm sat}=81$ kV/cm) | $u=0.04$ | 0.12 | 0.37 | 1.24 | 3.71 |
| one layer, $T_1$ | 0.7093 | 0.7053 | 0.6679 | 0.6713 | 0.7739 |
| two layers, $T_2$ | 0.4362 | 0.4341 | 0.4128 | 0.3778 | 0.4715 |
| $T_1^2$ (the naive estimate) | 0.5031 | 0.4975 | 0.4461 | 0.4506 | 0.5989 |
| $T_1^2/T_2$ | **1.153** | **1.146** | 1.081 | 1.193 | 1.270 |
| **$E_F=0.4$ eV** ($E_{\rm sat}=54$ kV/cm) | $u=0.06$ | 0.19 | 0.56 | 1.86 | 5.57 |
| one layer, $T_1$ | 0.8487 | 0.8422 | 0.7939 | 0.8358 | 0.8772 |
| two layers, $T_2$ | 0.6493 | 0.6434 | 0.5847 | 0.5894 | 0.7226 |
| $T_1^2$ | 0.7203 | 0.7093 | 0.6303 | 0.6986 | 0.7695 |
| $T_1^2/T_2$ | **1.109** | **1.103** | 1.078 | 1.185 | 1.065 |

Three things to read off it.

* **The measured $T\cdot T$ of a bilayer is an *over*estimate here, at every field and
  both dopings** — by 15 % at $E_F=0.6$ eV and 10 % at 0.4 eV in the linear regime.
  That is the reactive branch of Eq. (6a): a coherent doped sheet at 14 meV is an
  inductor, $\mathrm{Im}\,z\gg\mathrm{Re}\,z$ (band-averaged $z=0.017-1.195i$ at
  $E_F=0.6$ eV), and for $\mathrm{Im}\,z$ the second-order term raises $T_1^2/T_2$
  above 1, i.e. $\sqrt{T_2}<T_1$: inverting a measured bilayer transmission as
  $\sqrt{T_2}$ makes the monolayer look *darker* than it is and therefore
  **overstates** its conductance. A sample whose absorption is dominated by momentum
  relaxation sits on the dissipative branch, where the inequality reverses — see the
  propagated check in §4a.5.6 — so the direction has to be decided from $\sigma$'s
  phase, not assumed.
* **The single-frequency formula is not enough for this pulse.** Eq. (6a) with the
  band-averaged $z$ predicts $T_1^2/T_2=1.31$ where the propagated pair gives 1.153.
  Using the sheet's own frequency-resolved $\sigma(\omega)$ instead
  (`transmission.py --predict-layers 2`) gives 1.169 — within 1.5 % — because a
  single-cycle transient has a band as wide as its centre and $|\sigma|\propto1/\omega$
  across it. The prediction stays that good only while the response is linear: it is
  off by 8 % at $u=0.56$ and by 27 % at $u=5.6$, always in the same direction, because
  each layer of the real stack sees the field reduced by *both* currents and therefore
  sits at a smaller $u$, where it is more conductive — the true stack is darker than
  any scaling of one run.
* **The stacking itself is exact.** Bin by bin across the driven band the two-layer
  run returns exactly twice the one-layer $\sigma(\omega)$ (median $|\sigma_2/\sigma_1|
  = 1.997$ at $E_F=0.6$ eV, $2.007$ at 0.4 eV). The band-*averaged*
  $\mathrm{Re}\,\sigma$ ratio is only 1.83, which is a property of the diagnostic and
  not of the solver: the two-layer transmitted field is filtered more strongly at the
  low-frequency end of the band, where $|\sigma|$ is largest, so the average moves.

**Does the bilayer brighten at the sample's own doping?** Yes — and the monolayer at
the same doping essentially does not. $297^2$ (the first mesh that resolves
$k_F=0.0167$ a.u., $k_F/\Delta k=3.19$), $E_F=0.2$ eV, coherent, saturation at
27 kV/cm:

| $E_0$ [kV/cm] | 1 | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|---|
| $u=A_0/k_F$ | 0.04 | 0.37 | 1.11 | 2.23 | 3.71 | 11.1 |
| one layer, $T$ | 0.9502 | 0.9447 | 0.9438 | 0.9475 | 0.9438 | 0.9334 |
| two layers, $T$ | 0.8540 | 0.8445 | **0.8346** | 0.8616 | **0.8710** | 0.8539 |
| $|Z_0\sigma|$, one layer | 0.377 | 0.401 | 0.401 | 0.329 | 0.284 | 0.278 |

![one against two layers at the sample doping](figures/graphene_layers_ef02.png)

The bilayer rises **+4.4 %** from its minimum at 30 kV/cm (just past saturation) to
100 kV/cm; the monolayer moves by $0.4\,\%$ and then drifts down. Both sheets saturate
identically — $|Z_0\sigma|$ falls by $29\,\%$ from its peak in each — so the difference
is not in the physics but in how much transmission a given change in $\sigma$ can buy.
For $T=|2/(2+z)|^2$ with $z$ nearly imaginary, $T\simeq(1+y^2/4)^{-1}$ with
$y=|{\rm Im}\,z|$, so
$$
\frac{\delta T}{T}=-\frac{y/2}{1+y^2/4}\,\delta y \tag{4a.17}
$$
— the leverage grows with $y$, and doubling the layers doubles both $y$ and $\delta y$.
A sheet that already transmits 95 % has almost nothing to gain; one at 85 % has four
times as much. Meanwhile Landau–Zener absorption grows the same way in both ($A$ from
0.008 to 0.031 in one layer), and in the monolayer it cancels the gain outright.

**This is why the measurement needs the ring, not a bigger mesh.** The sample's
$Z_0\sigma=0.565$ is 1.5× the coherent monolayer's $0.377$ and its $\mathrm{Re}\,\sigma$
is six times larger, because the real conductance is $D\tau/\pi$ with a finite $\tau$
that a collisionless run does not have. At that conductance Eq. (4a.17) turns the same
$29\,\%$ saturation into roughly ten points of transmission — the measured
$0.68\to0.79$. The coherent run at $E_F=0.2$ eV cannot show the effect; the $E_F=0.6$
eV proxy of §4a.3 stands in for exactly that missing conductance.

**At the sample's own doping (coherent, against the intrinsic control).** The same scan at $E_F=0.2$ eV on the production
$147^2$ mesh is quieter in the raw transmission, because a $3\times10^{12}$ cm⁻²
sheet without dissipators is only a weak inductor ($T=0.968$ at 1 kV/cm, not 0.68 --
the coherent run has no momentum relaxation, so it cannot produce the measured
*absorption*, only the reactive screening). The doping's own signature is then read
against the intrinsic control:

![doped vs intrinsic sheet at the sample doping, 147x147](figures/graphene_doped_vs_intrinsic.png)

| $E_0$ [kV/cm] | 1 | 10 | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|---|---|
| $T$ doped ($E_F=0.2$ eV) | 0.9682 | 0.9666 | 0.9595 | 0.9629 | 0.9442 | 0.9090 |
| $T$ intrinsic | 0.999995 | 0.999993 | 0.999498 | 0.986694 | 0.933550 | 0.916273 |
| extinction added by the doping | 3.18 % | 3.34 % | **4.00 %** | 2.41 % | −1.14 % | 0.79 % |

The doping-induced extinction **peaks at 30 kV/cm** -- the predicted saturation field
for this doping is 27 kV/cm (Eq. 4a.2) -- and then collapses by a factor of five,
going briefly negative at 300 kV/cm, where the doped sheet transmits *better* than
the undoped one because its Drude response has saturated while its occupied states
Pauli-block part of the Landau-Zener pair creation. The onset therefore comes out at
the field Eq. (4a.2) predicts from $k_F$ alone, at the doping the transfer literature
[18] reports, without anything having been tuned. Reproduce both figures with
`bash samples/exercise_x14_graphene_self_induced_transparency/run_field_scan.sh`
(`NK`, `EF` as in §4a.0). The fitted $D$ sits below the analytic
   $E_F$ by the mesh factor of §4a.2 ($D_{\rm fit}/D_{\rm eq}=0.66$ / $0.89$ /
   $0.89$ for $E_F=0.2$ eV at $N=147$, $0.4$ eV at $N=147$, $0.6$ eV at $N=48$) --
   a scale error common to all fields, so the shape of $T(E_0)$ is intact.

The fit that produces those numbers is the run's own current obeying
$\dot J_s=(D/\pi)E_{\rm tot}-J_s/\tau$; least squares over the driven window
returns $D$ and $\tau$ separately (`drude_check.py`), so the three mechanisms above
are read off rather than assumed.

### 4a.4 Comparison with a measured sample
A monolayer on PET transmitting 60 % with the substrate included is, with
$n_{\rm PET}=1.65$ and two incoherent faces (bare substrate 88.3 %),
$Z_0\sigma=0.565$, i.e. $\sigma=24.7\,\sigma_{\rm univ}$; 70 % is
$14.3\,\sigma_{\rm univ}$ -- a 42 % fall. $E_F=0.2$ eV with $\tau=63$ fs
reproduces the low-field value exactly. Mechanism 1 can supply ~10 % of the 42 %; mechanism 3
alone produces $-56\,\%$ in the calculated scan above, and mechanism 2 adds to it
in the same direction. All three read "the carriers stop responding linearly",
not "the carriers disappear". Note that the onset
Eq. (7a) sits at 27 kV/cm for that doping -- inside the measured range, and the
approach to saturation is gradual ($\propto k_F/A_0$), which is why the
transmission creeps up by ten points instead of jumping.

### 4a.5 Why the transmission of a doped sheet falls and then rises

At normal incidence the transmission of a free-standing sheet is fixed by its sheet
conductance alone (Eq. 4),
$$
T=\Big|\frac{2}{2+Z_0\sigma}\Big|^{2},
$$
a monotone decreasing function of $|\sigma|$. The whole of the non-monotonic $T(E_0)$
therefore reduces to one question: how does the sheet conductance of a doped Dirac
sheet depend on the strength of the drive? Three contributions enter $\sigma$, and
they respond to the field in different ways and on different scales.

**The semiclassical picture first.** Everything in this section is one kinematic fact
about a cone, and it is worth having in front of the equations:

![semiclassical kinematics of a doped Dirac cone in a velocity-gauge field](figures/graphene_cone_kinematics.png)

*Where the field pushes, where the electron turns, when it can jump, and where Pauli
forbids it.* In the velocity gauge the canonical label $\mathbf k$ **does not move**;
the field enters only through $H_{\mathbf k}(\mathbf A)=v_F\boldsymbol\sigma\cdot
(\mathbf k+\mathbf A)$. So the occupied disc stays where it is and what travels
through it is the point of instantaneous degeneracy, $\mathbf k=-\mathbf A(t)$, the
only place the two branches touch.

* **(a), (b) — where the electron turns.** The velocity of an occupied conduction
  state is $v_F(\mathbf k+\mathbf A)/|\mathbf k+\mathbf A|$: it points radially *away
  from the moving degeneracy point*, and its modulus is $v_F$ whatever the field does.
  At $A=0$ the arrows point outward from the centre and cancel — no current. As $A$
  grows they all swing towards $+\mathbf A$, and by $u=2$ they are nearly parallel:
  the drift cannot exceed $v_F$, so the current saturates and the sheet brightens.
* **(a) versus (b) — where Pauli forbids the jump.** For $|A|<k_F$ the degeneracy sits
  *inside* the occupied disc: the upper-cone state there is already filled, so no
  interband transition is possible and the response is pure drift. For $|A|>k_F$ the
  degeneracy has left the disc, the upper cone is empty there, and pairs are created in
  a strip of half-width $\sqrt{E/\pi v_F}$ about its path.
* **(c) — the consequence.** $J/nev_F=G(u)$ rises linearly and then flattens at 1;
  the chord conductivity $G(u)/u$ is flat and then falls as $1/u$.
* **(d), (e) — the same thing in energy and in time.** The cone slides by $-A$ while
  the occupation stays put; the drive's own $A(t)/k_F$ shows when the excursion leaves
  the shaded disc, which is when panel (b) applies.

**One number does both.** $A_0=k_F$ ends the linear regime *and* opens pair creation —
27 kV/cm at $E_F=0.2$ eV, 81 kV/cm at 0.6 eV. Reproduce the figure with
`cone_kinematics.py --ef-ev 0.6 --drive DAST_E100kVcm.txt`.

#### 4a.5.1 The Fermi radius $k_F$: the yardstick the doping sets

Doping puts electrons into the conduction cone up to $E_F$, i.e. it fills a disc in
reciprocal space of radius
$$
k_F=\frac{E_F}{\hbar v_F}=\sqrt{\pi n_{2D}}\qquad(g=g_sg_v=4), \tag{4a.5}
$$
$k_F$ being the only reciprocal length the doping introduces ($n_{2D}=k_F^2/\pi$;
$E_F=0.2$ eV gives $k_F=0.0167\ a_0^{-1}$ and $3.2\times10^{12}$ cm⁻²). It is not a
fitted quantity and not a numerical parameter: it is the radius of the initial
occupation Eq. (4a.4) puts on the mesh, and `plot_occupation.py` prints and draws it.

The drive enters through the vector potential, and in the velocity gauge $\mathbf A$
is *literally a displacement in reciprocal space*: $H_{\mathbf k}(\mathbf A)$ has the
spectrum of $\mathbf k+\mathbf A$, so the occupied set is rigidly displaced by
$\mathbf A(t)$ while its canonical labels stay put. The relevant measure of the drive
is therefore the excursion
$$
A_0=\max_t|\mathbf A(t)|\ \ \big(=6.213\times10^{-4}\,a_0^{-1}\ \text{per kV/cm for
the scaled DAST transient}\big),
$$
and the single dimensionless control parameter of the problem is
$$
u=\frac{A_0}{k_F}=\frac{\text{displacement of the Fermi sea}}{\text{its own radius}}. \tag{4a.6}
$$
Everything below is a statement about $u$. (With scattering the sea only reaches
$A\simeq eE\tau/\hbar$ before it is randomised, so strictly $u=\min(A_0,eE\tau/\hbar)/k_F$;
at 3 THz with $\tau\approx60$ fs the two agree to tens of per cent, §4a.0.)

#### 4a.5.2 Regime I, $u\ll1$: the linear Drude sheet

The displaced disc carries $J=(D/\pi)A$ with $D=E_F$, so
$\sigma(\omega)=(D/\pi)/(\tau^{-1}-i\omega)$ and $T$ is *field-independent*. This is
the plateau at $T=0.7286$ (1 kV/cm) to $0.7257$ (10 kV/cm) in §4a.3, and it is where
the measured low-field conductance of §4a.4 is read.

#### 4a.5.3 Regime II, $u\gtrsim1$: drift saturation — the brightening

This is the mechanism behind the rise. It is a property of the Dirac cone, it is
derived here from the same Hamiltonian the solver propagates, and it involves no
dissipation, no heating and no change of carrier number.

**Model.** Near K the two-band EPM *is* the Dirac model: the momentum matrix is
$\mathbf p_{\mathbf k}=\partial H/\partial\mathbf k=v_F\boldsymbol\sigma$, so the
velocity-gauge Hamiltonian of §2 becomes exact,
$$
H_{\mathbf k}(\mathbf A)=\varepsilon_{\mathbf k}+\mathbf A\cdot\mathbf p_{\mathbf k}
=v_F\,\boldsymbol\sigma\cdot(\mathbf k+\mathbf A), \tag{4a.9}
$$
with eigenvalues $\pm v_F|\mathbf k+\mathbf A|$ and band velocity
$$
\mathbf v_\pm(\mathbf k,\mathbf A)=\pm v_F\frac{\mathbf k+\mathbf A}{|\mathbf k+\mathbf A|}
\qquad\Longrightarrow\qquad |\mathbf v_\pm|=v_F\ \ \text{always}. \tag{4a.10}
$$
Equation (4a.10) is the whole of the effect: the field can only **turn** the velocity,
never lengthen it.

**Occupations.** In the velocity gauge the canonical labels $\mathbf k$ do not move;
the field enters only through $H_{\mathbf k}(\mathbf A(t))$. While the evolution is
adiabatic the occupation numbers attached to the labels are constant, so the occupied
set stays the disc $|\mathbf k|\le k_F$ *in label space*, while every state's velocity
is that of the displaced point $\mathbf k+\mathbf A$. (Adiabaticity is not an
assumption here but a consequence: the only place it can fail is the instantaneous
degeneracy at $\mathbf k=-\mathbf A$, and for $A<k_F$ that point lies inside the disc,
where both branches are occupied and no transition is possible — §4a.5.4.)

**Current.** With $g=g_sg_v=4$ and $\hat e=\hat{\mathbf A}$, the valence sea cancelling
against the undoped pure-gauge reference of §6a,
$$
J(A)=\frac{g}{(2\pi)^2}\int_{|\mathbf k|\le k_F}v_F\,
\frac{(\mathbf k+\mathbf A)\cdot\hat e}{|\mathbf k+\mathbf A|}\,d^2k. \tag{4a.11}
$$
Scale $\mathbf k=k_F\mathbf x$ and put $u=A/k_F$. With $n=gk_F^2/4\pi=k_F^2/\pi$,
$$
J(A)=n\,e\,v_F\,G(u),\qquad
G(u)=\frac1\pi\int_{|\mathbf x|\le1}\frac{x_\parallel+u}{|\mathbf x+u\hat e|}\,d^2x. \tag{4a.12}
$$

**Reduction.** The integrand is $\partial|\mathbf x+u\hat e|/\partial u$, so
$\pi G=\partial_u\!\int_{|\mathbf x|\le1}|\mathbf x+u\hat e|\,d^2x$. More useful is the
shift $\mathbf y=\mathbf x+u\hat e$, which turns Eq. (4a.12) into the mean of
$\cos\varphi_y$ over a unit disc whose centre is at distance $u$ from the origin.
Integrating over $\varphi$ at fixed $|\mathbf y|=\rho$ (the disc boundary is
$\cos\varphi\ge(\rho^2+u^2-1)/2\rho u$, and shells that lie wholly inside contribute
zero) collapses the double integral to a single quadrature,
$$
G(u)=\frac{1}{\pi u}\int_{|1-u|}^{1+u}\!\!
\sqrt{\big[(\rho+u)^2-1\big]\big[1-(\rho-u)^2\big]}\;d\rho, \tag{4a.13}
$$
a complete elliptic integral (the radicand is a quartic with roots
$\rho=\pm1\pm u$), evaluated by quadrature in `drift_saturation.py`. Exact properties,
all verified against the two-dimensional integral to $10^{-8}$:
$$
G(1)=\frac{8}{3\pi}=0.848826\ldots,\qquad G'(1)=\frac{4}{3\pi}=\tfrac12G(1),
\qquad G(u)=u\,G(1/u), \tag{4a.14}
$$
$$
G(u)=u\Big(1-\frac{u^2}{8}+O(u^4)\Big),\qquad
G(u)=1-\frac{1}{8u^2}+O(u^{-4}). \tag{4a.15}
$$
The duality in Eq. (4a.14) maps the weak-field and strong-field branches onto each
other, and the same coefficient $1/8$ governs the leading correction on both.

**The two response ratios.** A field scan reports the current the *same* peak field
produces, i.e. the chord (secant) response, while a weak probe on top of a strong
pump would see the slope:
$$
\boxed{\ \frac{\sigma_{\rm eff}(E_0)}{\sigma_{\rm lin}}=\frac{G(u)}{u}
\ \xrightarrow[\ u\gg1\ ]{}\ \frac{1}{u}=\frac{k_F}{A_0}\propto\frac{1}{E_0}\ },
\qquad
\frac{\sigma_{\rm diff}}{\sigma_{\rm lin}}=G'(u)\ \xrightarrow[\ u\gg1\ ]{}\ \frac{1}{4u^3}. \tag{4a.16}
$$
Both at 300 K (thermal smearing changes them by less than $10^{-3}$ while
$E_F\gg k_BT$; the finite-$T$ version simply replaces the sharp disc in Eq. (4a.11)
by the Fermi weight $f_{\rm FD}(v_Fk;\mu,T)$):

| $u=A_0/k_F$ | 0.05 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| $G(u)/u$ — chord, what a scan measures | 1.000 | 0.995 | 0.980 | 0.953 | 0.912 | **0.849** | 0.627 | 0.484 | 0.329 | 0.199 |
| $G'(u)$ — differential | 0.999 | 0.985 | 0.938 | 0.853 | 0.713 | **0.424** | 0.085 | 0.033 | 0.010 | 0.002 |

**Limits recovered.** At $u\to0$, Eq. (4a.12) with Eq. (4a.15) gives
$J=nev_Fu=E_FA/\pi$, i.e. exactly the Drude form $J=(D/\pi)A$ with $D=E_F$ of §4a.3 —
the derivation reproduces the linear Drude weight rather than assuming it. At
$u\to\infty$, $J\to nev_F$: the drift velocity saturates at $v_F$ and every further
increment of field produces no further current.

**Why this is what the solver computes.** `calc_current_bloch` evaluates
$\mathrm{Tr}[(\mathbf p+\mathbf A)\rho]/V$ and subtracts the adiabatic ground-state
current of the *undoped* filled sea (§6a, Eq. 10). Writing
$\rho=\rho^{(0)}_{\rm undoped}+\delta\rho$, the subtraction removes the first term
identically and leaves $\sum_{\mathbf k}\delta f\,\langle p+A\rangle$, which is
Eq. (4a.11) plus the doped carriers' own uncancelled diamagnetic unit
($k_F/v_F=3.8\,\%$ at $E_F=0.2$ eV, 11 % at 0.6 eV, §4a.1). Evaluated on the actual
mesh this reproduces the simulated peak current to 3 % up to 300 kV/cm (§4a.5.4),
which is the direct check that Eq. (4a.12) is the solver's own physics and not a
side model.

**Consequences.** The saturation is *kinematic* (no dissipator required), *gradual*
($\propto1/E_0$, so the transmission creeps up over a decade of field rather than
switching), and its onset is set by the doping alone through $A_0=k_F$: 27 / 54 /
81 kV/cm at $E_F=0.2$ / 0.4 / 0.6 eV. With scattering the sea only reaches
$A\simeq eE\tau/\hbar$ before being randomised, so $u=\min(A_0,eE\tau/\hbar)/k_F$; at
3 THz with $\tau\approx60$ fs the two agree to tens of per cent and the onset is
unchanged.

#### 4a.5.4 Regime III, high field: interband pair creation — the darkening

The displaced-sea picture is adiabatic. It fails where the instantaneous gap
$2v_F|\mathbf k+\mathbf A|$ closes, i.e. within the Landau–Zener tube around
$\mathbf k=-\mathbf A$, and there the field creates real electron–hole pairs (§7).
Those pairs are *extra carriers*: they raise $\sigma$ and darken the sheet — this is
the whole of the intrinsic sheet's behaviour, $T:1.000\to0.90$ (§7.7 of the x14
README). In a doped sheet the same channel is **Pauli-suppressed as long as
$A_0<k_F$**, because the displaced Dirac point then lies inside the occupied disc and
the pair state at that $\mathbf k$ is already full; it switches on at essentially the
same $u\simeq1$ at which saturation begins.

Its size is measurable in the runs as the excess of the simulated current over the
adiabatic displaced-sea sum evaluated on the same mesh ($48^2$, $E_F=0.6$ eV):

| $E_0$ [kV/cm] | 1 | 10 | 30 | 60 | 100 | 200 | 300 | 500 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| $J_{\rm sim}/J_{\rm adiabatic}$ | 0.975 | 0.972 | 0.966 | 1.020 | 1.005 | 0.997 | 1.003 | 1.171 | 1.611 |

Up to 300 kV/cm the sheet is the adiabatically displaced Fermi sea to within 3 %;
only above that do created pairs add current. **So $T(E_0)$ of a doped sheet is the
competition of two channels with a common onset scale: brightening by drift
saturation, which starts at $A_0=k_F$ and grows as $1/E_0$, and darkening by pair
creation, which starts at the Landau–Zener threshold and grows as $E_0^{3/2}$.**
Saturation wins first because it acts on the pre-existing carriers, which at
$3\times10^{12}$–$3\times10^{13}$ cm⁻² outnumber the created ones until several
hundred kV/cm.

#### 4a.5.5 What the dissipators add, and what the mesh takes away

*Scattering.* $\sigma_{dc}=D\tau/\pi$, so any shortening of $\tau$ darkens the sheet
at fixed drive and brightens it as the field grows. Above the optical-phonon
thresholds (E$_{2g}$ 196 meV, A$_1'$ 160 meV) every emission randomises momentum, and
$v_FA_0=0.74$ eV at 100 kV/cm puts the whole distribution over them twice per cycle.
Switching the ring on changes the *character* of the sheet, and that part is robust.
At 100 kV/cm on $48^2$ the sheet goes from an inductive mirror ($R=0.33$, $A=0.007$,
$\tau=7.9\times10^3$ fs, i.e. collisionless) to a Drude absorber ($R=0.048$,
$A=0.231$) -- which is what a real sample is. The same happens at $72^2$:
$R=0.285\to0.020$, $A=0.010\to0.202$ at 10 kV/cm.

**The full dissipative scan** ($72^2$, $E_F=0.6$ eV, ring at a 300 K lattice, 100 fs
ring-down):

| $E_0$ [kV/cm] | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|
| $T$ | 0.7786 | 0.7789 | 0.7978 | 0.8279 | **0.8949** |
| $\mathrm{Re}\,\sigma/\sigma_{\rm univ}$ | 10.65 | 11.16 | 10.07 | 8.49 | **5.15** |
| $R$ | 0.020 | 0.026 | 0.024 | 0.020 | 0.016 |
| $A$ | 0.202 | 0.195 | 0.178 | 0.152 | 0.089 |

![T(E0) with and without the phonon ring](figures/graphene_T_of_field_ring.png)

*The dissipative scan (orange, $72^2$) against the best coherent one (red, $147^2$),
same doping and same transient.* The ring turns an inductive mirror into a Drude
absorber — $R$ falls from 0.28 to 0.02, $A$ rises from 0.01 to 0.20 — which moves the
whole curve up towards the measured band and, more importantly, removes the residual
darkening: the coherent curve still dips before the shaded saturation region, the
dissipative one is flat and then rises. The two are on different meshes because the
ring is $O(N_k^2)$ (see below), so the vertical offset between them mixes the mesh with
the physics; the *shapes* are the comparison.

Flat through 10–30 kV/cm, then a monotone rise to 0.895, with $\sigma$ down $54\,\%$
from its peak. **The mesh dip is gone**: the coherent curve on this same $72^2$ mesh
dips $7.2\,\%$ before rising (§4a.5.5 table above), the dissipative one does not dip at
all. Momentum relaxation broadens the Fermi surface over many mesh cells, which is
exactly what the discretization error needs.

*The 3 kV/cm point is omitted, and the reason is a caveat for every low-field
dissipative run.* Its residual current at the end of the record is
$3.225\times10^{-8}$ a.u. against $3.214\times10^{-8}$ at 10 kV/cm — **identical, hence
field-independent**: the ring's fixed point is not exactly the initial Fermi–Dirac
occupation, and on a finite mesh its relaxation is not isotropic, so the sheet carries
a current with no field at all. The zero-field `dark_diss` control on this mesh
measures it directly: $|J_{\rm dark}|_{\max}=2.42\times10^{-8}$ a.u., steady over the
whole 384 fs record ($2.19$, $2.31$, $1.93\times10^{-8}$ at 50, 200, 384 fs), sign
mixed, not growing. Being field-independent it matters in proportion to how weak the
drive is:

| $E_0$ [kV/cm] | 3 | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|---|
| $\|J_{\rm dark}\|/\|J\|_{\rm peak}$ | **36 %** | **14 %** | 3.3 % | 2.0 % | 1.3 % | 0.6 % |
| $T$ as reported | 0.7968 | 0.7786 | 0.7789 | 0.7978 | 0.8279 | 0.8949 |
| $T$ with the dark run subtracted | 0.8439 | 0.7826 | 0.7780 | 0.7973 | 0.8275 | 0.8947 |
| shift | **+0.047** | +0.004 | −0.001 | −0.000 | −0.000 | −0.000 |

**The rule.** Above $\sim10\,\%$ the point is unusable, and *not repairable by
subtraction either*: the dark current of a driven run is not the field-free one,
because the field changes the distribution the ring acts on, so the first-order
correction overshoots — at 3 kV/cm it lands on 0.844 where the linear-regime plateau
is 0.778. Below $\sim3\,\%$ the correction is under 0.001 and can be ignored. Between
them, subtract and quote both. Every dissipative scan therefore needs its `dark_*`
control, and `transmission.py --dark` prints the fraction and flags it past $10\,\%$.

**A number that did not survive the mesh.** On $48^2$ this scan gave $T=0.644\to0.721$
and a high-field $14.5\,\sigma_{\rm univ}$ against the measured $14.3$ — an almost exact
match. At $72^2$ the same calculation gives $0.779\to0.895$ and $5.15\,\sigma_{\rm
univ}$. The $48^2$ agreement was a coincidence of an under-resolved Fermi disc, in the
same family as the $D\to E_F$ agreement at `nstate = 2` (§6b). What is robust is the
*shape* — flat, then a rise of the right size once $A_0>k_F$ — not the absolute level,
which still carries the unconverged $\tau$ below.

**The absolute $\tau$, however, is not converged with the k-mesh** and must not be
quoted as a first-principles number:

| mesh | $E_0$ [kV/cm] | $\tau$ [fs] | mean free path [nm] | $\mathrm{Re}\,\sigma/\sigma_{\rm univ}$ | $T$ |
|---|---|---|---|---|---|
| $48^2$ | 10 | 141 | 136 | 18.1 | 0.644 |
| $48^2$ | 100 | 60 | 58 | 14.5 | 0.721 |
| $72^2$ | 10 | **33** | 31 | 10.6 | 0.779 |
| $72^2$ | 100 | **72** | 69 | 8.5 | 0.828 |

(the $72^2$ rows are the tail-inclusive repeat; the ring damps the current, so dropping
the tail moves $T$ by 0.0009 there against 0.025 in a coherent run.) A second,
independent route confirms the reversal: inverting $\sigma(\omega)$ bin by bin,
$\tau=-\mathrm{Im}\,\sigma/(\omega\,\mathrm{Re}\,\sigma)$, with no time-domain fit and no
driven window, gives $21.5$ fs at 10 kV/cm and $58.0$ fs at 100 — the same rise, from
data the fit never touches.

Two things are wrong with $\tau$ at these mesh densities, and they must both be said.
*Its magnitude is not converged*: the same field on a mesh with 2.25× the points gives
a $\tau$ four times shorter, and the e-ph energy-exchange rate of the channel ledger
doubles ($4.5\to9.2\times10^{-6}$ eV/cell/fs). *Its field dependence is not converged
either, and it reverses sign*: at $48^2$ $\tau$ falls by 2.3× from 10 to 100 kV/cm, at
$72^2$ it rises by 2.2× over the same interval. An earlier version of this page read the
$48^2$ fall as the optical-phonon thresholds switching on; the finer mesh gives the
opposite trend, so **that reading is withdrawn** — it was mesh noise, not threshold
physics. The cause is the same in both cases: the inter-$k$ golden-rule sum is
unresolved, the energy-matching window $\sigma_E=0.1$ eV containing too few final
states, so the rate samples the mesh rather than the phonon spectrum.

What *does* survive mesh refinement is the sign of every effect at fixed mesh when the
ring is switched on: scattering shortens $\tau$ against the collisionless run, converts
reflection into absorption, lowers $\mathrm{Re}\,\sigma$ and raises $T$.

**Why the dissipative curve cannot simply be moved to the converged mesh.** The
coherent scan can: it is $O(N_k)$, and $147^2$ costs minutes per field, which is how
the mesh series above was made. The ring is $O(N_k^2)$. Measured at $72^2$ on four
threads a dissipative field takes $\approx1.4$ h; at $111^2$, where the Fermi disc of
$E_F=0.6$ eV finally holds $3.6$ spacings, the same field is $2.38^2\times1.35\approx7.6$
times that, i.e. $\approx10$ h, and a six-point scan is close to three days. The two
halves of the comparison are therefore *not* on the same mesh, and that has to be
stated wherever they appear together: the coherent runs give the shape of $T(E_0)$ on
a mesh where the spurious bump is down to $+8.6\,\%$ and the dip in $T$ to $3.5\,\%$,
the dissipative runs give its character and absolute level on a mesh where the bump is
still $+19\,\%$.

One tempting shortcut does not work. The criterion Eq. (4a.3) can also be met by
*raising the doping* on the same mesh — $E_F=0.8$ eV puts $k_F/\Delta k$ at $3.09$ on
$72^2$ at no extra cost. But the doping is not free to choose: $E_F=0.6$ eV is what
puts the sheet conductance near the measured $24.7\,\sigma_{\rm univ}$ (§4a.0), and
$0.8$ eV would take it to $\approx35$, so the run would be better resolved and further
from the sample. Resolving the Fermi surface and matching the measurement are two
different constraints, and at fixed cost only one of them can be satisfied.

Accordingly the comparison of the *absolute* transmission with the measured
$0.68\to0.79$ is indicative, not quantitative.

*Discretization.* The one caveat. Equation (4a.12) is a continuum statement; a mesh
represents it only as well as it fills the Fermi disc. Evaluating the same adiabatic
sum on real meshes:

![drift saturation of a doped Dirac sheet and its representation on a k-mesh](figures/graphene_drift_saturation.png)

| $u=A_0/k_F$ | 0.05 | 0.2 | 0.4 | 0.6 | 1.0 | 2.0 | 5.0 |
|---|---|---|---|---|---|---|---|
| continuum | 1.000 | 0.995 | 0.979 | 0.952 | 0.851 | 0.487 | 0.200 |
| $300^2$, $E_F=0.2$ eV (140 partial pts) | 0.930 | 1.062 | 1.015 | 0.978 | 0.878 | 0.543 | 0.249 |
| $147^2$, $E_F=0.2$ eV (36) | 6.2 | 1.225 | 1.000 | 0.984 | 0.829 | 0.532 | 0.243 |
| $48^2$, $E_F=0.6$ eV (12) | 0.901 | 1.015 | 1.171 | 1.062 | 1.003 | 0.635 | 0.333 |

A disc holding one shell of mesh points reproduces the linear limit to $\sim10\,\%$
but develops a spurious **bump of 15–20 % near $u\simeq0.2$–$0.5$**, and on the
coarsest Fermi surfaces the small-$u$ limit itself breaks down (the $147^2$ entry at
$u=0.05$ is the K-point shell dominating a nine-point disc). That bump is exactly the
*dip* in the transmission scan of §4a.3, and the two scale together — the decisive
test, because a physical darkening would not care about the mesh:

Measured on the runs themselves rather than on the adiabatic sum, the same effect is
the spurious peak the run's own $\mathrm{Re}\,\sigma$ develops before saturation,
against its low-field plateau. The two numbers are not the same quantity — the table
above is the chord current of a *fixed* occupation, the one below the band-averaged
dynamical response of the driven sheet, which also carries the reactive part — but
they scale together, and the second is what a field scan actually shows. How much of
it reaches $T$ then depends on the sheet impedance, which is why the same artifact
makes a bigger dent at heavier doping. The control parameter is $k_F/\Delta k$ alone,
**not** the doping:

| run | $k_F/\Delta k$ | $\sigma$ plateau | $\sigma$ peak | **spurious bump** | $T$ plateau | $T$ minimum | dip in $T$ |
|---|---|---|---|---|---|---|---|
| $48^2$, $E_F=0.6$ eV | 1.54 | 23.6 | 31.1 | **+32 %** | 0.7275 | 0.6414 | −11.8 % |
| $72^2$, $E_F=0.4$ eV | 1.54 | 12.9 | 17.4 | **+35 %** | 0.8487 | 0.7939 | −6.5 % |
| $72^2$, $E_F=0.6$ eV | 2.32 | 25.9 | 30.8 | **+19 %** | 0.7081 | 0.6571 | −7.2 % |
| $111^2$, $E_F=0.6$ eV | 3.57 | 25.6 | 27.8 | **+8.6 %** | 0.7149 | 0.6900 | −3.5 % |
| $147^2$, $E_F=0.6$ eV | 4.73 | 27.0 | 29.0 | **+7.6 %** | 0.6997 | 0.6763 | −3.3 % |

The two runs at $k_F/\Delta k=1.54$ are different meshes at different dopings and give
the same bump to within the sampling of the field grid; filling the disc from 1.5 to
3.6 spacings takes the bump from a third to under a tenth, and the dip in $T$ with it.

**But it does not go to zero.** Between $111^2$ and $147^2$ — a 78 % increase in
k-points, taking the disc from 3.6 to 4.7 spacings — the bump moves only from $+8.6\,\%$
to $+7.6\,\%$ and the dip stays at $3.3\,\%$. The mesh-driven part of the darkening is
therefore essentially gone by $k_F\simeq3.5\,\Delta k$, and what remains is a residual
few-per-cent dip near $u\simeq0.7$–$1.2$. It is *not* the $O(N_k)$ sampling of the Fermi
disc, and **it is not the basis either**: the same $147^2$ scan at `nstate = 4` gives a
dip of $3.41\,\%$ against $3.39\,\%$ at `nstate = 2`, unchanged to two digits, even
though the basis moves the absolute transmission by $0.035$ and $\mathrm{Re}\,\sigma$ by
$11\,\%$ (§4a.2). So the two numerical suspects are both excluded and the dip is left
unexplained; the remaining candidate this page can name is the departure of the real
EPM band from an ideal cone over the excursion $A_0\simeq k_F$, which the continuum
derivation of §4a.5.3 assumes away. An earlier version of this page said the converged
curve is monotone, on the strength of the $48^2\to72^2$ trend alone; the finer meshes
do not support that and the claim is withdrawn — see the end of this section. The dip in $T$ follows the bump and
its position tracks the bump rather than any field scale of the physics.

As the disc is filled further the bump shrinks — the adiabatic sum is down to 1.06 at
its worst on $300^2$ — and the curve converges onto Eq. (4a.13).

Two deviations must not be confused with it. At **large $u$** the mesh curves lie
systematically above the continuum ($0.31$ against $0.199$ at $u=5$): there the
displacement is a sizeable fraction of the Brillouin zone, the EPM band is no longer a
linear cone, and the difference is trigonal warping — physics of the real band
structure, not discretization. At **small $u$** the doped carriers carry their own
uncancelled diamagnetic remainder $\eta k_F/v_F$ (§4a.1), $\sim1\,\%$ at $E_F=0.2$ eV
and $\sim3\,\%$ at 0.6 eV, which is a basis effect and shrinks with `nstate`.

**What the mesh series does and does not establish.** It establishes that most of the
darkening before the rise is numerical: three quarters of it disappears between
$k_F/\Delta k=1.5$ and $3.6$, and the disappearance is governed by that ratio alone, not
by the doping. It does **not** establish that the converged curve is monotone. A dip of
$3.3\,\%$ near $u\simeq0.7$–$1.2$ survives at both $111^2$ and $147^2$ and at both
`nstate = 2` and `nstate = 4`, unchanged, and the transmission then rises exactly as
Eq. (4a.13) requires. So the honest reading is:
flat while $A_0\lesssim0.5\,k_F$, a residual few-per-cent darkening around
$A_0\simeq k_F$ of unsettled origin, then the drift-saturation rise, and finally
Landau–Zener darkening at several hundred kV/cm. A scan meant to be compared with
experiment quantitatively must still satisfy Eq. (4a.3), $k_F\gtrsim3\Delta k$, because
below it the darkening is several times larger and is an artifact. Reproduce the figure and the table with
`python3 drift_saturation.py GSDIR/graphene_sit --ef-ev 0.2`.

#### 4a.5.6 The bilayer with the ring on: the $T\cdot T$ error changes sign

Everything above about layers was coherent, where the doped sheet is an inductor. With
the ring on it becomes a Drude absorber, and Eq. (6a) then predicts the opposite sign.
The propagated pair confirms it. Same $72^2$ mesh, $E_F=0.6$ eV, ring at a 300 K
lattice, 100 fs ring-down, `sbe_sheet_nlayers` 1 and 2:

| $E_0$ [kV/cm] | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|
| $u=A_0/k_F$ | 0.12 | 0.37 | 0.74 | 1.24 | 3.71 |
| one layer, $T$ | 0.7786 | 0.7789 | 0.7978 | 0.8278 | **0.8949** |
| two layers, $T$ | (0.6714) | 0.6164 | 0.6552 | 0.6936 | **0.7969** |
| $\mathrm{Re}\,\sigma_1/\sigma_{\rm univ}$ | 10.7 | 11.2 | 10.1 | 8.5 | 5.2 |
| $\mathrm{Re}\,\sigma_2/\sigma_{\rm univ}$ | (23.2) | 20.7 | 17.5 | 15.1 | 9.6 |
| $T_1^2/T_2$ | (0.903) | **0.984** | **0.971** | **0.988** | 1.005 |
| $\|A-A_E\|/A$, two layers | **0.81** | 0.07 | 0.02 | 0.01 | 0.00 |

The 10 kV/cm column is **quarantined**, hence the brackets, and the figure omits it.
Two independent diagnostics condemn it and agree with each other. Its fluence
absorption and its electron-energy ledger disagree by a factor 5 ($A=0.144$ against
$A_E=0.027$) where every other field agrees to $\le7\,\%$; and the two-layer `dark_diss`
control puts the field-independent ring current at $12.6\,\%$ of its peak there, against
$5.6$, $4.0$, $2.4$ and $1.3\,\%$ at 30, 60, 100 and 300 kV/cm — the only field over the
$10\,\%$ line of §4a.5.5. Subtracting the dark run moves it by only $+0.006$ in $T$ and
does not repair the ledger mismatch, exactly as the "not repairable above 10 %" rule
says. On the four clean fields the dark correction is $\le0.004$.

![one against two layers with the phonon ring](figures/graphene_layers_ring.png)

**The sign flips.** Coherent, the same doping gave $T_1^2/T_2 = 1.15$ at every field
(§4a.3); with the ring the trustworthy points give $0.97$–$0.99$ — below unity at every
one of them. The effect is *smaller* than the quarantined 10 kV/cm point suggested
($-1$ to $-3\,\%$, not $-10\,\%$), but its sign is unambiguous and it is the sign that
matters. The reason is the phase of $z$, exactly as
Eq. (6a) says: the coherent sheet has $\mathrm{Re}\,z/|z| = 0.014$ — an inductor — and
the dissipative one $0.899$ — a Drude resistor. Nothing else changed. So a real CVD
sample, which is resistive, sits on the branch where $T_1^2$ **under**estimates $T_2$,
and inverting a measured bilayer as $\sqrt{T_2}$ makes the monolayer look brighter than
it is, understating its conductance. That is the opposite of the coherent conclusion,
and it is the one that applies to a measurement.

**Both stacks brighten.** Over the clean points (30 to 300 kV/cm) one layer rises
$+14.9\,\%$ and two layers $+29.3\,\%$, both measured from their minimum at 30 kV/cm. The bilayer gains more, for the reason Eq. (4a.17) gives: it starts at
twice the conductance, where the same fractional saturation of $\sigma$ buys more
transmission. $\sigma_2/\sigma_1$ stays near 2 throughout (2.18, 1.85, 1.78, 1.86), so
the layers really are adding conductance and not something else.

#### 4a.5.7 Testing the drift law itself: the Drude weight, not the transmission

Everything above compared $T(E_0)$ with the law. That is what an experiment sees, but it
is a poor test of Eq. (4a.13), for two reasons that have nothing to do with the physics
of the cone. The sheet boundary condition $T=|2/(2+z)|^2$ is nonlinear in $z$, so the
same fractional change in the sheet response reads differently depending on where the
sheet sits; and it depends on the *phase* of $z$, which the drift law says nothing
about — Eq. (4a.13) predicts how the response shrinks, not whether it is reactive or
resistive. A test that removes both is to compare the quantity the law actually
predicts: the low-frequency weight.

The route is the fit-free inversion of §4a.5.5. The Drude form
$\sigma(\omega)=(D/\pi)/(1/\tau + i\omega)$ inverts bin by bin with nothing adjustable,

$$\tau(\omega) = -\frac{\mathrm{Im}\,\sigma}{\omega\,\mathrm{Re}\,\sigma},
\qquad
D(\omega) = -\pi\,\frac{\omega^2 + 1/\tau^2}{\omega}\,\mathrm{Im}\,\sigma ,
\tag{4a.18}$$

so each run returns its own $D$ *and* its own $\tau$, and the two can be read against
the field separately. That separation is the point: a sheet whose conductance falls
because the disc has drifted and a sheet whose conductance falls because the carriers
scatter more look identical in $T(E_0)$, and are told apart at a glance by whether
$\tau$ moved. $72^2$, $E_F = 0.6$ eV, ring on, anchored at 30 kV/cm (the lowest field
its dark control clears):

| $E_0$ [kV/cm] | 30 | 60 | 100 | 300 |
|---|---|---|---|---|
| $u = A_0/k_F$ | 0.371 | 0.742 | 1.237 | 3.711 |
| $D$ [eV] | 0.2581 | 0.2309 | 0.2010 | 0.1283 |
| $\tau$ [fs] | 59.1 | 53.9 | 58.0 | 62.7 |
| $D(u)/D(u_0)$ | 1.000 | **0.894** | **0.779** | 0.497 |
| $G(u)/u$, normalised the same way | 1.000 | **0.941** | **0.750** | 0.273 |
| residual | — | $-5.0\,\%$ | $+3.8\,\%$ | $+82\,\%$ (off the cone) |

![one layer with the ring against the drift law](figures/graphene_drift_fit_1layer.png)

**The law holds where it is supposed to, and $\tau$ says why.** Over $0.37 \le u \le
1.24$ the measured weight follows the parameter-free curve to $5\,\%$ and $4\,\%$, while
$\tau$ moves by $\pm8\,\%$ about 58 fs with no trend — the fall is weight, not
scattering, which is exactly the claim Eq. (4a.13) makes and the claim $T(E_0)$ alone
cannot establish. Past $u = 2$ the residual runs away to $+82\,\%$, as §4a.5.4 requires:
there the excursion is no longer small against the zone, the band is warped and pairs
are being made, and $G(u)$ is being quoted outside the picture it was derived in. The
figure shades that region rather than hiding it.

The curve is also carried *down* through the two fields the dark control threw out,
which turns out to be worth doing. At 10 kV/cm the quarantined point lands on the law to
$0.4\,\%$; at 3 kV/cm it sits $2.8\,\%$ above it. So the contamination is not a uniform
offset that could be calibrated away — it is invisible in $T$ at $14\,\%$ dark current and
plain at $36\,\%$ — which is the same conclusion §4a.5.5 reached from the energy ledger,
and the reason the rule there is a threshold and not a correction.

Two things this does **not** show, and both matter.

*It is two independent points.* Four fields survive the dark control, one of them is
the anchor, and one is off the cone. The agreement is a real test but a thin one; the
$297^2$, $E_F=0.2$ eV set (§4a.2) is where it should be repeated, because there $A_0=k_F$
falls at 27 kV/cm and the whole range $0 < u < 2$ is reachable at fields the dark
control clears.

*Only the shape is tested; the magnitude is off by a factor of two.* The same inversion
puts the ring run's absolute weight at $D/D_{\rm eq} = 0.43$ — the ring removes more
than half the low-frequency weight of the doping it was given. Two candidates are
excluded by measurement: carrier heating cannot do it (at $n=3\times10^{13}$ cm$^{-2}$
the equilibrium weight only falls to $0.87\,D$ even at 3000 K, and rises again above
that), and neither can the field-independent dark current (subtracting the `dark_diss`
trajectory from every run changes $D$ by $0.2\,\%$ and $\tau$ by $\le1\,\%$). Momentum
relaxation towards the lattice-frame distribution should broaden the Drude peak, not
drain it, so this is a genuine open item about the dissipator rather than about the
cone — and §4a.5.8, which reads the per-channel ledger, names e-ph as the suspect and
measures what it has already done to the sheet before the pulse arrives. It does not touch the field dependence — every entry in the table above is a
ratio taken within one run — but it does mean the ring's *absolute* low-field
conductance should not be compared with a measurement until it is understood.

**The collisionless sheet does not obey the law**, which is the opposite of what the
derivation would suggest, since Eq. (4a.13) is itself collisionless. The same inversion
on the coherent runs, both meshes, anchored the same way:

| $u$ | 0.012 | 0.371 | 0.742 | 1.237 | 3.711 |
|---|---|---|---|---|---|
| $D$ [eV], $147^2$ | 0.6080 | 0.6313 | 0.6503 | 0.6440 | 0.4536 |
| $D(u)/D(u_0)$, $147^2$ | 0.963 | 1.000 | 1.030 | 1.020 | 0.718 |
| $D(u)/D(u_0)$, $72^2$ | 0.896 | 1.000 | 1.055 | 1.006 | 0.698 |
| $G(u)/u$ | 1.018 | 1.000 | 0.941 | 0.750 | 0.273 |
| residual, $147^2$ | $-5\,\%$ | — | $+10\,\%$ | $+36\,\%$ | $+163\,\%$ |

Its absolute weight is right — $D/D_{\rm eq} = 1.01$ at $147^2$, so the doped ground
state and the f-sum restoration are doing their job — but its field dependence is
wrong: $D$ *rises* to a maximum near $u\approx0.74$ and does not begin to fall until
$u>2$, where the law wants it down by a quarter already at $u=1.24$.

That rise is the conductivity bump of §4a.5.5, arriving here by a second and independent
route. Measured from each run's own lowest field, $D$ peaks $+17.7\,\%$ above it at $72^2$
and $+7.0\,\%$ at $147^2$ — against $+19\,\%$ and $+7.6\,\%$ for the bump in
$\mathrm{Re}\,\sigma$ on the same two meshes. Two quantities extracted in different ways
(a band-averaged $\mathrm{Re}\,\sigma$ referred to the incident field, and $D$ inverted
from $\sigma(\omega)$ bin by bin) agreeing to under a per cent on both meshes is a
strong check that the bump is a property of the solution and not of either diagnostic.
Refinement halves it and does not remove it, exactly as §4a.5.5 found: partly the
discretization of the Fermi disc, partly still unexplained.

The reading this suggests — and it is a hypothesis, not a result — is that $G(u)/u$ is a
*quasi-static chord* response: it assumes the occupied disc is at the displaced position
belonging to the instantaneous $A(t)$, with nothing left over. A collisionless run has
no mechanism to enforce that, and keeps its weight up; a run whose momentum relaxation
(58 fs) is short against the drive period (300 fs) is held near the quasi-static
distribution, and the geometry shows through. The test that would settle it is cheap and
not yet done: vary $\tau$ through the lattice temperature and see whether the agreement
tracks $\tau/T_{\rm drive}$ rather than the presence of the ring.

Reproduce the figure and both tables with `drift_fit_plot.py` (x14 README §7.14).

#### 4a.5.8 The channel ledger: what e-ph does, and the one thing it does wrong

The solver keeps a cumulative per-cell ledger of every ring channel in
`*_sbe_channels.data` — $dN$, the conduction-population change, and $dE$ [Ha], the
eigenvalue-weighted energy the *electrons gained*, accumulated in `ring_ledger` as
$\sum_k \varepsilon_{nk}\,\delta f_{nk}/N_k$. It has to be read together with
`*_sbe_rt_energy.data`, whose `Eall` is **not** $\mathrm{Tr}\,\rho H$:
`realtime_ssbe.f90` accumulates `energy += (E_tot . -J) volume dt`, the work the local
field does on the sheet. So

$$\Delta E_{\rm all} = W_{\rm field}
\qquad\text{(an integration identity)},
\qquad
E_{\rm elec} = W_{\rm field} + \sum_{\rm ch} \Delta E_{\rm ch}.
\tag{4a.19}$$

Exactly one channel has a bath on the other side. $-\Delta E_{\rm eph}$ is what the
carriers hand to the phonons; Auger and impact ionization redistribute energy *inside*
the electron gas, so their $\Delta E$ must come out near zero while their $\Delta N$ does
not — which is the check that they are doing what they claim. $72^2$, $E_F=0.6$ eV, ring
on, dark run subtracted, per unit cell (the doping is 0.01512 carriers/cell):

| $E_0$ [kV/cm] | 3 | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|---|
| $W_{\rm field}$ [meV] | 0.0017 | 0.0203 | 0.1682 | 0.6335 | 1.5119 | 8.5315 |
| to the lattice, $-\Delta E_{\rm eph}$ [meV] | 0.0014 | 0.0173 | 0.0885 | 0.2466 | 0.6423 | 5.2112 |
| … as % of $W$ | 82 | 85 | 53 | 39 | 43 | 61 |
| $E_{\rm elec}$ left [meV] | 0.0003 | 0.0030 | 0.0797 | 0.3870 | 0.8694 | 3.3220 |
| $\Delta N_{\rm eph}$ /cell | 5.9e−11 | 4.8e−10 | 1.2e−08 | 2.4e−07 | 9.5e−07 | 4.0e−06 |
| $\Delta N_{\rm rana}$ /cell | −1.3e−08 | −1.5e−07 | −1.3e−06 | −6.1e−06 | −5.1e−05 | −2.3e−03 |
| $\Delta E_{\rm rana}$ [meV] | −4e−09 | −5e−08 | −1e−06 | −2e−07 | −2e−04 | 1.6e−03 |

**e-ph is an energy sink, not a population source.** It carries 39–61 % of the absorbed
work to the lattice inside the 400 fs window while making almost no pairs: $\Delta
N_{\rm eph}$ tops out at $4\times10^{-6}$ per cell at 300 kV/cm, $0.03\,\%$ of the doped
carriers. **Rana is the mirror image** — $\Delta N_{\rm rana} = -2.3\times10^{-3}$ per
cell, $15\,\%$ of the carriers, with $\Delta E_{\rm rana}$ four orders of magnitude below
$W$. That is exactly what Auger should look like: it moves carriers, not energy, and the
near-zero $\Delta E$ column is how one checks it. (The 3 and 10 kV/cm percentages are
ratios of numbers at the $10^{-3}$ meV level and mean nothing.)

**The zero-field pathology is e-ph, and only e-ph.** In the `dark_diss` control, with no
drive at all,

| channel | $\Delta N$ /cell | $\Delta E$ [meV/cell] |
|---|---|---|
| e-ph | $+2.8\times10^{-10}$ | $\mathbf{+2.7489}$ |
| Rana | $-4.6\times10^{-10}$ | $+0.0000$ |
| impact ionization, ring Auger | 0 | 0 |

The e-ph channel puts 2.75 meV per cell into the electron gas with the field switched
off — **182 meV per doped carrier**, of which 136 meV is already in by the time the pulse
peaks near 150 fs. The whole work done by the 100 kV/cm pulse is 1.51 meV/cell, so the
spurious pump is $1.8\times$ the entire signal at that field. This is the energetic face
of the dark current of §4a.5.5, and it names the channel: Rana at zero field contributes
$+0.0000$ meV and $-4.6\times10^{-10}$ carriers, i.e. nothing.

Two consequences follow, and they point in opposite directions.

*A clean dissipative sheet is available now.* Switching the phonons off and keeping the
Rana Auger channel gives a ring whose zero-field ledger is empty at this mesh. That is
**not** what §4a.5.5 found for the `nfs == 0` gain bug, where Auger alone still grew the
current — but that is a different failure, on a mesh that cannot represent the doping at
all, and the two should not be conflated. Where the doping is resolved, the phonon rates
are the thing that leaks.

*It is the leading suspect for the missing weight of §4a.5.7.* The pulse arrives at a
sheet that has already taken 136 meV per carrier from a bath it is supposed to be in
equilibrium with, so $D_{\rm eq}$ evaluated for $f_{\rm FD}(E_F, 300\,{\rm K})$ is the
wrong reference and the measured $D/D_{\rm eq}=0.43$ is in part a statement about that
reference. The k-resolved snapshots (§4a.5.9) show the pumped state is
approximately a *hot Fermi-Dirac* at $T_e\approx2260$ K, $\mu\approx0.42$ eV, which by
§4a.0's heating table costs only $-12\,\%$ in $D$ — so heating accounts for about a
quarter of the deficit and the rest is not yet explained. The repair belongs in the
detailed balance of the e-ph rates, not in the sheet, the mesh, or the drift law; §4a.5.9
is what it turned out to be.

Reproduce the ledger with `channel_budget.py` (x14 README §7.15); always pass `--dark`,
because without it every column carries the offset above.

#### 4a.5.9 What the zero-field pumping was: detailed balance on the wrong energy

A relaxation channel coupled to a bath at $kT$ has exactly one fixed point, the
Fermi-Dirac distribution at that $kT$. The e-ph ring did not have it, and §4a.5.8's
2.75 meV/cell of field-free heating is what that costs.

The Gaussian energy matching has width $\sigma$ (`sbe_search_sigma_e_ev`, 0.1 eV in
production), so a source is connected to partners at $|\Delta E|$ anywhere within a few
$\sigma$ of $\hbar\omega_p$ — but the emission/absorption split was taken **once per
mode**, from $N_B(\hbar\omega_p)$. A pair that actually transfers $\delta$ was then
weighted by $e^{\hbar\omega_p/kT}$ instead of $e^{\delta/kT}$, so its upward rate was too
large by $e^{(\delta-\hbar\omega_p)/kT}$. Graphene is where that runs away: the appended
acoustic mode is $\hbar\omega_{\rm ac} = 5.39$ meV against $\sigma = 0.1$ eV and carries
$98.2\,\%$ of the channel weight, the two optical modes (196 and 160 meV) having
$N_B\approx10^{-3}$ and weights 0.003 and 0.015.

Held at an *exact* $f_{\rm FD}(300\,\rm K)$ on a Dirac-like spectrum, `eph_interk_dpop`
returns, in meV per step:

| $\sigma$ [eV] | before | after | before, undoped ($\mu = 0$) |
|---|---|---|---|
| 0.100 | $+1.0418$ | $+1.8\times10^{-5}$ | $+0.0579$ |
| 0.020 | $+0.0118$ | $+5\times10^{-7}$ | $+1.25\times10^{-4}$ |
| 0.005 | $+7.9\times10^{-5}$ | $\sim0$ | $+1\times10^{-6}$ |
| 0.001 | 0 | 0 | 0 |

The leak is governed by $\sigma$, falling about $90\times$ per five-fold narrowing. The
Fermi level *amplifies* it $18\times$ — there are carriers with a sharp edge to smear —
but does not cause it. And the bath is innocent: it supplies the right $N_B$; the code
applied it to the wrong energy.

The dark run shows exactly that shape. Over 350 fs with no field the conduction
occupation goes from $f(0.50\text{--}0.58\ \rm eV) = 0.90$ to $0.39$ and from
$f(0.8\text{--}1.0\ \rm eV) = 3.5\times10^{-5}$ to $0.079$, the valence band moves by
$\le10^{-5}$, and the particle number is conserved to $10^{-11}$ — a Fermi edge
diffusing outward at fixed density. The tail fits a hot Fermi-Dirac at $T_e\approx2260$
K, $\mu\approx0.42$ eV, and $k_BT_e = 0.195$ eV is of the order of $\sigma = 0.1$ eV, not
of the bath's 25.9 meV. Its $\mathrm{Tr}\,\rho H$ rises by $+2.722$ meV/cell against the
ledger's $+2.749$, which is also the check that `Eall` really is the accumulated work and
the ledger really is the electronic energy (§4a.5.8).

**The correction.** Evaluate the split at the realized transfer. With $x = |\delta|/kT$,

$$f_{\rm emit} = \frac{N_B+1}{2N_B+1} = \frac{1}{1+e^{-x}},
\qquad
f_{\rm abs} = \frac{N_B}{2N_B+1} = \frac{1}{1+e^{x}},
\tag{4a.20}$$

a logistic in the transferred energy — one exponential, exact at both ends ($\tfrac12 :
\tfrac12$ for a degenerate pair, $1:0$ for a transfer far above the bath), and reducing
to the old expression when $\delta = \hbar\omega_p$. A second violation is removed with
it: the collision prefactor $\nu(\varepsilon)$ was read from the *source* alone, giving
one pair two different rates for its two directions (a factor 2 across 0.4–0.8 eV at
$\varepsilon_0 = 0.8$ eV); it is now the geometric mean over the pair, with $\nu$ at the
sink hoisted into a table so the pair loop — which runs $N_k^2 n_b^2 n_{\rm ph}$ times a
step — keeps one exponential rather than three.

**Scope.** Gated on the material (`mp%auger_2d_rana`): on for the gapless 2D Dirac
materials, off for Si, GaAs and CdS, which take the historical branch verbatim and stay
bit-identical to the validations published with them.

*The runs this section used to ask for have now been done, and they do not support
leaving it there.* An earlier revision reasoned from the optical modes alone
($\sigma/\hbar\omega \approx 1$) and called the 3D case a small correction. Two things
were missed. Silicon also gets an **appended acoustic mode** — the same structure that
makes graphene run away — at $\hbar\omega_{\rm ac} = 13.25$ meV on a $9^3$ mesh,
carrying **23.4 %** of the channel weight; and the production search width for the 3D
materials is the grid-matched default $\sigma = 0.2$ eV, not the 0.1 eV of the sheet
runs. Together those give $\sigma/\hbar\omega = 3.2 \dots 20$ across the silicon table
— the graphene condition, if anything more strongly.

Measured the way the ring is actually fed (the dressed-reference measure clamps the
filled valence sea to zero, so the sources are conduction electrons only) on a gapped
spectrum held at an exact $f_{\rm FD}(300$ K$)$, in meV per step:

| | historical split | realized transfer |
|---|---|---|
| Si, $\sigma = 0.2$ eV | **$+4.09\times10^{-2}$** | $+9.1\times10^{-7}$ |

Both conserve trace to $10^{-21}$: CPTP is intact, it is the energy balance that is
violated, and the sign is heating — a carrier gas already at the bath temperature is
warmed by it. The silicon mode energies and weights above are the ones the run's own
banner prints; repeating the probe with a GaAs-like table (Fröhlich LO at 36 meV plus
the five intervalley modes near 29 meV, weights assumed) gives $+1.39\times10^{-8}$ Ha
against $+1.2\times10^{-15}$, the same verdict — it is the mode *energies* against
$\sigma$ that decide this, not the detail of the weights. So the violation is **not** a
small correction for Si and GaAs; enabling the gate there is a physics decision about
re-validating those materials, not a question of whether the correction matters.

**Narrowing $\sigma$ is not an alternative.** The leak is governed by the search width,
so the obvious cheaper remedy is to shrink it rather than touch the split. Measured on
the real silicon spectrum, it does not work. The quantity that matters is the leak per
unit of *useful* relaxation — $\Sigma E\,\delta f$ at an exact $f_{\rm FD}(300$ K$)$,
where the channel must do nothing, over the same sum at $f_{\rm FD}(600$ K$)$, where it
must cool the carriers into the lattice:

| $\sigma$ [eV] | 0.400 | 0.200 | 0.100 | 0.050 | 0.030 | 0.020 |
|---|---|---|---|---|---|---|
| $\|{\rm leak/rate}\|$, historical | 9.7e-2 | 2.0e-2 | 3.0e-2 | 3.0e-2 | 8.6e-3 | 3.9e-6 |
| $\|{\rm leak/rate}\|$, realized | 2.1e-4 | 8.7e-10 | 1.3e-9 | 1.5e-10 | 1.5e-11 | 8.1e-12 |

The historical ratio is **flat at 2–3 %** from 0.2 down to 0.05 eV: narrowing removes
leak and useful work in equal measure. It improves only at 0.02–0.03 eV, and there the
rate has fallen four orders ($-3.6\times10^{-10} \to -6.2\times10^{-14}$) — the channel
is switched off, not corrected. There is no window in which the violation is small and
the channel still relaxes. (Note the grid guard is not the binding constraint here: it
compares $\sigma$ to the mean spacing of the *band ladder* at one $k$, 1.9 eV on this
mesh, whereas the density of final states within 0.5 eV of the CBM gives a real matching
floor near 10 meV.) The realized-transfer split reaches $8.7\times10^{-10}$ at the
production width with the rate intact, which the width alone never does.

**A dark control does show it — but only on a fine enough mesh, which is a trap.** The
reasoning that says it cannot is seductive and wrong: `dressed_ref` makes the ring read
the excess-carrier measure, that measure is zero for every band at $A = 0$, valence
included, so an undoped gapped material looks inert by construction. A $5^3$ silicon
dark run agrees — every channel reads exactly $0$, `nex` stays at $0$ for the whole run.
That is not health. The ring needs a source with $f > {\tt occ\_eps} = 10^{-12}$, and
the residue of the dressed projection only crosses that threshold once the mesh is fine
enough. Refining it, with the drive switched off throughout:

| mesh | $k$-points | $n_{\rm elec}$ [cm$^{-3}$] | $n_{\rm hole}$ [cm$^{-3}$] |
|---|---|---|---|
| $5^3$ | 125 | 0 | 0 |
| $7^3$ | 343 | $8.31\times10^{10}$ | $8.55\times10^{10}$ |
| $9^3$ | 729 | $1.32\times10^{12}$ | $1.32\times10^{12}$ |

(at $t \approx 31$ fs; the $9^3$ run is up to $2.21\times10^{13}$ by 50 fs and still
climbing linearly). Electrons and holes appear *together* to four
digits, so these are pairs promoted across the 1.07 eV gap — the Gaussian tail at
$5.3\sigma$, weighted for absorption of a 13 meV phonon. Meanwhile the trace holds at
8.000 and the current sits at $10^{-15}$: nothing else in the solver is disturbed. For
scale, $1.3\times10^{12}$ cm$^{-3}$ is some 100× silicon's intrinsic carrier density at
300 K, from a calculation with no field in it.

The lesson generalises past this bug: **a dark control run on a cheap mesh is not
evidence of health**, because the very threshold that makes it cheap also makes it
silent. Run the control on the production mesh.

*How to know these three rows are signal and not roundoff.* The columns of
`_sbe_nex.data` are $(\mathrm{tr}\,\rho - \mathrm{tr}_{\rm vb}\,\rho)/V$ and
$(n_{\rm elec} - \mathrm{tr}_{\rm vb}\,\rho)/V$, so their difference is the drift of
the total trace — the solver's own noise floor, free with every run. Here the $7^3$ row
misses by $2.4\times10^{9}$ cm$^{-3}$ against a signal of $8.6\times10^{10}$ (36×) and
the $9^3$ row agrees to three digits, so both are real pairs across the gap. That test
is not optional: in a *driven* 6600-step Si run the same floor reaches
$\sim10^{11}$ cm$^{-3}$ and swallows the whole signal. See `wiki/06`, addendum
2026-09-15.

**What the widened gate does to it.** The same $9^3$ control, same input, rebuilt with
`eph_db_realized = .true.`:

| | $n_{\rm elec}$ at 10 fs | at 20 fs | at 30 fs | ledger, whole run |
|---|---|---|---|---|
| historical split | $4.39\times10^{11}$ | $8.78\times10^{11}$ | $1.32\times10^{12}$ | $\delta N_{\rm eph} = 5.44\times10^{-11}$, $\delta E_{\rm eph} = +3.25\times10^{-12}$ Ha |
| realized transfer | 0 | 0 | 0 | identically 0, every channel, every row |

Not reduced — gone. $\sum|{\rm ledger}|$ over the entire run is exactly zero, so e-ph,
impact ionization and the rest all stand down when there is nothing to act on, which is
what a dissipator with no carriers and no field is supposed to do. (`nhole` in the
corrected run carries a static $-1.14\times10^{10}$ cm$^{-3}$ offset, bit-identical from
the first output to the last; it is reference bookkeeping, not dynamics, and it does not
grow.) The full test suite is 32/32 on the widened gate, and the 2D Dirac materials are
bit-identical since they already took this branch.

`tests/test_eph_detailed_balance.f90` holds a Dirac spectrum at an exact
$f_{\rm FD}(T_{\rm bath})$ and requires the channel to leave it alone; checks 5–6 repeat
the criterion on the gapped 3D spectrum with the silicon table at $\sigma = 0.2$ eV.

**The rescan.** Repeating the whole dissipative monolayer scan with the corrected channel —
same $72^2$, same $E_F$, same pulse and time grid as §4a.5.5, so the two are comparable
point by point — settles what the correction was worth.

The zero-field pump is gone rather than reduced: the ring current falls from
$|J|_{\max} = 2.4\times10^{-8}$ to $5.3\times10^{-19}$ and the energy it injects from
$+2.749$ to $+7\times10^{-5}$ meV/cell, which is $800\times$ below the 0.055 meV that
heating this sheet from 0 to 300 K costs. The dark fraction is $0.0\,\%$ at every field
against $35.8$ and $14.0\,\%$ at 3 and 10 kV/cm, and $\Delta E_{\rm eph}$ is now negative
everywhere — the carriers cool into the lattice instead of being warmed by it.

The drift law, anchored at 30 kV/cm as before:

| $u$ | 0.371 | 0.742 | 1.237 | 3.711 |
|---|---|---|---|---|
| $D(u)/D(u_0)$ | 1.000 | 0.957 | 0.770 | 0.461 |
| $G(u)/u$ | 1.000 | 0.941 | 0.750 | 0.273 |
| residual, corrected | — | $+1.7\,\%$ | $+2.7\,\%$ | $+68.8\,\%$ |
| residual, before | — | $-5.0\,\%$ | $+3.8\,\%$ | $+81.8\,\%$ |

The worst on-cone deviation drops from $5.0$ to $2.7\,\%$, and the absolute weight recovers
to $D/D_{\rm eq} = 0.69$ from $0.43$ — so more than half of the deficit §4a.5.7 could not
explain was the pump destroying the weight. The sheet brightens $+26\,\%$ from its minimum
(against $+14.9\,\%$), and the fluence absorption and the electron-energy ledger, which
disagreed by a factor 5 on the worst pre-fix point, now agree to $0.1$–$7\,\%$.

![one layer with the corrected ring against the drift law](figures/graphene_drift_fit_1layer.png)

**A second fault is now visible underneath, and it is not this one.** At 3 and 10 kV/cm the
runs pass the dark control and fail the energy ledger: $A = -0.538$ against $A_E = -0.906$
at 3 kV/cm, $T + R = 1.54$, and with the drive exactly zero over the last 60 fs the current
sits at its maximum and is still growing. The phonon bath is not the source — $\Delta
E_{\rm eph}$ is negative there too. Two suspects were tested and both cleared: switching the Rana Auger channel off
reproduces the numbers to five figures, and the self-consistent sheet field is on in the
run that behaves. **It is e-ph.** With `yn_sbe_eph = 'n'` at 3 kV/cm the absorption turns
positive ($A = +0.0095$ against $-0.538$), the balance closes ($T+R = 0.990$ against
$1.538$) and the current decays to $14\,\%$ of its peak instead of standing at its maximum
with the drive at zero. So the low-field gain lives in the same channel the
detailed-balance correction repaired — a second, smaller fault that the zero-field control
cannot see, because it wakes only once the field displaces the distribution. The next
discriminator is the parameter behind the first fault: repeat at
`sbe_search_sigma_e_ev = 0.005` and see whether the gain follows $\sigma$. The fields from 30 kV/cm up are unaffected,
which is why the drift-law test above stands; the two low points are drawn hollow and not
fitted. The existing CPTP
tests could not have caught this — trace was conserved, populations stayed in range,
transfers went to energy-matched partners. A dissipator can be a perfectly valid CPTP map
and still have the wrong fixed point, and only a stationarity test says which.

## 5. Sheet electrodynamics

### 5.1 Boundary condition
At normal incidence the tangential $E$ is continuous across a current sheet and $H$ jumps by the sheet current $J_s$; with $Z_0=4\pi/c$ and a substrate of index $n_s$ behind the sheet,
$$
E_t=\frac{2E_{\rm inc}-Z_0J_s}{1+n_s},\qquad E_r=E_t-E_{\rm inc},\qquad
J_s=-J_m L_z, \tag{4}
$$
where $J_m$ is the electron current per cell volume written by the solver (its energy ledger is $dW=-\mathbf E\!\cdot\!\mathbf J_m\,V\,dt$, so the charge current is $-J_m$). Fluence-integrated,
$$
T=n_s\frac{\int E_t^2}{\int E_{\rm inc}^2},\quad R=\frac{\int E_r^2}{\int E_{\rm inc}^2},\quad A=1-T-R,\qquad
\frac{c}{4\pi}\big(E_{\rm inc}^2-E_t^2-E_r^2\big)=E_tJ_s\ \ (n_s=1) \tag{5}
$$
pointwise: $A$ is the Joule absorption of the sheet in the *local* field. For the universal conductance $\sigma=e^2/4\hbar$: $T=(1+\pi/2c)^{-2}=0.97746$, $A=\pi\alpha/(1+\pi\alpha/2)^2=0.02241$ (`tests/test_sheet_transmission.py`).

### 5.2 Radiation reaction in the velocity gauge
If the solver is driven by $E_{\rm inc}$ alone, its absorption $A_E=\int E_{\rm inc}J_s/F$ exceeds the sheet's by exactly
$$
S_{rr}=\frac{Z_0}{2}\frac{\int J_s^2\,dt}{F_{\rm inc}},\qquad A=A_E-S_{rr}, \tag{6}
$$
which is $O((Z_0\sigma)^2/2)\approx3\times10^{-4}$ on a mesh that resolves the response but becomes comparable to $A$ when a few discrete near-resonant k-points carry a large reactive current, and is not small at all for a THz-driven plasma ($Z_0\sigma/2\sim0.1$). The single-cell driver therefore propagates in the **local** field (`yn_sbe_sheet_field`): with $E=-\dot A$,
$$
\frac{dA_{\rm ind}}{dt}=-\frac{2\pi}{c}L_z\,J_m(t),\qquad A_{\rm tot}=A_{\rm ext}+A_{\rm ind},\qquad E_{\rm tot}=E_{\rm ext}+\frac{2\pi}{c}L_zJ_m, \tag{7}
$$
integrated explicitly with the current of the previous step (lag error $O(\Delta t\,\omega\,Z_0\sigma/2)$). `Ac_tot/E_tot` in `*_sbe_rt.data` are then the transmitted field, `E_ext` the incident one, the energy ledger uses $E_{\rm tot}$, and $A_{\rm ind}$ is part of the checkpoint state. The post-processor `transmission.py` detects this mode, uses $E_{\rm tot}$ directly and reports the deviation from the boundary-condition reconstruction (Eq. 4) as a consistency number.

## 5a. The driving transient

Every field number in this page refers to the same waveform: the maintainer's DAST
single-cycle terahertz transient, read by the solver as `ae_shape1 = 'input'` (a table
of $t$, $A_x$, $A_y$, $A_z$) and rescaled by `make_inputs.py` to an exact peak $|E|$,
with the offset removed and 10 fs cos² windows at both ends so that $A(\pm\infty)=0$.

![the DAST single-cycle transient: field, vector potential and spectrum](figures/graphene_thz_drive.png)

*Left* — the waveform $E(t)=-\dot A(t)$, one clean cycle of 284 fs. *Middle* — the
vector potential, which in the velocity gauge **is** the displacement of the Fermi sea
in reciprocal space (§4a.5.1): its peak is
$$
A_0 = 6.213\times10^{-4}\ a_0^{-1}\ \text{per kV/cm},\qquad
A_0 = 0.0621\ a_0^{-1}\ \text{at 100 kV/cm} = 3.71\,k_F\ (E_F=0.2\ \text{eV}),
$$
drawn against $k_F$ so the saturation criterion $A_0=k_F$ can be read off the picture.
*Right* — the intensity spectrum, peaking at **14.5 meV (3.52 THz)** with a centroid at
15.8 meV and a FWHM band of 7.3–22.7 meV, which is the band `transmission.py`
integrates $\mathrm{Re}\,\sigma$ over.

A single-cycle transient has no carrier frequency worth the name: the band is as wide
as its centre. That is why the drive is measured by $A_0$ rather than by $E_0/\omega$ of
a nominal carrier — the two differ by a factor 1.6 here — and why the transmission is
reported as a fluence ratio rather than at a single frequency. The figure is produced
by `plot_drive.py`, which prints all of these numbers for any drive file.

## 6. Numerical realization

1. Ground state: in-SALMON EPM, 43 plane waves, $N\times N\times1$ half-shifted MP mesh; **K is on the mesh only for odd multiples of 3** ($(2i-N-1)/2N$ contains $2/3$ iff $N$ is odd) — 147 or 153 put the Dirac point on the mesh, 12/24/150 straddle it by half a spacing. With K on the mesh the two levels at K are exactly degenerate and the ground-state occupation there must be the **group average** (1, 1 per spin-summed level), not the integer filling (2, 0): the latter is LAPACK's arbitrary basis inside the degenerate pair, a broken-symmetry state with a velocity expectation of order $v_F$ that radiates a field-independent current $\sim2v_F/N_k$ per valley and, with the self-consistent sheet field, pins the local field to zero at low fields (x14 README §7.3). `gs_info_ssbe` now averages every degenerate partially filled group (inert for gapped materials).
2. Per time step: CF4 unitary half-steps around the Strang-split dissipators; one Houston pass + gather per ring step; channels in order e-ph (with the phonon-line memory), Coulomb (with the plasmon-line source filter and $T_e$); sheet field update by Eq. (7); outputs.
3. Cost: the graphene ring is e-ph + Rana only — **$O(N_k^2)$**, exponential-bound (no $O(N_k^3)$ impact-ionization kernel). Measured $\approx40$ ms/step at $24^2$ on 4 threads $\Rightarrow\approx1.3$ s/step at $147^2$ on 48 threads. The 285-fs single-cycle THz transient + 100 fs tail at $\Delta t=0.1$ fs is 3844 steps: $\approx1.4$ h per dissipative run, minutes per coherent run.
4. Mesh rules. *Near-IR*: the resonance shell $k_{\rm res}=\hbar\omega/2\hbar v_F$ needs $\gtrsim3$ mesh points per radius ($N\ge150$ at 0.8 eV). *THz*: the relevant scale is the Zener excursion $A_0=E_0/\omega$ ($0.062$ a.u. at 100 kV/cm, 3.36 THz) against the spacing $|\mathbf b|/N$ ($0.0106$ at $N=147$) and the Landau–Zener tube $k_\perp^{LZ}=\sqrt{E/\pi v_F}=3.7\times10^{-3}$ a.u.; $A_0\ge2\,|\mathbf b|/N$ holds for $E_0\gtrsim30$ kV/cm at $N=147$; below that the pair creation is analytic (Eq. 8) rather than mesh-resolved.
5. Time step: the CF4 exponential is exact per step; $\Delta t$ is set by the field interpolation, the dissipator splitting and — with many bands — the stiffness of the high bands in the S4 composition: 0.1 fs is clean for $n_b\le4$ (bands $\le39$ eV), $n_b\ge8$ (bands to 90 eV) needs 0.05 fs (§6a).

## 6a. The velocity-gauge f-sum rule, the diamagnetic sheet current and its parameter-free removal (2026-09-04)

**Observation.** The first self-consistent THz run (24², $n_b=2$, 100 kV/cm) gave $T=0.15$, $R=0.85$: the empty sheet behaved as a plasma mirror. The sheet current was perfectly anti-correlated with the vector potential, $\mathrm{corr}(J_s,A_{\rm tot})=-1.000$ and uncorrelated with $E$, i.e. purely reactive, $J_s\simeq-\eta\,(N_e/A_{2D})\,A_{\rm tot}$ with $\eta=0.30$.

**Origin.** In the velocity gauge the current $J=\mathrm{Tr}[(\mathbf p+\mathbf A)\rho]/V$ carries the diamagnetic term $\mathbf A N_e/V$ for *every* electron of the filled π band. A uniform static $\mathbf A$ is a pure gauge, so this term must be cancelled exactly by the paramagnetic (interband) response — which happens only if the basis is complete: per band $n$ and direction $a$,
$$
S_n^a(\mathbf k)=\sum_{m\ne n}\frac{2|p^a_{nm}(\mathbf k)|^2}{\varepsilon_m-\varepsilon_n}=1-\frac{\partial^2\varepsilon_n}{\partial k_a^2},\qquad \big\langle S_n^a\big\rangle_{\rm full\ band}=1, \tag{9}
$$
and with $n_b$ bands the occupied-state average $\langle S^{(n_b)}\rangle<1$. The uncancelled fraction $\eta_a=1-\langle S^{a,(n_b)}\rangle$ multiplies $A=E/\omega$: harmless in the near-IR, decisive at 3 THz where $A_0$ is 60 times larger for the same field. This is the `wiki/06` basis-sufficiency issue in its most violent form — a 2D sheet whose whole valence band responds as if it were free.

**It is a property of the truncated basis, not of the mesh, the time step or adiabaticity.** This was checked directly (24², $n_b=2$, 1 kV/cm, no sheet field): the SBE current reproduces the *exact adiabatic response of the same truncated Hamiltonian* $H_{\mathbf k}(A)=\varepsilon_{\mathbf k}+A\,p^x_{\mathbf k}$ (ground state of the 2×2 problem at the instantaneous $A(t)$, Hellmann–Feynman current) to $1\times10^{-4}$ at every instant, both at 1 and at 100 kV/cm; the in-phase coefficient is $\eta_x=0.2995$ for $\Delta t=0.1,\,0.05,\,0.02$ fs alike, equal to Eq. (9) evaluated on the ground-state data; and driving at 0.1 and 0.4 eV instead of 14 meV moves it to 0.2985 and 0.2808, exactly the dispersive value $1-\big\langle\sum_m 2|p_{nm}|^2\Delta\varepsilon/(\Delta\varepsilon^2-\omega^2)\big\rangle$ (0.2985, 0.2806). The solver integrates its truncated model correctly; the model's ground state carries a current a complete basis would not.

**Convergence with $n_b$** (24², DAST 100 kV/cm, coherent, self-consistent sheet, no correction; static $\eta$ from the ground-state data, `sumrule_check.py`):

| $n_b$ | $\eta_x$ static | $A_{\rm tot}/A_{\rm ext}$ | $T$ | $R$ | $A$ |
|---|---|---|---|---|---|
| 2 | 0.300 | 0.20 | 0.147 | 0.851 | 0.002 |
| 4 | 0.097 | 0.53 | 0.646 | 0.321 | 0.033 |
| 8 | 0.036 | 0.90 | 0.937 | 0.013 | 0.049 |
| 16 | 0.030 | 0.92 | — | — | — |

The residual mirror at $n_b=16$ ($R\approx(Z_0\sigma/2)^2\approx(9.6\,\eta)^2\approx8\%$ at 3.36 THz) is still larger than the physical signal, so basis enlargement alone does not converge fast enough in the THz regime; the missing weight sits in bands $\gtrsim10$ eV up.

**The linear static correction is not admissible.** A first remedy, $J_a\to J_a-\eta_a(N_e/V)A_a(t)$ with the static $\eta$, removes the artifact at small $A$ but *over*-corrects at 100 kV/cm: the adiabatic ground-state current of the truncated model is **non-linear** in $A$ once $A$ is comparable to the k-distance of the mesh points nearest to K ($\partial J_{gs}/\partial A$ falls from $0.2995\,N_e/V$ to $0.268$ at $A=0.03$ and $0.276$ at $A_0=0.062$ a.u. on the 24² mesh — the level repulsion of the near-K pairs grows as $A(p_{cc}-p_{vv})$ closes their gap). An over-corrected sheet has a *negative* kinetic inductance, $J_{\rm phys}=+\kappa A$; with the self-consistent field this mode is unstable — after the pulse $A_{\rm ind}$ and $J$ grow together and the ledger shows gain ($T=1.21$, $R=0.50$, $A=-0.71$ at $n_b=2$; the same runaway appeared at $n_b=8$ in the dissipative runs). A fitted coefficient also sits badly in a first-principles code. It was removed.

**Remedy adopted: pure-gauge restoration, no parameter.** `yn_sbe_vg_sumrule='y'` now subtracts, inside `calc_current_bloch`, the adiabatic ground-state current of the *same* truncated Hamiltonian at the instantaneous vector potential,
$$
J_a(t)\ \to\ J_a(t)-J^{\rm gs}_a\big(\mathbf A(t)\big),\qquad
J^{\rm gs}_a(\mathbf A)=\frac{1}{V}\sum_{\mathbf k}w_{\mathbf k}\sum_n f^{(0)}_{n}\,\big\langle\phi_{n\mathbf k}(\mathbf A)\big|\,p_a+A_a\,\big|\phi_{n\mathbf k}(\mathbf A)\big\rangle
=\frac{\partial}{\partial A_a}\frac1V\sum_{\mathbf k}w_{\mathbf k}\sum_n f^{(0)}_n\Big[E_{n\mathbf k}(\mathbf A)+\tfrac12A^2\Big], \tag{10}
$$
where $\phi_{n\mathbf k}(\mathbf A)$, $E_{n\mathbf k}(\mathbf A)$ are the eigenvectors and eigenvalues of $H_{\mathbf k}(\mathbf A)=\varepsilon_{\mathbf k}+\mathbf A\cdot\mathbf p_{\mathbf k}$ (the propagator's own coupling, including any nonlocal or coset correction) and $f^{(0)}_n$ the ground-state occupations in energy order (adiabatic continuation). This is the statement that a uniform $\mathbf A$ is a pure gauge, imposed within the truncated space: in a complete basis $E_{n\mathbf k}(\mathbf A)=\varepsilon_n(\mathbf k+\mathbf A)$ and the BZ sum of Eq. (10) vanishes identically for every $\mathbf A$, so the correction is zero there; in the truncated basis its linear limit is the static $\eta_aN_eA_a/V$ and beyond it the exact non-linear ground-state current. Because the k-trace of $\rho_{\mathbf k}$ is conserved, the artifact resides in the ground-state sum alone — the subtraction is exact for **any** population (real pairs, their Drude current and the polarization currents are untouched: $J-J^{\rm gs}=\sum_{\mathbf k}\sum_n\big(f_n-f^{(0)}_n\big)\langle p+A\rangle_n+$ coherences). Cost: one $n_b\times n_b$ diagonalization per k per step, beside the six of the S4 propagator. The energy ledger, the sheet self-field (Eq. 7) and the outputs all use the restored current. Unit test `tests/test_vg_sumrule.f90` (Hellmann–Feynman identity at $A=10^{-5},0.05,0.3$, linear limit $=\eta$, non-linear departure); `sumrule_check.py` recomputes Eq. (10) from the ground-state files and reports the residual A-projection of any run.

**With the restoration** (24², DAST, coherent, self-consistent sheet): the post-pulse field is stationary ($A_{\rm ind}\to$ const, $J\to0$), $A_{\rm tot}/A_{\rm ext}=1.000$, and the residual A-projection of the sheet current is $<10^{-4}$ at all $n_b$:

| $n_b$ | $E_0$ [kV/cm] | $T$ | $R$ | $A$ | $\eta_{\rm dyn}$ residual |
|---|---|---|---|---|---|
| 2 | 1 | 1.000000 | 6×10⁻⁸ | 0.0000 | 0.0000 |
| 2 | 100 | 0.999996 | 1.3×10⁻⁶ | 0.0000 | 0.0000 |
| 4 | 100 | 0.999996 | 1.4×10⁻⁶ | 0.0000 | 0.0000 |
| 8 | 100 | 0.972 (Δt = 0.1 fs: stiff-band leakage, see text) | 3.4×10⁻⁴ | 0.027 | 0.0001 |
| 8 | 100 | 0.999985 (Δt = 0.05 fs) | 1.4×10⁻⁶ | 1.4×10⁻⁵ | 0.0000 |

The 24² mesh (smallest π–π* gap 0.78 eV, K half a spacing off the mesh) has no interband channel at 14 meV and no k-point inside the Landau–Zener tube, so the physical answer at this resolution is $T\simeq1$, which $n_b=2,3,4$ give identically ($n_b=3$: $T=0.999996$ as well). The earlier $A=0.049$ at $n_b=8$ without restoration was part of the artifact, not absorption; the remaining 2.7 % of $n_b=8$ *with* restoration at $\Delta t=0.1$ fs is a **time-step effect of the stiff high bands**: the deposited energy sits in bands 4–8 (22–58 eV above π, $5.6\times10^{-6}$ per cell carrying the whole ledger), which a 14 meV field cannot populate; at $\Delta t=0.05$ fs the same run gives $T=0.999985$, ledger $1.1\times10^{-7}$ eV per cell (2400× less). The S4/CF4 exponential is exact per step, but the composition with a backward sub-step applied to levels 34–90 eV up (13.7 rad per step) leaks population; $n_b\le4$ ($\le39$ eV) is clean at 0.1 fs. **Production recipe: $n_b=2$ with the restoration** (the restoration makes $n_b=2,3,4$ agree to $10^{-6}$ in $T$, the THz physics lives on the cone, and the ring cost is lowest); $n_b\ge8$ only with $\Delta t\le0.05$ fs. The absorption physics (Eq. 8, §7) lives on the 147² mesh with K on it (§8 and x14 README §7).

**Remark for the experiment.** Real CVD samples *are* strongly absorbing at THz — but through the Drude conductance of doping-induced carriers ($\sigma_{dc}\sim20$–$50\,\sigma_{\rm univ}$, sheet transmission 0.5–0.7 already in the linear regime), not through the artifact above; reproducing that requires the FD$(E_F,T)$ initial state (§9). A measured transmission below 90 % "without subtracting the substrate" also contains the substrate's own Fresnel loss ($\approx11\%$ per face for $n\approx2$); the sheet contribution is the ratio to the bare-substrate reference, Eq. (4) with $n_s$. Two electronically decoupled layers (a large-angle twisted or incoherently stacked bilayer) at $d\ll\lambda$ sit in the same local field and add their sheet currents: `sbe_sheet_nlayers = 2` (Eq. 7 with $2L_zJ_m$). The naive estimate $T_2\approx T_1^2$ is then only second-order accurate, and **its sign depends on whether the sheet is resistive or reactive**. With $z=Z_0\sigma$ complex, Eq. (4) at $n_s=1$ gives
$$
\frac{T_1^2}{T_2}=\frac{16\,|1+z|^2}{|2+z|^4}=1+\tfrac12\big[(\mathrm{Im}\,z)^2-(\mathrm{Re}\,z)^2\big]+\mathcal O(z^3). \tag{6a}
$$
For a purely **dissipative** sheet ($z$ real) $T_1^2<T_2$: squaring *under*estimates the bilayer transmission — at the sample's $z=0.565$ by $-9.6\,\%$ ($T_1^2/T_2=0.904$). For a purely **reactive** (inductive) sheet ($z$ imaginary), the same $|z|$ gives $T_1^2/T_2=1.133$, a $+13\,\%$ *over*estimate. A doped THz sheet is a mixture of the two, so the error can have either sign and $|z|$ alone does not bound it; only $\lesssim1\,\%$ for $|z|\le0.15$ ($T_1\ge0.87$) is safe. Both estimates also ignore that each layer sees the field reduced by *both* currents, which is what makes the non-linear (field-dependent) part genuinely non-multiplicative.

## 6b. How large must `nstate` be?

`wiki/03` states the general rule — *keep `nstate` large for the basis and pay only
for the window you dissipate*, because a strong field pushes population up through the
high bands and brings it back, so velocity-gauge basis sufficiency has to be preserved
even where nothing is dissipated. x14 runs production at `nstate = 2`, which is a
departure, and this section says on what grounds and where the grounds run out. Two
different requirements are involved and they bind at opposite ends of the field range.

### 6b.1 The ponderomotive requirement, and why it is not the binding one here

In this solver the velocity-gauge Hamiltonian is built as
$H_{nm}(\mathbf A)=\varepsilon_n\delta_{nm}+\mathbf A\cdot\mathbf p_{nm}$
(`build_HVG`): **the $A^2/2$ term is not there**. It is band-uniform, so it is a global
phase and drops out of every commutator; the same cancellation is why the scattering
thresholds must not restore it (`wiki/00`, 2026-07-12). The ponderomotive energy
therefore never shifts a level in this code, and "the top level must lie above $U_p$"
is a statement about whether the basis can represent the *dressed* states, whose scale
is the off-diagonal $\mathbf A\cdot\mathbf p$. On the Dirac cone
$\mathbf p=v_F\boldsymbol\sigma$ exactly, so that scale is $v_FA_0$ — the same energy
the drift picture of §4a.5.3 uses, and larger than $U_p$ everywhere below
$A_0=2v_F=0.88$ a.u.

| $E_0$ [kV/cm] | 1 | 10 | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|---|---|
| $A_0$ [a.u.] | 0.0006 | 0.0062 | 0.0186 | 0.0621 | 0.1864 | 0.6213 |
| $U_p=A_0^2/2$ [eV] | $7\times10^{-6}$ | 0.001 | 0.005 | 0.053 | 0.473 | **5.25** |
| $v_FA_0$ [eV] | 0.007 | 0.074 | 0.223 | 0.742 | 2.23 | **7.42** |
| $A_0/|\mathbf b|$ | 0.0004 | 0.004 | 0.012 | 0.040 | 0.119 | **0.398** |

Against the `nstate = 2` basis of this exercise, whose two bands are the full π
($-7.79$ to $0$ eV) and π\* ($0$ to $+19.68$ eV) pairs over the whole zone, with the
third band starting at $+13.13$ eV above the Dirac point:

* **the ponderomotive criterion is satisfied at every field in the scan** — the basis
  top is $19.68$ eV against $U_p\le5.25$ eV, a factor 3.7 even at 1000 kV/cm;
* the *dressing* scale is the tighter one, $v_FA_0=7.42$ eV at 1000 kV/cm against the
  $13.13$ eV to the third band — a factor 1.8, i.e. **marginal**. At 300 kV/cm and
  below it is a factor 6 or more and the two-band basis is uncontroversial.

So the answer to "should `nstate` be raised until the top level clears $U_p$?" is that
on this system it already does, and by the tighter $v_FA_0$ measure the two-band basis
is comfortable to $\sim300$ kV/cm and marginal at 1000.

### 6b.2 The requirement that does bind: f-sum completeness of the doped current

The constraint that actually limits `nstate` here is field-*independent*. A truncated
basis captures only the fraction $S=\langle\sum_m2|p_{nm}|^2/\Delta\varepsilon\rangle$
of the f-sum strength — 0.70 for two bands — and the missing part appears as an
uncancelled diamagnetic current $\propto A$ (§6a). The pure-gauge restoration removes
it *exactly*, but only for the occupation it is evaluated with, which is the undoped
reference: for the intrinsic sheet $n_b=2,3,4$ then agree to $10^{-6}$ in $T$ (§6a).
The carriers **added** on top of that reference get no such subtraction, and they do
not converge as fast. Same mesh, same doping, 1 kV/cm — where $U_p=7\times10^{-6}$ eV
and no ponderomotive argument can apply at all:

| $147^2$, $E_F=0.6$ eV, 1 kV/cm | `nstate = 2` | `nstate = 4` |
|---|---|---|
| $T$ | 0.70004 | 0.73477 |
| $\mathrm{Re}\,\sigma/\sigma_{\rm univ}$ | 26.92 | 23.77 |
| $D_{\rm spec}$ [eV] (analytic $E_F=0.600$) | 0.608 | 0.548 |

An $11\,\%$ change in the sheet response from the basis alone, at a field whose vector
potential ($6\times10^{-4}$ a.u.) can mix nothing. The natural reading is the
second-order repulsion of the π\* level by the bands above it,
$\Delta E\propto A^2$ — precisely a Drude-weight shift (§4a.2). Consequently:

**The control that separates the two.** Repeating the *same* mesh, pulse and field
with the intrinsic filling isolates it completely:

| $72^2$, 1 kV/cm | $n_b=2$ | $n_b=4$ |
|---|---|---|
| intrinsic, $T$ | 0.999996 | 0.999996 |
| doped $E_F=0.6$ eV, $T$ | 0.70961 | 0.74452 |

The undoped sheet is basis-independent to six digits, exactly as §6a promises; the
doped one moves by 3.5 points. So this is not a general basis insufficiency — it is
specifically the carriers the doping adds on top of the reference the restoration
subtracts.

**How many bands the doped sheet needs.** A scan at $72^2$, $E_F=0.6$ eV, 1 kV/cm:

| `nstate` | 2 | 3 | 4 | 6 | 8 ($\Delta t=0.05$ fs) |
|---|---|---|---|---|---|
| $T$ | 0.7096 | 0.7443 | 0.7445 | 0.7298 | 0.7473 |
| $\mathrm{Re}\,\sigma/\sigma_{\rm univ}$ | **25.81** | 22.62 | 22.60 | 23.11 | 22.34 |
| $D_{\rm spec}$ [eV] | 0.586 | 0.526 | 0.525 | 0.528 | 0.521 |
| $D_{\rm spec}/E_F$ | **0.976** | 0.876 | 0.875 | 0.881 | 0.868 |

The whole error is in the $2\to3$ step. Everything from 3 to 8 agrees to $\pm1.7\,\%$ in
$\sigma$ and $\pm0.7\,\%$ in $D$, while `nstate = 2` sits $14\,\%$ above all of them:
it is the outlier, not a converged choice. **Production is therefore `nstate = 4`** —
within $1.2\,\%$ of the $n_b=8$ answer and the largest basis still clean at
$\Delta t=0.1$ fs (§6a: $n_b\ge8$ needs $0.05$ fs) — and `make_inputs.py` defaults to
it.

Two consequences for what has been quoted from $n_b=2$ runs:

* the converged Drude weight is $\approx0.87\,E_F$, not $E_F$. The apparent agreement
  with the analytic Dirac value at $n_b=2$ (§4a.2) was **coincidental**: near K the
  $2\times2$ EPM block *is* the Dirac Hamiltonian, so it returns the ideal-cone answer
  by construction, and adding the bands above it brings in the real band structure's
  second-order response.
* ratios, shapes and field dependences at fixed `nstate` are unaffected — the
  transmission dip of §4a.5.5 is $3.39\,\%$ at $n_b=2$ and $3.41\,\%$ at $n_b=4$ — so
  every conclusion of §4a.5 about the *shape* of $T(E_0)$ stands. Only absolute
  conductances change.

## 7. Expected physics (analytic anchors for the validation)

*Near-IR, 1–100 kV/cm.* Perturbative interband regime: $A_0\ll|\mathbf b|/N$; pulse area $\theta\approx A_0v_F\tau_p/2=0.25$ rad at 100 kV/cm (8 cycles), so the coherent bleaching of the resonant shell is $\Delta A/A\approx-\theta^2/12\approx-0.5\%$ and, equivalently, the Pauli factor of the shell within the pulse bandwidth, $n/N_{\rm shell}\approx3.6\times10^{10}/6\times10^{12}$, gives $\Delta A/A\approx-1\%$: the transmission changes by $\sim10^{-4}$ absolute — the sheet stays at the universal $T=0.9775$.

*THz (DAST, 3.36 THz), 1–100 kV/cm, intrinsic sheet.* The field sweeps $\mathbf k+\mathbf A(t)$ through the Dirac point; pair creation is the massless Landau–Zener/Schwinger process with $P(k_\perp)=\exp(-\pi v_Fk_\perp^2/E)$ [10,11], i.e. a rate per area
$$
\Gamma=\frac{g}{4\pi^2}\,\frac{E^{3/2}}{v_F^{1/2}}\quad(g=4), \tag{8}
$$
which for the single-cycle transient gives $n\approx(1/\pi^2)A_0\sqrt{E_{\rm peak}/v_F}\approx1.5\times10^{12}$ cm⁻² per passage at 100 kV/cm (two passages per cycle, Stückelberg interference neglected), scaling as $E^{3/2}/\omega$: $\sim5\times10^{10}$ at 10 kV/cm. The created carriers absorb: coherently only their creation energy ($\sim2v_Fk_\perp^{LZ}\approx0.1$ eV per pair, a few per cent of the pulse energy at 100 kV/cm), with phonon scattering also the intraband (Drude) energy they acquire while accelerated to $v_FA_0\approx0.7$ eV. **For the intrinsic sheet the model therefore predicts THz-induced *absorption* growing with the field** (as observed in undoped graphene [12]); the self-induced *transparency* of doped CVD graphene [13–15] is the Drude-weight reduction of pre-existing carriers by heating and requires a doped / finite-temperature initial occupation, which the present solver does not have (§9).

*Size of the coherent effect.* In the pure Landau–Zener picture the pairs are created at the tube energy $2v_Fk_\perp^{LZ}=2\sqrt{v_FE/\pi}$, so the energy taken per passage is $\Gamma\tau\cdot2\sqrt{v_FE/\pi}\propto E^2$ — the same scaling as the fluence: the coherent absorbed *fraction* is field-independent, $A_{LZ}\approx5\,\%$ for the DAST transient at any amplitude (28 meV per pair at 10 kV/cm, 90 meV at 100 kV/cm; two passages), while the *number* of pairs scales as $E^{3/2}$. The transmission of the intrinsic coherent sheet is therefore expected to stay within a few per cent of unity over 1–100 kV/cm, and its field dependence is set by how the mesh resolves the tube: below $\sim30$ kV/cm the tube is thinner than one mesh cell and the mesh gives 0 (K averaged) to 0.7 % (K off, the driven near-K points), at 100 kV/cm the resolved tube gives 1.0 % (147², K averaged), 2.1 % (147², y polarisation), 3.7 % (150²) and 2.7 % (300²) — the x14 README §7.3 table; the spread is the transverse sampling of the strip (0.7 cells wide at 147², 1.4 at 300²). Dissipation adds the Drude heating of the created carriers (the `diss`/`mem` production runs).

## 8. Validation

| item | test / run | result |
|---|---|---|
| Dirac levels | `test_graphene_dirac_levels` | 43 PW: gap $6.5\times10^{-6}$ eV, $v_F=0.960\times10^6$ m/s, e–h asymmetry 0.8 %, linear to 0.5 eV; 7 PW: 0.2125 eV spurious gap (Fortran bandpath: 0.2125 / 0.0000 eV) |
| Rana balance | `test_rana_saturation` | $n_0=n_i(T)$ to $10^{-4}$, $T^2$ law, two-sided monotone CPTP saturation |
| plasmon line + filter | `test_colmem_2d` | limits, fixed point $10^{-12}$, $|R(2\omega)|$ transmission, Rana overrides |
| two-temperature fit | `test_dirac_te_fit` | mesh moments 0.5 %, $T_e$ to $10^{-4}$, degenerate limit, fallbacks |
| sheet BC | `test_sheet_transmission` | universal sheet $T=0.97746$, $A=0.02241$, energy identity, Fresnel |
| pipeline, 24² | x14 smoke (17 runs) | electrons = 2.000 every step; 9.5×10¹⁰ pairs/cm² at 100 kV/cm ⇔ ledger 5.7 %; Rana sign flips at $n_i$; Markovian ring +5 % "dephasing ionization", removed by the memory analog |
| VG f-sum rule / pure gauge | `test_vg_sumrule`; §6a tables | SBE = exact adiabatic response of the truncated $H(A)$ to $10^{-4}$ ($\Delta t$, $\omega$ checks); plasma mirror at $n_b=2$; linear static correction over-corrects and runs away with the sheet field; Eq. (10) restoration: residual $<10^{-4}$, sheet stationary, $A_{\rm tot}/A_{\rm ext}=1.000$ |
| calc-level THz, 147²/150²/300², $n_b=2$, pure gauge, coherent, sheet | x14 README §7.3 | $T$: 1.000 (1–10 kV/cm) → 0.9995 (30) → 0.97 ± 0.01 (100: 0.987 on 147² K-averaged, 0.974 for the y polarisation, 0.956 on 150², 0.969 on 300²) — induced *absorption* of 1–4 % at 100 kV/cm from the Landau–Zener pairs ($4$–$11\times10^{11}$ cm⁻² after the pulse vs $1.5\times10^{12}$ per passage from Eq. 8; the spread is the transverse sampling of the 0.7–1.4-cell-wide LZ strip); reflection $\le0.7$ %; the pair Drude weight screens $A_{\rm tot}$ by 4–6 %; ledger = fluence deficit to 3 digits; sheet passive after the pulse |
| Dirac point on the mesh | x14 README §7.3 | integer filling at a degenerate K = arbitrary broken-symmetry state → field-independent current $2v_F/N_k$ per valley, relay screening at 1 kV/cm ($R>1$); group-average occupation restores a current-free, continuous reference |
| time step vs. stiff bands | x14 README §7.2 | $n_b=8$ at $\Delta t=0.1$ fs leaks 2.7 % of the fluence into the 22–58 eV bands; 0.05 fs: $10^{-5}$; $n_b\le4$ clean at 0.1 fs |
| polarisation / bilayer | x14 README §7.4 | x vs y at 147², 100 kV/cm: $A=1.0$ vs 2.1 % — the transverse sampling of the Landau–Zener strip (0.7 cells wide), not crystal anisotropy (C₆ᵥ: isotropic through χ⁽³⁾, warping at χ⁽⁵⁾ $\lesssim10^{-4}$); two decoupled layers (`sbe_sheet_nlayers=2`): intrinsic $T=0.9585$ vs $T_1^2=0.9736$; doped ($72^2$, $E_F=0.6$ eV, 3 kV/cm) $T_2=0.4362$ vs $T_1^2=0.5031$ — $T\cdot T$ high by 15 %, the reactive branch of Eq. (6a); bin-by-bin $\sigma_2/\sigma_1=1.997$ |
| near-IR πα, ring with $T_e$ | x14 README §7.5–7.6 | 0.8 eV, 147²: $A=0.80\,\pi\alpha$ (3.2 points per shell radius), coherent bleaching $\Delta A/A=-0.4\,\%$ ($\theta^2/12$), $\Delta T=8\times10^{-5}$; 24² ring at 100 kV/cm: Markovian `diss` fabricates $A=0.12\,\%$ from the dressing, the memory analog (`mem`) gives the coherent $3\times10^{-6}$ (ring-visible density $10^5\times$ smaller); $T_e$ held at the lattice value while $n+p<10^{-3}n_i(T_L)$ |

## 9. Limitations and outlook

1. **Initial state at $T=0$, undoped.** `gs%occup` is integer filling; the thermal/doped Drude background of real samples and the bleaching physics of doped graphene need an FD$(E_F,T)$ initial occupation (which also generalizes the dressed-reference formula from full/empty bands to fractional ones).
2. **Mesh at THz.** Below $\sim30$ kV/cm the Landau–Zener tube is thinner than the $147^2$ spacing; the pair creation is then given by Eq. (8) analytically rather than by the mesh. A K-refined (non-uniform) mesh would need the momentum-map machinery of the acoustic mode to be generalized.
3. **Hot phonons.** The lattice is a fixed-temperature bath; the optical-phonon bottleneck of graphene (hot $A_1'$/$E_{2g}$ populations) is not included.
4. **Quasi-equilibrium form of the Coulomb sector.** The R07 rates and the two-temperature fit assume thermal branch distributions; the early, strongly non-thermal shell is mapped onto an effective $T_e$.
5. The sheet field is free-standing (vacuum both sides); a substrate index enters Eq. (4) trivially and can be added to the driver.

## References

1. S. Piscanec, M. Lazzeri, F. Mauri, A. C. Ferrari, J. Robertson, Phys. Rev. Lett. **93**, 185503 (2004).
2. M. Lazzeri, C. Attaccalite, L. Wirtz, F. Mauri, Phys. Rev. B **78**, 081406(R) (2008).
3. E. H. Hwang, S. Das Sarma, Phys. Rev. B **77**, 195412 (2008).
4. F. Rana, Phys. Rev. B **76**, 155431 (2007).
5. H. Haug, A.-P. Jauho, *Quantum Kinetics in Transport and Optics of Semiconductors* (Springer), build-up of screening.
6. E. H. Hwang, S. Das Sarma, Phys. Rev. B **75**, 205418 (2007).
7. L. A. Falkovsky, A. A. Varlamov, Eur. Phys. J. B **56**, 281 (2007).
8. T. Winzer, A. Knorr, E. Malic, Nano Lett. **10**, 4839 (2010).
9. D. Brida *et al.*, Nat. Commun. **4**, 1987 (2013).
10. D. Allor, T. D. Cohen, D. A. McGady, Phys. Rev. D **78**, 096009 (2008).
11. B. Dóra, R. Moessner, Phys. Rev. B **81**, 165431 (2010).
12. S. Tani, F. Blanchard, K. Tanaka, Phys. Rev. Lett. **109**, 166603 (2012).
13. H. Y. Hwang *et al.*, Phys. Rev. B **87**, 115413 (2013).
14. M. J. Paul *et al.*, New J. Phys. **15**, 085019 (2013).
15. Z. Mics *et al.*, Nat. Commun. **6**, 7655 (2015).
16. N. Boroumand *et al.*, Rep. Prog. Phys. **88**, 070501 (2025) — the non-Markovian framing (`wiki/10` §6).
17. R. Ramanujam, M.S. thesis, Arizona State University (2015) — the π-EPM form factors.
18. M. Hassanpour Amiri, J. Heidler, A. Hasnain, S. Anwar, H. Lu, K. Müllen, K. Asadi,
    *Doping free transfer of graphene using aqueous ammonia flow*, RSC Adv. **10**,
    1127–1131 (2020), doi:10.1039/C9RA06738H (open access) — residual ionic dopants of
    typical density $4\times10^{12}$ cm⁻² on transferred CVD graphene, removed by an
    ammonia wash; the reference point for the doping of §4a.0.
