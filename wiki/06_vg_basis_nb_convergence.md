# VG Basis Sufficiency & N_b Convergence

Standalone specification for checking that the number of bands carried into the velocity-gauge (VG) dynamics is sufficient. This is a **separate** correctness axis from the plane-wave cutoff and is **not** fixed by rotating into the Houston basis. Belongs in the long-term reference because the band budget must be re-verified for every new material and every new driver wavelength.

> **⚠️ Read §6 first (measured, 2026-07-23).** At sub-gap THz the real-carrier "over-generation" is a **time-step artifact**, not a band-count problem: with `dt` converged (≤ 0.05 fs for Si) the clean VG reproduces the < 10¹⁶ theory bound and is **flat from N_b ≈ 8**. Always converge `dt` on the **non-adiabatic** measure (`nex_proj`/`nex_dref`) before running an N_b study — an unconverged `dt` fakes a basis-insufficiency.

> **Implementation in this fork.**
> - Primitives (pure, unit-tested): `vg_eta_admixture`, `vg_trunc_shift2`, `vg_conv_error`, `vg_ptop_exceeds` in [`../src/ssbe/sbe_superres_ssbe.f90`](../src/ssbe/sbe_superres_ssbe.f90).
> - Test: [`../tests/test_vg_basis_nb.f90`](../tests/test_vg_basis_nb.f90) verifies the three criteria plus the Hylleraas-Undheim-MacDonald interlacing/upper-bound theorem and the 2nd-order truncation-shift formula (with a self-contained Jacobi eigensolver — no LAPACK).
> - Runtime hook (criterion (a)): the real-time solver computes `P_top = max_k ρ̃_{N_b,N_b}(k,t)` at the `out_projection_k_step` cadence and, if it exceeds `1e-3`, writes a WARNING to the **error channel** and **continues** the run (`src/ssbe/realtime_ssbe.f90`). Criteria (b)/(c) are operator procedures (run two N_b; a priori η estimate), not automated.

---

## 1. Two independent truncations (the core point)

There are two distinct cutoffs; conflating them is the trap.

**(A) Plane-wave truncation N_PW (the EPM cutoff, e.g. 11 Ry).** Sets the accuracy of the EPM bands E_n(k) and Bloch functions |u_nk> themselves. Static, diagonalized once. 11 Ry for GaAs (~59 G-vectors/spin) converges the low conduction and valence bands to ~meV. Not the subject of this test.

**(B) Band-count truncation N_b carried into the dynamics (e.g. 32 spinor bands).** The density matrix rho(k) is 2N_b x 2N_b, and the VG Hamiltonian

> H_VG = H_0(k) + A(t) * pi,   pi_mn = <m|p|n>

couples bands through the interband momentum matrix element pi_mn. The field admixes high bands into low ones. Truncating at N_b discards the couplings pi_{m,N_b+1}, pi_{m,N_b+2}, ... — and **no subsequent diagonalization of H_VG in the Houston basis can restore those discarded columns.** This test targets (B).

---

## 2. Why the Houston basis does NOT cure an insufficient N_b

The Houston (adiabatic) basis is U(t) from diagonalizing the **truncated** VG Hamiltonian H_VG^(N_b) = P_{N_b} H_VG P_{N_b}, where P_{N_b} projects onto the retained N_b bands:

> H_VG^(N_b) U = U diag(eps_a),   eps_a = eps_a[ P_{N_b} H_VG P_{N_b} ].

These are eigenvalues of the **projected** operator. By the Hylleraas-Undheim-MacDonald variational/interlacing theorem, the truncated eps_a are upper bounds to the true levels, and the truncation error of level a is

> delta eps_a = Sum_{c > N_b} |<a| A pi |c>|^2 / (eps_a - eps_c) + O(A^4)
>             = Sum_{c > N_b} A^2 |pi_ac|^2 / (eps_a - eps_c) + ...

So the Houston basis **inherits** the error of the truncated VG Hamiltonian — it diagonalizes an already-corrupted matrix. Changing representation (Bloch -> adiabatic) changes neither the retained subspace nor the discarded coupling. The missing band space is absent in both representations.

**What an insufficient N_b corrupts, concretely:**
- adiabatic energies eps_a(k,t) -> shifts the impact-ionization gate eps_kin = E_h(k+A) - E_CBM and all energy bins;
- group velocities V_a = (U^dag pi U)_aa + A -> distorts X_a and the Kuhn-Zurek dephasing;
- current j = Tr[(pi + A) rho] -> distorts the HHG spectrum, most sensitively the high harmonics (which live on the high bands). The HHG plateau cuts off where bands run out — the high-harmonic spectrum is the most sensitive indicator.

---

## 3. Three convergence criteria (increasing rigor)

### Criterion (a): top-band occupation P_top — cheapest, run every production job
Track the adiabatic population of the highest retained band:

> P_top(t) = max_k rho~_{N_b,N_b}(k,t)   (rho~ = U^dag rho U)

Practical threshold: **P_top < 1e-3** at all times. If P_top reaches ~1e-4 to 1e-3, the field is pushing population to the basis edge -> bands above the cutoff would also have been populated -> enlarge N_b. Necessary condition: a populated top band means the discarded bands above it mattered. *(Automated here: a WARNING on the error channel; the run continues.)*

### Criterion (b): N_b convergence study — the gold standard
Run with N_b and N_b + Delta (e.g. 32 and 48) and compare an observable O (current, HHG spectrum, carrier number):

> eps_conv = || O_{N_b+Delta} - O_{N_b} || / || O_{N_b+Delta} ||

Band-count convergence is the only rigorous check. For HHG, compare the spectrum **on a log scale up to the harmonic order of interest** — the plateau must coincide. This is mandatory validation, not optional.

### Criterion (c): adiabatic admixture parameter eta — a priori estimate
Dimensionless strength of admixing the first discarded band:

> eta_ac = A_max |pi_ac| / |eps_a - eps_c|

The basis is sufficient for level a if eta_{a,N_b} << 1 for the coupling to the **first discarded** band. When eta >~ 1 perturbation theory fails, admixture is nonperturbative, and the band cannot be discarded. This is the bridge to the field estimates below.

---

## 4. Worked field estimate: 1 MV/cm THz (GaAs)

**Field and A_max.** E_0 = 1 MV/cm = 1e8 V/m. With E_a.u. = 5.142e11 V/m:
> E_0 = 1e8 / 5.142e11 = 1.945e-4 a.u.

For a monochromatic drive A_max = E_0/omega. Chefonov THz, f = 1.5 THz:
> hbar*omega = 2 pi hbar f = 6.2e-3 eV = 2.28e-4 a.u.
> A_max = E_0/omega = 1.945e-4 / 2.28e-4 = 0.85 a.u.

**This A_max is of order the Brillouin-zone size.** GaAs 2pi/a = 2pi/10.68 = 0.588 a.u., so
> A_max / (2pi/a) ~ 0.85 / 0.588 ~ 1.45

— the field sweeps the electron across **more than a full BZ** per half cycle. This is why the adiabatic/Houston language is mandatory; but the band budget still must be checked.

**Per-level admixture.** GaAs interband momentum element from Kane's E_P = 2|pi_cv|^2/m ~ 25.7 eV -> |pi_cv| ~ 0.97 a.u. Characteristic interband coupling energy:
> A_max |pi_cv| ~ 0.85 * 0.97 ~ 0.82 a.u. = 22.4 eV.

Compare to gaps:
- **Low bands (gap, lowest conduction), gaps ~1-3 eV:** eta = 22.4/2 ~ 11 >> 1 -> strongly nonperturbative admixture. This does NOT mean "basis too small" — it means the lowest ~10-15 bands are strongly mixed and must ALL be retained (as they are). The adiabatic basis diagonalizes them correctly **if they are in the basis**.
- **Basis edge (N_b = 32 -> 16 orbital bands):** the 16th GaAs band sits ~20-30 eV above the valence-band top. The naive A_max|pi| ~ 22 eV looks alarming, but the **formal perturbative estimate is useless here because everything is nonperturbative**. Real band occupation is cut off **energetically**: the populated region extends to eps_kin ~ a few eV above CBM (Chefonov: "up to several eV"), not to 22 eV, because the oscillating THz field returns the packet.

**THz-specific regime caveat.** At hbar*omega = 6 meV << gap 1.42 eV the Keldysh parameter gamma_K = omega sqrt(2 m E_g)/(e E_0) is small -> **tunneling / quasi-static regime**. High-band admixture is then predominantly **virtual** (adiabatic following); real interband population is exponentially suppressed by the Zener factor. So despite A_max ~ 1.45 x BZ, the number of **really populated** bands is modest, and **N_b = 32 (16 orbital) is likely converged with margin for 1 MV/cm THz** — but this must be **confirmed** by a 32-vs-48 run, not assumed.

**Bottom line for 1 MV/cm THz:** the formal-perturbative estimate is uninformative (all nonperturbative); the only reliable criteria are **(b) N_b convergence + (a) P_top monitoring**. Physically the eps_kin ~ 5 eV cap corresponds to roughly GaAs bands #10-14 above the valence top, so 32 spinor bands should have headroom — verify, do not assume. For shorter-wavelength drivers (IR/visible: A_max = E_0/omega is smaller at fixed E_0, but the photon is larger and the energy reach is higher) the band budget must be re-checked separately.

---

## 5. Procedure (what to actually do)

1. **Every production run:** log P_top(t) = max_k rho~_{N_b,N_b}(k,t). Flag if it exceeds 1e-3. *(Done automatically: a warning on the error channel, run continues.)*
2. **Once per (material, driver) setup:** run N_b vs N_b+16 (e.g. 32 vs 48), compute eps_conv on the current/HHG spectrum/carrier number; require the HHG plateau to coincide on a log scale to the harmonic order of interest. Treat a failed match as "increase N_b and repeat."
3. **A priori sizing:** estimate eta_{a,N_b} = A_max |pi| / gap-to-first-discarded-band; if >~ 1 near the top of the energetically populated region, add bands before running.
4. **Re-verify on any change** of material or driver wavelength — the converged N_b does not transfer.

---

## 6. MEASURED CASE STUDY — the sub-gap-THz "over-generation" is a **dt artifact**, not a band-count problem (2026-07-23)

A clean-velocity-gauge convergence study on **Si primitive, driven by the maintainer's
DAST optical-rectification THz transient** (peak E ≈ 100 kV/cm, 3.3 THz, Keldysh
γ_K ≈ 5.7 — deep sub-gap, where theory expects real carriers **< 10¹⁶ cm⁻³**),
`5×5×5`, no dissipation. This settled a live puzzle: the frozen-window VG appeared to
over-generate real carriers by ~10³ at this working point. **It does not — the coherent
kernel reproduces the theory bound once `dt` is converged.**

**Method note.** The intended knob `nstate_sbe < nstate` for shrinking the coherent basis
is **currently broken** (heap overflow: `dt_evolve_bloch_cf4` copies the `(nstate,nstate)`
`gs%p_tm_matrix` into a `(nstate_sbe,nstate_sbe)` buffer — see the decisions log). So each
N_b point here is a **separate ground state** regenerated at that band count (the EPM's
lowest-N_b eigenpairs are truncation-invariant, so this is exactly criterion (b)).

### (i) N_b convergence is meaningless at an unconverged dt

![N_b convergence, dt=0.25 vs dt=0.05](figures/vg_nb_convergence_dt.png)

At **dt = 0.25 fs** the non-adiabatic real-carrier density `nex_proj` **climbs** with N_b
(×2 → ×80 theory) and "converges" to a **wrong, dt-inflated** value — the classic symptom
that would send you to add ever more bands (the §4 THz worry). At **dt = 0.05 fs** it is
**flat from N_b ≈ 8** (×1.4) — the added bands are not needed. The dt-error was filling
each newly-available band, *faking* a basis-insufficiency.

### (ii) Refine dt at fixed N_b: the real measure collapses, the dressing does not

| dt [fs] | `nex_proj` (real, non-adiabatic) | × theory | diabatic `nelec` (dressing) |
|---|---|---|---|
| 0.25 | 5.5×10¹⁷ | **55** | 2.554×10²¹ |
| 0.10 | 1.5×10¹⁶ | 1.5 | 2.556×10²¹ |
| 0.05 | 1.3×10¹⁶ | 1.3 | 2.555×10²¹ |
| 0.02 | 1.4×10¹⁶ | 1.4 | 2.555×10²¹ |

![dt convergence and time series](figures/vg_dt_convergence.png)

`nex_proj` falls **×40** and converges to ≈ 1.3–1.4×10¹⁶ ≈ the theory bound (residual ×1.3
is the coarse `5³` grid / finite N_b / genuine small multiphoton). The **diabatic `nelec`
is dt-flat to 4 digits** — it is the reversible A²(t) dressing (~2.5×10²¹), 700× larger and
dt-insensitive. The right panel shows the mechanism: at dt = 0.25 the spurious population
**accumulates over the pulse**; at dt ≤ 0.1 it tracks the theory bound.

### Three lessons (add to the criteria in §3)

1. **Judge convergence on the non-adiabatic real-carrier measure** (`nex_proj`/`nex_dref`
   in `_sbe_nex_nonad.data`), **never on the diabatic `nelec`/`nhole`** — the latter is
   dominated by the reversible dressing and is dt-flat, so it will falsely certify
   "dt-converged" while the real excitation is off by ×40.
2. **Criterion (b) N_b-convergence MUST be run at a converged `dt` first.** An
   under-resolved `dt` pumps the non-adiabatic sector into every band you add, so `nex`
   rises with N_b and mimics a band-budget problem. Converge `dt` on the real-carrier
   measure, *then* converge N_b.
3. **Time-step for the real measure:** for Si's ~14 eV band spread use **`dt ≤ 0.05 fs`**
   (`0.1` is ~10 % high; `0.25` is ×40). The CF4 stays **unitary** — electrons = 8.000 at
   every `dt` — so the failure is invisible in the trace and in the diabatic density; it is
   a phase-accuracy failure of the fast interband coherences (5.3 rad/step at dt = 0.25 fs
   × 14 eV). This is a *separate* axis from `dt` for the absorbed **energy** (wiki/11 §3c).

> **Reproduce:** `samples/exercise_x08_Si_primitive_hhg_basis` GS at several `nstate`, clean
> `&sbe`, a sub-gap THz field, and scan `dt` — read `_sbe_nex_nonad.data` col 2/3, not
> `_sbe_nex.data`.

### (iv) With dissipators ON, the over-generation is a SEPARATE, dt-DIVERGENT pathology

The clean-kernel `dt` cure above does **not** carry over to the dissipative run — the
opposite is true. Same field/material, `4³`, `nstate=16`, frozen window, full ring
(e-ph + acoustic + II + Auger), to 1000 fs:

| `dt` [fs] | `nex_proj` (diss ON) | cumulative ring-Auger [e⁻/cell] |
|---|---|---|
| 0.25 | 1.3×10²² (×10⁶) | −59 |
| 0.05 | 3.6×10²² (×10⁶) | **−24 900** |

![dissipator dt-divergence](figures/vg_dissipator_dt_divergence.png)

The green curve is the **clean VG at the same dt=0.05** (~10¹⁶ = physical). Turning the
dissipators on pumps `nex_proj` to **~10²² at any `dt`, and refining `dt` makes it worse**
(nex ~×2, Auger churn **×420**); electrons = 8.000 throughout (trace is conserved — this is
*not* a trace leak). **Mechanism:** the ring/frozen-sector CP-decoherence realifies the
reversible A²(t) dressing **per scattering event** (wiki/00, the "collision-assisted
generation" note), and that per-step realification is **not rate-normalized (∝ dt)** — so
more steps (smaller `dt`) ⇒ more spurious real carriers ⇒ more Auger (∝ n²).

**Consequence.** Two independent over-generation layers with *opposite* `dt` behaviour:
the **coherent** one is cured by `dt ≤ 0.05 fs`; the **dissipative** one is a separate bug
that `dt` cannot fix (it worsens it). The real fix must exclude the reversible dressing from
the collision source (the Option-A direction, `yn_sbe_dressed_ref`) **and/or** rate-normalize
the per-step realification so the ledger converges as `dt → 0`. Until then, sub-gap dissipative
absolute yields are unreliable regardless of `dt` — see the wiki/04 flag box.

---

## 7. Band budget by field strength (1 MV/cm and above)

Two independent reasons a stronger field needs more `nstate`, and both now bite
the **same** knob because the dressed projection is full-basis (`yn_sbe_full_dressed`):

1. **VG unitary sufficiency** (§2–3): the field shifts `k → k + A(t)`; the excursion
   is `A_max / (2π/a)` of the BZ. Population must not reach the top band (`P_top`).
2. **Dressed projection** (§6): the dissipators/measure diagonalise H_VG — with the
   fix that uses **all `nstate` bands**, so `nstate` must also span the states the
   field actually dresses. A too-small `nstate` now over-generates at the source,
   not just at the readout.

**The excursion scales linearly with E** (single-cycle THz: read `A_max` off the
field file directly). For the Si primitive (`2π/a = 0.612 a.u.`) driven by the
DAST 3.3 THz transient, scaling the measured `A_max = 0.070 a.u.` (at ~100 kV/cm):

| peak E | `A_max/BZ` | regime | γ_K (Si, E_g 3.34 eV) | **start `nstate`** (Si, 4 val) |
|---|---|---|---|---|
| ~100 kV/cm | 0.12 | perturbative, sub-cycle | ≈ 6 | 16–20 (x08 headroom) |
| **1 MV/cm** | **1.15** | sweeps a **full BZ**/half-cycle; ε_kin ~ few eV [Chefonov] | ≈ 0.6 (tunnelling onset) | **24–32** |
| 3 MV/cm | 3.4 | multi-BZ, hot tail toward X | < 0.3 (tunnelling) | 40–48 |
| ≥ 10 MV/cm | ≥ 11 | strongly non-perturbative | ≪ 1 | 48–64+, re-verify hard |

These are **starting points, not converged values** — always confirm with the
procedure below. The numbers assume a THz driver (large `A_max = E/ω`); a mid-IR/
optical driver at the same E has a **smaller** `A_max` (higher ω) but reaches
**higher energy** per photon, so the band budget must be re-checked separately (§3).

**Procedure at high field (mandatory):**
1. **Converge `dt` first** (§6 lesson): on the **non-adiabatic** measure
   (`nex_proj`/`nex_dref`), never the diabatic `nelec`. At γ_K ≲ 1 the coherences
   are faster — expect `dt ≤ 0.02–0.05 fs` (tighter than the ~0.05 fs of the weak
   field; `wiki/11 §3c`).
2. **`P_top < 1e-3`** at all times (criterion a). At 1 MV/cm+ the warning fires
   readily — raise `nstate` until it clears.
3. **`nstate` convergence at the converged `dt`** (criterion b): run e.g. 24/32/48,
   require `nex_proj` to plateau (§6 shows it is *flat from N_b ≈ 8 only at weak
   field*; a strong field needs many more before the plateau).
4. Keep `yn_sbe_full_dressed = 'y'` — a narrowed frozen window truncates the
   projection and re-introduces the over-generation exactly where the strong field
   populates the high bands.

**Cost note:** with the full-basis projection the ring scales `O(nk²·nstate²)`, so
doubling `nstate` for a strong field is ~4× on the ring (`wiki/11 §3d`). Budget the
grid/`dt`/channel set accordingly, or use the **cost-preserving source mask**
(landed): keep the full dressed projection but restrict the dissipator *sources*
to the `[fermi+core, fermi+free]` window (`frozen_*_threshold_ev` +
`yn_sbe_full_dressed='y'`). It covers the e-ph **and** the 2-particle II/Auger/Rana
kernels. **The window width follows the same band budget as this section:** a
measured Si 5³/`nstate=40` II run under a strong THz drive kept `dN_ii` *identical*
to all-active with a **±10 eV** window (1.45× faster) but lost **−14.5 %** with a
too-tight **±6 eV** window that froze the `E−E_F ≳ 6 eV` bands the multi-BZ field
actually populates (`wiki/11 §3d` for the table). Widen the window for stronger
fields, exactly as you widen `nstate`.

---

## 8. References
- Variational upper-bound / interlacing of truncated eigenvalues: E. A. Hylleraas & B. Undheim, Z. Phys. 65, 759 (1930); J. K. L. MacDonald, Phys. Rev. 43, 830 (1933).
- Velocity-gauge band-coupling and the need for many bands / gauge care: M. S. Wismer & V. S. Yakovlev, Phys. Rev. B 97, 144302 (2018); L. Yue & M. B. Gaarde, J. Opt. Soc. Am. B 39, 535 (2022).
- Kane momentum matrix element / E_P for GaAs: E. O. Kane, J. Phys. Chem. Solids 1, 249 (1957); E_P ~ 25.7 eV, I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001).
- Keldysh parameter / tunneling-vs-multiphoton crossover: L. V. Keldysh, Sov. Phys. JETP 20, 1307 (1965).
- THz carrier energies reaching several eV in n-Si (populated-band reach): O. V. Chefonov et al., Phys. Rev. B 98, 165206 (2018).


---

## Addendum (2026-09-04): 2D sheets at THz — the uncancelled diamagnetic current

The most violent form of the basis-sufficiency problem: in the velocity gauge the
filled band's diamagnetic current $A N_e/V$ is cancelled by the interband
response only for a complete basis; with $n_b$ bands the fraction
$\eta = 1 - \langle\sum_{m} 2|p_{nm}|^2/(\varepsilon_m-\varepsilon_n)\rangle$ survives as a
reactive current $\propto A = E/\omega$. Graphene, $n_b = 2$: $\eta = 0.29$ — at 3 THz and
100 kV/cm the sheet reflects 85 % (a plasma mirror); $n_b = 8$: 0.036; 16: 0.030. The
solver prints $S,\eta$ at start-up. A linear static subtraction $-\eta N_e A/V$ was
tried first and **withdrawn**: the truncated ground-state current is non-linear in $A$
once $A$ reaches the k-distance to the band-touching points, the static form then
over-corrects, and an over-corrected (anti-inductive) sheet is dynamically unstable
under the self-consistent field. `yn_sbe_vg_sumrule='y'` now performs a **pure-gauge
restoration**: it subtracts the adiabatic ground-state current of the same truncated
$H_{\mathbf k}(\mathbf A(t))$ — identically zero in a complete basis, exact at every
$A$ and for any population, no adjustable quantity. Full account, tables and recipe:
[`wiki/12`](12_graphene_sheet_solver.md) §6a.

## Addendum (2026-09-07): the same effect in bulk silicon, and what it costs

The 2D sheet made this dramatic, but nothing about it is two-dimensional. Bulk Si on a
$9^3$ mesh, EPM, driven by the measured DAST THz transient — peak $|E_{\rm tot}|$ of
1043 kV/cm inside the 600 fs window these runs cover — dissipators
off so the only thing under test is the basis:

| `nstate` | captured strength $S$ | $\eta$ (banner) | $\eta$ from $J/A$ at low field |
|---|---|---|---|
| 8 | 0.9023 | 9.77 % | 9.75 % |
| 12 | 0.9773 | 2.27 % | 2.27 % |
| 20 | 0.9974 | **0.26 %** | — |

The start-up banner and a direct measurement of $|J/A|/n_e$ in the leading tail of the
pulse agree to three digits, so the printed $\eta$ can be trusted as the real thing and
not just an estimate. Two consequences worth stating separately.

**What the residue costs.** At $n_b = 8$, the same run with and without
`yn_sbe_vg_sumrule='y'` over the full 600 fs:

| | sum rule off | sum rule on | |
|---|---|---|---|
| $\eta$ ($t < 150$ fs) | 9.77 % | 0.089 % | ÷110 |
| $\|J\|_{\max}$ [a.u.] | $1.02\times10^{-2}$ | $3.20\times10^{-5}$ | ÷319 |
| absorbed work [eV/cell] | $3.546\times10^{-1}$ | $1.939\times10^{-4}$ | **÷1828** |

So the absorbed energy — the quantity one actually plots against field — is 99.95 %
artifact at that basis size. This is the reason a THz absorption curve computed without
the sum rule cannot be compared to experiment at face value.

**It cannot be repaired afterwards.** The tempting shortcut is to subtract
$\eta n_e A$ from an existing $J(t)$ using the banner's $\eta$, and so rescue archived
runs without recomputing them. It does not work, and it does not become workable at a
larger basis where $\eta$ is small:

| $n_b$ | $\eta$ | raw | $-\eta n_e A$ post-hoc | sum rule (exact) | post-hoc miss |
|---|---|---|---|---|---|
| 8 | 9.77 % | $3.546\times10^{-1}$ | $8.647\times10^{-2}$ | $1.939\times10^{-4}$ | 446× |
| 12 | 2.27 % | $9.773\times10^{-2}$ | $3.544\times10^{-2}$ | $2.401\times10^{-4}$ | 148× |
| 20 | 0.26 % | $2.953\times10^{-2}$ | $2.240\times10^{-2}$ | $4.189\times10^{-4}$ | **53×** |

(eV/cell, full 600 fs). The reason is that $\eta$ is not a constant — it grows with the
field, because real carriers join the truncation residue. Along this pulse $n_b = 8$
runs 9.75 → 12.1 → 26.4 → 47.9 % at 60, 180, 300 and 420 fs, and even $n_b = 20$ runs
0.23 → 1.13 → 1.73 → 7.35 %. A single multiplier undercorrects exactly where the
absorbed energy accumulates, which is why the miss only falls from 446× to 53× as
$\eta$ falls 38-fold. `tests/test_vg_sumrule` shows the same from the other side: the
pure-gauge current gives $\eta = 0.740$ at $A = 0$ but $0.483$ at $A = 0.3$, against a
linear $\eta_{\rm lin} = 0.740$ throughout. Only the adiabatic subtraction tracks it,
and it has to run inside the propagation.

**The sum rule is basis-independent; the physics it uncovers is not.** With
`yn_sbe_vg_sumrule='y'` the residual $\eta$ is the same at every basis size — 0.0895,
0.0894 and 0.0894 % at $n_b$ = 8, 12, 20 in the leading tail, and $\le 0.15$ % anywhere
along the pulse — so the correction itself needs no convergence study. What *does* still
need one is the response underneath: the corrected absorbed work climbs
$1.94 \to 2.40 \to 4.19 \times 10^{-4}$ eV/cell over the same $n_b$, monotonically and
without flattening. At 1043 kV/cm the ponderomotive reach is large and each added
conduction band opens real absorption, so **$n_b = 20$ is not converged for silicon at
this field even with the sum rule on** — the two requirements are independent, and
satisfying one does not excuse the other. (For the band budget by field strength, §7.)

## Addendum (2026-09-08): the same scan for GaAs, which behaves differently

Silicon's corrected absorbed work was still climbing at the largest basis tried, so it
is worth recording that this is a property of the material and the field, not a general
verdict on the method. Repeating the scan for GaAs — its production ground state
(9³ mesh, $a = 10.683$ bohr, `epm_pw_cutoff_ry = 12`), the pure-gauge restoration on,
dissipators off, the same measured DAST transient:

| $n_b$ | $\eta$ (banner) | absorbed work [eV/cell] | step |
|---|---|---|---|
| 24 | 0.65 % | $2.0976\times10^{-4}$ | |
| 32 | 0.44 % | $2.2559\times10^{-4}$ | +7.5 % |
| 40 | 0.37 % | $2.3600\times10^{-4}$ | +4.6 % |

Each step is *smaller* than the last (ratio 0.61), where silicon's grew (+24 % then
+74 %, ratio 3.1). Summing the geometric tail puts $n_b = 40$ within about 7 % of the
converged value and the production $n_b = 32$ within about 12 %. The residual $\eta$
after the restoration is 0.0164–0.0167 % at all three, i.e. basis-independent, as it is
for silicon. **So GaAs needs the sum rule and nothing else; silicon needs the sum rule
*and* a larger basis.** The two must be decided per material and per field, not once.

*Scope.* The absorbed-work column was measured on a $5^3$ mesh over a 491.7 fs window
(which contains the pulse peak at 441.8 fs); the $\eta$ column is the production $9^3$.
The 9³ runs at $n_b = 24$ and 32 were abandoned at 13 h and 42 h of wall time. Whether a
sequence converges should not depend on the $k$-mesh, but that is an assumption here
rather than something these runs measured.

Two other things worth knowing about GaAs, both from runs with no field in them. Its
$\eta$ **plateaus** in the band count (1.61, 0.65, 0.44, 0.37 % at $n_b$ = 16, 24, 32,
40) but responds to the plane-wave cutoff: at $n_b = 32$, raising `epm_pw_cutoff_ry`
from 12 to 18 takes $\eta$ from 0.44 % to 0.21 % while the direct gap holds at 1.423 eV
against the measured 1.42. The cutoff sizes only the ground-state basis, so this costs
nothing at propagation time. Past 18 Ry the captured strength overshoots ($S = 1.005$ at
21 Ry) and $\eta$ changes sign, so 18 is the useful setting rather than "as high as
possible". And its dark control is **silent both before and after** the e-ph gate fix,
unlike silicon's — with no carriers at zero field the ring has no sources, and GaAs's
1.42 eV gap puts an across-gap transfer at $7.1\sigma$ rather than silicon's $5.3\sigma$.
That is not a clean bill of health for GaAs: it means the dark control cannot test it,
and the thermal-gas criterion (wiki/12) is what does.

## Addendum (2026-09-15): read the plateau, and check it against the solver's own floor

The **absorbed-work** columns of the two 2026-09 addenda above, and the $n_b$ and mesh
steps quoted from them in `samples/exercise_x15_.../README.md`, report series that refuse
to settle. They share one measurement recipe, and that recipe is wrong twice over. Both
faults are identified below; neither is a convergence failure.

This does **not** touch §6. That study converged a *carrier density* (`nex_proj`, falling
×40 with `dt` to $1.3$–$1.4\times10^{16}$ cm$^{-3}$ against a theory bound of
$10^{16}$), which is five orders above the floor of Fault 2 and is not an integral over
the drive at all, so Fault 1 cannot reach it either. §6's conclusion — that the apparent
over-generation was a `dt` artifact and that $n_b \gtrsim 8$ suffices at 100 kV/cm —
stands as measured. The contrast is the point: **the same solver gave a clean,
convergent series the moment the observable was one it could actually resolve.**

### Fault 1: the absorbed work was read before the drive stopped

The work done on the crystal,

$$W(t) \;=\; -\int_0^{t}\! \mathbf{E}(t')\cdot\mathbf{J}(t')\,V_{\rm cell}\,{\rm d}t',$$

is the absorbed energy only once $\mathbf{E}$ has switched off. While the pulse is on,
$W(t)$ is dominated by the polarisation the field has *lent* the crystal and not yet
taken back — energy that is on its way out, not in. Those runs used the measured DAST
field file, whose support runs to 3274 fs, inside windows that ended at 129 fs and
491.7 fs. Both are mid-pulse. On the 129 fs window the instantaneous power
$-\mathbf{E}\cdot\mathbf{J}$ changed sign **1554 times** and $W$ was still rising at the
edge.

A difference of two such numbers is a difference of two quantities that are not yet
defined, so it has no reason to converge in anything. Refining `dt`, adding bands or
refining the mesh each shifts *where in its swing* the integral is truncated, which is
why the steps changed sign and size at random and why refining one knob appeared to make
another knob worse.

**Rule.** Use a drive with compact support, run past its end, and read the plateau. The
analytic `Acos2` pulse of `samples/exercise_x15_.../` has $\mathbf{E} \equiv 0$ for
$t > $ `tw1` exactly, to the last bit; $W(t)$ then goes flat and *is* the absorbed
energy. A field *file* generally does not have this property — check its support before
you trust any energy taken from it.

### Fault 2: the residue was below the propagator's own noise floor

`_sbe_nex.data` writes two columns that are equal when the propagator is exact:

| column | expression | |
|---|---|---|
| 2 | $(\mathrm{tr}\,\rho - \mathrm{tr}_{\rm vb}\,\rho)/V$ | conduction population |
| 3 | $(n_{\rm elec} - \mathrm{tr}_{\rm vb}\,\rho)/V$ | valence depletion |

Their **difference is $(\mathrm{tr}\,\rho - n_{\rm elec})/V$, the drift of the total
trace** — and it is free. It is the error bar on either column, so quoting column 2
alone hides it. A run whose excitation is smaller than this difference has measured
nothing.

Measured on Si, $5^3$, $n_b = 28$, `dt` = 0.05 fs, 6600 steps, 1000 kV/cm single cycle,
all dissipators off: the drift grows with $|A(t)|$, peaks at $4.56\times10^{-12}$
electrons per cell and freezes at $3.71\times10^{-12}$ (i.e. $4.6\times10^{-13}$ of
$n_{\rm elec} = 8$) the moment the field stops. That is the double-precision floor of a
matrix-exponential chain — $\approx 2000\,\varepsilon$ over 6600 steps — not a leak:
with the ring off every step is the exponential of an anti-Hermitian matrix and
conserves the trace exactly in exact arithmetic.

$$\boxed{\text{floor} \;\approx\; 1\times10^{11}\ {\rm cm^{-3}}\ \text{in carrier
density},\quad \approx 1\times10^{-12}\ {\rm eV/cell}\ \text{in } W_{\rm plateau}}$$

for a run of this length. Nothing below it is resolvable, at any mesh or band count.

### The mesh question, answered on the plateau

Si, single-cycle 1000 kV/cm, `tw` = 273 fs, run to 330 fs, coherent, sum rule on:

| mesh | $k$-points | $W_{\rm plateau}$ [eV/cell] | tail drift | $n_{\rm ex}$ post-pulse [cm$^{-3}$] | own floor [cm$^{-3}$] |
|---|---|---|---|---|---|
| $5^3$ | 125 | $1.7795\times10^{-12}$ | 0.00e+00 | $4.88\times10^{8}$ | $9.28\times10^{10}$ |
| $7^3$ | 343 | $5.8159\times10^{-13}$ | 0.00e+00 | $2.66\times10^{8}$ | $6.90\times10^{10}$ |

Both sit ~200× **below** their own floor. Refining the mesh does not raise the absorbed
energy; it lowers it, and lowers the floor with it. The $k$-resolved check agrees: at
$t = 300$ fs the diabatic conduction population summed over all 125 points of the $5^3$
mesh is $3.0\times10^{-12}$ electrons, max $1.7\times10^{-13}$ at any single point — the
trace drift and nothing else.

`_sbe_nex.data` is the *diabatic* measure and carries the reversible $A^2(t)$ dressing
(it peaks at $7.7\times10^{23}$ cm$^{-3}$ mid-pulse and comes back down), so the real
carriers must be read from `_sbe_nex_nonad.data` — column 3, `nex_dref`, is the
Option-A dressed-reference density **the ring dissipators actually see**. Post-pulse it
gives $1.21\times10^{9}$ cm$^{-3}$ at $5^3$ and $9.69\times10^{8}$ at $7^3$: the same
verdict, one order higher, still two orders under the floor.

**Independent check that this null is physics and not a dead solver.** At this working
point — $E_{\rm peak}$ = 1000 kV/cm, $\hbar\omega = \hbar\pi/t_w$ = 7.57 meV — the
Keldysh parameter is $\gamma_K \approx 0.13$–$0.20$, i.e. deep *tunnelling*, so the
Zener rate applies:
$G = \frac{e^2E^2\sqrt{m_r}}{18\pi\hbar^2\sqrt{E_g}}
\exp\!\big[-\tfrac{\pi\sqrt{m_r}E_g^{3/2}}{2\sqrt{2}\,e\hbar E}\big]$.
With $E_g$ = 1.07 eV this gives $\sim\!9\times10^{10}$ cm$^{-3}$ over the pulse for
$m_r = 0.2\,m_e$ and $\sim\!1\times10^{6}$ for $m_r = 0.5\,m_e$ — the exponent is
$\approx 20$–$32$, so the estimate is order-of-magnitude at best. But it brackets the
measured $10^{9}$ and, decisively, **it lands on the floor itself.** The signal this
scan was trying to converge is genuinely of the same size as the arithmetic noise. No
mesh and no band count can fix that; only a longer-lived observable or a stronger field
can.

*Noted in passing, not yet diagnosed:* the $7^3$ `nex_dref` trace carries a single
out-of-family sample at $t$ = 180 fs ($2.8\times10^{13}$ cm$^{-3}$, eight orders above
its neighbours at 160 and 200 fs). That column is exactly what the ring dissipators
read, so with dissipators **on** such a spike would inject real carriers. It does not
affect anything here (this scan is coherent) but it should be chased before the
production runs.

So **the mesh was never the problem** *at this band count*. In the coherent below-gap
regime silicon absorbs nothing measurable at this field, the two meshes agree on that
zero, and every series built on top of it was a series of ratios of noise. (At
$n_b$ = 36 the $5^3$ run does leave the floor — but only at `dt` = 0.05 fs, and Fault 3
below shows that is the step error, not the mesh.) Note that the formal step here is
$-67\,\%$ — as uncitable as the $+33\,\%$ it replaces, and for the same reason. **When
both endpoints are at the floor, quote the absolute values and the floor, never the
percentage.**

This floor is a property of the run length, not of the physics: with dissipators on, the
carrier densities the induced-transparency experiment is about are 6–7 orders above it.
The floor only ever obstructed what these scans were doing — comparing zeros.

### Fault 3: at 1000 kV/cm, `dt` = 0.05 fs invents carriers once the basis is large

The two rows above are both `nstate` = 28. Raising the basis to 36 at the *same* mesh
and the *same* step changes the answer by a factor of 7000 — and then halving the step
takes it all back:

| $n_b$ | `dt` = 0.05 fs | `dt` = 0.025 fs | step in `dt` |
|---|---|---|---|
| 28 | $1.7795\times10^{-12}$ | $1.6945\times10^{-12}$ | $-4.8\,\%$ |
| 36 | $\mathbf{1.2306\times10^{-8}}$ | $1.9897\times10^{-12}$ | $\mathbf{-99.98\,\%}$ |
| step in $n_b$ | $\mathbf{\times 6914}$ | $+17\,\%$ | |

(Si, $5^3$, coherent, $W_{\rm plateau}$ in eV/cell.) At `dt` = 0.05 fs and $n_b$ = 36 the
solver reports $4.53\times10^{12}$ cm$^{-3}$ of carriers; at 0.025 fs it reports
$8.4\times10^{8}$ with $n_{\rm hole}$ back to $-2.5\times10^{10}$, i.e. the floor that
$n_b$ = 28 was already on. **99.98 % of that population was the step error.**

This is §6's mechanism — "the dt-error was filling each newly-available band, *faking* a
basis-insufficiency" — at 10× the field §6 used, where it amplifies ×6185 instead of
×40. The step a run needs falls as the band ceiling and the field rise, so a `dt` that
was adequate at 100 kV/cm and $n_b$ = 8 is not adequate at 1000 kV/cm and $n_b$ = 36.

**Why the series looked divergent: the artifact needs BOTH knobs at once.** Three of the
four corners agree at $1.7$–$2.0\times10^{-12}$ eV/cell — the floor. Only $(n_b = 36,\,
{\rm d}t = 0.05)$ escapes it. So a scan along either edge of that table is misleading on
its own: sweeping $n_b$ at $\rm{d}t = 0.05$ walks into the bad corner and reports a
band-count dependence that never flattens (this is the withdrawn "~10 % per 8 bands"),
while sweeping `dt` at $n_b = 28$ walks along the safe edge and reports $-4.8\,\%$,
i.e. "`dt` is already converged". Both readings are real and both are wrong, because
each holds the *other* knob where the artifact is dormant. **One-knob-at-a-time
convergence testing cannot see this class of error.** At ${\rm d}t = 0.025$ the
band-count dependence is gone: $+17\,\%$ between two numbers that are both floor.

**The physics conclusion survives, now by two independent routes.** $(n_b = 28,\,
{\rm d}t = 0.05)$ and $(n_b = 36,\, {\rm d}t = 0.025)$ give $1.78\times10^{-12}$ and
$1.99\times10^{-12}$ eV/cell — the same floor, from opposite corners. Coherent Si at
1 MV/cm and 1.8 THz promotes nothing above $\sim10^{11}$ cm$^{-3}$, as the Zener
estimate said. And the `nstate` dependence that wrecked every earlier scan was never a
basis insufficiency; it was this.

### The e/h test is necessary but NOT sufficient — three checks, all required

The fake $4.5\times10^{12}$ cm$^{-3}$ came with $n_{\rm elec} = 4.527\times10^{12}$ and
$n_{\rm hole} = 4.509\times10^{12}$ — **balanced to three digits**, because a step error
drives a coherent valence→conduction transfer exactly as a real excitation does. Matched
pairs prove the transfer is not trace drift; they do not prove it is physical. Before
believing any residue:

1. $n_{\rm elec} \approx n_{\rm hole}$ — separates real transfer from trace drift;
2. the residue exceeds the trace drift (column 2 − column 3) — Fault 2;
3. **it survives halving `dt`** — Fault 3, and the only check that catches this one.

Criterion 3 is not optional and not expensive relative to being wrong by 6185×.

### Does this retract the dark-control mesh table (wiki/12)?

No, and the check is the one above. The signature of a noise-limited reading is
$n_{\rm elec} \neq n_{\rm hole}$. In that table the $7^3$ row reads $8.31$ vs
$8.55\times10^{10}$ cm$^{-3}$ — a mismatch 36× smaller than the signal — and the $9^3$
row agrees to three digits at $1.32\times10^{12}$. Those are real pairs across the gap.
The table stands.

### Procedure, amended

Add to §5, before anything else:

0. **Check the support of the drive.** If it is a file, find where it actually ends.
   Read energies only after it does, and confirm the tail drift of $W$ is zero.
1. **Take the floor from column 2 minus column 3** of `_sbe_nex.data`. If the residue
   you are about to converge is not several times that difference, stop: the series will
   be noise, and refining anything will make it look worse or better at random.
