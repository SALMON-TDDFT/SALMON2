# exercise x14 — graphene under the DAST THz transient (1–100 kV/cm, normal incidence): field before / field after, Zener pair creation, two-temperature Coulomb sector

**What this exercise delivers.** The cluster-ready study the maintainer asked for
before launching the graphene self-induced-transparency (SIT) runs: the **dense
production grid**, the **self-consistent sheet field** (so the solver's total field
*is* the transmitted field), the **DAST single-cycle THz drive** rescaled to
1–100 kV/cm, the **two-temperature model** of the Coulomb (Auger / carrier-
multiplication) sector with cooling through the phonon channel, the **level
check** that exposed a spurious 0.21 eV Dirac-point gap in the old basis, and
the **population-saturation check** (Auger and impact ionization balance at
n_i(T), Rana 2007). The methods write-up is
[`wiki/12_graphene_sheet_solver.md`](../../wiki/12_graphene_sheet_solver.md);
the collisional-memory design is `wiki/10` §8.11.

Geometry: THz beam at **normal incidence** on a free-standing monolayer, E in the
sheet plane (x); everything in `A_eV_fs` units.

---

## 1. Physics and the honest expectation (read before running)

**The drive.** `DAST_singlecycle_100kV.txt` (exercise x12): a Gaussian A(t), so
E = −dA/dt is one clean cycle at **3.36 THz** (ħω = 13.9 meV, period 298 fs),
285 fs long. `make_inputs.py` rescales it to the exact peak field (1, 3, 10, 30,
100 kV/cm), removes the initial offset and cos²-windows the ends (the SALMON
file reader returns 0 outside the file — an offset would be a δ-spike in E).
Peak k-excursion `A₀ = E₀/ω`: **0.062 a.u. at 100 kV/cm** = 6 spacings of the
147² mesh (0.0106 a.u.); 6×10⁻⁴ a.u. at 1 kV/cm.

**What the field does to the cone.** At THz the photon is far below any
interband transition the mesh can resolve; the physics is the field *sweeping*
k + A(t) through the Dirac point: massless **Landau–Zener / Schwinger pair
creation** with P(k⊥) = exp(−π v_F k⊥²/E) [Allor–Cohen–McGady 2008; Dóra–Moessner
2010], rate per area Γ = (g/4π²) E^{3/2}/v_F^{1/2} (g = 4) — see `wiki/12` §7:

| E₀ [kV/cm] | A₀ [a.u.] | A₀ / Δk (147²) | LZ tube k⊥ = √(E/πv_F) [a.u.] | pairs per passage (analytic) |
|---|---|---|---|---|
| 1 | 6.2×10⁻⁴ | 0.06 | 3.7×10⁻⁴ | ~1.5×10⁹ cm⁻² |
| 10 | 6.2×10⁻³ | 0.6 | 1.2×10⁻³ | ~5×10¹⁰ |
| 30 | 1.9×10⁻² | 1.8 | 2.0×10⁻³ | ~2.5×10¹¹ |
| 100 | 6.2×10⁻² | 5.9 | 3.7×10⁻³ | ~1.5×10¹² |

(two passages per cycle, Stückelberg interference neglected). The mesh
*resolves* this process when A₀ spans ≥ 2 spacings — **E₀ ≳ 30 kV/cm at 147²**;
below that the created density is the analytic one and the mesh only samples
the K point itself (147 = odd×3 ⇒ K exactly on the half-shifted MP mesh).

**Expected transmission range.** The sheet starts *empty* (T = 0 filling, no
thermal or doping carriers — §9), so at low field it is transparent at THz:
T → 1 (the interband σ_univ response sits at k = ħω/2ħv_F ≈ 6×10⁻⁴ a.u., far
inside one mesh cell; the universal 2.3 % is *not* available to the mesh at
3.4 THz). As the field grows, the created pairs absorb: coherently only their
creation energy (~2v_F k⊥ ≈ 0.1 eV per pair → a few per cent of the pulse energy
at 100 kV/cm), with phonon scattering also the intraband (Drude) energy they
acquire while accelerated to v_F A₀ ≈ 0.7 eV. **For the intrinsic sheet the
model predicts THz-induced *absorption* that grows with the field** — ΔT from
≈ 0 at 1–10 kV/cm to −(several %) coherent and up to −(tens of %) with
dissipation at 100 kV/cm — as observed in undoped graphene [Tani et al. 2012].
The self-induced *transparency* (bleaching) reported for doped CVD graphene
[Hwang 2013; Paul 2013; Mics 2015] is the Drude-weight reduction of
*pre-existing* carriers by heating and requires a doped / finite-temperature
initial occupation, which this solver does not have yet (§9). The numbers from
the local runs are in §7.

For comparison, the near-IR interband regime (0.8 eV, `--field acos2`) at the
same fields is perturbative: pulse area 0.25 rad at 100 kV/cm ⇒ ΔA/A ≈ −θ²/12
≈ −0.5 %, ΔT ~ 10⁻⁴.

## 2. Readiness matrix (graphene, in-plane field, 1–100 kV/cm)

| effect | status | note |
|---|---|---|
| EPM ground state (π/π* Dirac pair) | ✅ **fixed** | 43-PW basis (`epm_pw_cutoff_ry = 29.4`); the old 7-PW basis is gapped at K by 0.21 eV (§4) |
| coherent velocity-gauge SBE (CF4 + Houston) | ✅ | 2 bands; the VG "basis-edge" monitor fires for ANY excitation (top band = the conduction band) — expected |
| **2D-sheet self-consistent field** (`yn_sbe_sheet_field`) | ✅ **NEW** | radiation reaction in the single-cell driver: `E_tot` = transmitted field, ledger in the local field, checkpointed (§3) |
| e-ph optical E2g/A1′ + acoustic (inter-k ring) | ✅ | Piscanec/Lazzeri; Hwang–Das Sarma acoustic, TF-screened |
| 2D Rana Auger / carrier multiplication (ring) | ✅ | R−G on the gathered n, p |
| **two-temperature Coulomb sector** (`yn_sbe_rana_te`) | ✅ **NEW** | T_e and quasi-Fermi levels from the cone moments; R−G, Q_TF, plasmon line at T_e; lattice = e-ph bath (cooling via phonons); `*_sbe_te.data` (§5) |
| **2D collisional-memory analog** (`yn_sbe_colmem`, `_pop`) | ✅ | phonon lines for e-ph, 2D Dirac-plasmon line for the Rana source (`wiki/10` §8.11) |
| Option A dressed reference (`yn_sbe_dressed_ref`) | ✅ | removes the Dirac-point rotation background |
| 2D-sheet Σ^HF (`yn_sbe_coulomb`) | ✅ (off by default) | renormalizes v_F; enable deliberately |
| impact ionization / eeh / Kuhn–Zurek | 🚫 physics | gapless cone (error stop) |
| doped / finite-T initial occupation; hot phonons; substrate in the sheet field | — limits | §9 |

## 3. Field before → field after → transmission

The driver writes `E_ext` (incident) and, with `yn_sbe_sheet_field = 'y'`,
propagates in the **local** field of the sheet (Hartree a.u., Z₀ = 4π/c):

    E_t = E_inc − (Z₀/2) J_s,   J_s = −Jm·L_z,   dA_ind/dt = −(2π/c) L_z Jm,
    A_tot = A_ext + A_ind,      E_tot = E_ext + (2π/c) L_z Jm

so the `E_tot`/`Ac_tot` columns of `*_sbe_rt.data` **are** the transmitted field
and the energy ledger is the work of the local field (`wiki/12` §5). Then

    T = ∫E_t²/∫E_inc²,  R = ∫E_r²/∫E_inc² (E_r = E_t − E_inc),  A = 1 − T − R,

fluence-integrated (Parseval-exact; a single-FFT-bin "T at the carrier" is not
bounded by 1 for a reshaped pulse and is not used). `transmission.py` detects
the self-consistent mode (`SC`), prints the deviation of `E_tot` from the
boundary-condition reconstruction on `Jm` as a consistency number (`dEt`, the
explicit-Euler lag), and — for runs without the flag (`pert`) — the
radiation-reaction term `S_rr = (Z₀/2)∫J_s²/F = A_E − A` that tells whether the
perturbative estimate is trustworthy. Linear universal sheet: T = 0.97746,
A = 0.02241 (`tests/test_sheet_transmission.py`).

## 3a. Basis: the velocity-gauge f-sum rule and the pure-gauge restoration

In the velocity gauge every electron of the filled π band carries the diamagnetic
current A·N_e/V; a *complete* basis cancels it exactly (a uniform A is a pure
gauge), an `nstate`-band basis only to the fraction S = ⟨Σ_m 2|p_nm|²/Δε⟩ (0.70 for
2 bands, 0.90 for 4, 0.964 for 8, 0.970 for 16 — the rest sits > 10 eV up). The
remainder η N_e A/V is a reactive current ∝ E/ω: negligible in the near-IR,
decisive at 3 THz, where the bare 2-band sheet reflects 85 % (a plasma mirror) and
even 16 bands leave R ≈ 8 %. The solver prints S and η at start-up. With
`yn_sbe_vg_sumrule = 'y'` it subtracts, at every step, the **adiabatic ground-state
current of the same truncated H_k(A(t))** (one ZHEEV per k): identically zero in a
complete basis, exact at every A and for any population, no fitted quantity
(wiki/12 §6a, `tests/test_vg_sumrule.f90`). A linear static form (−η N_e A) was
tried first and withdrawn — it over-corrects at 100 kV/cm and the anti-inductive
sheet runs away under the self-consistent field (T > 1). `sumrule_check.py` reports
S, η and the residual A-projection of any run (must be ≈ 0). With the
restoration the transmission is independent of `nstate` (2/3/4/8 agree at 24²,
§7.1–7.2), so production runs on the cheapest basis, `nstate = 2`.

## 3b. Doped / finite-temperature initial state (a real sample)

`sbe_ef_ev` (Fermi level from the Dirac point, eV) and `sbe_temp_init_k` set the
initial occupation to occ_max·f_FD(ε; μ, T) instead of integer filling. That is
what turns the sheet from an intrinsic semimetal into a *metal*: a partially
filled band, a Drude weight D = 2k_BT ln[2cosh(μ/2k_BT)] → E_F, an intraband
conductivity, and the 60–70 % THz transmission real CVD graphene has. Three
things go with it:

- **The pure-gauge reference stays undoped.** The f-sum-rule restoration of §3a
  subtracts the adiabatic ground-state current of the truncated basis. Evaluated
  on the *doped* occupation it would also subtract the physical intraband
  current, which is the whole point of a doped run. The solver therefore keeps
  `gs%occup_ref` (the undoped, T→0 filling) for that subtraction alone. The
  doped carriers' own truncation error is smaller than their Drude weight by
  η/⟨∂²ε/∂k²⟩ ≈ 0.28/13.5 ≈ 2 %.
- **The mesh has to resolve the Fermi surface.** k_F = E_F/ħv_F must exceed a few
  mesh spacings, else the density is set by a handful of points. At E_F = 0.2 eV
  (k_F = 0.0167 a.u.) that means N ≳ 147 (9 points per valley inside the Fermi
  disc, 3 % density error); N = 24 gives *zero* and 1.8×10¹⁰ cm⁻² instead of
  3.4×10¹². The solver counts the partially occupied k-points at start-up and
  warns below 20; `plot_occupation.py` draws the initial occupation of the doped
  and the undoped sheet on the actual mesh and prints the sheet density, k_F and
  that count — run it before committing to a mesh. Cheap alternative for a smoke test: raise E_F (E_F = 0.6 eV is
  resolved at N = 48).
- **The dressed reference and the basis-edge monitor follow the doped baseline.**
  `dressed_ref_delta` now takes the initial occupation vector (adiabatic rotation
  of ρ₀ minus ρ₀, still trace-neutral), and the VG basis-edge warning measures
  the *excess* over the initial occupation — a doped metal has its top band full
  at every k inside the Fermi surface, which is not a basis failure.

## 4. Level check — the Dirac cone the SBE runs on

`python3 tests/test_graphene_dirac_levels.py`:

| basis (`epm_pw_cutoff_ry`) | PW | gap at K | v_F | Γ bottom / M dip |
|---|---|---|---|---|
| 2.94 (old x11 input) | 7 | **0.2125 eV — spurious** (Python = Fortran bandpath) | — | −8.5 / −2.8 eV |
| **29.4 (x11 + x14 now)** | 43 | 6.5×10⁻⁶ eV (Python), 0.0000 (Fortran) | 0.960×10⁶ m/s | −7.78 / −2.70 eV |

The 7-vector set is not closed under the little group C₃ᵥ of K, so the
truncation breaks the symmetry protection of the Dirac degeneracy. Also: the
half-shifted Monkhorst–Pack mesh (2i−N−1)/2N contains K = (2/3, 1/3) **only for
odd multiples of 3** (…, 147, 153); 12, 24, 150 straddle K by half a step —
hence nk = 147 here.

## 5. Population saturation and the two-temperature model

`tests/test_rana_saturation.f90`: the Coulomb balance R = G holds at the intrinsic
density n_i(T) = (π/6)(k_BT/ħv_F)² = **8.08×10¹⁰ cm⁻² at 300 K** (∝ T²); the
CPTP channel relaxes to it monotonically from above (Auger) and below (carrier
multiplication). Evaluated at the lattice temperature this under-estimates the
plateau of a hot plasma by (T_e/T_L)². With `yn_sbe_rana_te = 'y'` the common
carrier temperature and the quasi-Fermi levels are **read from the first two
moments (n, p, energy) of the gathered cone populations** at every ring step
(`dirac_fit_te`; `tests/test_dirac_te_fit.f90` recovers T to 10⁻⁴ from
explicit-mesh moments); R−G, Q_TF and the plasmon line run at T_e, the phonon
Bose factors at the lattice T — the two-temperature model with cooling through
the e-ph channel, without a separate T_e rate equation. `*_sbe_te.data`
records t, T_e, μ_c, μ_h, n, p, n_i(T_e), T_bath; `saturation_check.py` plots it.

## 6. Run

**The production scans, shipped ready to run on the cluster.** Three input sets sit
in the exercise; they are the server runs this exercise exists for.

| set | mesh | occupation | variants | fields [kV/cm] | what it gives | cost |
|---|---|---|---|---|---|---|
| `prod_nk297_doped/` | 297² | E_F = 0.2 eV, 300 K | coh | 1…1000 (7) | **the converged T(E₀) curve — the headline result** | O(N_k): minutes/field per node, few core-hours total |
| `prod_nk297_intrinsic/` | 297² | undoped | coh | 1…1000 (7) | the control the doping is measured against | same |
| `prod_nk147/` | 147² | E_F = 0.2 eV, 300 K | diss, mem | 1…300 (5) + dark | τ, the mean free path, T_e, the absolute absorption | ring is O(N_k²): ≈7.5 h/run on 48 threads |
| `prod_nk297_ef02_layers/L1`, `L2` | 297² | E_F = 0.2 eV, 300 K | ring **commented out** | 1…300 (7) + dark | **one against two layers at the sample doping** — the T(E₀) rise, and whether a bilayer brightens where a monolayer does not | coherent as shipped: minutes/field. With the ring uncommented: O(N_k²), a cluster job |

**`prod_nk297_ef02_layers` is the set to run on the cluster.** 297² is the first mesh
that resolves the 0.2 eV Fermi disc (k_F/Δk = 3.19; 147² gives 1.58 and 72² only 0.77 —
below one spacing, where the doping stops being representable at all, §7.12b). It ships
**coherent**: every ring switch is present but commented, so the cheap scan runs as-is
and the expensive one is a one-character edit per line:

```bash
sed -i 's/^  !\(yn_sbe_\(superres\|eph\|eph_acoustic\|auger\)\|sbe_\(eph_temperature_k\|coulomb_epsilon\|search_sigma_e_ev\)\)/  \1/' prod_nk297_ef02_layers/L*/rt_*.inp
```

Note what is *not* commented: `yn_sbe_vg_sumrule` and `yn_sbe_sheet_field` stay on in
both modes — they are the pure-gauge restoration and the sheet boundary condition, not
dissipation. Run `rt_dark_diss.inp` whenever the ring is on: §7.11 shows the zero-field
control is mandatory at the low-field end.

All four sets are `nstate = 4` since 2026-09-05 (§7.13): the pure-gauge restoration
makes an *undoped* sheet basis-independent, but it subtracts the reference occupation
only, so a **doped** run at `nstate = 2` carries a 12 % error in σ. And the mesh is
297², not 300²: 297 = 3 × 99 is an odd multiple of 3, so the Dirac point sits **on**
the half-shifted MP mesh, which is the grid rule this exercise states for 147² and
which the old 300² set quietly broke (k_F/Δk is 3.19 at 297² against 3.22 at 300², so
nothing is lost).

```bash
cd samples/exercise_x14_graphene_self_induced_transparency
# 1) the converged transmission curve (cheap, this is the main result)
for d in prod_nk297_doped prod_nk297_intrinsic; do
  ( cd $d && cp ../run_scan.sh . && OMP_NUM_THREADS=48 SALMON=../../../build/salmon bash run_scan.sh )
done
python3 field_scan_plot.py --doped 'prod_nk297_doped/runs/*/graphene_sit_sbe_rt.data'                            --intrinsic 'prod_nk297_intrinsic/runs/*/graphene_sit_sbe_rt.data'                            --t-meas 0.60 0.70 --n-sub 1.65 --continuum --out T_of_field_nk297.png
python3 plot_occupation.py prod_nk297_doped/graphene_sit --ef-ev 0.2   # pre-flight: 140 partial k-points
python3 drift_saturation.py prod_nk297_doped/graphene_sit --ef-ev 0.2  # the saturation curve vs continuum

# 2) the dissipative half (tau and the absolute absorption)
cd prod_nk147 && cp ../run_scan.sh . && OMP_NUM_THREADS=48 SALMON=../../../build/salmon bash run_scan.sh
python3 ../drude_check.py 'runs/*/graphene_sit_sbe_rt.data' --t-meas 0.60 0.70 --n-sub 1.65
python3 ../saturation_check.py runs/E100kVcm_mem/graphene_sit runs/dark_mem/graphene_sit --plot
```
Everything can also be driven end to end by `run_field_scan.sh` with `NK`, `EF`,
`FIELDS` and `VARIANTS`. A 24² smoke set is in `smoke_nk24/`.

**Watching a long scan.** `run_scan.sh` appends every finished field to
`runs/scan_progress.log` as well as printing it. Use the file, not the driver's
stdout: a driver that pipes the script through `grep` or `tee` block-buffers it, so a
watcher tailing the driver log sees nothing at all until the whole scan ends (this bit
us on a 14-run 297² scan). `tail -f runs/scan_progress.log`, or poll the row count —
a finished run has `nt + 7` lines in `graphene_sit_sbe_rt.data`. A watcher should also
alarm on "no `salmon` process while the scan is unfinished", since a crash and a
still-running job look identical to a completion-only filter.

**Why 297² for the transmission and 147² for the dissipators.** The Fermi surface
must be resolved (§3b: k_F ≳ 3 mesh spacings means N ≳ 280 at E_F = 0.2 eV), and at
~300² the drift-saturation curve is converged onto the continuum (§7.12). The unitary
propagation is O(N_k), so that mesh is cheap. The graphene ring is **O(N_k²)**, so
the same mesh costs (88209/21609)² ≈ 17× the 147² dissipative run; the dissipative
production therefore runs at 147², where the Fermi surface is marginal but the
scattering physics is not mesh-critical. A converged dissipative 297² scan needs MPI
over k across nodes (`wiki/11`).
Variants: `coh` (no dissipation), `diss` (ring e-ph + acoustic + Rana at the lattice
T, Markovian — as x11), `mem` (diss + 2D colmem analog + dressed reference + T_e).
All carry the sheet field (`--no-sheet` to switch it off) and the velocity-gauge
**pure-gauge restoration** `yn_sbe_vg_sumrule = 'y'` (`--no-sumrule` to switch it
off; §3a). `--field acos2 --hw-ev 0.8 --cycles 8` gives the near-IR pulse of the
first x14 version. Other switches: `--nstate 2` (default; ≤ 4 at dt = 0.1 fs, ≥ 8 needs dt ≤ 0.05 fs — §7.2),
`--pol x|y` (in-plane polarisation; the beam is at normal incidence, k ∥ z, so E is
always in the sheet plane), `--n-layers 2` (two electronically decoupled sheets in
the same local field — incoherent/large-angle-twisted bilayer; `sbe_sheet_nlayers`),
`--snap-fs 50` (k-resolved level-population snapshots for `plot_levels.py`),
`--ef-ev 0.2 --temp-init-k 300` (a doped sheet with Drude carriers, §3b).

**The headline figure in one command** (doped + intrinsic scans on the same mesh,
then the T(E₀) / σ(E₀) plot of §7.9):

```bash
NK=147 EF=0.2 OMP_NUM_THREADS=48 SALMON=../../build/salmon bash run_field_scan.sh
# -> scan_nk147_ef0.2/doped_vs_intrinsic.png
```
`NK` must resolve the doping: k_F = E_F/ħv_F needs ≈ 3 mesh spacings, so NK ≳ 280
at E_F = 0.2 eV, 140 at 0.4, 93 at 0.6 (wiki/12 §4a.0). NK = 147 at 0.2 eV is the
cheapest setting that resolves the Fermi surface at all and is what the shipped
figure uses; the Drude weight is then ≈ 30 % low, a field-independent scale error
that leaves the shape of T(E₀) intact. Coherent runs only, so the whole pair of
scans is minutes per field.

```bash
python3 ../plot_levels.py runs/E100kVcm_mem/graphene_sit --times 100,150,200,300   # band populations vs t + k-maps around K
python3 ../sumrule_check.py graphene_sit --runs 'runs/*/graphene_sit_sbe_rt.data'  # basis f-sum rule + residual reactive current
python3 ../drude_check.py runs/*/graphene_sit_sbe_rt.data --t-meas 0.60 0.70 --n-sub 1.65  # doped sheet: D, tau, mean free path
```

**Cost (measured 2026-09-04, 4 threads).** Coherent runs: seconds to minutes at
147². Dissipative runs: the graphene ring (e-ph inter-k + Rana) is O(N_k²) and
exponential-bound — ≈ 60 ms/step at 24² with the THz drive (all cone states are
sources) ⇒ ≈ 84 s/step at 147² on 4 threads ⇒ **≈ 7 s/step on 48 threads ⇒ ≈ 7.5 h
per 3844-step run**. Levers: `--dt-fs 0.2` (halves it; the CF4 unitary is exact
per step, check `nex` against 0.1 fs), MPI ranks over k (the ring gather is
MPI-parallel, `wiki/11`). Essential set = 100 kV/cm × {coh, diss, mem} + `dark_mem`.

**Do not economise on `--tail-fs`.** The 100 fs ring-down after the 284 fs transient
looks like 25 % of the cost for nothing, but the sheet current is still ringing when
the drive ends, and every reported number is a fluence integral. Measured at 72²,
E_F = 0.6 eV, coherent, with and without it:

| E₀ [kV/cm] | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|
| T, `--tail-fs 0` | 0.6846 | 0.6807 | 0.6439 | 0.6454 | 0.7495 |
| T, `--tail-fs 100` | 0.7093 | 0.7053 | 0.6679 | 0.6713 | 0.7739 |
| A, `--tail-fs 0` | 0.0589 | 0.0593 | 0.0592 | 0.0554 | 0.0775 |
| A, `--tail-fs 100` | 0.0095 | 0.0100 | 0.0112 | 0.0035 | 0.0287 |

T comes out 0.024 low at every field (the shape survives, the absolute value does
not) and the absorption of a *coherent* run — which has no dissipation at all and
must be near zero away from pair creation — is inflated six-fold, because the
transmitted fluence that is still in flight at the last step is simply not counted.

**Dissipative runs are the exception, and for the reason you would guess:** the ring
damps the current, so by the end of the drive there is little left in flight. Same
mesh and doping, `diss`:

| | T, no tail | T, +100 fs | A, no tail | A, +100 fs |
|---|---|---|---|---|
| 10 kV/cm | 0.77773 | 0.77864 | 0.20363 | 0.20182 |
| 100 kV/cm | 0.82761 | 0.82785 | 0.15297 | 0.15248 |

— a shift of 0.0009 and 0.0002 against the coherent run's 0.0246. So the tail is
cheap insurance in general and indispensable for `coh`; if the ring is on and the
budget is tight, it is the one place the 25 % can honestly be saved.

`transmission.py` also refuses to analyse a record shorter than the run's own `nt`
(`--allow-partial` overrides), so a job that is still going cannot quietly contribute
a number.

## 7. Local validation (2026-09-04, 4 threads, container)

All runs: DAST single-cycle 3.36 THz proxy, x-polarised, normal incidence,
sheet field ON, `coh` unless stated; T, R, A fluence-integrated from
`transmission.py`; "resid." = residual A-projection of the sheet current in
units of N_e A/A_2D (`sumrule_check.py`, must be ≈ 0).

**7.1 The basis artifact and its removal (24², 100 kV/cm).**

| nstate | restoration | T | R | A | A_tot/A_ext | after the pulse | resid. |
|---|---|---|---|---|---|---|---|
| 2 | off | 0.147 | 0.851 | 0.002 | 0.20 | A_ind relaxes | 0.2935 |
| 4 | off | 0.646 | 0.321 | 0.033 | 0.53 | relaxes | — |
| 8 | off | 0.937 | 0.013 | 0.049 | 0.90 | relaxes | — |
| 16 | off | 0.945 | 0.005 | 0.050 | 0.92 | relaxes | — |
| 2 | linear static −ηN_eA (withdrawn) | 1.213 | 0.499 | −0.712 | 1.71 | **runaway** (A_ind, J grow together) | −0.023 |
| 8 | linear static (withdrawn), `diss` | — | — | — | 1.40 | **runaway** | — |
| 2 | pure gauge (Eq. 10) | 0.999996 | 1.3e-6 | 0.000002 | 1.000 | stationary | 0.0000 |
| 3 | pure gauge | 0.999996 | 1.4e-6 | 0.000003 | 1.000 | stationary | 0.0000 |
| 4 | pure gauge | 0.999996 | 1.4e-6 | 0.000003 | 1.000 | stationary | 0.0000 |
| 8 | pure gauge | 0.972 | 3.4e-4 | 0.027 | 0.986 | slow DC drift | 0.0001 |

The 24² mesh has no interband channel at 14 meV (smallest π–π* gap 0.78 eV, K
half a spacing off the mesh, no k-point in the Landau–Zener tube), so the
physical answer at this resolution is T ≃ 1 — which nstate = 2, 3 and 4 now
give identically. The 2.7 % of nstate = 8 is population leaked into bands
4–8 (22–58 eV above π; 5.6×10⁻⁶ per cell carrying 2.5×10⁻⁴ eV per cell, the
whole ledger) — a 14 meV field cannot do that; see 7.2 for the time-step test.
Linear checks of the mechanism (24², nstate = 2, 1 kV/cm, no sheet): the SBE
current equals the exact adiabatic response of the truncated H(A) to 1×10⁻⁴
(dt 0.1/0.05/0.02 fs), in-phase coefficient η = 0.2995 = Eq. (9); at 0.1 and
0.4 eV drive: 0.2985 / 0.2808 vs the dispersive prediction 0.2985 / 0.2806.

**7.2 Time step at nstate = 8 (24², 100 kV/cm).** With dt = 0.05 fs the nstate = 8
run gives T = 0.999985, A = 1.4×10⁻⁵, ledger 1.1×10⁻⁷ eV per cell — the same
answer as nstate = 2/3/4, and 2400× less deposited energy than at dt = 0.1 fs.
The 2.7 % at dt = 0.1 fs was the S4/CF4 step applied to bands 34–90 eV above
the cone (13.7 rad of phase per step): population leaks into them. Rule: at
dt = 0.1 fs use nstate ≤ 4; nstate ≥ 8 needs dt ≤ 0.05 fs. **Production recipe
therefore: nstate = 2 + pure-gauge restoration** (`make_inputs.py` default) —
the restoration makes 2, 3 and 4 bands agree to 10⁻⁶ in T, the THz physics
lives on the cone, and the ring cost is lowest.

**7.3 147² (K on the mesh), nstate = 2, pure gauge, coherent — the field scan.**
Static η = 0.2777 (x = y: the 147² mesh restores the hexagonal isotropy the 24²
mesh breaks at the 10 % level) — removed by the restoration; what remains in the
A-projection is the *physical* Drude response of the pairs the field creates.

| E₀ [kV/cm] | T | R | A | ledger A_E | pairs after the pulse [cm⁻²] | peak diabatic n_c (per cell) | A_tot/A_ext | Drude resid. |
|---|---|---|---|---|---|---|---|---|
| 1 (K arbitrary, see below) | 0.816 | 1.82 | −1.63 | −2.08 | 1.7×10³ | 8×10⁻¹¹ | 0.002 | 60.6 |
| 10 (K arbitrary) | 0.9569 | 0.0348 | 0.0083 | 0.0083 | 2.1×10¹⁰ (4.2×10⁻⁵ per cell) | 2.2×10⁻⁴ | 0.818 | 0.0100 |
| 30 | 0.9936 | 0.0039 | 0.0025 | 0.0025 | 3.5×10¹⁰ | 2.4×10⁻³ | 0.938 | 0.0029 |
| 100 | 0.9824 | 0.0043 | 0.0134 | 0.0134 | 6.7×10¹¹ (3.5×10⁻⁴ per cell) | 2.0×10⁻² (99 % virtual, returns) | 0.941 | 0.0048 |
| 1 (K averaged, rerun) | 0.999995 | 1.2×10⁻⁶ | 4×10⁻⁶ | 4×10⁻⁶ | 0 net | 9.3×10⁻⁵ (= baseline) | 1.000 | 0.0000 |
| 10 (K averaged, rerun) | 0.999993 | 2.4×10⁻⁶ | 4×10⁻⁶ | 4×10⁻⁶ | 0 net (returns to the 0.9×10⁻⁴ baseline) | 3.8×10⁻⁴ | 1.000 | −0.0001 |
| 30 (K averaged, rerun) | 0.9995 | 2.4×10⁻⁴ | 2.6×10⁻⁴ | 2.6×10⁻⁴ | 2×10⁹ net (1×10⁻⁶ per cell) | 2.6×10⁻³ | 0.995 | 0.0005 |
| 100 (K averaged, rerun) | 0.9867 | 0.0030 | 0.0103 | 0.0103 | 4.0×10¹¹ net (3.0−0.9×10⁻⁴ per cell; 0.9×10⁻⁴ = the two half-filled K points) | 2.0×10⁻² | 0.958 | 0.0041 |
| 10, 150² (K off mesh) | 0.9925 | 0.0004 | 0.0071 | 0.0059 | 1×10¹⁰ (5×10⁻⁶ per cell) | 3.8×10⁻⁴ | 1.001 | −0.0003 |
| 100, 150² (K off mesh) | 0.9563 | 0.0067 | 0.0370 | 0.0365 | 1.1×10¹² (5.6×10⁻⁴ per cell) | 2.0×10⁻² | 0.956 | 0.0035 |
| 100, 300² (K off mesh) | 0.9688 | 0.0046 | 0.0266 | 0.0264 | 7.7×10¹¹ (4.0×10⁻⁴ per cell) | 2.0×10⁻² | 0.950 | 0.0041 |

**The Dirac point on the mesh (found in this scan).** With K exactly on the
mesh the two π/π* levels at K are degenerate and the integer ground-state
filling (band 1 = 2, band 2 = 0) picks LAPACK's arbitrary basis inside the
degenerate pair: a broken-symmetry state with a velocity expectation of order
v_F at that k-point — a spurious current ~ 2v_F/N_k per valley that does not
depend on the field. At 100 kV/cm it is buried under the physical response;
at 10 kV/cm it is the whole "Drude residual" (0.010) and the 3.5 % reflection;
at 1 kV/cm, together with the sheet field, it acts as a relay (the pure-gauge
reference jumps by ±2v_F w_K at A = 0 and the sheet pins A_tot to zero: the
nonsense row above, R > 1). Fix (`gs_info_ssbe.f90`, all materials): a
degenerate, partially filled level group is given its **group-average
occupation** — the T → 0 density matrix, the identity on the block, current-free
and rotation-invariant; the pure-gauge reference becomes continuous. The rows
marked "K averaged" are the rerun with this fix; the 150² rows have K half a
spacing off the mesh (no degeneracy at all) and bracket the result from the
other side. Below ~30 kV/cm the Landau–Zener tube is thinner than one mesh
cell (0.01 / 0.05 / 0.3 / 1.6 / 9.6 cells at 1 / 3 / 10 / 30 / 100 kV/cm on
147²), so the pair creation there is Eq. (8) of wiki/12 analytically —
1.5×10⁹ / 4.7×10¹⁰ / 2.5×10¹¹ / 1.5×10¹² cm⁻² per passage at 1 / 10 / 30 / 100
kV/cm — rather than mesh-resolved; the 30 and 100 kV/cm rows are resolved.

At 100 kV/cm the sheet current is 96 % anti-correlated with A_tot *after* the
restoration: that is the kinetic inductance of the ~10¹² cm⁻² pairs created by
the Landau–Zener sweep through K (Eq. 8 of wiki/12 predicts 1.5×10¹² per
passage; two passages with Stückelberg interference give the 6.7×10¹¹ left),
their creation energy is the 1.3 % absorbed (ledger = fluence to 3 digits), and
the induced field relaxes after the pulse (A_ind −6.1 → −4.0×10⁻³ a.u. over the
tail, J decaying) — a stable, passive sheet. **Expected range of the
transmission change for the intrinsic sheet, coherent limit: T ≈ 1.000 up to
30 kV/cm (induced absorption < 0.1 %), T = 0.97 ± 0.01 at 100 kV/cm (induced
absorption 1–4 % across the 147²/150²/300² meshes and the two polarisations,
2.7 % on 300²; reflection ≤ 0.7 %)** — i.e. the intrinsic sheet *darkens* by a
few per cent at the top of the DAST range and does not bleach. With phonon
dissipation the carriers are also heated while accelerated, and the
`diss`/`mem` production runs will raise A (README §8). The mesh spread is the
transverse sampling of the Landau–Zener strip (0.7 cells at 147², 1.4 at 300²);
convergence to ≲ 0.5 % needs N ≳ 600 or a K-refined mesh.

**7.4 Polarisation and two decoupled layers (147², 100 kV/cm).** The `input`
field file carries its own three Cartesian components — `epdir_re1` is *not*
applied to it (a first "y" run silently reproduced the x result to six digits;
`make_inputs.py --pol y` now writes the waveform into the A_y column). True y
run (147², K averaged, 100 kV/cm): T = 0.9735, R = 0.0056, A = 0.0209 against
x: T = 0.9867, R = 0.0030, A = 0.0103. This factor of 2 is **not** the crystal's
anisotropy (see below) but the mesh's: the Landau–Zener tube is a strip
2k_⊥ = 7.5×10⁻³ a.u. wide (0.7 mesh cells at 147²) and 2A₀ = 0.12 long (12
cells), and how many mesh points fall inside it depends on its orientation on
the hexagonal mesh — the same sampling spread as 147² vs 150² for one
polarisation (1.0 % vs 3.7 %). The transverse tube width sets the convergence
requirement: ≳ 3 cells across the tube means N ≳ 600 at 100 kV/cm (or a
K-refined mesh); the 300² run below is the first step of that series.

Two electronically decoupled sheets in the same local field
(`sbe_sheet_nlayers = 2`, 147², K averaged, 100 kV/cm): T = 0.9585, R = 0.0112,
A = 0.0304 against the single layer's T = 0.9867, R = 0.0030, A = 0.0103. The
naive estimate T₁² = 0.9736 is 1.5 % (absolute) too high: in the coherent
non-linear regime the two layers do not simply multiply — they share one local
field (A_tot/A_ext = 0.92 instead of 0.96) and the re-radiated field reshapes
the waveform each layer sees. For the *linear* regime of real CVD samples the
sum-of-conductances formula applies (d ≪ λ): with z = Z₀σ complex,

    T₁ = |2/(2+z)|²,   T₂ = |2/(2+2z)|²,   T₁²/T₂ = 16|1+z|²/|2+z|⁴
                                                  = 1 + ½[(Im z)² − (Re z)²] + O(z³)

**The sign of the T·T error depends on what kind of sheet it is.** For a purely
dissipative (real σ) sheet T₁² is LOW: 0.03 % at T₁ = 0.98, 0.6 % at 0.89,
9.6 % at the sample's z = 0.565 (T₁²/T₂ = 0.904), 25 % at T₁ = 0.40. For a
purely reactive (inductive, imaginary σ) sheet — which is what the coherent
doped runs here are — T₁² is HIGH: the same |z| = 0.565 gives T₁²/T₂ = 1.133,
+13 %. A doped THz sheet is a mixture, so |z| alone does not even fix the sign;
only |z| ≤ 0.15 (T₁ ≥ 0.87) keeps the error under 1 % either way. The
substrate's own Fresnel factor must also be divided out *before* squaring. For Bernal (AB) bilayer or small-angle moiré the electronic structure
itself changes and neither estimate applies.

**The same comparison on a DOPED sheet (72², +100 fs tail, coherent), which is what
a real CVD bilayer is** — figure `wiki/figures/graphene_layers_1_vs_2.png`, rebuilt by

```bash
python3 layers_plot.py --set 0.6 'ef0.6_L1/runs/*/*_rt.data' 'ef0.6_L2/runs/*/*_rt.data' \
                       --set 0.4 'ef0.4_L1/runs/*/*_rt.data' 'ef0.4_L2/runs/*/*_rt.data' \
                       --t-meas 0.60 0.70 --out layers.png
```


| E₀ [kV/cm] | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|
| E_F = 0.6 eV: T₁ | 0.7093 | 0.7053 | 0.6679 | 0.6713 | 0.7739 |
| T₂ | 0.4362 | 0.4341 | 0.4128 | 0.3778 | 0.4715 |
| T₁²/T₂ | 1.153 | 1.146 | 1.081 | 1.193 | 1.270 |
| E_F = 0.4 eV: T₁ | 0.8487 | 0.8422 | 0.7939 | 0.8358 | 0.8772 |
| T₂ | 0.6493 | 0.6434 | 0.5847 | 0.5894 | 0.7226 |
| T₁²/T₂ | 1.109 | 1.103 | 1.078 | 1.185 | 1.065 |

T₁² is HIGH at every field and both dopings — the reactive branch, because a coherent
doped sheet at 14 meV is an inductor (band-averaged z = 0.017 − 1.195i at 0.6 eV).
√T₂ is then SMALLER than T₁ (0.660 against 0.709 at 3 kV/cm), so inverting a measured
bilayer as √T₂ makes the monolayer look darker than it is and **overstates** its
conductance. With the ring on the inequality reverses — see §7.11b. Two further checks: `--predict-layers 2`, which uses the run's
own σ(ω), reproduces T₂ to 1.5 % in the linear regime (against 14 % for the
single-frequency formula) and degrades to 27 % at u = 5.6, always predicting the
stack too bright; and bin by bin across the driven band the two-layer run returns
exactly twice the one-layer σ(ω) (median ratio 1.997 and 2.007), which is the check
that `sbe_sheet_nlayers` itself is right. The band-averaged Re σ ratio is 1.83 rather
than 2 because the two-layer transmitted field is filtered more strongly where |σ| is
largest — a property of that diagnostic, not of the solver.

**Answer to "does the direction matter?"** For the linear and the χ⁽³⁾ response
the hexagonal sheet is isotropic in-plane (any 2nd- and 4th-rank tensor of C₆ᵥ
is); anisotropy (zigzag vs armchair) enters at χ⁽⁵⁾ through trigonal warping,
i.e. at relative order (A₀/|K−M|)² ≈ (0.06/0.78)² ≈ 0.6 % of the already-small
non-linear part — below 10⁻⁴ in T at 100 kV/cm. On a mesh the x/y difference is
a resolution artifact: the 24² mesh breaks C₆ (S_x = 0.700 vs S_y = 0.728), the
147² mesh does not (0.7223 both).

**7.5 Near-IR πα check (147², acos2 0.8 eV, 8 cycles, nstate = 2, pure gauge, sheet).**
1 kV/cm: T = 0.98128, R = 6.8×10⁻⁴, A = 0.01804 (ledger 0.01804); 100 kV/cm:
T = 0.98136, A = 0.01796. The sheet absorbs 0.80 of the universal πα value
(A = 0.0224; Re σ/σ_univ = 0.90 over the pulse band) — the 3.2 mesh points per
resonance-shell radius and the finite-basis matrix elements at k_res, not a
solver error (the reactive residual is −0.013, i.e. the interband *capacitive*
response above resonance, the correct sign). The field dependence is the
predicted coherent bleaching: ΔA/A = −0.4 % (θ²/12 ≈ −0.5 %), ΔT = +8×10⁻⁵ —
near-IR self-induced transparency at 1–100 kV/cm is a 10⁻⁴ effect.

**7.6 Ring with T_e (24², nstate = 2, 100 kV/cm, `diss` / `mem`).** The 24² mesh
has no real pair channel at this field (coherent: A = 2×10⁻⁶), so whatever the
dissipative runs absorb is what the ring does to the *virtual* dressing (peak
diabatic n_c = 2.3×10⁻² per cell, all of it returning in the coherent run).
`diss` (Markovian ring, lattice T): T = 0.9988, A = 0.0012; ring-visible density
7.6×10⁹ cm⁻² generated out of the dressing, Rana ledger +1.9×10⁹ (net
multiplication, n < n_i) — the "dephasing ionization" of `wiki/10` §8.7 in its
graphene form. `mem` (2D colmem analog + dressed reference + T_e): T = 0.999996,
A = 3×10⁻⁶ — identical to the coherent run; ring-visible density 5×10⁴ cm⁻²
(10⁵ times less), Rana ledger +3.5×10⁶: the memory filters remove the
fabricated generation completely, as they do for Si/GaAs/CdS (`wiki/10` §8.11).
The two-temperature fit is meaningless when there are no carriers (it returned
T_e ≈ 4×10⁴ K for 10⁸ cm⁻²); the solver now holds T_e at the lattice value while
n + p < 10⁻³ n_i(T_lattice), so `*_sbe_te.data` reads 300 K there. Real carrier
heating (T_e in the 10³ K range, cooling on the phonon timescale) needs the
147² production runs, where pairs are actually created.

**7.7 Above 100 kV/cm — the intrinsic sheet keeps darkening (147² K-averaged and
150², nstate = 2, pure gauge, self-consistent sheet, coherent).**

| E₀ [kV/cm] | T (147²) | R | A | Re σ/σ_univ | T (150²) | A (150²) |
|---|---|---|---|---|---|---|
| 1 | 0.999995 | 1.2e-6 | 4e-6 | 0.000 | — | — |
| 10 | 0.999993 | 2.4e-6 | 4e-6 | 0.000 | 0.9925 | 0.0070 |
| 30 | 0.9995 | 2.4e-4 | 2.6e-4 | 0.03 | — | — |
| 100 | 0.9867 | 0.0030 | 0.0103 | 0.78 | 0.9563 | 0.0370 |
| 200 | 0.9505 | 0.0099 | 0.0396 | 2.83 | 0.9515 | 0.0389 |
| 300 | 0.9336 | 0.0178 | 0.0486 | 4.02 | 0.9404 | 0.0459 |
| 500 | 0.9369 | 0.0200 | 0.0430 | 3.94 | 0.9230 | 0.0535 |
| 1000 | 0.9163 | 0.0335 | 0.0503 | 5.52 | 0.9103 | 0.0546 |

Three things settle here. **(i)** The transmission of the *intrinsic* sheet falls
monotonically over three decades of field and never rises: there is no
self-induced transparency of undoped graphene in this model. **(ii)** The
absorbed fraction saturates near 5 % above ~300 kV/cm, which is the analytic
Landau-Zener value (wiki/12 §7: the pairs created per passage scale as E^{3/2}
and their creation energy as √E, so the absorbed *energy* scales as the fluence
and the *fraction* is field-independent). **(iii)** The two meshes, which
differ by a factor 3.7 in A at 100 kV/cm, agree to 1–2 % above 200 kV/cm: once
the Landau-Zener strip is several mesh cells wide the sampling ambiguity of
§7.3 is gone. The growing reflection (3.4 % at 1000 kV/cm) and Re σ = 5.5 σ_univ
are the Drude response of the created plasma, not an artifact — the residual
A-projection stays below 0.005 throughout.

**7.8 A doped sheet — what the measured 60 % → 70 % needs.** A sample on PET
transmitting 60 % (substrate included) is not the intrinsic sheet of §7.7. With
n_PET = 1.65 and two incoherent faces the bare substrate passes 88.3 %, so the
graphene itself accounts for 0.60/0.883 = 0.68 and 0.70/0.883 = 0.79, i.e. a
sheet conductance Z₀σ of 0.565 → 0.327, **σ = 24.7 → 14.3 σ_univ (−42 %)**.
That is a Drude conductance: σ_dc = Dτ/π with the Dirac-cone Drude weight
D = 2k_BT ln[2cosh(μ/2k_BT)] → E_F. In transport units it is R_s = 666 Ω/sq at low
field and 1153 Ω/sq at high field — an ordinary as-transferred CVD monolayer.
Splitting it across dopings, each with the τ that reproduces the same σ:

| E_F [eV] | n_2D [cm⁻²] | τ [fs] | mobility [cm²/V s] | mean free path [nm] | onset A₀ = k_F [kV/cm] |
|---|---|---|---|---|---|
| 0.1 | 8.0×10¹¹ | 127 | 11700 | 122 | 13 |
| **0.2** | **3.2×10¹²** | **64** | **2900** | **61** | **27** |
| 0.3 | 7.2×10¹² | 43 | 1300 | 41 | 40 |
| 0.4 | 1.3×10¹³ | 32 | 730 | 31 | 54 |
| 0.6 | 2.9×10¹³ | 21 | 330 | 20 | 81 |

E_F = 0.2–0.4 eV is where an ordinary sample sits; 0.1 eV would need a mobility
CVD-on-polymer does not have, and 0.6 eV a mobility that is low even for
chemically doped film. Independent corroboration: Hassanpour Amiri *et al.*,
*Doping free transfer of graphene using aqueous ammonia flow*, RSC Adv. **10**,
1127 (2020), attribute the unintentional doping of transferred CVD graphene to
ionic etch residue of **typical density 4×10¹² cm⁻²** — on the Dirac cone that is
E_F = 0.224 eV and a saturation onset of 30 kV/cm, within 12 % of what the
transmission alone gives. Their ammonia-washed (doping-free) films should follow
the *intrinsic* curve instead and darken with field: a cheap test of the mechanism. The field at which the transmission starts to rise picks
between them (last column). The E_F = 0.6 eV used in the 48² runs below is a
mesh-affordable proxy, not the sample.

`sbe_ef_ev` / `sbe_temp_init_k` now put that state into the solver (§3b). Which
of the two factors, D or τ, can give −42 %?

- **D cannot.** At *fixed* carrier density, heating lowers μ but adds thermal
  carriers; D passes a shallow minimum and comes back up. Over 300–3000 K the
  deepest excursion is 0.88 at n = 10¹³ cm⁻² and 0.90 at 3×10¹² (table in
  `drude_check.py`, pinned by `tests/test_doped_drude.py`). Carrier heating on
  its own is worth ~10 %, not 42 %.
- **τ can, and so can the drift itself.** Two field scales control it. The
  Fermi sea is displaced by the vector potential, A₀ = 6.2×10⁻⁴ a.u. per kV/cm
  for this transient, against k_F = E_F/ħv_F = 0.0167 a.u.: they are equal at
  **27 kV/cm**. Above that the displacement exceeds the Fermi radius every
  half-cycle, the drift velocity saturates at v_F, and the differential
  conductivity falls like k_F/A₀. The collisional excursion eEτ/ħ crosses k_F at
  the same place for τ = 63 fs. Simultaneously the carriers are pushed past the
  optical-phonon threshold (E₂g 196 meV, A₁′ 160 meV; v_F A₀ = 0.74 eV at
  100 kV/cm), and every emission randomises momentum — τ drops.

So the mechanism the measurement is showing is the one you suspected: the
conductivity falls because the carriers stop responding linearly (mobility /
mean free path), not because the Drude weight disappears. The onset the model
puts at ~27 kV/cm for E_F = 0.2 eV is inside the measured range, and the
predicted saturation is *sub-linear in the drift*, which is why the transmission
creeps up by ten points rather than jumping.

**7.10 At the sample's own doping (147², E_F = 0.2 eV) — the signature is the ratio
to the intrinsic control.** A 3×10¹² cm⁻² sheet without dissipators is only a weak
inductor, so the raw transmission is quiet (T = 0.968 at 1 kV/cm, not the measured
0.68 — a coherent run has no momentum relaxation and cannot produce absorption,
only reactive screening). Divide by the intrinsic control on the same mesh and the
doping carriers stand alone:

| E₀ [kV/cm] | 1 | 10 | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|---|---|
| T doped (E_F = 0.2 eV) | 0.9682 | 0.9666 | 0.9595 | 0.9629 | 0.9442 | 0.9090 |
| T intrinsic | 0.999995 | 0.999993 | 0.999498 | 0.986694 | 0.933550 | 0.916273 |
| extinction added by the doping | 3.18 % | 3.34 % | **4.00 %** | 2.41 % | −1.14 % | 0.79 % |
| Re σ/σ_univ (doped) | 2.73 | 2.86 | 3.45 | 2.81 | 3.73 | 6.11 |

The doping-induced extinction **peaks at 30 kV/cm**, against the 27 kV/cm that
A₀ = k_F predicts from k_F alone (k_F = E_F/ħv_F = √(πn) is the radius of the
occupied disc, and in the velocity gauge A is literally its displacement in
reciprocal space — wiki/12 §4a.5.1), then collapses by a factor of five. At 300 kV/cm
it is *negative*: the doped sheet transmits better than the undoped one, because its
Drude response has saturated while its occupied states Pauli-block part of the
Landau–Zener pair creation. Figure: `doped_vs_intrinsic.png` (three panels: T(E₀),
the extinction ratio, σ(E₀)); the 48²/E_F = 0.6 eV proxy of §7.9 is the same physics
with a sheet conductance large enough to move the raw transmission.

**7.9 The doped sheet, calculated: transmission rises above the saturation
onset.** 48² mesh, E_F = 0.6 eV (n = 3.0×10¹³ cm⁻², 8 mesh points per valley
inside the Fermi disc), T_init = 300 K, nstate = 2, pure gauge, self-consistent
sheet, **coherent** (no dissipators at all — so everything below is mechanism 3
alone). D_fit and τ_fit come from fitting the run's own current to
dJ_s/dt = (D/π)E_tot − J_s/τ.

| E₀ [kV/cm] | T | R | A | D_fit [eV] | Re σ/σ_univ |
|---|---|---|---|---|---|
| 1 | 0.7286 | 0.2497 | 0.022 | 0.533 | 23.6 |
| 3 | 0.7284 | 0.2498 | 0.022 | 0.533 | 23.6 |
| 10 | 0.7257 | 0.2519 | 0.022 | 0.537 | 23.8 |
| 30 | 0.6990 | 0.2737 | 0.027 | 0.580 | 25.8 |
| 100 | 0.6622 | 0.3304 | 0.007 | 0.664 | 30.3 |
| 200 | 0.7202 | 0.2736 | 0.006 | 0.546 | 25.4 |
| 300 | 0.7716 | 0.2066 | 0.022 | 0.424 | 20.5 |
| 500 | 0.8122 | 0.1455 | 0.042 | 0.325 | 15.8 |
| 1000 | 0.8263 | 0.1124 | 0.061 | 0.290 | 13.4 |

The **intrinsic control on the same mesh**, same pulse, same everything except the
initial occupation:

| E₀ [kV/cm] | 1 | 10 | 30 | 60 | 100 | 200 | 300 | 500 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| T, doped (E_F = 0.6 eV) | 0.7286 | 0.7257 | 0.6990 | 0.6414 | 0.6622 | 0.7202 | 0.7716 | 0.8122 | 0.8263 |
| T, intrinsic | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9970 | 0.9772 | 0.9519 | 0.9230 | 0.9026 |

(figure: `doped_vs_intrinsic.png`). The two curves run in opposite directions and
cross near 700 kV/cm.

Read from the bottom up, this is the measurement. **T is flat at 0.726 through
the linear regime (1–10 kV/cm), dips to 0.66 near 100 kV/cm, then rises
monotonically to 0.83.** The turn is where predicted: A₀ = k_F at 81 kV/cm for
E_F = 0.6 eV. Above it the fitted Drude weight collapses, 0.66 → 0.29 eV, and
with it the sheet conductivity, **30.3 → 13.4 σ_univ (−56 %)** — against the
−42 % (24.7 → 14.3 σ_univ) the PET measurement implies. Same sign, same
mechanism, same order. **The dip before the turn is not physics.** §7.12 shows it is the
discretization bump of a Fermi disc holding only 12 mesh points: the continuum
drift-saturation curve is monotone, and the T minimum at 60 kV/cm falls exactly on
u = A₀/k_F = 0.386, where the 48² mesh curve peaks at 1.171. Converged, the doped
sheet is flat and then brightens, with no darkening in between.

Note what is **not** in this table: no phonons, no Auger, no heating — the
dissipators are off. Mechanism 3 (current saturation on the Dirac cone) accounts
for the whole effect on its own. Mechanisms 1 and 2 add to it: heating supplies
≈10 % more, and optical-phonon emission above 196 meV shortens τ, which lowers σ
further in the same direction. So the measured rise does not need all three —
but none of them makes the sheet *darker*, and the model has no mechanism that
would.

*Drude-weight accuracy — measured convergence.* The fitted D is below the analytic
D = E_F because the mesh under-resolves ∂²ε/∂k² ∝ 1/k near the Dirac point. At
1 kV/cm (linear regime):

| mesh | E_F [eV] | partially occupied k-points | n_2D vs analytic | D_fit/E_F |
|---|---|---|---|---|
| 147² | 0.2 | 36 | 3.27 vs 3.36×10¹² (−2.6 %) | 0.659 |
| 300² | 0.2 | 116 | 3.345 vs 3.36×10¹² (−0.5 %) | **0.930** |
| 147² | 0.4 | 72 | 1.284×10¹³ | 0.894 |
| 48² | 0.6 | 8 | 3.04 vs 2.89×10¹³ (+5 %) | 0.888 |

The density needs one mesh shell inside the Fermi circle, the Drude weight three
or four. The deficit is a scale error common to all fields, so the *shape* of the
T(E₀) curve above is unaffected; absolute conductivities want k_F ≳ 3 mesh
spacings (N ≳ 280 at E_F = 0.2 eV).

**7.11 With the phonon ring on: the scattering time and the mean free path.**
Same 48² mesh and doping, `diss` (e-ph ring at a 300 K lattice), against the
coherent runs:

| E₀ [kV/cm] | τ [fs] | mean free path [nm] | T_e | Re σ/σ_univ | T | R | A |
|---|---|---|---|---|---|---|---|
| 10 | 141 | 136 | 361 K | 18.1 | 0.6442 | 0.061 | 0.294 |
| 100 | 60 | 58 | 2050 K | 14.5 | 0.7213 | 0.048 | 0.231 |
| 100, coherent | 7930 | 7613 | — | 30.3 | 0.6622 | 0.330 | 0.007 |

Two things are robust here. **(i)** Scattering turns the sheet from an inductive
mirror (R = 0.33, A = 0.007) into a Drude absorber (R = 0.048, A = 0.231) — a real
sample, not a mirror. **(ii)** T rises with field.

**The full scan on the better mesh** (72², E_F = 0.6 eV, ring at 300 K, 100 fs tail):

| E₀ [kV/cm] | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|
| T | 0.7786 | 0.7789 | 0.7978 | 0.8279 | **0.8949** |
| Re σ/σ_univ | 10.65 | 11.16 | 10.07 | 8.49 | **5.15** |
| A | 0.202 | 0.195 | 0.178 | 0.152 | 0.089 |

Flat through 10–30, then monotone to 0.895, σ down 54 % from its peak — and **no mesh
dip**, where the coherent curve on the same 72² dips 7.2 %. Momentum relaxation
broadens the Fermi surface over many mesh cells, which is what the discretization
error needs.

**The 3 kV/cm point is omitted; read this before running any low-field dissipative
scan.** A dissipative run started from a Fermi–Dirac occupation carries a current at
**zero field**: the ring's fixed point is not exactly that occupation, and on a finite
mesh its relaxation is not isotropic. The `dark_diss` control measures it —
|J_dark|max = 2.42×10⁻⁸ a.u., steady over the whole 384 fs record. Being
field-independent it matters in proportion to how weak the drive is:

| E₀ [kV/cm] | 3 | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|---|
| dark / peak current | **36 %** | **14 %** | 3.3 % | 2.0 % | 1.3 % | 0.6 % |
| T as reported | 0.7968 | 0.7786 | 0.7789 | 0.7978 | 0.8279 | 0.8949 |
| T, dark subtracted | 0.8439 | 0.7826 | 0.7780 | 0.7973 | 0.8275 | 0.8947 |
| shift | **+0.047** | +0.004 | −0.001 | −0.000 | −0.000 | −0.000 |

Above ~10 % the point is unusable and **not repairable by subtraction**: the dark
current of a driven run is not the field-free one (the field changes the distribution
the ring acts on), so the first-order correction overshoots — at 3 kV/cm it lands on
0.844 where the linear plateau is 0.778. Below ~3 % the correction is under 0.001.
Between them, subtract and quote both. Run `dark_*` with every dissipative scan;
`transmission.py --dark runs/dark_diss/graphene_sit_sbe_rt.data` prints the fraction
and flags anything past 10 %.

**A number that did not survive the mesh.** On 48² this scan gave T = 0.644 → 0.721
and 14.5 σ_univ at high field against the measured 14.3 — an almost exact match. At
72² the same calculation gives 0.779 → 0.895 and 5.15 σ_univ. That agreement was a
coincidence of an under-resolved Fermi disc, in the same family as the D → E_F
agreement at nstate = 2 (§7.13). The *shape* is robust; the absolute level carries the
unconverged τ below.

**Caveat: neither the absolute τ nor its field dependence is mesh-converged.**
Repeating the scan on 72² (2.25× the k-points):

| mesh | E₀ [kV/cm] | τ [fs] | mean free path [nm] | Re σ/σ_univ | T |
|---|---|---|---|---|---|
| 48² | 10 | 141 | 136 | 18.1 | 0.644 |
| 48² | 100 | 60 | 58 | 14.5 | 0.721 |
| 72² | 10 | 33 | 31 | 10.6 | 0.779 |
| 72² | 100 | 72 | 69 | 8.5 | 0.828 |

(72² rows: the tail-inclusive repeat. A second, independent route agrees — inverting
σ(ω) bin by bin, τ = −Im σ/(ω Re σ), with no time-domain fit and no driven window,
gives 21.5 fs at 10 kV/cm and 58.0 fs at 100. Same rise, from data the fit never
touches, so the reversal against 48² is not an artifact of the fitting window.)

The magnitude moves by 4× at fixed field, and the *trend* reverses: at 48² τ falls
2.3× from 10 to 100 kV/cm, at 72² it rises by 2.2× over the same interval. An earlier
version of this file read the 48² fall as the carriers crossing the optical-phonon
thresholds (v_F A₀ = 0.74 eV at 100 kV/cm against E₂g 196 meV, A₁′ 160 meV); the
finer mesh does not reproduce it, so **that reading is withdrawn** — it was mesh
noise, not a mobility drop that the calculation can claim. The e-ph ledger rate
doubles between the two meshes for the same reason: the inter-k golden-rule sum is
unresolved, the σ_E = 0.1 eV energy window holding too few final states, so the rate
samples the mesh. What is robust is the sign
of every effect at fixed mesh when the ring is switched on: scattering shortens τ
against the collisionless run, converts reflection into absorption, lowers σ and
raises T. Quoting an absolute τ, or comparing the absolute transmission with the
measurement quantitatively, needs the ring on a Fermi-surface-resolving mesh: the
O(N_k²) production run of §6.

The carrier temperature comes from the total electronic energy (361 K and 2050 K),
*not* from the per-channel `dE_eph` column, which for a doped initial state is a
gross exchange counter: it grows at 7×10⁻⁶ eV/cell/fs from t = 0 regardless of
field, while the electronic energy of the doped sea moves by only 2.5×10⁻⁵ eV/cell
over 384 fs at 10 kV/cm and the carrier number is conserved to 0.14 %. The doped
Fermi sea is stationary under the dissipators; the ledger column is not a net loss.

*Still to run (cluster).* The same with `--ef-ev 0.2 --temp-init-k 300 --variants
diss,mem` on a Fermi-surface-resolving mesh (147² or denser, §3b), which replaces
the E_F = 0.6 eV proxy by the sample's own doping. Cost as in §6.

**7.12 Drift saturation, and separating it from the mesh.** The semiclassical picture
behind all of this is drawn by `cone_kinematics.py` (figure
`wiki/figures/graphene_cone_kinematics.png`):

```bash
python3 cone_kinematics.py --ef-ev 0.6 --drive prod_nk147/DAST_E100kVcm.txt
```

In the velocity gauge the label k does not move; the field enters only through
H_k(A) = v_F σ·(k + A). So the occupied disc stays put and what travels through it is
the point of instantaneous degeneracy, k = −A. An occupied conduction state has
velocity v_F (k+A)/|k+A| — it points radially *away* from that moving point, with
modulus v_F always — so at A = 0 the arrows cancel and by A ≫ k_F they are parallel
and the current saturates. Whether the electron can jump to the upper cone is pure
Pauli: for |A| < k_F the degeneracy sits inside the occupied disc and the target state
is already filled; for |A| > k_F it has left the disc and pairs are created in a strip
of half-width √(E/πv_F). The same number A₀ = k_F therefore ends the linear regime and
opens pair creation — 27 kV/cm at E_F = 0.2 eV, 81 kV/cm at 0.6 eV.

The rise of T is the
kinematics of the Dirac cone: the band velocity has fixed modulus v_F, so a Fermi
disc displaced by A saturates at J = n e v_F, and the differential conductivity is

    sigma_eff/sigma_lin = G(u)/u,  u = A_0/k_F,   G(u)/u -> k_F/A_0  for u >> 1

with G(u) the exact displaced-disc integral (wiki/12 Eqs. 4a.12-4a.13). This is the
CHORD response -- the current the same peak field produces, which is what a field
scan measures. The differential response a weak probe on a strong pump would see is
G'(u), and it falls faster, as 1/(4u^3) (wiki/12 Eq. 4a.16); the two must not be
swapped. At 300 K the continuum values are

| u = A₀/k_F | 0.05 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| σ_eff/σ_lin (continuum) | 1.000 | 0.995 | 0.979 | 0.952 | 0.912 | 0.851 | 0.630 | 0.487 | 0.332 | 0.200 |
| 300², E_F = 0.2 eV (140 partial pts) | 0.930 | 1.062 | 1.015 | 0.978 | — | 0.878 | — | 0.543 | — | 0.249 |
| 147², E_F = 0.2 eV (36) | 6.2 | 1.225 | 1.000 | 0.984 | — | 0.829 | — | 0.532 | — | 0.243 |
| 48², E_F = 0.6 eV (12) | 0.901 | 1.015 | 1.171 | 1.062 | — | 1.003 | — | 0.635 | — | 0.333 |

The continuum curve is **monotone** — no darkening at any field. A Fermi disc holding
one shell of mesh points reproduces the linear limit to ~10 % but develops a spurious
15–20 % bump near u ≈ 0.2–0.5, and at very small u on the coarsest discs the ratio
breaks down entirely (the 147² entry at u = 0.05 is nine points, two of them the
degenerate K pair). That bump is the dip of §7.9: the 48² transmission minimum sits
at 60 kV/cm = u 0.386, exactly the peak of the 48² mesh curve. It shrinks as the disc
fills (300²: 1.06 at worst).

Measured on the runs instead of on the adiabatic sum, the artifact is the spurious
peak Re σ develops before saturation, and it is governed by k_F/Δk alone — not by the
doping (same transient, +100 fs tail, coherent, E_F = 0.6 eV unless noted):

| run | k_F/Δk | σ plateau | σ peak | spurious bump | T plateau | T minimum | dip in T |
|---|---|---|---|---|---|---|---|
| 48², E_F = 0.6 eV | 1.54 | 23.6 | 31.1 | **+32 %** | 0.7275 | 0.6414 | −11.8 % |
| 72², E_F = 0.4 eV | 1.54 | 12.9 | 17.4 | **+35 %** | 0.8487 | 0.7939 | −6.5 % |
| 72², E_F = 0.6 eV | 2.32 | 25.9 | 30.8 | **+19 %** | 0.7081 | 0.6571 | −7.2 % |
| 111², E_F = 0.6 eV | 3.57 | 25.6 | 27.8 | **+8.6 %** | 0.7149 | 0.6900 | −3.5 % |
| 147², E_F = 0.6 eV | 4.73 | 27.0 | 29.0 | **+7.6 %** | 0.6997 | 0.6763 | −3.3 % |

Two different meshes at two different dopings but the same k_F/Δk give the same bump;
filling the disc from 1.5 to 3.6 spacings takes it from a third to under a tenth —
and then it stops, 111² and 147² differing by 78 % in k-points and by nothing in the
dip. The four-mesh figure is `wiki/figures/graphene_T_of_field_mesh.png`; rebuild it
with

```bash
V=<scan root>
python3 field_scan_plot.py \
  --doped "$V/n147/runs/*/graphene_sit_sbe_rt.data" --doped-label '147² ($k_F/\Delta k$ = 4.7)' \
  --series "111² (3.6):$V/n111/runs/*/graphene_sit_sbe_rt.data" \
           "72² (2.3):$V/n72/runs/*/graphene_sit_sbe_rt.data" \
           "48² (1.5):$V/n48/runs/*/graphene_sit_sbe_rt.data" \
  --continuum --t-meas 0.60 0.70 --n-sub 1.65 --out T_of_field_mesh.png
```

The four curves agree above A₀ = k_F to better than 0.02 and disagree below it, which
is what makes the disagreement a discretization error rather than a field scale. How
much of the bump reaches T depends on the sheet impedance, which is why the same
artifact makes a 11.8 % dent at 0.6 eV and 6.5 % at 0.4 eV — the dip depth is not the
artifact, the conductivity bump is. Reproduce the comparison with
`field_scan_plot.py --continuum`, which draws the parameter-free continuum curve
beside the data and prints the residual.

**Most of the darkening before the rise is numerical, but not all of it.** The mesh
series settles the first half: three quarters of the dip disappears between
k_F/Δk = 1.5 and 3.6, and the disappearance is governed by that ratio alone, not by
the doping. It does not settle the second: the dip stops shrinking after that —
−3.5 % at 111² and −3.3 % at 147², with 78 % more k-points between them — so a
residual few-per-cent darkening near A₀ ≈ k_F survives mesh refinement. It is not the
basis either: the same 147² scan at nstate = 4 dips by 3.41 % against 3.39 % at
nstate = 2, unchanged to two digits, though the basis moves T by 0.035 and Re σ by
11 % (§7.13). Both numerical suspects are therefore excluded and the dip is left
unexplained; the remaining candidate is the departure of the real EPM band from an
ideal cone over the excursion A₀ ≈ k_F, which the continuum derivation assumes away. The honest shape is therefore: flat while A₀ ≲ 0.5 k_F, a residual few-per-cent
dip around A₀ ≈ k_F, the drift-saturation rise as 1/E₀, and Landau–Zener darkening at
several hundred kV/cm. An under-resolved Fermi surface multiplies that residual dip by
three or four, which is why Eq. (4a.3) still has to be satisfied.
`drift_saturation.py` produces the table and `graphene_drift_saturation.png`.

**7.12b A doping the mesh cannot represent produces GAIN — check `nfs` before
believing any dissipative number.** Compressing the mesh to exercise the ring (32² and
33², E_F = 0.2 eV, two layers, `diss`) turned this up. 33 = 3×11 puts the Dirac point
*on* the half-shifted mesh; at that size the Fermi disc (k_F = 0.0167) is smaller than
one spacing (0.0473), so **no** k-point is partially occupied and the doping charge
lands entirely on fully filled or empty levels:

| | 32² (K off mesh) | 33² (K **on** mesh) |
|---|---|---|
| partially occupied k-points | 2 | **0** |
| A at 100 kV/cm | +0.168 | **−0.240** |
| Eall − Eall0 | +8.1×10⁻⁴ eV/cell | −1.2×10⁻³ eV/cell |
| \|J\|last/\|J\|max | 0.36 | 0.90 |

Energy is conserved — the electrons *lose* and the field *gains*. That is gain from an
inverted population. The **zero-field control proves the inversion is in the initial
state, not made by the drive**: with `ae_shape1 = 'none'` the 33² sheet still has a
current that grows monotonically to the end of the record (|J| = 1.32×10⁻⁹ a.u., its
own maximum, still rising) while the electronic energy falls, against 1.2×10⁻¹¹ and
flat at 32². The dark current alone is far too small to explain A = −0.24, so the
large negative absorption is the field *amplifying* a state that was already inverted
at t = 0.

**Switching channels off does not fix it — only the mesh does.** The same zero-field
33² state with the ring configured three ways:

| channels | \|J\|max | still growing at t_end | Eall − Eall0 [eV/cell] |
|---|---|---|---|
| e-ph ON, Auger ON | 1.32×10⁻⁹ | yes | −8.2×10⁻¹⁰ |
| e-ph OFF, Auger ON | 6.4×10⁻¹⁰ | **yes** | −5.0×10⁻¹¹ |
| both OFF (coherent) | 1.3×10⁻¹⁶ | no | −3.4×10⁻²³ |

Dropping the phonons and keeping Auger halves the amplitude but the current still grows
— the instability is fed by *whichever* channel is on, because the pathological initial
state is not a fixed point of any of them and on this mesh the relaxation runs towards
more inversion, not less. Only disabling every dissipator removes it, and that is not a
fix, it is deleting the physics. The cause is upstream of the ring: a Fermi disc the
mesh cannot hold. Raise `num_kgrid`.

Two guards now exist. `gs_info_ssbe` separates `nfs == 0` from `nfs < 20` and calls the
first an error: the doping cannot be represented, the run can develop gain, do not use
it. `transmission.py` refuses to let A < 0 pass silently. This is the same family as
the degenerate-K bug of §7.1, which group-averaging the occupation fixed — that fix
does not cover a Fermi disc containing no mesh point at all.

**7.11b The bilayer with the ring on — the T·T error changes sign.** Everything in
§7.4 about layers was coherent, where the doped sheet is an inductor. With the ring it
is a Drude absorber, and Eq. (6a) predicts the opposite sign. Same 72² mesh,
E_F = 0.6 eV, ring at 300 K, +100 fs tail, `sbe_sheet_nlayers` 1 and 2:

| E₀ [kV/cm] | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|
| u = A₀/k_F | 0.12 | 0.37 | 0.74 | 1.24 | 3.71 |
| 1 layer, T | 0.7786 | 0.7789 | 0.7978 | 0.8278 | **0.8949** |
| 2 layers, T | (0.6714) | 0.6164 | 0.6552 | 0.6936 | **0.7969** |
| Re σ₁/σ_univ | 10.7 | 11.2 | 10.1 | 8.5 | 5.2 |
| Re σ₂/σ_univ | (23.2) | 20.7 | 17.5 | 15.1 | 9.6 |
| T₁²/T₂ | (0.903) | **0.984** | **0.971** | **0.988** | 1.005 |
| \|A−A_E\|/A, 2 layers | **0.81** | 0.07 | 0.02 | 0.01 | 0.00 |

The 10 kV/cm column is quarantined and the figure omits it. Two independent
diagnostics condemn it and agree: the fluence absorption and the electron-energy ledger
disagree by 5× (A = 0.144 vs A_E = 0.027, others ≤7 %), and the two-layer `dark_diss`
control puts the field-independent ring current at 12.6 % of peak there against 5.6 /
4.0 / 2.4 / 1.3 % at 30 / 60 / 100 / 300 — the only field past the 10 % line of §7.11.
Subtracting the dark run moves T by +0.006 and does not repair the ledger mismatch,
just as the "not repairable above 10 %" rule predicts. On the four clean fields the
dark correction is ≤0.004.

Coherent, the same doping gave T₁²/T₂ = 1.15 at every field; with the ring the
trustworthy points give 0.97–0.99 — below unity at all of them. The effect is smaller
than the quarantined point suggested (−1 to −3 %, not −10 %) but the sign is
unambiguous, and the sign is what matters. The cause is the phase of z: Re z/|z| = 0.014 coherent (inductor) against
0.899 dissipative (Drude resistor). A real CVD sample is resistive, so it sits on the
branch where T₁² **under**estimates T₂ and √T₂ understates the monolayer conductance —
the opposite of the coherent conclusion, and the one that applies to a measurement.

Both stacks brighten, measured from their minimum at 30 kV/cm over the clean points:
+14.9 % for one layer, +29.3 % for two. The bilayer gains more for the
reason Eq. (4a.17) gives — it starts at twice the conductance, where the same
fractional saturation buys more transmission. σ₂/σ₁ stays near 2 throughout (2.18,
1.85, 1.78, 1.86). Figure: `wiki/figures/graphene_layers_ring.png`, built by
`layers_plot.py`, which reads the ring flag from the runs and labels the branch from
the measured phase rather than assuming it.

**7.14 Testing the drift law itself, not the transmission it produces.** $T(E_0)$ is a
compressed view of Eq. (4a.13): the sheet BC $T=|2/(2+z)|^2$ is nonlinear in $z$ and
depends on the *phase* of $z$, which the drift law says nothing about. The test that
removes both is the fit-free inversion of §7.11 — $\tau(\omega)=-\mathrm{Im}\,\sigma/(\omega\,
\mathrm{Re}\,\sigma)$, $D(\omega)=-\pi(\omega^2+1/\tau^2)\mathrm{Im}\,\sigma/\omega$ — because it
returns $D$ and $\tau$ separately, and a conductance that falls because the disc drifted
looks identical in $T$ to one that falls because the carriers scatter more. 72², E_F = 0.6 eV,
ring on, anchored at 30 kV/cm (the lowest field its dark control clears):

| E₀ [kV/cm] | 30 | 60 | 100 | 300 |
|---|---|---|---|---|
| u = A₀/k_F | 0.371 | 0.742 | 1.237 | 3.711 |
| D [eV] | 0.2581 | 0.2309 | 0.2010 | 0.1283 |
| τ [fs] | 59.1 | 53.9 | 58.0 | 62.7 |
| D(u)/D(u₀) | 1.000 | **0.894** | **0.779** | 0.497 |
| G(u)/u, same normalisation | 1.000 | **0.941** | **0.750** | 0.273 |
| residual | — | −5.0 % | +3.8 % | +82 % (off the cone) |

The law holds to 5 % and 4 % where it applies, and τ moves by ±8 % about 58 fs with no
trend — so the fall is weight, not scattering, which is the claim T(E₀) alone cannot
establish. Past u = 2 the residual runs away, as §7.12 requires. Build the figure
(`wiki/figures/graphene_drift_fit_1layer.png`) with

```bash
V=<scan root>
python3 drift_fit_plot.py \
  --series "one layer, ring on, 72\$^2\$:$V/runs/E{30,60,100,300}kVcm_diss/*_rt.data" \
  --excluded "$V/runs/E3kVcm_diss/*_rt.data" "$V/runs/E10kVcm_diss/*_rt.data" \
  --anchor-kvcm 30 --u-max 4.0 --e0-max 1000
```

`--e0-max` carries the law past the last run so its whole bend is visible; the curve is
also drawn *down* through the quarantined fields, which is worth doing: at 10 kV/cm the
condemned point lands on the law to 0.4 %, at 3 kV/cm it is 2.8 % off. The contamination
is not a uniform offset that could be calibrated away — invisible in T at 14 % dark
current, plain at 36 % — which is why §7.11 gives a threshold and not a correction.

Two caveats. It is **two independent points** (four fields clear the dark control, one is
the anchor, one is off the cone); the `prod_nk297_ef02_layers` set is where to repeat it,
because at E_F = 0.2 eV the saturation field is 27 kV/cm and the whole range 0 < u < 2 is
reachable at fields the dark control clears. And **only the shape is tested**: the same
inversion puts the ring run's absolute weight at D/D_eq = 0.43 — see §7.15.

**The collisionless sheet does not obey the law**, which is the opposite of what the
derivation suggests, Eq. (4a.13) being collisionless itself. Same inversion, coherent
runs, both meshes, same anchor:

| u | 0.012 | 0.371 | 0.742 | 1.237 | 3.711 |
|---|---|---|---|---|---|
| D(u)/D(u₀), 147² | 0.963 | 1.000 | 1.030 | 1.020 | 0.718 |
| D(u)/D(u₀), 72² | 0.896 | 1.000 | 1.055 | 1.006 | 0.698 |
| G(u)/u | 1.018 | 1.000 | 0.941 | 0.750 | 0.273 |
| residual, 147² | −5 % | — | +10 % | +36 % | +163 % |

Its absolute weight is right (D/D_eq = 1.01 at 147², so the doped ground state and the
f-sum restoration are doing their job) but its field dependence is not: D *rises* to a
maximum near u ≈ 0.74 and does not start falling until u > 2. That rise is the
conductivity bump of §7.11, reached by a second route: measured from each run's own
lowest field it is +17.7 % at 72² and +7.0 % at 147², against +19 % and +7.6 % for the
bump in Re σ on the same meshes. Two differently-extracted quantities agreeing to under a
per cent on both meshes says the bump belongs to the solution, not to either diagnostic.
Refinement halves it and does not remove it, exactly as §7.11 found.

The reading this suggests — a hypothesis, not a result — is that G(u)/u is a *quasi-static
chord* response, assuming the occupied disc sits at the position belonging to the
instantaneous A(t). A collisionless run has nothing to enforce that and keeps its weight
up; a run whose momentum relaxation (58 fs) is short against the drive period (300 fs) is
held near the quasi-static distribution and the geometry shows through. The test that
would settle it is cheap and not done: vary τ through the lattice temperature and see
whether the agreement tracks τ/T_drive rather than the presence of the ring.

**7.15 The channel ledger: what e-ph does, and the one thing it does wrong.**
`*_sbe_channels.data` carries a cumulative per-cell ledger — dN (conduction-population
change) and dE [Ha] (the eigenvalue-weighted energy the electrons *gained*) for each ring
channel. Read it together with `*_sbe_rt_energy.data`, whose Eall is **not** Tr(ρH):
`realtime_ssbe.f90` accumulates `energy += (E_tot·−J)·volume·dt`, the work the local field
does on the sheet. So Eall = W_field is an integration identity, and the energy left in
the electron gas is E_elec = W_field + Σ_ch dE_ch. Exactly one channel has a bath on the
other side — e-ph — so −dE_eph is what the carriers hand to the phonons; Auger and impact
ionization redistribute *within* the electron gas and their dE must come out near zero
while their dN does not. `channel_budget.py` does this, with the dark run subtracted:

```bash
python3 channel_budget.py "$V/runs/E*kVcm_diss/*_rt.data" --dark "$V/runs/dark_diss/*_rt.data"
```

72², E_F = 0.6 eV, ring on, dark-subtracted, per unit cell (0.01512 carriers/cell):

| E₀ [kV/cm] | 3 | 10 | 30 | 60 | 100 | 300 |
|---|---|---|---|---|---|---|
| W_field [meV] | 0.0017 | 0.0203 | 0.1682 | 0.6335 | 1.5119 | 8.5315 |
| to the lattice, −dE_eph [meV] | 0.0014 | 0.0173 | 0.0885 | 0.2466 | 0.6423 | 5.2112 |
| … as % of W | 82 % | 85 % | 53 % | 39 % | 43 % | 61 % |
| E_elec left [meV] | 0.0003 | 0.0030 | 0.0797 | 0.3870 | 0.8694 | 3.3220 |
| dN_eph /cell | 5.9e−11 | 4.8e−10 | 1.2e−08 | 2.4e−07 | 9.5e−07 | 4.0e−06 |
| dN_rana /cell | −1.3e−08 | −1.5e−07 | −1.3e−06 | −6.1e−06 | −5.1e−05 | −2.3e−03 |
| dE_rana [meV] | −4e−09 | −5e−08 | −1e−06 | −2e−07 | −2e−04 | 1.6e−03 |

**e-ph is an energy sink, not a population source.** It removes 39–61 % of the absorbed
work to the lattice inside the 400 fs window while creating almost no pairs: dN_eph tops
out at 4.0e−6 per cell at 300 kV/cm, 0.03 % of the doped carriers. **Rana is the
opposite** — dN_rana reaches −2.3e−3 per cell, 15 % of the carriers, with dE_rana four
orders of magnitude below W. That is exactly right for Auger: it moves carriers, not
energy, and the near-zero dE column is the check that it does. The 3 and 10 kV/cm
percentages are ratios of numbers at the 1e−3 meV level and are not meaningful.

**The zero-field pathology is e-ph, alone.** In the `dark_diss` control — no drive at all —
the ledger reads

| channel | dN /cell | dE [meV/cell] |
|---|---|---|
| e-ph | +2.8e−10 | **+2.7489** |
| Rana | −4.6e−10 | +0.0000 |
| impact ionization, ring Auger | 0 | 0 |

e-ph pumps 2.75 meV per cell into the electrons with the field switched off — **182 meV
per doped carrier**, of which 136 meV is already in by the time the pulse peaks at
~150 fs. For scale, the entire work done by the 100 kV/cm pulse is 1.51 meV/cell, so the
spurious pump is 1.8× the whole signal at that field. This is the energetic face of the
dark current of §7.11, and it names the channel: at this mesh Rana contributes
+0.0000 meV and −4.6e−10 carriers at zero field, i.e. nothing. **So switching the phonons
off and keeping Auger is a viable way to run a clean dissipative sheet at this mesh** —
which was not true of the `nfs == 0` gain bug of §7.12b, a different failure on a mesh
that cannot represent the doping at all, where Auger alone still grew the current.

It is also the leading suspect for the missing Drude weight of §7.14. The pulse arrives
at a sheet that has already absorbed 136 meV per carrier from a bath it should be in
equilibrium with, so the D_eq computed for FD(E_F, 300 K) is the wrong reference and the
measured D/D_eq = 0.43 is partly a statement about the reference. It cannot be the whole
story — a *thermal* distribution carrying that much excess energy would still keep
≥0.87 D_eq (§7.9's heating table) — so what the pump makes is non-thermal. The repair
belongs in the e-ph rates' detailed balance, not in the sheet or the mesh.

**7.16 The zero-field pumping was a broken detailed balance in the e-ph ring.**
§7.15 measured it and named the channel; this is what it turned out to be. The bug is
not in the bath temperature and not in the doping — it is in the energy-matching width.

The Gaussian matching has width `sbe_search_sigma_e_ev` (0.1 eV in production), so a
source is connected to partners at |ΔE| anywhere within a few σ of ħω_p, while the
emission/absorption split was taken **once per mode** from N_B(ħω_p). A pair that
actually transfers δ was therefore weighted by exp(ħω_p/kT) instead of exp(δ/kT), and
for δ ≫ ħω_p the upward rate was too large by exp((δ−ħω_p)/kT). Graphene is where that
becomes a runaway: its appended acoustic mode is **5.39 meV against a 0.1 eV search
width** and carries **98.2 %** of the channel weight (the two optical modes, 196 and
160 meV, have N_B ≈ 10⁻³ and weights 0.003 and 0.015 — they can barely absorb at 300 K).

Held at an *exact* FD(300 K) on a Dirac-like spectrum — where a channel in equilibrium
with its bath must do nothing — `eph_interk_dpop` gives, in meV per step:

| σ [eV] | before | after | before, **undoped** (μ = 0) |
|---|---|---|---|
| 0.100 | **+1.0418** | +0.000018 | +0.0579 |
| 0.020 | +0.0118 | +0.0000005 | +0.000125 |
| 0.005 | +0.000079 | ~0 | +0.000001 |
| 0.001 | 0 | 0 | 0 |

Three things to read off. The leak is **governed by σ**, falling ~90× per 5× narrowing.
The Fermi level **amplifies** it 18× — there are carriers with a sharp edge to smear —
but does not cause it: an undoped sheet leaks too. And the bath is innocent: it supplies
the right N_B, the code applied it to the wrong energy.

That is what the k-resolved dark run shows directly (§7.15): 350 fs with no field moves
the conduction occupation from f(0.50–0.58 eV) = 0.90 to 0.39 and from f(0.8–1.0 eV) =
3.5×10⁻⁵ to 0.079, valence untouched, particle number conserved to 10⁻¹¹ — a Fermi edge
diffusing outward. Fitting the tail gives a hot Fermi–Dirac at **T_e ≈ 2260 K,
μ ≈ 0.42 eV**, and k_BT_e = 0.195 eV is of the order of σ = 0.1 eV, not of the bath's
25.9 meV. (The distribution is approximately *thermal*, which corrects the reading in
an earlier revision of §7.15.)

**The fix**, `eph_thermal_split_de`: evaluate the split at the realized transfer. With
x = |δ|/kT the pair (N+1, N)/(2N+1) collapses to a logistic, f_emit = 1/(1+e^{−x}),
f_abs = 1/(1+e^{x}) — one exponential, exact at both ends (½:½ for a degenerate pair,
1:0 far above the bath). Second correction: the collision prefactor ν(ε) was taken from
the *source* alone, so one pair had two different rates for its two directions (a factor
2 across 0.4–0.8 eV at ε₀ = 0.8 eV); it is now the geometric mean over the pair, with
ν at the sink hoisted into a table so the pair loop keeps one exponential, not three.

**Scope.** Gated on the material (`mp%auger_2d_rana`): on for the gapless 2D Dirac
materials, off for Si, GaAs and CdS, which take the historical branch verbatim and stay
bit-identical to the validations published with them. The same violation exists there —
it is not a graphene-specific bug — but their optical modes make σ/ħω ≈ 1 rather than
≈ 20, so it is a small correction, not a runaway. Whether to enable it for them needs
its own runs and is deliberately not decided here.

`tests/test_eph_detailed_balance.f90` holds a Dirac spectrum at an exact FD(T_bath) and
requires the channel to leave it alone, checks that the residual shrinks with σ, and
checks that the gate still reproduces the historical path.

**The rescan.** The whole dissipative monolayer scan was repeated with the corrected
channel — same 72², same E_F = 0.6 eV, same pulse and time grid as §7.11, so the two are
directly comparable.

*The pump is gone, not reduced.* In the zero-field control the ring current fell from
|J|max = 2.4×10⁻⁸ to **5.3×10⁻¹⁹** (eleven orders) and the energy it injected from
+2.74889 to **+0.00007 meV/cell** — 40 800× down, and 800× below the 0.055 meV it costs to
heat this sheet from 0 to 300 K. The dark fraction is **0.0 % at every field**, against
35.8 % and 14.0 % at 3 and 10 kV/cm before. The sign of the heat flow also flipped to the
physical one: dE_eph is now negative at every field (−0.015, −0.061, −0.213 meV/cell at 3,
10, 30 kV/cm) — the carriers cool into the lattice instead of being warmed by it.

*The law holds better.* Anchored the same way, at 30 kV/cm:

| u = A₀/k_F | 0.371 | 0.742 | 1.237 | 3.711 |
|---|---|---|---|---|
| D(u)/D(u₀) | 1.000 | 0.957 | 0.770 | 0.461 |
| G(u)/u | 1.000 | 0.941 | 0.750 | 0.273 |
| residual, corrected | — | **+1.7 %** | **+2.7 %** | +68.8 % (off the cone) |
| residual, before | — | −5.0 % | +3.8 % | +81.8 % |

The worst on-cone deviation falls from 5.0 % to 2.7 %. And the absolute weight recovers:
**D/D_eq = 0.69 against 0.43** before, so more than half of the missing Drude weight of
§7.14 was the spurious pump destroying it. The remaining 31 % is still unexplained.

*Transmission and the two absorption measures now agree.*

| E₀ [kV/cm] | 30 | 60 | 100 | 300 |
|---|---|---|---|---|
| T | 0.679 | 0.717 | 0.750 | **0.856** |
| A, fluence | 0.253 | 0.247 | 0.221 | 0.133 |
| A_E, energy ledger | 0.236 | 0.243 | 0.219 | 0.132 |
| disagreement | 7 % | 2 % | 1 % | **0.1 %** |

The sheet brightens by +26 % from its minimum (against +14.9 % before), and the fluence and
the electron-energy ledger — which disagreed by a factor 5 on the worst pre-fix point —
now agree to between 0.1 and 7 %.

**What the rescan did NOT fix, and it is a different fault.** At 3 and 10 kV/cm the runs
now pass the dark control and fail the *energy ledger*: A = −0.538 against A_E = −0.906 at
3 kV/cm, T + R = 1.54, and with the drive exactly zero over the last 60 fs the current sits
at its maximum and is still growing. The phonon bath is not the source — dE_eph is negative
there too. The suspect is the Rana Auger channel: dN_rana < 0 at every field, i.e. it is
*recombining*, and a doped gapless sheet has no holes for its Drude carriers to recombine
with. Those two points are drawn hollow in the figure and are not fitted. The fields from
30 kV/cm up are unaffected, which is why the drift-law test above stands.

**The Auger suspicion is wrong.** Repeating 3 and 10 kV/cm with `yn_sbe_auger = 'n'` —
everything else identical — reproduces the Auger-on numbers to five figures: T = 0.930077
against 0.930049, A = −0.538351 against −0.538301, A_E = −0.905860 against −0.905792 at
3 kV/cm, and the same at 10. The Rana channel is not the source of the low-field gain.

What that leaves: e-ph and the self-consistent sheet field. The dark control is clean
(|J|max = 5.3×10⁻¹⁹), so whatever it is needs the field to seed it.

**It is e-ph.** Repeating 3 kV/cm with `yn_sbe_eph = 'n'` and the sheet field still on:

| 3 kV/cm | e-ph on | e-ph off |
|---|---|---|
| T | 0.930 | 0.709 |
| R | 0.608 | 0.281 |
| A | **−0.538** | **+0.0095** |
| A_E | −0.906 | +0.0044 |
| T + R | **1.538** | **0.990** |
| \|J\|end / \|J\|max | **1.00** | **0.14** |

Absorption turns positive, the energy balance closes, and the current decays to 14 % of
its peak instead of sitting at its maximum with the drive at zero. So the sheet self-field
is exonerated along with Auger, and the low-field gain belongs to the same channel the
detailed-balance fix repaired — a second, smaller fault in it that the zero-field control
cannot see, because it only wakes once the field displaces the distribution.

(The residual A = 0.0095 against A_E = 0.0044 flags because the test is relative; in
absolute terms it is 0.005 on a 1 % absorption, with 14 % of the current still un-decayed
at the end of the window.)

The next discriminator is the parameter that produced the first fault: repeat 3 kV/cm with
`sbe_search_sigma_e_ev = 0.005` instead of 0.1, i.e. a matching width of order the acoustic
phonon itself. If the gain falls with σ it is the same disease, curable by narrowing the
window or by matching energies more carefully.

**The law extrapolated to the sample's own doping.** `drift_fit_plot.py --predict-ef 0.2`
draws T(E₀) for another doping with nothing fitted: the equilibrium Drude weight
D_eq(E_F, T) fixes the linear sheet response, G(u)/u with that doping's own k_F gives the
field dependence, and the sheet boundary condition turns the pair into a transmission.
At E_F = 0.2 eV the saturation field moves to **27 kV/cm** from 81 — linearly in E_F — so
the brightening starts three times earlier and goes further: the law gives T ≈ 0.99 at
300 kV/cm against 0.86 at 0.6 eV, because the sheet is loaded half as heavily (|z| = 0.32
against 0.96) while saturating at the same u. `--predict-dscale 0.69 1.0` shades the band
between the equilibrium weight and the 0.69 D_eq the 0.6 eV runs actually carry, which is
the honest width of what is not yet known; the shape inside the band is the same.

Two tooling repairs came out of running this on a box that restarts: `read_rt` now drops
the duplicated rows a checkpoint resume leaves at the seam (19 and 24 rows in two of these
fields) and `read_energy_delta` accepts an energy file with no header, which is what a run
whose FIRST launch was a resume produces — the accumulator is restored from the checkpoint,
so the last row is still right and only the earlier rows are missing. The existing CPTP tests could
not have caught this: trace was always conserved, populations always stayed in range,
and transfers always went to energy-matched partners. A dissipator can be a perfectly
valid CPTP map and still have the wrong fixed point.

**7.13 A doped sheet is not nstate-converged the way an intrinsic one is.** The
pure-gauge restoration (§7.1) makes the *undoped* filling basis-independent to 10⁻⁶ in
T, because it subtracts the adiabatic ground-state current of exactly the same
truncated Hamiltonian. Carriers added on top of that reference get no such
subtraction. Repeating the 147², E_F = 0.6 eV scan at nstate = 4 (dt = 0.1 fs, energy
ledger 5×10⁻⁹ eV/cell, nex baseline identical to nstate = 2 — no leakage):

| E₀ [kV/cm] | u | T (nb = 2) | T (nb = 4) | Re σ/σ_univ (nb = 2) | (nb = 4) |
|---|---|---|---|---|---|
| 1 | 0.01 | 0.70004 | 0.73477 | 26.92 | 23.77 |
| 60 | 0.74 | 0.67629 | 0.70975 | 29.00 | 25.95 |
| 100 | 1.24 | 0.68124 | 0.71716 | 28.69 | 25.46 |

An 11 % change in the sheet response from the basis alone, at 1 kV/cm, where the
vector potential (6×10⁻⁴ a.u.) cannot mix anything — so this is not field-induced.
The natural reading is the second-order repulsion of the π* level by the bands above
it, ΔE ∝ A², which is precisely a shift of the Drude weight.

**The control that pins it.** Same mesh, same pulse, same field, intrinsic filling:

| 72², 1 kV/cm | nb = 2 | nb = 4 |
|---|---|---|
| intrinsic, T | 0.999996 | 0.999996 |
| doped E_F = 0.6 eV, T | 0.70961 | 0.74452 |

The undoped sheet is basis-independent to six digits, exactly as §7.1 promises; the
doped one moves by 3.5 points. It is not a general basis insufficiency — it is
specifically the carriers the doping adds on top of the reference that the
restoration subtracts.

**How many bands are enough** (72², E_F = 0.6 eV, 1 kV/cm):

| nstate | 2 | 3 | 4 | 6 | 8 (dt = 0.05 fs) |
|---|---|---|---|---|---|
| T | 0.7096 | 0.7443 | 0.7445 | 0.7298 | 0.7473 |
| Re σ/σ_univ | **25.81** | 22.62 | 22.60 | 23.11 | 22.34 |
| D_spec [eV] | 0.586 | 0.526 | 0.525 | 0.528 | 0.521 |
| D_spec/E_F | **0.976** | 0.876 | 0.875 | 0.881 | 0.868 |

The whole error is in the 2 → 3 step. Everything from 3 to 8 agrees to ±1.7 % in σ and
±0.7 % in D, while nstate = 2 sits 14 % above all of them.
**`nstate = 2` is the outlier, not a converged choice, and production is now
`nstate = 4`** — within 1.2 % of the nb = 8 answer and the largest basis still clean
at dt = 0.1 fs (§7.2). `make_inputs.py` defaults to it and the four shipped sets were
regenerated on 2026-09-05.

Two consequences for numbers quoted from earlier nb = 2 runs. **(i)** The converged
Drude weight is ≈0.87 E_F, not E_F: the agreement with the analytic Dirac value at
nb = 2 (§7.12) was coincidental, because near K the 2×2 EPM block *is* the Dirac
Hamiltonian and returns the ideal-cone answer by construction. **(ii)** Ratios,
shapes and field dependences at fixed nstate are unaffected — the dip of §7.12 is
3.41 % at nb = 4 against 3.39 % at nb = 2 — so the conclusions about the *shape* of
T(E₀) stand; only absolute conductances move.

## 8. What to look for at production (147², DAST)

- `coh` 1 kV/cm: T ≈ 1 (empty sheet, THz far below the resolvable interband
  window); the pair density after the pulse follows the analytic Γ ∝ E^{3/2}
  table of §1 only from ~30 kV/cm up (mesh-resolved LZ tube).
- `coh` 100 kV/cm: pairs ~10¹² cm⁻² after the pulse; absorption a few per cent
  (creation energy only); `R` small.
- `diss`/`mem` 100 kV/cm: the created carriers are accelerated to ~0.5 eV and
  cooled by optical-phonon emission — the ledger's absorption grows to tens of
  per cent; T_e in `*_sbe_te.data` rises to thousands of K during the cycle and
  decays toward 300 K on the phonon timescale (tens of fs); the Rana channel
  *multiplies* carriers while n < n_i(T_e) and recombines once T_e has dropped.
  `mem` differs from `diss` by the removed dephasing-ionization share (`wiki/10`
  §8.7 logic) and by the T_e-consistent balance.
- Reflection stays ≪ absorption (R ~ (Z₀σ/2)² — a few 10⁻³ even for σ ~ 5 σ_univ).

## 8a. Why the energy in this exercise is read where it is (2026-09-15)

`wiki/06` (addendum 2026-09-15) records three faults that made every Si convergence
series in this repo diverge. Two of them are about *where* and *whether* an absorbed
energy can be read at all, so they are worth stating here — this exercise is safe from
both, but by construction rather than by luck.

**Compact support, checked.** $W(t) = -\int \mathbf{E}\cdot\mathbf{J}\,V\,{\rm d}t$ is
the absorbed energy only after the drive stops; before that it is mostly polarisation the
field has lent the crystal and will take back. The DAST proxy used here
(`DAST_E*kVcm.txt`) has support **0 → 284.4 fs** (|E| under 1e-4 of peak past 282.6 fs)
against a run window of **0 → 384.4 fs** (`nt = 3844`, `dt = 0.1`). That leaves a 100 fs
field-free tail, `transmission.py` takes the **final** `Eall - Eall0` from
`*_sbe_rt_energy.data`, and its fluence integral and energy ledger are cross-checked to
25 %. All correct.

*What would break it:* regenerating or extending the field file. The Si exercise (x15)
was fighting a DAST file whose support ran to **3274 fs** inside windows that ended at
129 and 491.7 fs — every energy taken there was mid-pulse and none of the resulting
series could converge. If you swap the drive, check its support against `nt*dt` first.

**The resolution floor.** In `*_sbe_nex.data`, column 2 minus column 3 is the drift of
the total trace — the solver's own noise floor, free with every run. On a driven
6600-step Si run it reaches $\sim10^{11}$ cm$^{-3}$. Any residue below it has measured
nothing. Use it on the `dark` control in particular.

**The `dt` × basis artifact does not apply here — but only because `nstate` is 4.**
Measured on Si: at `nstate` = 36 a 0.05 fs step manufactures $4.5\times10^{12}$ cm$^{-3}$
of carriers, and halving the step removes 99.98 % of them, while at `nstate` = 28 the
same `dt` change moves nothing ($-4.8\,\%$). The artifact needs a large basis *and* a
coarse step together, which is why a one-knob-at-a-time scan cannot see it. The 4-band
Dirac basis here has no stiff high manifold to leak into. **If you raise `nstate`, redo
the `dt` check** — and redo it at the highest field, not the lowest, since the required
step falls as the field rises. `dt` = 0.1 fs here has not been verified against 0.05 at
100 kV/cm.

## 9. Limits recorded (not blockers for this study)

1. **Initial state at T = 0, undoped** (`gs%occup` = integer filling): no thermal
   / doping Drude background, hence no *bleaching* channel; the FD(E_F, T)
   initial occupation is the next increment (it also generalizes the dressed
   reference to fractional filling).
2. Below ~30 kV/cm at 147² the LZ tube is thinner than the mesh spacing — the
   pair creation is then analytic rather than mesh-resolved.
3. Hot phonons are not included (fixed-temperature bath).
4. The Coulomb sector and the T_e fit assume quasi-thermal branch distributions.
5. The sheet field is free-standing; a substrate index enters the boundary
   condition trivially (`wiki/12` Eq. 4) but is not yet a driver option.

## 10. Tests added by this exercise

`test_graphene_dirac_levels.py`, `test_sheet_transmission.py`,
`test_rana_saturation.f90`, `test_colmem_2d.f90`, `test_dirac_te_fit.f90`,
`test_vg_sumrule.f90` (velocity-gauge f-sum rule and the pure-gauge restoration,
§3a), `test_doped_drude.py` (Fermi-Dirac occupation on a k-mesh and its
resolution requirement, Dirac-cone Drude weight and its weak temperature
dependence at fixed density, measured transmission → sheet conductance,
current-saturation field scale, §3b/§7.9) — `python3 tests/run_all.py`, 31/31.

Analysis scripts: `transmission.py` (T/R/A, sheet BC, energy ledger, Re σ),
`saturation_check.py` (populations, Rana ledger, T_e), `sumrule_check.py`
(basis f-sum rule, residual reactive current), `plot_levels.py` (level
populations vs t and k-maps), `drude_check.py` (doped sheet: D, τ, mean free
path, experiment conversion), `field_scan_plot.py` (the T(E₀)/σ(E₀) figure),
`run_field_scan.sh` (the whole scan-and-plot pipeline), `plot_occupation.py`
(the initial level occupation of a doped vs an undoped sheet on the actual mesh:
sheet density, k_F, partially occupied k-points — the pre-flight check for a doped
run), `drift_saturation.py` (the universal drift-saturation curve σ_eff/σ_lin
against A₀/k_F, continuum against any mesh — the diagnostic that separates the
physical brightening from the discretization bump, wiki/12 §4a.5),
`drift_fit_plot.py` (one layer against the analytic law, with the Drude weight and τ
separated so drift saturation and scattering can be told apart, §7.14),
`channel_budget.py` (the per-channel energy ledger: how much of the absorbed work goes
to the lattice, and which channel is pumping at zero field, §7.15).
