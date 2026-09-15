# Exercise x15 — Si under a single-cycle THz transient: induced transparency

The bulk counterpart of `exercise_x14` (graphene). One silicon EPM→SBE pipeline per
field strength, driven by the **measured DAST single-cycle terahertz transient —
reproduced analytically, with no field file** — and dissipated by the production
channel set. The question is the experimental one: *does silicon bleach or darken as
the THz field is raised, and by how much?*

---

## 1. Run it

```bash
python3 make_runs.py --dt-scan          # writes runs/<tag>/ , nothing else needed
cd runs/E1000kVcm
salmon < Si_prim_epm_gs.inp             # step 1: ground state, writes Si_prim_{k,eigen,tm}.data
salmon < Si_prim_sbe_rt.inp             # step 2: the real-time SBE
```

**Every run directory is self-contained and independent.** That is not cosmetic:
the ground state writes `Si_prim_k.data`, `Si_prim_eigen.data`, `Si_prim_tm.data`
into the *current* directory and the SBE step reads them from there, then writes
`Si_prim_sbe_*.data` back. Two runs launched in one directory silently overwrite
each other's ground state **and** each other's output. Keep one run per directory
and they can all be submitted at once.

`make_runs.py --help` lists the knobs (`--fields`, `--nk`, `--nstate`, `--dt-fs`,
`--t-end-fs`, `--no-ring`, `--dt-scan`).

**If a job is killed, restart with `Si_prim_sbe_rt_resume.inp`, not the original.**
Every run checkpoints each 200 steps, but the plain input reopens the outputs with
`status='replace'`, so restarting with it begins again from t = 0 — the checkpoint
gets written and never read, which is the worst of both worlds. The resume twin adds
`yn_sbe_checkpoint_restart='y'` and continues; `run_lomonosov.sbatch` picks it
automatically when output and a checkpoint are both already present.

## 2. The drive: no field file

The DAST transient used throughout this fork is a single cycle of E at ≈3.5 THz,
which is the time derivative of a single *hump* in A. That is exactly what SALMON's
analytic `Acos2` shape produces when the carrier-envelope phase is set so the
carrier is a cosine and the carrier period equals twice the envelope support:

```
A(t) = -(f0/w) cos^2(pi*tt/tw) sin(w*tt + 2*pi*phi_CEP),   |tt| < tw/2

  phi_cep1 = 0.75          ->  sin(...) = -cos(w*tt)
  omega1   = hbar*pi/tw    ->  the cosine runs 0 -> 1 -> 0 across the whole support
  tw1      = 273 fs
```

Measured against `x12/DAST_singlecycle_100kV.txt`:

| | analytic | DAST file |
|---|---|---|
| corr(A) | **+0.9994** (1.8 % rms) | — |
| corr(E) | **+0.9952** (6.3 % rms) | — |
| spectral centroid | 3.81 THz | 3.83 THz |
| spectral peak | 3.50 THz | 3.38 THz |
| peak \|A\| | low by 2.6 % | — |

The envelope makes `peak|E| = (2/sqrt3) f0` exactly, so the generator sets
`E_amplitude1 = E0 * sqrt(3)/2` and the requested peak field is hit to 5 digits
(verified in the solver: 1000 kV/cm requested → 1000.0 kV/cm in `*_sbe_rt.data`).

The residual 6.3 % in E is the measured waveform's asymmetry, which the symmetric
analytic form cannot carry. If you need the measured asymmetry, point
`ae_shape1='input'` at the x12 file instead; nothing else changes.

## 3. Two settings that are not optional

**`yn_sbe_vg_sumrule = 'y'`.** In the velocity gauge the current carries a
diamagnetic term for every electron of the filled band, cancelled by the interband
response only in a complete basis. What survives is `eta * N_e * A / V`, and since
`A = E/omega` it is negligible in the near-IR and **decisive at THz**. For Si on
9³ at nstate = 20 the raw absorbed energy is 70× the restored value — 98.6 %
artifact. It **cannot be repaired afterwards**: eta is not a constant, it runs
0.23 % → 7.35 % along this pulse, so subtracting a fixed `eta*n_e*A` from an
archived J(t) still leaves the absorbed work 53× too large. See
`wiki/06` §"Addendum (2026-09-07)".

**The e-ph detailed-balance split** now acts for every material. The emission/
absorption split must be taken at the energy a pair actually transfers, not at the
phonon energy, whenever the Gaussian matching width is not small against the phonon
energies — and for silicon `sigma/hbar*omega` runs 3.2…20. With the historical
split a Si dark control on 9³ generates 1.32e12 cm⁻³ of electron–hole pairs across
the gap **with the drive switched off**; with the corrected split the same control
is identically zero. See `wiki/12` §6a scope note.

## 4. Converge dt before you converge nstate

`--dt-scan` emits `dtscan_{0.1,0.05,0.025,0.0125}_E1000`. Run those **first**.

This ordering is not pedantry. `wiki/06` §6(i) measured, on this same material and
this same drive, that at an unconverged dt the carrier density *climbs* with the
band count and "converges" to a wrong, dt-inflated value — "the dt-error was filling
each newly-available band, **faking** a basis-insufficiency". At nstate = 28…36 the
Si levels reach 54 eV, which is 8.2 rad of phase per 0.1 fs step, and the S4
composition leaks there. The default here is therefore `dt = 0.05 fs`, and `wiki/06`
§7 asks for ≤0.02–0.05 fs at this field strength.

**Measured, and it is worse than §7 suggests.** Si $5^3$, 1000 kV/cm, coherent,
post-pulse plateau:

| nstate | dt = 0.05 fs | dt = 0.025 fs | step in dt |
|---|---|---|---|
| 28 | 1.7795e-12 eV/cell | 1.6945e-12 | −4.8 % |
| 36 | **1.2306e-08** | **1.9897e-12** | **−99.98 %** |
| step in nstate | **×6914** | +17 % | |

Three corners agree at the noise floor; only (36, 0.05) escapes it. Sweep nstate at
dt = 0.05 and you walk into that corner and see a band-count dependence that never
flattens; sweep dt at nstate = 28 and you walk along the safe edge and conclude dt is
already converged. **Neither one-knob scan can see this**, which is why every earlier
series here refused to settle. At dt = 0.025 the nstate dependence is gone (+17 %,
floor to floor).

At nstate = 36 a 0.05 fs step manufactures 4.5e12 cm⁻³ of carriers out of nothing;
halving the step removes 99.98 % of them and lands on the same floor nstate = 28
already sat on. The default is therefore `dt = 0.025 fs`, not 0.05.

The fake carriers are **electron–hole balanced to three digits** (4.527e12 vs
4.509e12), because the step error drives a coherent valence→conduction transfer like a
real one. So e/h agreement does not prove a residue is physical. Three checks, all
required:

1. `nelec` ≈ `nhole` — separates population transfer from trace drift;
2. the residue exceeds the trace drift (column 2 − column 3 of `_sbe_nex.data`);
3. **it survives halving `dt`** — the only one that catches this failure.

### Read the absorbed work only AFTER the drive stops

The earlier version of this section reported that the absorbed work grows ~10 % per
8 bands with no flattening, and that refining 5³ → 7³ makes that step larger. **Both
numbers are withdrawn.** They were measured on the DAST field file, whose support
runs to 3274 fs, inside a window that ended at 129 fs — i.e. *mid-pulse*. Before the
drive ends, W(t) = -∫E·J dt is not absorbed energy at all: it is dominated by the
reversible polarization the field has lent the crystal and not yet taken back. On
that window the instantaneous power changed sign 1554 times and W was still rising
at the edge. Differences of such a quantity do not converge in dt, in nstate, or in
the mesh, because the quantity itself is not yet defined. That, and not a basis or
mesh insufficiency, is why every series in the old text refused to settle.

The analytic drive used here has **compact support**: E ≡ 0 for t > `tw1`, exactly,
to the last bit. W(t) therefore goes flat and its plateau *is* the absorbed energy.
Run to `tw1` + a field-free tail (the defaults give 273 fs + 57 fs) and read the
plateau. Only then are dt, nstate and mesh comparable.

**And check the plateau against the solver's own floor before you converge it.** In
`*_sbe_nex.data`, column 2 is (tr - tr_vb)/V and column 3 is (nelec - tr_vb)/V, so
their *difference* is the drift of the total trace. That difference is the error bar on
either column. Measured here — Si, 5³, nstate 28, dt 0.05, 6600 steps, 1000 kV/cm,
dissipators off — it is 3.7e-12 electrons/cell, i.e. ~1e11 cm^-3, while the post-pulse
excitation is 4.9e8 cm^-3. Coherent below-gap Si at this field absorbs nothing this run
length can resolve, on 5³ *or* 7³ (W_plateau 1.78e-12 vs 5.82e-13 eV/cell, both ~200x
under the floor). Turning the dissipators on puts the densities this exercise is about
6-7 orders above the floor, so this is a caveat for the convergence scans, not for the
experiment. `wiki/06`, addendum 2026-09-15, has the numbers.

## 5. What to measure

`../exercise_x11_full_dissipation_showcase/thz_permittivity.py` reads
`*_sbe_rt.data` and reports the transmission and the absorbed energy; point it at a
run directory. The `dark` run is the control — subtract it, and if it is not
identically zero, stop and find out why before reading anything else.

| directory | what it is for |
|---|---|
| `E100kVcm` … `E3000kVcm` | the field scan: T(E0), the experiment's observable |
| `dark` | zero-field control; every channel must read exactly 0 |
| `dtscan_*_E1000` | the mandatory dt convergence, §4 |
