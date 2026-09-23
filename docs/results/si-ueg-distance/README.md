# Si: distance to a homogeneous electron-gas reference

Offline evaluation of the user-approved **unclipped** model

\[
 d(t)=\|P(t)-F(t)\|_{\rm HS}^2,\qquad D(t)=d(t)/d(0),\qquad \alpha(t)=0.2D(t).
\]

P is the propagated one-spin occupied density matrix; F is a specified homogeneous free-electron ensemble. Equal spin multiplicity cancels in D. This model is a hypothesis, not a dielectric-response derivation. No production Fortran or feedback was changed. The analyzed trajectories remain the previous alpha0=1, K0=0, instantaneous bare-K polarization-closure runs, not trajectories propagated with this proposed alpha.

## Reference and implementation

Use the existing Si8 12³ real-space grid and shifted4³ k mesh, 16 occupied orbitals at each k, through1600 RT steps (3.09617 fs). Evaluate at steps0,10,...,1600: 0.0193511 fs between snapshots. No k convergence scan. The supercell has110592 plane-wave modes and1024 one-spin occupied states. Fill the free-electron energies |q|²/2 in ascending order, with equal fractional occupancy of a degenerate final shell. The96 partially occupied states give rest-reference Tr(F²)/Tr(F)=0.977213542.

The primary reference matches the SALMON electron-number current: target UEG canonical drift is J/n - a_total, n=32/L³ and a_total=Ac_tot+Ac_xc. The ordinary rt.data Ac_tot is classical only. Both arrays are sampled at the same checkpoint time. At the initial state the macroscopic current is set to zero. This is a **current-equivalent UEG drift**, not necessarily the spectral mean momentum of the SALMON orbitals: SALMON includes finite-difference and nonlocal-pseudopotential velocity terms.

On this finite momentum grid, implement a noninteger drift as a convex mixture of the8 surrounding integer-shifted Fermi seas. Electron number and mean momentum are conserved (maximum drift error below6e-16 au); occupations remain in[0,1]. Momentum wrapping is rejected. This is a legitimate uniform ensemble, but not an exact continuously boosted zero-temperature Fermi sea or a fixed-flow energy-minimum state. Its purity depends on the fractional shift, which affects the distance. There is no added electronic-temperature parameter.

Two controls use the actual FFT mean canonical momentum or a resting reference. Each control is normalized by its own initial distance, so all begin at alpha0.2. The initial spectral momentum has a small finite-grid asymmetry (~3.35e-4 au), while physical current is effectively zero. This motivates retaining the actual-current and spectral conventions separately. A common electromagnetic gauge is used throughout; arbitrary continuous-boost invariance of the finite-grid ensemble interpolation is not claimed.

With Fourier coefficients C_knG, P_k=C_k C_k†. Compute

\[
 d=\frac1{N_k}\sum_k\left[\mathrm{Tr}(C_k^\dagger C_k)^2+\sum_G f_{kG}^2-2\sum_G f_{kG}\sum_n|C_{knG}|^2\right].
\]

Thus no full real-space density matrix or new MLWF optimization is required. The result is invariant under occupied-band unitary rotations, so computing it in a Wannier gauge gives the same value. All483 snapshots took about26 s for the final run; this is an observed local postprocessing time, not an integrated runtime benchmark.

## Results

The primary initial squared distance per primitive cell and spin is7.0499377113.

| Case | Minimum alpha | Maximum alpha | Final alpha | Final D |
|---|---:|---:|---:|---:|
| No pump |0.199999998|0.200000000|0.199999998|0.999999988|
| Weak, 1e8 W/cm² |0.19991734|0.20000000|0.19994890|0.99974452|
| Strong, 1e13 W/cm² |0.18771176|0.21855607|0.20244860|1.01224302|

Strong minimum occurs at0.73534 fs; maximum at1.87705 fs. At the maximum, mean MLWF spread is5.54526 Å² compared with2.05768 initially; invariant spread is4.65632 compared with1.88561. At the final time, these are4.93676 and4.35097 Å². **D>1 is not evidence of increased localization.** The total density's relative spatial RMS also decreases from0.66830 to0.64647 while this distance rises slightly.

For the strong trajectory, final alpha is0.20144370 using spectral-momentum matching and0.21650900 with the resting reference. Maximum primary-versus-spectral alpha difference is0.00324159; primary-versus-rest difference is0.02224524. Reference choice has material effects.

The primary initial→final distance decomposition is:

| Term | Change |
|---|---:|
| Occupied purity Tr(P²) |−3.3e−8|
| Reference purity Tr(F²) |−0.539248|
| Minus twice overlap −2Tr(PF) |+0.625560|
| Total distance |+0.086313|

The decreased overlap with the cold-gas reference dominates the reference-purity reduction. The resting-reference distance also increases, so the increase is not solely an interpolation artifact. The metric detects differences in the electronic state, including excitation-induced changes in momentum occupation, rather than isolating real-space localization or screening. Alpha above0.2 is an allowed outcome of the specified model; these data do not establish it as a physical enhancement of excitonic attraction.

## An obstruction to reaching the zero endpoint

The time-evolved Slater determinant retains rank16 at every k, while the reference electron counts vary from14.5 to17 per k and include fractional occupations. Even the closest rank16 projector has a positive distance to this particular reference. For reference eigenvalues sorted descending,

\[
 d_{\min}=\frac1{N_k}\sum_k\left[16+\sum_G f_{kG}^2-2\sum_{i=1}^{16}f_{k,i}^{\downarrow}\right].
\]

The resting-reference lower bound is0.46875, equivalent to alpha>=0.013298 with this normalization. The current-matched final reference gives alpha>=0.013948. This is not a clipping rule; it is a representability bound for the chosen trajectory class and reference. Algebraically P=F still gives alpha=0 (tested), but the specified finite-grid pure occupied subspace cannot reach that mixed reference. Consequently, this experiment does not demonstrate a physically accessible metal limit alpha→0.

## Decision and reproduction

The unclipped distance model was implemented and evaluated correctly. Weak excitation remains near0.2. A unique relationship between this distance, localization and screening was **not** established, so production feedback remains unchanged. Refining the reference/metric to separate excitation from spatial nonuniformity is a subsequent modeling question.

Seven tests cover Fermi count/symmetry, convex boost number/current/positivity, rejected boundary crossings, explicit-matrix HS agreement, occupied-gauge invariance, FFT normalization and momentum indexing, pure/mixed UEG zero endpoints, fixed-rank bound, integer-boost invariance and alpha>0.2. Independent review checked the actual full momentum mesh and all483 three-term decompositions (maximum residual7.1e-15).

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s docs/results/si-ueg-distance -p 'test_*.py'
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 docs/results/si-ueg-distance/analyze.py
```

Requires the ground-state restart and `/private/tmp/salmon-si-time-wannier-dense/{none,weak,strong}` wavefunction checkpoints, plus the stored MLWF and current data in the adjacent `si-time-wannier` directory. The previous Wannier result record contains checkpoint regeneration instructions. JSON records retain all sampled distances, alpha values, reference purities, overlap terms, lower bounds and matched-flow errors. No new TDDFT runs were needed.
