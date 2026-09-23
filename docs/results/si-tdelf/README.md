# Current-corrected time-dependent ELF in Si

ELF was evaluated every10 RT steps on the existing no-pump, weak and strong Si trajectories (4³ k mesh, 12³ real-space grid, 16 occupied orbitals per spin,1600 steps=3.09617 fs). This is postprocessing of the previous bare-K TDCDFT trajectories. No alpha-from-ELF mapping, feedback, new dynamics or k convergence scan was introduced.

## Definition

Use [Burnus, Marques and Gross, Time-dependent electron localization function](https://arxiv.org/abs/physics/0404126), equations1–4. For one spin:

\[
 D_s=\tau_s-\frac{|\nabla n_s|^2}{4n_s}-\frac{|\mathbf j_{p,s}|^2}{n_s},\quad
 D_{s,0}=\frac35(6\pi^2)^{2/3}n_s^{5/3},\quad
 \mathrm{ELF}=\frac1{1+(D_s/D_{s,0})^2}.
\]

Here tau has no factor1/2. n, tau, and orbital paramagnetic current use the equal k weights1/64 and one-spin occupancy1. Wavefunction gradients include the Bloch term (grad+ik)u; grad(n)=2Re sum(u* grad(u)) is evaluated consistently. Current is summed over occupied orbitals and k **before** squaring. A common vector potential cancels from tau-j²/n if included consistently, so the calculation uses canonical orbital gradients. The macroscopic SALMON current, which contains nonlocal-pseudopotential velocity contributions, must not replace this local orbital current.

ELF measures short-range same-spin exclusion/localization of the Slater determinant. It is invariant under occupied-unitary rotations and needs no MLWF minimization. It can distinguish changes in bonding structure, including the bonding-to-antibonding example in the cited TDELF paper. It is not a transport coefficient or dielectric response. ELF=0.5 refers to the cold uniform free-electron reference; hot or nonequilibrium uniform distributions need not have the same kinetic-energy density and need not yield0.5. A spatial mean of0.5 can also arise from nonuniform values on both sides of0.5.

The finite-mesh cold gas check gives ELF=0.50010526 without rescaling; see `reference.json`. This small departure is the Fermi-sphere quadrature error on the same mesh, not a fitted correction.

## Grid derivative and existing SALMON output

The main analysis uses Fourier derivatives with the odd-derivative Nyquist component set to zero on the even12-point grid, while retaining ik. This is a symmetric discrete derivative convention, not a claim of exact continuum differentiation. An alternative that retains the signed complex-FFT Nyquist component changes both tau and current. Endpoint comparison is in `derivative_check.json`:

| Quantity | Symmetric Nyquist convention (main) | Signed complex-FFT convention |
|---|---:|---:|
| Initial electron-weighted ELF |0.59969302|0.59894949|
| Final strong electron-weighted ELF |0.57410710|0.57302222|
| Initial bond-region ELF |0.92504592|0.92490378|
| Final strong bond-region ELF |0.89390216|0.89370829|

The excitation trend survives this derivative check. The main convention yields zero current correction in the initial time-reversal-symmetric state; the alternative leaves a maximum0.001186 ELF correction there. This is a coarse-grid interpolation ambiguity, not evidence of equilibrium currents. These are valence pseudowavefunction diagnostics, with no all-electron core reconstruction.

The existing `src/io/write_field.f90:write_elf` was inspected but not used or modified. In the current checkout its orbital sums omit occupation/k weighting and the explicit ik term, and its current expression sums individual squared orbital currents rather than squaring the summed current. It therefore cannot be used directly as the periodic current-corrected reference for this comparison. This issue is separate from the Nyquist convention used in the new postprocessor.

## Results at3.09617 fs

| Density-weighted diagnostic | Initial | Weak | Strong |
|---|---:|---:|---:|
| Overall mean ELF |0.599693|0.599693|0.574107|
| Bond-region mean ELF |0.925046|0.925046|0.893902|
| Outside-bond mean ELF |0.528572|0.528572|0.507099|
| Electron fraction0.45<ELF<0.55 |17.91%|17.91%|22.43%|
| Electron fractionELF>0.75 |44.64%|44.64%|26.96%|
| Mean squared deviation from0.5 |0.074899|0.074899|0.066173|

All483 snapshots completed. No-pump maximum field change is8.2e-9. The strong-pulse current correction reaches0.00897 locally and0.00435 in the electron-weighted mean, so it is retained. The signed complex-FFT derivative control gives20.77% instead of22.43% in the narrow near-half bin at the final strong time; smooth averages have much smaller convention sensitivity. No precision claim based on the bin fractions is intended.

## Spatial diagnostics and interpretation

All averages and distribution fractions are weighted by the instantaneous electron density. Bond regions are the union of fixed radius1-bohr spheres around the16 initial Si bond centers, rather than dynamically fitted ELF basins. The complement is reported separately. Histograms use bins of width0.02; the displayed near-gas interval0.45<ELF<0.55 and high-localization intervalELF>0.75 are descriptive thresholds, not metallicity criteria. On a coarse spatial grid, threshold fractions are more sensitive than smooth averages.

The strong pulse lowers both the bond-region ELF and the overall mean, while appreciable high ELF remains near bonds. Weak excitation has almost no change in these averaged measures. This supports reduced bonding localization; it does not demonstrate a fully metallic transient state or determine alpha. Real-space maps and distributions are essential because a single mean can conceal spatially distinct behavior.

Outputs: `metrics.json`, all161 records per case in JSON, all ELF and one-spin density fields plus electron-weighted histograms in `{case}_fields.npz`. `dynamics.png` shows time traces, `maps.png` the x=y section at initial/1.548/3.096 fs, and `histogram.png` the corresponding distributions. Shared color scales are used.

## Verification and reproduction

Five tests check the single-orbital ELF=1 limit, occupied-unitary invariance, consistent constant-boost invariance, analytic plane-wave derivatives, a uniform-gas ELF=.5 benchmark, and the symmetric Nyquist convention. Independent review verified spin factors, the Bloch term, summed-current correction, vector-potential cancellation and the derivative sensitivity interpretation. No-pump time traces test stability under field-free orbital phase evolution.

From repository root with NumPy/Matplotlib:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s docs/results/si-tdelf -p 'test_*.py'
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 docs/results/si-tdelf/reference.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 docs/results/si-tdelf/check_derivative.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 docs/results/si-tdelf/analyze.py
```

Requires the ground-state restart and existing `/private/tmp/salmon-si-time-wannier-dense` checkpoints. Their regeneration instructions are in the neighboring `si-time-wannier` record. All quantities are instantaneous; no temporal averaging is introduced.
