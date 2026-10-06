# Shared EXX MLWF controls

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

HSE06 and PBEh40 use the same localization controls in `&functional`:

```fortran
  exx_mlwf_interval=5
  exx_mlwf_maxiter=100
  exx_mlwf_tolerance=1d-7
  exx_mlwf_radius=0d0
  exx_local_fft='auto'
```

The unchanged defaults are interval 10, maximum iterations 200, tolerance
1e-6 and radius 0. Interval counts changed-source exchange refreshes, not ionic
or electronic time steps. Tolerance is always in atomic units (bohr squared).
The three old names `hse_mlwf_interval`, `hse_mlwf_maxiter`, and
`hse_mlwf_tolerance` remain compatible inputs. Supplying both names with different
values is an error; matching values are accepted. Logs use the canonical names.

## Meaning of the support radius

`exx_mlwf_radius` uses the length unit selected by `unit_system` (bohr in `a.u.`,
angstrom in `A_eV_fs`); it is logged in bohr. A positive **R takes precedence**
over `exx_mlwf_norm_fraction`. The fraction is then only a diagnostic target;
0 or omission uses 0.999 as the warning target. If any source retains less than
the target, `WARNING EXX fixed radius retains less than target` is printed once
per source team, with R, the minimum retained squared norm and the target.
The code does not enlarge R or stop solely because this target is missed.

```fortran
  exx_mlwf_radius=5d0           ! R in the selected input length unit
  exx_mlwf_norm_fraction=.999d0 ! warning threshold when R > 0
```

With R=0, a positive fraction selects the adaptive norm-based radius; with both
zero, sources remain full support. A positive radius applies the existing LCFO periodic sphere-mask geometry to each
occupation-weighted Wannier source `Q=Psi sqrt(f/2) U`, about its periodic center
in the Born–von Karman supercell (fragment supercell for DC). Both source
appearances in the Fock action use this mask. Targets are not masked and sources
are not renormalized. Exactly zero pair densities are skipped, and compact
support accelerates pair formation and output accumulation through the existing
Wannier backend. With `exx_local_fft='auto'` (default), compact sources use a padded local-box
FFT whenever its volume is smaller than the global grid. `'off'` retains the
full-grid FFT for comparison. The local kernel is sampled from the inverse FFT
of the original global multiplier, including G=0; it is not a new local periodic
Coulomb model. Padding at least 2*m-1 along every occupied box axis makes the
restricted convolution equal to the global discrete action, up to roundoff.
The source mask remains the only additional approximation in this comparison.

Periodic boxes are unwrapped across the largest empty gap on each axis. Broad
sources fall back to the global FFT. The local path currently runs serially per
pair and caches one box shape; the global path keeps its OpenMP/batched FFTs.
`EXX_FFT` logs local/global pair counts and actual/full pair-grid points. Those
counts exclude kernel/plan setup, source gathering and other SCF work. Global
source/target arrays and an O(G) periodic kernel are still stored. Dense target
states can still overlap every source: this does not establish overall linear
scaling, nor does it enable the missing PBEh optical/DC-MD adapters.

Centers use circular density moments. As in the existing LCFO path, a source
with normalized moment below 0.1 on any axis stays uncut because its center is
ambiguous. Zero-norm sources are also protected. `EXX_MLWF` output reports the
radius, protected count, total discarded norm fraction, and maximum discarded
fraction of any source at localization updates. A large radius reproduces full
support. Compare several radii against radius 0 for energy, density and spectrum;
a small norm loss alone is not an exchange-energy error bound.

This is an explicit, gauge-dependent source approximation. It preserves the
Hermiticity of the exchange action and reuses ACE, but it does not establish a
variational energy functional or consistent ionic forces. Positive radius is
currently limited to static DFT (including explicitly enabled DC MLWF), with MD, relaxation,
real-time propagation, restart and legacy snapshot export rejected. Any force
values printed by a static finite-radius run are not validated for use in MD or
optimization. Unconverged localization is reported during iteration and does not silently
disable the requested mask. If a source is actually truncated and the last
spread minimization was not converged, the final SCF result is rejected, even
when the density criterion passed. Increase localization iterations or use full
support; density convergence alone cannot certify a gauge-dependent cutoff.
Finite-radius initialization uses physical k indices rather than k-rank layout
for random seeds. This removes an avoidable MPI dependence of the initial state,
but does not guarantee a unique localization minimum or radius convergence.

The separate `pbeh_coulomb_radius` truncates the Coulomb **interaction kernel**,
not orbitals. Changing it is a different approximation. The existing LCFO-RT
`hse_lcfo_wf_radius` remains a separate control in the input length unit (internally bohr); it is not an alias for
`exx_mlwf_radius` and is unaffected by this SCF extension.

## Adaptive support on the spatial mesh

For Gamma-point native mesh calculations, use:

```fortran
&functional
  exx_mlwf_norm_fraction=0.999d0
  exx_mlwf_radius=0d0
  exx_mlwf_interval=5
  exx_mlwf_maxiter=100
  exx_mlwf_tolerance=1d-7
  exx_local_fft='auto'
/
```

`exx_mlwf_norm_fraction` defaults to 0 (disabled); 0.999 retains at least
99.9% of each occupation-weighted source's squared norm. Each periodic sphere
has its own radius, determined from local grid rows and small reductions.
Sources with ambiguous circular centers remain uncut. There is no source
renormalization and no cutoff on arbitrary target wavefunctions. A value of
1 is the full-support reference through the same spatial backend. A positive
fixed `exx_mlwf_radius` overrides the adaptive radius; the fraction then serves
only as a retained-norm warning threshold.

The masked source appears in both factors of the exchange action. Compact
convolution uses the same discrete global kernel (including G=0), but stores
only its distributed grid rows and gathers compact displacement/source/target
tiles. Target work is batched in four columns; no full wavefunction grid is
gathered. Broad sources fall back to pencil FFTs. `exx_local_fft='off'` is the
full-grid convolution reference for the **same masked sources**.
`EXX_ADAPTIVE` reports retained fraction, maximum radius and norm loss, and
local/global pair counts with local FFT point counts for orbital group zero.
Only exactly zero pair densities are skipped. Dense ACE/gauge algebra and
source-target pair work remain; this is not a claim of linear overall scaling.

Without PBE preconvergence (see below), in SCF full support is used until the selected convergence residual is below
`sqrt(threshold)`. Masking then requires converged MLWF localization. This
allows occupied-only calculations without adding empty states just to measure
a gap. Mode changes reset the density mixing history and invalidate that
iteration's convergence result. Failed gauge transport or localization restores
full support; a result without established adaptive localization is rejected.
With explicit `yn_exx_dc_mlwf='y'`, in DC all source geometry and exchange convolutions use the buffered fragment
cell and its communicators. SCF switching readiness is synchronized across
fragments so global mixing history is reset consistently.

Adaptive support is also admitted for fixed-ion native mesh RT after DC
initialization, or after conventional GS for PBE0/PBEh40/PBEh40+rVV10,
with the existing Taylor4/ACE predictor-corrector. See
[conventional GS to RT](conventional-hybrid-rt.md) for checkpoint requirements. RT starts localization
immediately, without a gap or electronic-temperature criterion, and recomputes
support on exchange refreshes. This is distinct from propagating in an LCFO
basis. Moving ions, ionic relaxation, retained-LCFO RT, restart and snapshot
export remain rejected for adaptive support. The fixed-radius legacy option
keeps its previous static-only contract. As with existing finite support,
energy/force variational consistency and long-time optical accuracy are not
established by norm retention alone.


## HSE pair screening (experimental, opt-in)

```fortran
&functional
 xc='hse06'
 exx_mlwf_norm_fraction=0.999d0
 exx_pair_screening='diagnose'
 exx_pair_tolerance=1d-6
/
```

`exx_pair_screening` is `off` (default), `diagnose` (evaluate candidates but retain
all pairs), or `on` (attempt omission). `exx_pair_tolerance` defaults to zero and
is always in atomic units, independent of `unit_system`. It budgets the grid-
weighted Frobenius norm of the **raw exchange action on the current occupied
columns**, before the exchange mixing fraction. It is not an energy tolerance,
a force tolerance, or a bound on arbitrary-target ACE propagation error.

Supports HSE06 (positive omega), PBEh(40) and PBEh(40)+rVV10 with
`exx_mlwf_radius=0` and `exx_mlwf_norm_fraction > 0`; use fraction 1 for
a full-source-support reference. Existing Gamma, native-mesh, static-SCF/fixed-
ion-RT and restart/snapshot restrictions apply. Canonical DC fragment SCF does
not enter this localized pair path. No memory ceiling is introduced.

The estimator uses the actual discrete functional multiplier (including G=0), source
amplitudes and localized pair densities. It does not use erfc of MLWF center
separation. At most half the budget is assigned to omitted pairs; a measured
Hermitian-completion correction must fit the remaining budget. Strict ACE
Hermitian/positive metric checks remain unchanged. If correction or ACE checks
fail, that refresh is recomputed without pair omission.

`EXX_PAIR mode/candidates/skipped/action bound/max rank CPU seconds` reports
attempted screening. CPU time measures only the estimator, not rotation,
completion, ACE, FFTs, MPI wall time or total speedup. The following
`EXX_PAIR unscreened ACE fallback` line says whether the attempt was discarded
and gives the accepted action bound (zero for diagnosis or fallback). Bounds
are relative to **the same source mask** and subject to floating-point rounding.
A finite WF mask introduces a separate approximation.

Before grid products, a sparse catalogue of target block-amplitude maxima is
queried over each source's nonzero support bounding box. A pair absent from the
query has `||q C(q* t)||₂ <= max|q| max|K(G)| ||q||₂ max_support|t|` within its
allocated action budget. The block width is an indexing choice, not a physical
cutoff. Periodic boundary boxes can be conservative (retain extra candidates).
Targets are rotated to the transported MLWF gauge but are not additionally
truncated. The remaining tighter pair norms use nonzero source grid rows only.
There is no dense source-by-target candidate table. Delocalized states may still
produce a dense catalogue; sparsity is measured, not assumed.

`EXX_PAIR generated grid products/catalogue entries` reports expensive pair
norm evaluations and stored sparse envelope entries (rank-boundary duplicates
can occur). `EXX_PAIR evaluated product points` reports actual source-supported
points used in those pair norms. These are attempted-screening diagnostics;
the following fallback line must be checked before claiming accepted savings.
Compact FFT target/action work contains surviving columns only, with one pair
per available spatial worker per batch. Native RT wavefunctions and ACE factors
are still dense mesh arrays: this change does not make total storage linear.

In RT, an already accepted MLWF gauge is now transported and retained if later spread
minimization fails. The failed attempt still appears in the status log, with a
separate `retained accepted transported gauge` line. For fraction < 1, RT stops
if it cannot obtain an accepted initial/transported gauge; it no longer
silently switches to a full-support operator. This does not freeze the adaptive
radius or remove errors from moving mask boundaries.

The Gamma-only spatial and spatial/orbital EXX paths share the existing Gamma
Jacobi localizer. It consumes their six periodic overlap links directly, avoiding
an extra copy and the long-axis convergence stalls of the generic gradient
minimizer. The MLWF tolerance retains its gradient-norm meaning; gauge transport
and accepted-gauge retention follow the same rules above. Canonical DC sources
do not invoke this localizer.


## PBE preconvergence before hybrid SCF

```fortran
&functional
  xc='pbeh40_rvv10'  ! also hse06 or pbeh40
  exx_pre_scf_threshold=1d-4
  exx_pre_scf_steps=3
  exx_mlwf_norm_fraction=0.999d0
/
```

`exx_pre_scf_threshold=0` (default) disables this stage and retains the existing
SCF route. A positive value is a threshold in the **selected density residual**
(`rho_dne`, `norm_rho`, or `norm_rho_dng`), not an energy tolerance; the example
is not a universal recommended accuracy. `exx_pre_scf_steps` defaults to 3 and
requires that many consecutive residuals below the threshold. A failed residual
resets the counter. Switching does not use a gap criterion.

PBE exchange and correlation are used from the first potential evaluation;
EXX, ACE, MLWF and rVV10 are inactive. At readiness, the requested hybrid
functional (including rVV10 when requested) is restored, the local potential and
first EXX/ACE are built, and mixing history is restarted. No pre-stage ACE/gauge
exists to reuse. Subsequent hybrid MLWF updates reuse the transported U normally.
For whole-system SCF or explicitly enabled DC MLWF, adaptive support becomes
eligible at this transition, still requiring successful
localization; otherwise the established full-support fallback applies. The final
`threshold` belongs to the target hybrid SCF. `nscf` counts both stages, and an
unfinished PBE stage or unconverged final hybrid stage is rejected before final
GS/restart/LCFO export. Logs explicitly mark the stage and transition.

In DC, PBE and hybrid operators act on each buffered periodic fragment. The
existing global Hartree update is preserved; the transition readiness reduction
synchronizes all fragments so density mixing resets on the same iteration.
Neither stage introduces a global EXX or global rVV10 evaluation. DC staged SCF
requires a positive specified electronic temperature (use `temperature_k=300d0`
in `&system` for the water workflow). Both stages use that same Fermi–Dirac
occupation rule and a common chemical potential satisfying the total core-weighted
charge. Current Ritz states and eigenvalues are updated before assigning
occupations in both stages. No electronic temperature is imposed on laser RT.

This option currently requires a fresh static SCF with density mixing. Restart,
checkpoints/shutdown timers, ionic optimization/MD, RT, and eigen/MLWF diagnostic
snapshots are rejected because they do not encode the temporary stage. Normal
final GS restart output is supported. For subsequent native mesh RT, use that
completed hybrid GS and omit `exx_pre_scf_threshold` from the RT input.

`exx_mlwf_tolerance` remains independent: it controls the localization gradient
norm, not the density residual or physical observable error. Its existing default
is **1d-6**; examples/benchmarks explicitly using **1d-7** do not impose it as a
requirement. This change does not normalize or loosen that criterion automatically.

## DC fragment SCF: full support without MLWF (default)

`yn_exx_dc_mlwf='n'` is now the default. DC already limits the exchange domain
to each buffered periodic fragment. Its SCF therefore skips localization,
gauge seeding/transport, and orbital masks; ACE and the existing spatial/orbital
MPI distribution remain active. The canonical sources are weighted by
`sqrt(occupation/2)`, preserving the finite-temperature density matrix exactly.
The spatial route allocates neither a gauge matrix nor a previous-gauge grid.
The full-k route retains its Fourier transform and periodic translation sum,
with no gauge rotation. Native whole-system SCF/RT is unaffected.

Use `exx_mlwf_norm_fraction=0` (default) or 1 and `exx_mlwf_radius=0` for this
DC route. Fractional masks, positive fixed radii, pair screening and Wannier
snapshots are rejected rather than silently ignored. Set `yn_exx_dc_mlwf='y'`
only to opt into the former DC localization route for comparisons. A 300 K
DC test with empty states and 0.999 support did not converge on that old route;
this does not limit the new full-fragment route. See `docs/pbe-pre-scf-ja.md`.


Fixed-radius static SCF also supports the distributed spatial/orbital backend.
Its support is enabled after PBE preconvergence (when selected), or the existing
SCF residual warmup, and after successful localization. `EXX_FIXED` reports R,
minimum retained squared norm and maximum loss. A radius enclosing the complete
periodic cell is exact and does not require converged localization. The existing
ambiguous-center protection still leaves such sources uncut. This correction
does not add fixed-radius RT/MD or fixed-radius pair-screening support.

## 局所支持軌道で構築する ACE（検証用・固定イオン RT）

`exx_ace_support='source'` は、切り詰め済み MLWF `S` と交換作用 `K_S S` で ACE を構築する。
既定値は `'occupied'`（従来方式）。`exx_mlwf_norm_fraction=0.999d0` と組み合わせる。
支持が重ならない対を誤差ゼロで除外するが、ACE の構築空間を変える近似は別に存在する。
99.9% の規格化条件を交換作用・電流・エネルギーの誤差上限とは解釈しない。

時間発展する軌道は実空間メッシュのままで、元の占有軌道に ACE を作用させて交換エネルギーも計算する。
非直交支持軌道の ACE 計量に従来と同じ Hermitian・正定値・条件数検査を適用し、
構築できなければ従来の占有軌道 ACE に戻る。各更新の `EXX_SUPPORT_ACE accepted` で確認できる。
全支持（fraction=1）では従来 ACE と一致する。

初期実装は完全占有の native 固定イオン RT に限定する。DC-SCF、LCFO 基底 RT、MD では使用しない。
`exx_pair_screening='off'` が必要（内部では支持の厳密なゼロだけを利用する）。
source モードでは交換作用 W の厳密な非ゼロ要素と計量の因子 A を保持し、
`-W A A† W†`（内積には格子体積を含む）として作用させる。密な ACE 因子は展開しない。
計量の逆行列を明示的に作らず、従来の条件数制限を維持する。
実空間軌道とゲージ運搬用軌道は密なまま、計量因子は軌道数の二乗の格納を要するため、
総メモリの線形スケーリングを保証しない。
