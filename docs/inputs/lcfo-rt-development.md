# LCFO RT development status

The DC input guard remains unchanged: `yn_dc=y` currently accepts `theory=dft`
only. An experimental fixed-LCFO-subspace RT adapter now reuses native
`initialization_rt`, density, Hartree, semilocal XC, pseudopotential, Taylor4
predictor/corrector and current routines. No additional SCF is performed.
Si128/16-fragment self-consistent RT has passed a four-step integration check;
there is still no converged dielectric spectrum.

## Experimental native path

Set environment `SALMON_LCFO_RT=1` and use `theory='tddft_response'`,
`yn_dc='n'`, `yn_conventional_from_dcdft='y'` in `&calculation`.
Use Gamma, unpolarized HSE06, `yn_hse_wannier='n'`, default `hse_taylor4`,
one real-space MPI rank per fragment, and `nproc_k=1`. Orbital groups are
allowed: e.g. eight fragments use `nproc_rgrid=8,1,1`, `nproc_ob=2`, MPI16.
This support applies to LCFO RT; the DC-HSE GS orbital-MPI restriction is unchanged.
Each rank's grid must coincide exactly with its fragment core. The existing
`./data_dcdft/fragments` LCFO records provide both the initial orbitals and fixed
orthonormal core bases. Native `hpsi` is projected after all Hamiltonian terms;
Hartree and density are evaluated on the whole physical system by existing code.

The current implementation requires occupied-only fixed occupations (for Si128,
256 states, 512 electrons, omit temperature). It loads the first256 saved LCFO
states without reconverging the accepted DC-to-LCFO density difference.
Restart input/output, checkpoints and time_shutdown are rejected until LCFO
basis and exchange state metadata are supported.

`hse_lcfo_rt.f90` reconstructs periodic fragment+buffer bases and density factors,
uses the existing screened FFT exchange kernel, core-weights and Hermitianizes
the projected exchange, and sums contributions once across spatial ranks.
The global trace energy is passed to the existing native energy bookkeeping.
Existing ACE is built in coefficient space; a rejected indefinite/singular
metric falls back explicitly to the full projected operator without clipping.
The predictor/corrector averages endpoint operators. By default, full support
uses density eigenfactors (not MLWFs). Its relative density eigenvalue threshold
1e-14 only removes numerical null modes, with discarded trace logged.

### Opt-in transported MLWF sources

Set `SALMON_LCFO_RT_MLWF=1`. Initial localization uses global occupied states in
the distributed fixed LCFO basis, not independently diagonalized fragment density
factors. A Gamma-specific SU(2) Jacobi optimizer minimizes the same six-link MV
functional as the existing gauge code. `hse_mlwf_maxiter` is the maximum sweep
count and `hse_mlwf_tolerance` the gradient tolerance; Si128 uses200 and1e-6.
The original general-k optimizer is unchanged. A coefficient-space pivoted trial
seed avoids gathering the complete real-space wavefunctions on one rank.

After this one localization, polar overlap transport aligns each new occupied
frame with the previous accepted WFs. This removes band dynamical phases;
merely holding the numerical U fixed would spread WFs even at stationary density.
Predictor and corrected endpoint frames share the accepted step-start reference;
trial updates are rolled back before accepting the corrected frame. Cache keys
include both coefficients and occupations. No repeated MV minimization occurs.

`SALMON_LCFO_RT_RADIUS=0` (default) retains full support. A positive value is a
three-dimensional sphere radius in bohr around fixed INITIAL WF centers.
Each xyz displacement uses the minimum image of the orthorhombic global cell.
A radius at least half the box diagonal retains full support.
Distances use the whole Si128 cell, not the shorter fragment period. Centers with
circular reliability below 0.1 on ANY axis are protected at full support. No source
renormalization, occupation cutoff or density/Hartree truncation is introduced.
Initial localization must converge before a truncated-support test proceeds.
Historical benchmark tables below used the former x-only halfwidth; they do not
measure the current spherical support.

Diagnostics `lcfo_mlwf_initial.bin` and `lcfo_mlwf_links.bin` preserve the initial
occupied coefficients, unitary/centers and optimizer inputs. They are not restart
files. The initial dump is now version 2 with centers(3,no) in Fortran order;
version 1 stored centers(no) for x only. The links dump remains version 1.
The active fragment-source count and global discarded source norm fraction
are logged. Existing zero-source skipping avoids unnecessary FFTs after masking;
no estimated speedup is inferred from source counts alone.

### Optical and continuity interpretation

Source-only masks are explicitly nonvariational. Hermitian assembly and ACE do
not by themselves restore local charge continuity. Set
`SALMON_LCFO_RT_CONTINUITY=1` to evaluate the exchange density source before LCFO
output projection, including both core-weighted adjoint halves and summing the
overlapping fragment buffers. Logs contain the signed integral and L1 integral
(electrons per atomic time) and maximum density-rate magnitude. Full-support Fock
cancellation should hold to roundoff; a masked case can conserve total charge
while producing a nonzero local defect. The additional diagnostic FFT work must
be held equal when comparing runtimes.

Native kinetic-plus-pseudopotential current remains the recorded observable.
Its masked transformed spectrum is an approximate-model diagnostic, not yet a
certified dielectric function: cutoff acceptance requires assessing the
continuity defect/field-current consistency as well as current and spectral
errors. Do not add an unwrapped-position exchange commutator in this periodic
fragment model as an ad hoc correction.

The environment switch is a development opt-in, not a new production input.
Without it the ordinary reconstruction and native RT paths are unchanged.


`src/rt/lcfo_rt_core.f90` supplies a fixed orthonormal complex basis kernel:

- `lcfo_cayley_step(H,C,dt,next,status)` solves
  `(I+i dt H/2) next=(I-i dt H/2) C`. It does not diagonalize H or normalize C.
- `lcfo_density(C,f,P,status)` forms `P=C diag(f) C†`, retaining complex coherence
  and fractional occupations. f is the physical occupation; no implicit factor2.
- `lcfo_midpoint_step` iterates a Hamiltonian callback on `(P_start+P_end)/2`
  at `time+dt/2` to the specified relative Frobenius density residual. Each trial
  propagates from the unchanged starting C. It is a small reference solver;
  fixed-point convergence is not guaranteed for arbitrary dt or functionals.
- `lcfo_grid_density` reconstructs density on a core's grid from its basis and
  the corresponding coefficient rows of ALL global occupied states.
- `lcfo_project_potential` projects one real grid potential using grid quadrature.

Outputs use separate storage from inputs. Failure retains the input C in `next`
when its shape permits; status1 invalid input,2 non-Hermitian H,3 linear solve
failure,4 fixed point not converged,5 callback failure. A caller must reject the
step on any nonzero status. No silent symmetrization of caller H is performed.

The independent Python reader `samples/dc_hse/lcfo_rt_reference.py` accepts only
Gamma, unpolarized V1 complex LCFO data. It verifies completion/footer lengths,
run identity, geometry, periodic halo mapping, orthonormal core bases and saved
eigenvectors. It reconstructs the native symmetrized matrix, retaining a separate
raw-directed anti-Hermitian diagnostic. Basis arrays use Fortran grid ordering.
The `lcfo_frozen_probe.f90` helper tests real saved LCFO eigenstates against exact
Cayley phases and time reversal; this is deliberately NOT a TDHSE driver.

Still required: longer-time/dt/LCFO-basis convergence, zero-field subtraction,
field-current assessment of the source-mask approximation, and dielectric comparisons.
The accepted initial density mismatch is retained; stationarity is not a gate
for starting RT.
A Hermitian matrix alone establishes neither energy-functional consistency nor
optical-response accuracy. In particular, the static source-cutoff scripts are
not silently promoted into a variational time-dependent functional.

Tests: CTest `lcfo_rt_core`; `OPENBLAS_NUM_THREADS=1 python3
testsuites/unit_lcfo_rt/test_reference.py`. The standalone check.py compiles
against the local Homebrew BLAS for the current development machine; CTest uses
the build's selected BLAS/LAPACK and is the portable verification path.

Native integration regression: `python3 testsuites/unit_lcfo_rt/test_native.py
--binary /absolute/path/to/salmon --pseudo /absolute/path/to/H_rps.dat`.
It runs MPI2 jobs sequentially in a fresh temporary directory: Gamma DC-SCF,
LCFO RT, half-dt RT, and rejection of unsupported restart. It checks normal
completion, finite output, endpoint-current agreement and post-impulse energy
width. This small integration test is not a production convergence criterion.

Additional regressions: `testsuites/unit_hse_wannier/test_gamma.py` (known
anisotropic Gamma solution and noncommuting links), and
`testsuites/unit_lcfo_rt/test_transport.py` (stationary-density phase invariance,
predictor rollback and corrected-state cache acceptance). These standalone tests
currently use the local Homebrew BLAS path. `test_native.py` also checks full
source parity, large-radius identity, initial-U identity across support cases,
and full/masked exchange continuity diagnostics with MPI2 jobs run sequentially.

### ACE reuse across physical time steps (2026-09-26)

`SALMON_LCFO_RT_ACE_INTERVAL=N` accepts a positive integer (default1).
For an impulse, step1 always rebuilds ACE **before the first predictor**, as
well as at predicted/corrected endpoints; later refresh steps are1+N,1+2N,….
For a smooth laser starting at zero field, initial ACE is retained and endpoint
refreshes occur at N,2N,… . The present LCFO path constructs initial ACE from
loaded GS orbitals; it does not deserialize fragment GS ACE factors.

Between refresh steps, only the nonlocal exchange operator is frozen. Density,
Hartree and semilocal XC continue to update through native routines. MLWF frames
still undergo polar transport for each predictor/corrected state. An invalid ACE
or changed occupations bypass reuse. Unchanged endpoint operators retain a
single set of ACE factors instead of doubling their rank by averaging duplicates.
Exact-coefficient cache keys are invalidated when a retained step transports WFs.

The skipped-refresh exchange energy is evaluated as0.5 sum f<C|K_ACE|C> at current
orbitals, explicitly logged as a frozen-operator trace diagnostic. It is neither
the instantaneous self-consistent HSE energy nor a conserved frozen-Hamiltonian
energy; compare its width only as a diagnostic. Local-field updates are not skipped.
The existing continuity diagnostic evaluates rebuilt exchange only and cannot
certify continuity of a retained ACE at evolving orbitals.

MPI2 regression checks default/interval1 identity, fewer builds at2/4, identical
initial MLWF data, finite output, invalid-interval rejection, impulse pre-predictor
rebuild and smooth-laser first-step reuse. Si128 wall-time/current comparisons
are recorded below once measured; a short pilot does not certify dielectric
accuracy or select a production refresh interval.

<!-- ACE_REUSE_20260926 -->
## ACEを複数ステップで再利用した測定（2026-09-26）

実装 `a9e3716`。Si128・16フラグメント、MLWF全範囲、MPI16/OpenMP1/BLAS1、16ステップ、dt=0.16 a.u.（合計0.06192 fs）。同じ初期状態・U・励起1e-4 a.u.、1計算ずつ実行。連続性診断は全条件で無効とし、診断コストを含めずに比較。

impulseでは第1ステップの最初の予測伝播前に必ずACEを再構築し、そのステップの予測・修正状態でも更新する。以後は1+N、1+2N…ステップで更新。再利用中も密度・Hartree・半局所XCとMLWFの極分解による位相追跡は継続する。滑らかなレーザーでは初段更新を省き、GS軌道から初期化時に構築したACEを再利用できる（保存済みGSのACE因子の読み込みではない）。

| 更新間隔 N | ACE構築回数 | 全実測時間 s | 全時間の高速化率 | RT部分 s | RT部分の高速化率 | E∞：全時刻の最大差 / 基準ピーク | 診断エネルギー幅 Ha |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 35 | 336.73 | 1.000倍 | 324.19 | 1.000倍 | 0.00000% | 3.73290e-08 |
| 2 | 19 | 219.16 | 1.536倍 | 203.07 | 1.596倍 | 0.31990% | 1.67213e-04 |
| 4 | 11 | 153.20 | 2.198倍 | 136.51 | 2.375倍 | 0.96453% | 4.29557e-04 |
| 8 | 7 | 119.28 | 2.823倍 | 103.02 | 3.147倍 | 2.26090% | 7.39916e-04 |

旧版と更新間隔1の最初の8ステップの電流差は最大 5.477e-18 a.u.。小規模試験で更新間隔1と既定動作の一致、impulseの予測前更新、レーザー初段再利用、無効間隔の拒否を確認。既存回帰7件とMLWF位相追跡試験も通過。

再利用中の交換エネルギーは、現在の軌道と保持したACEによる半トレースの診断値。瞬時の自己無撞着HSEエネルギーや保存すべき凍結Hamiltonianエネルギーではないため、この幅を物理的なエネルギー保存誤差と解釈しない。

電流差は無励起差引き前の同一初期状態の軌道比較。短時間・各条件1回の測定であり、速度の統計誤差、長時間の誘電関数の誤差、採用可能な更新間隔は未確定。既存の局所連続性診断は交換を再構築した時点を測るものであり、再利用中ACEの連続性を検証したことにはならない。長時間の範囲比較は更新間隔1を基準として維持する。

<!-- CURRENT_ERROR_DEFINITION_20260926 -->
## 電流誤差の定義の統一（2026-09-26訂正）

旧版の最初の範囲比較は最終時刻の相対差 `(J(T)-Jref(T))/Jref(T)`、ACE比較は最大差を基準ピークで規格化した量で、指標が不統一だった。最初の表を生データから再計算し、主指標を全表で次に統一した。

- `Jpeak = max_t |Jref(t)|`
- **主指標 `E∞ (%) = 100 max_t |J(t)-Jref(t)| / Jpeak`**（非負）
- 補助指標 `D_T (%) = 100 [J(T)-Jref(T)] / Jpeak`（符号付き）

E∞は「ピーク値同士の差」でも「最終時刻だけの誤差」でもない。D_Tも同じ分母を使い、Jref(T)で割らない。対象はx電流・無励起差引き前である。

各比較は同じ時刻列の全範囲・ACE間隔1を参照する。最初の範囲比較は8ステップ（0.03096 fs）、ACE単独・併用試験は16ステップ（0.06192 fs）。窓が異なるE∞を直接順位比較しない。以下の補助表とcurrent-error-metrics.jsonに評価窓・基準ピーク・最大差の発生時刻も記録した。

| 試験 | 条件 | ステップ数 | 実測秒 | 高速化率 | E∞ (%) | D_T (%) |
|---|---|---:|---:|---:|---:|---:|
| 範囲比較・8ステップ | R=full | 8 | 227.62 | 1.000倍 | 0.00000 | +0.00000 |
| 範囲比較・8ステップ | R=9 | 8 | 101.25 | 2.248倍 | 0.19196 | +0.19196 |
| 範囲比較・8ステップ | R=8 | 8 | 90.77 | 2.508倍 | 0.46714 | +0.46714 |
| 範囲比較・8ステップ | R=7 | 8 | 87.31 | 2.607倍 | 0.82320 | -0.82320 |
| ACE比較・16ステップ | 全範囲・N=1 | 16 | 336.73 | 1.000倍 | 0.00000 | +0.00000 |
| ACE比較・16ステップ | 全範囲・N=2 | 16 | 219.16 | 1.536倍 | 0.31990 | -0.31990 |
| ACE比較・16ステップ | 全範囲・N=4 | 16 | 153.20 | 2.198倍 | 0.96453 | -0.96453 |
| ACE比較・16ステップ | 全範囲・N=8 | 16 | 119.28 | 2.823倍 | 2.26090 | -2.26090 |

高速化率は各試験群の全範囲・N=1を基準とする。8ステップ群は227.62秒（連続性診断あり）、16ステップ群は336.73秒（同診断なし）。異なる群の時間・高速化率を直接比較しない。

<!-- ACE_LOCAL_COMBINED_20260926 -->
## ACE再利用とMLWF局所範囲を併用した測定（2026-09-26）

Si128・16フラグメント、MPI16/OpenMP1/BLAS1。16ステップ、dt=0.16 a.u.、計0.06192 fs。前回と同じバイナリ・初期状態・初期U・励起強度を使用し、連続性診断を無効にして全9条件を逐次実行した。impulseの最初の予測伝播前には全条件でACEを再構築する。

局所範囲は初期WF中心からの全体系の周期的x距離±R bohr。y,z方向は保持する。交換のソースWFのみを切り、フラグメント全体の周期的FFTカーネルは維持する。FFT格子自体を縮小した試験ではない。密度・Hartree・半局所XC・MLWF位相追跡は継続する。

基準は前回の全範囲・毎ステップ更新（336.73秒、RT部分324.19秒）。全範囲で同じ更新間隔の測定も再利用し、二つの効果を分ける。各条件1回で統計誤差は未測定。

| x半幅 R (bohr) | ACE間隔 N | 全実測秒 | RT部分秒 | 全範囲N=1比の高速化 | 同じNで局所化を加えた高速化 | E∞ (16ステップ) | D_T (16ステップ終端) | 同じRのN=1からの最大差 / 全範囲ピーク |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 1 | 124.05 | 115.57 | 2.715倍 | 2.715倍 | 0.3718% | +0.3718% | 0.0000% |
| 9 | 4 | 81.81 | 72.86 | 4.116倍 | 1.873倍 | 0.5722% | -0.5722% | 0.9439% |
| 9 | 8 | 76.07 | 66.49 | 4.426倍 | 1.568倍 | 1.8435% | -1.8435% | 2.2152% |
| 8 | 1 | 128.83 | 119.27 | 2.614倍 | 2.614倍 | 1.0062% | +1.0062% | 0.0000% |
| 8 | 4 | 81.26 | 72.31 | 4.144倍 | 1.885倍 | 0.2366% | +0.0477% | 0.9584% |
| 8 | 8 | 73.62 | 64.39 | 4.574倍 | 1.620倍 | 1.2545% | -1.2545% | 2.2606% |
| 7 | 1 | 121.96 | 113.18 | 2.761倍 | 2.761倍 | 3.4325% | -3.4325% | 0.0000% |
| 7 | 4 | 79.95 | 70.73 | 4.212倍 | 1.916倍 | 4.3657% | -4.3657% | 0.9332% |
| 7 | 8 | 72.23 | 63.47 | 4.662倍 | 1.651倍 | 5.6223% | -5.6223% | 2.1898% |

電流差はx電流の全時刻における最大絶対差を、全範囲N=1のピーク電流で規格化したもの。無励起差引き前の比較であり、誘電関数誤差ではない。局所範囲とACE保持の誤差は符号が逆なら相殺し得るため、合計誤差が小さいだけで精度が向上したとは判断しない。

ソースマスクは非変分的な近似で、保持したACEと局所電流の整合性は未検証。今回は速度の実測と短時間の軌道比較に限定し、採用可能な積分範囲・更新間隔は未確定。再利用中の半トレースエネルギーは診断値であり、保存エネルギーの誤差とは解釈しない。

元の1024ステップ計算はプロセスを一時停止して保存し、今回の測定終了後に同じプロセスを再開する。この長時間計算の経過時間には停止時間が含まれるため、そのまま速度比較には使わない。

相殺の実例：8 bohr・N=4の最終時刻では、局所範囲だけによる電流差が +2.783117e-08 a.u.、ACE保持による追加差が -2.651099e-08 a.u.、合計が +1.320178e-09 a.u.となる。これは各時刻の最大差を示す上表とは別の量である。

### 同一定義での評価窓と最終時刻の補助指標

主指標E∞は全時刻の最大差を基準ピークで規格化した量。最終時刻の補助指標D_Tも同じ分母を使う（Jref(T)では割らない）。E∞の8ステップ列は初期範囲比較と評価窓を揃えたものである。

| R (bohr) | N | 実測秒 | 高速化率（16ステップ） | E∞・最初の8ステップ (%) | D_T・16ステップ終端 (%) |
|---:|---:|---:|---:|---:|---:|
| 9 | 1 | 124.05 | 2.715倍 | 0.19196 | +0.37177 |
| 9 | 4 | 81.81 | 4.116倍 | 0.27617 | -0.57215 |
| 9 | 8 | 76.07 | 4.426倍 | 0.96085 | -1.84346 |
| 8 | 1 | 128.83 | 2.614倍 | 0.46714 | +1.00616 |
| 8 | 4 | 81.26 | 4.144倍 | 0.18516 | +0.04773 |
| 8 | 8 | 73.62 | 4.574倍 | 0.68886 | -1.25446 |
| 7 | 1 | 121.96 | 2.761倍 | 0.82320 | -3.43248 |
| 7 | 4 | 79.95 | 4.212倍 | 1.30278 | -4.36568 |
| 7 | 8 | 72.23 | 4.662倍 | 1.97511 | -5.62233 |

<!-- LCFO_RT_OPTIMIZATION_20260926 -->
## WF更新・分散ACE・行列演算の効率化（2026-09-26）

同じSi128・16フラグメント・MPI16/OpenMP1/BLAS1、R=9 bohr、16ステップ、dt=0.16 a.u.で新旧バイナリを逐次比較した。今回は両方とも既存の `&analysis out_rt_energy_step=10` を設定した。密度・Hartree・XC・電流の更新間隔は変更していない。GNU Fortran15/AArch64では従来どおりループベクトル化を無効化した。

変更は、(1) ACE保持中はMLWF係数の位相追跡のみ行い未使用の実空間ソースを作らない、(2) ACEの重なりをMPI集約し担当コアの出力行だけ計算する、(3) 局所Hamiltonianの射影と交換作用の加算を係数空間でまとめて1回だけ再構築する、の3点。B†との積には複素共役転置をBLASへ直接指定し、大きな転置一時配列を避けた。物理近似、初期U、impulse第1ステップでのACE再構築、更新スケジュールは維持した。

| ACE間隔 N | 実装 | 全実測秒 | RT部分秒 | 同一入力・旧実装比の高速化 | E∞ (%) | D_T (%) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | 旧 | 125.55 | 116.71 | 1.000倍 | 0.37177 | +0.37177 |
| 1 | 最適化後 | 117.94 | 108.59 | 1.065倍 | 0.37177 | +0.37177 |
| 4 | 旧 | 67.49 | 59.44 | 1.000倍 | 0.57215 | -0.57215 |
| 4 | 最適化後 | 58.20 | 49.41 | 1.160倍 | 0.57215 | -0.57215 |

E∞とD_Tは既存表と同じ16ステップ窓・全範囲N=1のピークで規格化した無励起差引き前のx電流。高速化率は今回の同じN・同じエネルギー間隔の旧実装を基準（旧時間／新時間）としており、上の全範囲N=1基準の倍率とは異なる。各条件1回の測定で、ばらつきは未評価。

N=1の新旧電流最大差は 9.190e-18 a.u.（基準ピーク比 3.322e-10%）。RT部分のみの高速化は1.075倍。

N=4の新旧電流最大差は 7.970e-18 a.u.（基準ピーク比 2.881e-10%）。RT部分のみの高速化は1.203倍。

N=4のHamiltonian全体の最大ランク時間は24.488秒→16.155秒（1.516倍）。この計時には交換作用・射影に加えて差分演算と擬ポテンシャルも含まれる。

エネルギー計算間隔だけを変えた旧実装との比較：

| N | 旧実装・間隔1（前回）秒 | 旧実装・間隔10（今回）秒 | 参考倍率 |
|---:|---:|---:|---:|
| 1 | 124.05 | 125.55 | 0.988倍 |
| 4 | 81.81 | 67.49 | 1.212倍 |

上の参考倍率は別時点の単回測定であり、エネルギー間隔の効果の厳密な統計評価ではない。今回の最適化速度にはこの間隔変更を混ぜていない。

ACE適用1回のMPI集約要素数は旧1024×256複素数から、通常256×256、中点512×256へ減少（それぞれ1/4、1/2）。各ランクが出力する係数行も1024行から担当64行へ限定した。これは配列寸法からの計算で、通信時間の個別プロファイルではない。再構築・エネルギー評価・位相追跡時の全係数集約は残る。

検証：不均等2+3行分割の複素数ACE試験（ACEランク2/4/6）は従来の密行列作用＋射影と約4e-16で一致。MLWFの位相不変性、予測子巻き戻し、保持中のフレーム追跡、native LCFO RTの間隔1/2/4、impulse/laser再構築規則の回帰試験を通過した。短時間の演算同等性と速度を検証したもので、長時間の誘電関数精度の検証は別途継続する。

追加のCTest 7件（LCFO core、複素DC-LCFO、DC-HSE）も全件成功。コードレビューで重大な指摘はなく、初期化直後の不要な位相輸送は除去した。測定と検証の終了後、停止していた長時間スペクトル計算を元のプロセス・元のバイナリで再開した。長時間計算には今回の最適化はまだ適用していない。

<!-- MLWF_SPHERICAL_SUPPORT_20260926 -->
## WF積分範囲を三次元球状カットへ変更（2026-09-26）

`SALMON_LCFO_RT_RADIUS=R` は、ここから初期WF中心からの三次元周期距離の半径R（bohr）を表す。直交セルで各軸の最小像変位を求め、dx²+dy²+dz² > R² のソースWFをゼロにする。y,z方向も判定に含む。半径0は全範囲、半径がセルの半対角長以上でも全範囲になる。中心はxyzの周期的モーメントから求め、時間発展中は初期中心に固定する。いずれかの軸で中心信頼度が0.1未満なら、そのWF全体を切らずに保持する。

ソースWFと破棄ノルムの診断には同じ三次元判定を使う。密度・Hartreeの範囲、交換カーネル、初期Uの局在化と位相追跡、ACE更新スケジュールは従来どおり。初期診断ファイル `lcfo_mlwf_initial.bin` はversion 2へ更新し、centers(3,no)をFortran配列順で記録する。旧version 1はx中心のみ。linksファイルはversion 1のまま。

検証：異方的な格子間隔を持つセルで、周期境界をまたぐWF、R=1.51、x半セル長を超えるR=4.1、R=0と十分大きなR=9、弱いy中心の保護を確認した。xyz中心と保護フラグのバイナリ保存も検証した。既存の位相不変性、予測子巻き戻し、ACE保持中のWF追跡試験は成功。

MPI2のnative LCFO RT回帰試験も成功。全範囲MLWF経路と基準の最大差は約9.0e-21、十分大きな半径は全範囲と一致し、有限半径では異なる有限応答を確認した。ACE間隔1/2/4、impulse第1予測子前の再構築、レーザー初期ACE再利用も検証した。この試験の16×8×8 bohrセルの半対角長は約9.80 bohrのため、全範囲試験の半径は旧8から10へ変更した。

**これより前のSi128の9/8/7 bohrの表はx方向だけのカットであり、三次元球状カットの速度・精度を表さない。** 三次元版のSi128の速度・電流差・誘電関数はまだ再測定していない。既存の長時間計算は旧バイナリによるx方向カットの比較として保持し、試験終了後に元のプロセスを再開した。

<!-- LCFO_TRACE_OPENMP_20260926 -->
## 分散ACE診断とWFマスクのOpenMP化（2026-09-26）

保持中ACEの診断エネルギーを、全係数のACE作用Wを各MPIランクで再構築する方法から、各ランクのF_local†C_localを集約してそのノルムから求める方法へ変更した。E = -dv/2 Σ_n f_n ||Σ_rank F_local†C_local||² は従来の半トレースと同じ量であり、追加近似はない。診断値は引き続き各refresh時点で計算する。out_rt_energy_stepでこの診断自体を間引いたわけではなく、同じ値の計算方法を軽くした。

三次元マスクをWF列ごとの連続メモリアクセスへ並べ替え、WF間をOpenMP並列化した。格子位置はWFループの外で一度作り、破棄ノルムはWFごとに集計して固定順序で合算する。MPI通信は並列領域の外で行う。局在範囲、球状マスク、保護WF、ACE更新時刻は維持している。

Si128・16フラグメント・三次元球R=9 bohr・ACE間隔4、16ステップ、dt=0.16 a.u.、energy間隔10。同じGSデータ、初期Uと中心、MPI16/OMP1/BLAS1で新旧を逐次比較した。

| 実装 | 全実測秒 | RT部分秒 | 旧実装比の全時間高速化 | E∞ (%) | D_T (%) |
|---|---:|---:|---:|---:|---:|
| 変更前 | 63.63 | 54.85 | 1.000倍 | 0.31977 | +0.31977 |
| 変更後 | 63.52 | 53.78 | 1.002倍 | 0.31977 | +0.31977 |

E∞とD_Tは既存の16ステップ窓・全範囲N=1のx電流ピークで規格化した指標。今回は三次元球の結果であり、以前のx方向カットと区別する。新旧の出力電流と出力エネルギーの最大差はともに0（ファイルの出力精度内）。RT部分は約1.020倍だが、全時間は約1.002倍で、単回測定から明確な高速化とは判断できない。局所範囲とACE保持の誤差は相殺し得るので、この短時間のE∞だけでは誘電関数の精度を保証しない。

このマシンは物理18コア。16フラグメントに必要なMPI16を維持してOMP2にすると32スレッドになるため、Si128の比較ではOMP1とした。追加したOpenMPのスケーリングは、この大系では未測定。BLASも1スレッドに固定した。

検証：分散半トレースは不均等MPI分割・複素数・dv=0.7・分数占有で従来のACE作用と最大2.22e-16の差。球状マスクと2WFの時間追跡をOMP1/2/4で検証し、native MPI2の回帰試験はOMP1/2で成功。同じGSデータを再利用する三次元有限半径のMPI2時間発展では、OMP1対2/4で初期診断ファイルがバイト一致し、出力電流差も0だった。

別々のOMP設定でGSから作り直す比較では、中心の微小差で鋭い半径境界のマスクが変わり、有限半径の電流に最大1.19e-7 a.u.の差が生じた。同じGSを用いた比較では解消したため、今後もスレッド数や速度の比較には同じGS・初期U・中心を使う。マスクは現時点で鋭い切断のままである。

短時間試験終了後、元の長時間計算を旧バイナリのまま再開した。今回の計算は逐次実行し、数値ジョブの同時実行は行っていない。


<!-- LCFO_TWO_LEVEL_MPI_20260926 -->
## フラグメント×軌道の2段階MPI（2026-09-26）

LCFO RTで既存の `icomm_r`（同じ軌道群を持つフラグメント間）と `icomm_o`（同じフラグメント内の軌道群間）を利用する。native波動関数、Hamiltonian作用、時間発展は `io_s:io_e` の担当軌道だけを保持・処理し、密度・Hartree・半局所XC・電流の更新は既存ルーチンを使う。軌道グループ0が各フラグメントの全軌道係数を受け取り、MLWF・交換行列・ACEを一度構築する。構築時だけ交換行列とACE因子を同じフラグメントの他軌道群へ配布する。保持ステップでは因子を再配布しない。

Gamma点・非スピン分極・直交セル・1空間ランク/フラグメントという条件は維持する。全体の交換行列とACE係数因子は各ランクに複製され、再構築担当には全ソースWFが残る。したがって全交換処理の完全分散ではない。DC-HSE基底状態計算のフラグメント内部の軌道MPI制限は今回変更していない。

### Si64・8×1×1比較

通常セル8個を直列配置したSi64、82.08×10.26×10.26 bohr、格子128×16×16、8フラグメント、buffer=8,0,0、LCFO基底64/フラグメント、占有128軌道。共通のDC-HSE基底状態は1148回、密度差9.9821446e-8で収束し、複素LCFO固有方程式残差は2.1988e-16。共通の保存LCFOデータを使い、追加SCFなしでRTを開始した。

三次元球R=9 bohr、ACE間隔4、impulse=1e-4 a.u.、16ステップ、dt=0.16 a.u.、energy出力間隔10。両ケースOMP1/BLAS1、GNU Fortran15/AArch64のループベクトル化無効。数値ジョブは常に1本ずつ実行した。以前キャンセルした長時間Si128計算は再開していない。

MPI8は `nproc_rgrid=8,1,1; nproc_ob=1`、MPI16は `nproc_rgrid=8,1,1; nproc_ob=2`、両方 `nproc_k=1`。MPI16でも `num_fragment` と `nproc_rgrid_tot` は8,1,1のままとする。

| フラグメント×軌道 | MPI | 全実測秒 | RT部分秒 | 全時間高速化 | 最大ランクRSS MiB | ランク別ピーク平均 MiB | 同時合計RSS最大 MiB | E∞ (%) | D_T (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8×1 | 8 | 24.46 | 21.545 | 1.000倍 | 451.4 | 444.7 | 3557.5 | 0 | 0 |
| 8×2 | 16 | 26.93 | 23.602 | 0.908倍 | 402.9 | 284.6 | 4553.2 | 1.734e-10 | +1.734e-10 |

ここでのE∞ = 100 max_t|Jx_MPI16−Jx_MPI8|/max_t|Jx_MPI8|、D_T = 100(Jx_MPI16(T)−Jx_MPI8(T))/max_t|Jx_MPI8|。**比較基準は同じR=9・ACE4のMPI8であり、全範囲ACE1に対する物理近似誤差ではない。** 電流3成分の最大差7.943e-18 a.u.、出力時刻のエネルギー差0（出力精度内）、16ステップ目の密度最大差2.102e-15 a.u.。初期U・中心・占有係数の診断ファイルはバイト一致し、ACE構築11回と全更新時刻も一致した。impulse第1予測子前の再構築を維持している。

RSSは約0.5秒間隔の観測値で、厳密な瞬間ピークではない。ランクごとに複製されるライブラリ・格子配列なども含む。MPI16では8ランクが約397–403 MiB、残る8ランクが約169–171 MiBのピークだったが、PIDだけからランク役割を同定してはいない。最大ランクは約11%減、ピークのランク平均は約36%減、同時合計は約28%増。波動関数1配列のコア部分だけなら、4096点×128軌道×16 byte=8 MiBから4 MiBへ半減する（haloや作業配列を除く）。

この系では全時間は約10%増え、高速化は得られなかった。Hamiltonian全体の最大ランク時間は5.7582→4.0236秒（約1.43倍）へ改善した一方、RT全体は21.545→23.602秒。伝播作用の軌道分割には効果があるが、ソース・ACE再構築を8つの軌道グループ0で行う部分と、新たな群間通信は残る。再構築と通信の寄与の個別計測は未実施で、今回の全時間差の原因割合は断定しない。各条件1回の短時間測定であり、大系・長時間のスケーリングや誘電関数精度は未検証。

### 回帰試験

2フラグメントのMPI2対MPI4を、ACE間隔1/4、有限球R=3、`process_allocation='orbital_sequential'`、滑らかなレーザーで比較した。電流・エネルギー・最終密度・更新時刻を確認し、行数と時刻も照合した。3占有軌道を2+1へ不均等分割する試験では、共通の保存基底を両側で同じように再占有して演算同等性を確認した（この試験はGSの物理精度の評価ではない）。最大差は電流9.632e-15、エネルギー1.022e-14、密度4.902e-14。既存の半時間刻み、全範囲・大半径一致、有限半径応答、連続性診断、restart拒否、impulse/laser更新規則も成功した。2段階MPIでのACE棄却→全交換行列へのフォールバックを意図的に誘発する統合試験は未実施。

### Distributed exchange construction (2026-09-26)

The LCFO coefficient rows now remain on their spatial/core owner. Transported
WF frames, W=Hx C, and ACE factors have the same row distribution. Each fragment
requests only the coefficient rows touched by its core+buffer basis; the adjoint
exchange sums overlapping fragment contributions back onto their owners. The
exchange matrix is retained as a Hermitian fragment block, including both
midpoint endpoints. Rejected ACE therefore uses the same distributed fragment
operator instead of collecting a dense global Hx. Every orbital group owns a
matching halo plan for its target columns. Exact-cache decisions are collective
across spatial ranks.

Polar overlap and ACE metrics are sums of local-row products. Small SVD/EVDs are
computed once on spatial rank zero. A build with `USE_SCALAPACK=ON` uses a BLACS
process grid and block-cyclic PZGESVD/PZHEEV for >=128 occupied states and multiple
spatial ranks. This LCFO backend selection is automatic at compile-time/size;
`yn_scalapack` still controls the original SALMON eigensolver and does not toggle
this new backend. Without ScaLAPACK the root LAPACK path remains available.

Remaining limits are explicit: orbital-space metric/U/transform arrays are
replicated; initial MLWF gauge optimization and its coefficient gather run once
on the root. Per-fragment exchange construction still runs on orbital group zero,
while target application is shared across orbital groups. Real-space WF arrays
still contain all occupied source columns locally. Thus this removes global-row
replication and redundant dense decompositions, but does not establish linear
scaling or fully distribute all orbital-space storage.

`LCFO distributed storage rank/local/global/halo rows` reports the actual retained
coefficient rows and requested halo rows. `LCFO timing pack/source/exchange/ACE`
reports the spatial maximum per phase (seconds); maxima for separate phases may
belong to different ranks and their sum is not the total wall time. The source
phase includes gauge transport and reconstruction; retained-step source time
also includes the frozen exchange trace. Initial MLWF work belongs to the first
source phase. These diagnostics supplement the native `rt iterations` timer.

The standalone `testsuites/unit_lcfo_rt/test_distributed_build.py` accepts
`--scalapack` and `--ranks 2|4`. It covers unequal local rows, arbitrary halo row
order, the sum of overlapping Hermitian fragment operators versus a dense
reference, root-only gather, dense-unitary polar transport, local ACE factors and
action, globally/locally zero exchange, and positive/singular metric rejection.
The large case uses129 states to exercise partial block-cyclic tiles, including
a2x2 BLACS grid with MPI4. Native tests also assert local storage and compare
fragment-only versus fragment x orbital parallel trajectories.

Validation/measurement: all seven sequential standalone/native suites passed. MPI4/129-state ACE action error was6.36e-16. Same-GS Si64/128 R9 ACE4 MPI8/16 OMP1 BLAS1 comparisons against b1773c78 preserved the initial dump byte-for-byte, with current differences <=7.26e-18, density <=2.51e-15, and no difference in printed total energy. RT times old→new:20.780→19.847s and48.339→43.238s (single runs). Maximum sampled rank RSS:448.1→438.3MiB and662.0→652.1MiB. The change improves RT time by1.047x/1.118x here; it does not remove the inverse weak-scaling trend. Native forced ACE-invalid fallback was not added; direct fragment action and rejection are covered at the algebra level.

### Select active WF columns before reconstruction (2026-09-26)

`lcfo_wf_support` caches the original minimum-image point mask because the initial
centers, radius, grid and protected flags remain fixed during polar transport.
Only columns with at least one retained point are reconstructed. Protected WFs
remain unmasked. Compact sources preserve their original order; the HSE source
sum is unchanged and an empty source list produces zero exchange. The core
loss diagnostic uses `trace(F^H (B^H B) F) - ||(BF)_retained||^2`; the basis Gram
matrix is cached once. This identity also holds for a nonorthogonal complex B.
Near-zero loss can suffer subtraction roundoff and is not a meaningful measure
below floating-point accuracy. This diagnostic does not modify exchange sources.

The optimization adds no mask approximation and leaves polar transport, initial
MLWF optimization, ACE scheduling and native fields unchanged. Dense full-column
coefficient frames and halo communication remain; orbital metrics/U and initial
root localization still limit large systems. Support preparation uses a temporary
point-by-all-WF boolean mask, while only active masks are retained afterward.

Tests: `test_active_wf.py` verifies dense reconstruction/mask equivalence for a
complex nonorthogonal basis, periodic edges, protected WFs, empty/full support,
and the norm identity. Existing transport and sphere fixtures pass at OMP1/2/4;
native MPI2/MPI4 tests pass including unequal orbital blocks and impulse/laser
ACE schedules. Static review found no blocking defect.

Diamond weak scaling at33dbea5a: C64/8 MPI versus C128/16 MPI, core16^3,
buffer8,0,0, dt0.02, nt16, ACE1, OMP1/BLAS1 on one18-core host. Full RT:
104.530→304.010s (34.38% weak efficiency). R6 RT:57.684→100.370s (57.47%).
Both R6 cases retain60 active source WFs per fragment. Native target orbital
columns still increase128→256 per rank because nproc_ob=1.

Same-GS R6 active-WF optimization: repeated WF phase2.6216→1.2369s for C64 and
7.8446→2.7070s for C128 (2.12x and2.90x). RT57.684→56.595s and100.370→94.138s
(1.019x and1.066x); new R6 weak efficiency60.12%. Reconstructed source/core
columns are60/44 for both sizes. Initial dumps match byte-for-byte; max current
difference3.90e-18, density2.00e-15; printed energy and norm-loss logs agree.
These are single-run measurements and short-trajectory implementation parity,
not a long-time dielectric accuracy certificate. Raw benchmark data and the
expanded Japanese notebook include both measured versions separately.

## Selected WF column halo (2026-09-26)

The fixed source support now seeds a cached column-request plan on the row halo. Each owner packs the requester’s ordered WF columns; the compact frame feeds reconstruction directly. Zero requests and zero owned rows are supported. The row topology, support columns and occupied dimension must remain fixed for the cache lifetime. Dense U transport and core Gram diagnostics remain unchanged; no physical approximation is added.

Validation: standalone MPI2/4 includes a row-owning rank with no local column requests that still sends to peers, duplicate/reordered requests, empty rows, full columns and repeated gets. Compact reconstruction, transport/sphere OMP1/2/4, and native fragment x orbital MPI2/4 tests pass. Static review found no blocking issue.

Paired sequential diamond R6 benchmarks, same GS, core16^3, dt=.02, nt=16, ACE1, OMP1/BLAS1, GNU loop vectorization disabled, baseline c9abacb9:

| Case | Halo receive complex values/rank before→after | Repeated WF seconds before→after | RT seconds before→after |
|---|---:|---:|---:|
| C64 / MPI8 | 24576→11520 | 1.206884→1.223684 | 56.062→57.873 |
| C128 / MPI16 | 49152→11520 | 2.637532→2.594449 | 100.290→99.992 |

Receive counts include self transfers, not measured network-only bytes. Current max difference 4.31e-18, density 2.00e-15, printed energy/loss identical, initial dump byte-identical. No clear speed benefit on this shared-memory host in one sample per case; exchange remains dominant. Weak efficiency changes 55.90%→57.88%, not statistically established improvement. The objective achieved is bounded WF halo payload for fixed-radius fixed-size fragments; dense U and initial root localization still scale with total occupied states. Raw results are in the diamond notebook output wf-column-summary.json.
