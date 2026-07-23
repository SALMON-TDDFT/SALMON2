# 新セッション指示書・レビューポイント

## 新セッションへの最初の指示

次の一文から開始すること。

> この指示書と参照計画書を読み、`executing-plans` を使って occupied-W bootstrap 計画の Task 6 から開始してください。

Task 1–5を最初からやり直さない。Task 5までの実装をfocused verificationで再確認し、Task 6のdensity-carrying solveとSi64 physical gateを進める。

## 作業場所

```text
worktree: /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
branch:   codex/singlescale-vortex-observables
HEAD:     daa982838f8972ab93c1787f46a96c6473818728
```

これはローカルworktreeである。Codexの環境表示に現れるOneDrive側の
`SALMON-v.2.2.2`を今回の作業場所として使わない。

## 最初に読む文書

1. `docs/plans/2026-07-20-wpw-occupied-w-bootstrap.md`
2. `docs/plans/2026-07-19-wpw-density-carrying-seed-design.md`
3. `docs/plans/2026-07-19-wpw-density-carrying-seed.md`
4. この指示書

最初の計画レビューでは、occupied-W bootstrap計画のTask 6を現在の開始点とする。
古いdensity-carrying seed計画にもTask 5/6という番号があるため、どの文書のTask番号かを必ず明記する。

## dirty worktree保護

現在のworktreeには、Task 5までの意図した未commit変更と未追跡の実行結果がある。
すべてユーザー資産として保存する。

- `checkout`、`reset`、`clean`、stash、破壊的復元を行わない。
- 無関係な差分を編集・整形しない。
- `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_runs/`を削除・上書きしない。
- 新しい計算は必ず一意な新規run directoryで行う。
- commit、push、PR作成は、ユーザーが新セッションで明示的に指示しない限り行わない。

開始確認:

```bash
cd /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG
pwd
git branch --show-current
git rev-parse HEAD
git status --short
git diff --check
```

HEADが上記と異なる場合でも即座に戻さず、履歴と差分を調べてユーザーへ報告する。

## 引き継ぐ実装状態

occupied-W bootstrap計画のTask 1–5は実装済みである。

- W非依存のcore-grid/P bootstrap
- deterministicなoccupied projected-Wannier descriptor、stable ID、16 W/fragment・128 W/global契約
- descriptorから構築するW owner/support layoutとtail halo
- tail-carrying実空間WW/WP/PP metric
- occupied-Wのkinetic/local-potential、stable projector-channelによるcomplex nonlocal、canonical-face WW/WP項
- dynamic local-potential refreshのoccupied-W support-tail経路
- checkpointのstable W IDおよびWP Z support-tail routing
- fixed-H H0/interface分離、frozen snapshot、fingerprint、rollback

Task 5の無条件停止と、WPW production経路からのlegacy
`import_wpw_lcfo_ww_components`呼出しは除去済みである。到達可能なproduction経路は
legacy `wpw_ww_components`に依存しない。legacy helper自体はdead codeとして残るが、
Task 6開始前の整理対象にはしない。

最終コードレビューではP0–P3 findingなし。直前の修正では次を確認した。

- real-space WW metricをproduction readinessとして受理する。
- final WWのmetric/H0/interface/projectorをfingerprintする。
- checkpoint WP Zをsupport-tailまでstable IDでrouteし、exact coverageを検査する。
- MPI collective/reductionの順序を全rankで一致させる。
- checkpoint provenanceはatomic publish後のoccupied-W descriptor fingerprintと
  frozen production snapshotから取る。

## 直前にPASSした検証

- `python3 tests/dg/run_dg_wpw_production_operator_mpi.py`
- `python3 tests/dg/test_dg_wpw_checkpoint_roundtrip.py`
- `python3 tests/dg/test_dg_dc_rt_handoff.py`
- `python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- `python3 tests/dg/check_dg_wpw_checkpoint_publication.py`
- `cmake --build build-mpi-eigenexa -j 4`
- `git diff --check`

新セッションでは少なくとも次を再実行してから物理runへ進む。

```bash
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/check_dg_wpw_checkpoint_publication.py
python3 tests/dg/check_dg_wpw_quadrature_assembler.py
python3 tests/dg/run_dg_wpw_candidate_halo_mpi.py
python3 tests/dg/run_dg_wpw_rank_local_quadrature_mpi.py
python3 tests/dg/run_dg_wpw_production_operator_mpi.py
python3 tests/dg/test_dg_wpw_checkpoint_roundtrip.py
python3 tests/dg/test_dg_dc_rt_handoff.py
cmake --build build-mpi-eigenexa -j 4
git diff --check
```

MPIテストがlocal socket制限で失敗した場合は、数値失敗と決めつけず、同じコマンドを必要な権限で再実行する。

## 次の作業: occupied-W bootstrap計画 Task 6

1. density seedがbootstrap時にpublishした単一のoccupied-W descriptorを、source ensembleとW basis provenanceの双方に使うことを確認する。第二のprojected-Wannier構築や別gaugeを許さない。
2. W/P両方のoverlapから同じcoupled metric方程式 `S C = B` を解く。
3. metric solver未収束時のbest-iterate continuationは、既存の明示的なnonpublishable diagnostic flagの下だけで許可する。
4. focused testsとfull MPI buildの後、過去結果を上書きしない新規directoryでSi64を実行する。
5. physical gateを満たさなければfixed-H、checkpoint/manifest publication、RTへ進まない。

このTaskでは「コード経路が到達可能」であることと「Si64で物理的に妥当」であることを区別する。
Task 5 stopを外した後のSi64 production physical runはまだ実行していない。

### 2026-07-21 Task 6実行結果

Task 6のSi64 production physical runを
`stage2d_wpw_runs/20260721_task6_realspace_fix_gate`で実行した。固定6点bufferと
片側4点stencilから、DG項の共通value/gradient領域をcore外側2点までとし、準備済み
fragment grid外への再帰的buffer要求は除去した。

runは128 global occupied W、routed/direct capture identity、assembled/real-space
WW/WP/PP identityを通過した。その途中で、real-space WW metricが有効なのに未allocateの
legacy `ww_*`をbounded operatorがコピーする不具合を再現テスト付きで修正した。

最終的にはpre-fixed-H physical gateでcollective stopした。主な値は
`captured_norm=9.6927E-01`、density projection residual `1.7052E-01`、projected
charge error `9.8323E-01`、DC-density residual `4.3791E-02`であり、設定値
`1.0E-08`を満たさない。checkpoint/manifestはpublishされていない。したがって次は
buffer閾値を緩和するのではなく、約3.1%のmissing occupied-subspace captureの由来を
調査する。fixed-Hへは進まない。

### 2026-07-22 B=10 image-periodization再実行結果

ユーザー指定bufferを収束パラメータとして扱い、片側4点stencilに対して
`D = 1-B+r:n+B-r`とする実装で、B=10のSi64 runを
`stage2d_wpw_runs/20260722_task6_b10_image_periodized/run_retry.log`まで完走させた。
同一canonical点へ写る異なるfragment-box imageは値・勾配を先に加算し、二次形式は
加算後の場から構築した。これにより以前の`tail_route`停止を越え、128 occupied W、
projection/Gram/chargeの代数的identityを通過した。

ただしpre-fixed-H physical density gateは通過していない。B=10では
`captured_norm=1.0000E+00`、projection residual `1.2430E-10`、projected charge
error `1.9826E-11`である一方、DC-density residualは`9.3116E-01`だった。
比較対象のB=6 common-source runはそれぞれ`1.0000E+00`、`7.3312E-11`、
`3.7168E-11`、`5.0195E-01`である。source chargeもB=6の`3.2217E+01`から
B=10の`1.0386E+01`へ低下し、post-periodization W normの最大値は
`1.2968`から`0.75210`へ変化した。したがってbuffer増加による密度収束は確認できず、
むしろ悪化した。runは`fixed_h_stage=density_seed`でcollective fatal stopし、
fixed-H、checkpoint/manifest publication、RTには進んでいない。

またB=10ではfragment box 32点がtotal cell 24点を越えるため、現行の
`outer_shell_norm=0`診断は有効な収束指標になっていない。次は許容値を変更せず、
source chargeとW normがBに強く依存するデータ経路、およびouter-shell診断の
box>cell条件をroot-cause tracingする。

## Si64 physical gate

### 2026-07-22 unwrapped-P修正後のfresh buffer系列

周期境界をfragment storage生成時に完了し、cyclic storageをunwrapped連続Pへ
並べ替えてからDを作る修正後バイナリで、過去結果を上書きせずB=6/10/12を実行した。
canonical aliasの後段加算は行っていない。

- B=6: `20260722_task6_unwrapped_b6/run.log`。DC-SCFは87反復で
  `diff=9.8498E-10`、charge 256を保った。occupied-W buffer gateは
  `outer_shell_norm=9.3038E-01`、`outer_ratio=7.2686E-03`でcollective stop。
- B=10: `20260722_task6_unwrapped_b10/run.log`。DC-SCFは85反復で
  `diff=9.4464E-10`、charge 256を保った。buffer gateは
  `outer_shell_norm=3.6267E-01`、`outer_ratio=2.8333E-03`でcollective stop。
- B=12（初回、無効）: `20260722_task6_unwrapped_b12/run.log`。36点coverageにより
  fragment当たり864電子となる一方、`nstate_frag=400`は800電子分しか保持できず、
  DC-SCF 88時点で`integral(rho_tot)=236.86594`（target 256）だったため手動停止した。
  これはbuffer実装の失敗ではなく、状態数不足を検出した診断runである。
- B=12（修正版）: `20260722_task6_unwrapped_b12_n500/run.log`。`nstate_frag=500`へ
  増やすとDC-SCFは88反復、`diff=8.1943E-10`で収束し、charge 256を保った。
  buffer gateは`outer_shell_norm=1.3883E+00`、`outer_ratio=1.0846E-02`で
  collective stopした。

B=6からB=10ではouter ratioが減少したが、有効なB=12点では
`1.0846E-02`へ再増加し、いずれも基準`1E-8`には遠い。したがってfixed-H、
checkpoint publication、RTへは進まない。状態数不足を除去しても非単調性が残るため、
次の問題は後段periodizationやDC chargeではなく、outer-shell normをbuffer収束の
必須ゲートとする物理的定義と、B依存で選択されるoccupied subspace/Wannier表現である。

### 2026-07-22 occupied-W spread診断

Wannier90のGamma点離散per-W spread規約を実装し、Pから一つの24点canonical cellを
選んで128個のruntime projected-Wを測定した。B=6 fresh run
`20260722_task7_spread_b6/run.log`では全128中心が有効で、`sqrt(Omega)`は
最小`1.26877 A`、平均`1.28540 A`、median/p90/max`1.29094 A`だった。以前の
約`1.2 A`よりわずかに大きいが、異常な非局在化は認められない。

B=10 fresh run `20260722_task7_spread_b10/run.log`はDC-SCF 85反復、
`diff=9.4464E-10`、charge 256で収束した後、`occupied_w_link_assembly`で停止した。
32点Pには同一canonical点の複数preimageが含まれるため、duplicate-image不一致と
整合するが、現ログは抽出とlink組立てを同じstageにしており、失敗座標と最大差を
まだ直接出力していない。aliasを加算も任意選択もしていない。
したがって現時点では最局在化を実装せず、抽出/link境界を計測してB>6のfragment
projected-W Pが物理24点周期像の同一性を満たすかを次のroot-cause対象とする。outer-shell failureを
Wannier幅の過大さで説明する根拠は得られなかった。

instrumented B=10 run `20260722_task8_image_diag_b10/run.log`で抽出stageを
link組立てから分離した。全8 fragmentで、P座標`[1,16,19]`と`[25,16,19]`
（x方向に物理周期24点差、同一canonical点）のW値が一致せず、最大絶対差は共通して
`3.92919E-02`だった。runは`occupied_w_canonical_cell_extraction`でcollective
stopした。したがってraw B=10 fragment projected-W Pは物理24点並進の周期像同一性を
満たしていないことが直接確認された。次はこのraw boxを後段で任意選択するのではなく、
fragment arrayのbuffer生成時に一意なphysical-cell representativeから周期像を構築する
設計を行う。

静的追跡では`src/gs/dc/dcdft.f90:init_fragment`がfragment計算セルとgridを
`domain+2B`へ置換している。B=10のraw fragment orbitalは32点fragment-cell周期で
あり、`dc%jxyz_tot`だけがその32 storage点を24点physical cellへfoldする。したがって
storage 15/23が同じglobal点へ写っても値は一致しない。発生源はWannier変換やP reorder
ではなく、raw fragment wavefunctionをphysical gridへ対応付ける境界である。修正では
raw fragment orbitalへ24点周期性を要求せず、`dc%jxyz_tot`ごとに一意なrepresentativeを
選び、その値からWPW用Pの全周期像を微分前に構築する必要がある。

### 2026-07-22 physical-periodic PとWannier規格化

WPW用occupied-orbital Pを`unwrapped_buffer_order`直後にphysical-periodic化した。
各canonical点についてclosed core boxへ最も近いpreimageを一つ選び、P内の全aliasへ
同じ値をcopyする。加算・平均は行わず、raw fragment `spsi`も変更しない。このPから
projected Wと解析微分の双方を作る。実装commitは`e15f7de`、production placementは
`1d7c22a`である。

B=10でcanonical mismatchは解消したが、32点fragment cell規格化から24点physical
representativeを選ぶことで初期W normが`0.987447--0.990242`となることを計測した。
各Wをprepared P内のunique canonical点上で規格化し、同じ列scaleを
`polar_transform`、core値、P値へ一度だけ適用した。解析微分はscaled transformから
生成する。全rankがnormとlocal fieldをcollective preflightした後にのみmutationし、
補正後normを全Bで再測定する。helperは`e612c94`、canonical norm assemblerは
`7625f45`、production integrationは`ae4b754`である。

fresh結果は次の通り。

- B=5: `20260722_task11_normalized_b5/run.log`。Si64 fragment SCFがiteration 75以降
  不安定化し、iteration 167で`diff=1.2466E-01`のため停止した。半径4に対してB=5は
  形式上有効だが、この物質ではupstream fragment SCFが未収束で、WPW bootstrapの
  end-to-end検証は未完了である。partial-Pがbuffer拡張を要求しない契約はfocused
  helper/source-contract testで検証しているが、このrunの実測結果ではない。
- B=6: `20260722_task11_normalized_b6/run.log`。SCF 87反復、
  `diff=9.8498E-10`。128 Wの最大norm誤差`1.11022E-15`、mismatchなし、width
  `1.26877--1.29094 A`、mean `1.28540 A`。
- B=10: `20260722_task11_normalized_b10/run.log`。SCF 85反復、
  `diff=9.4464E-10`。128 Wの最大norm誤差`1.15463E-14`、mismatchなし、width
  `1.23411--1.26600 A`、mean `1.25543 A`。

B=6/B=10の最大center component差は出力精度で`7.531E-03 A`、最大Omega差
`8.676E-02 A^2`、最大width差`3.466E-02 A`であり、局在幅は従来想定の約1.2 Aに近い。
一方、未変更outer-shell gateはB=6で`7.2686E-03`、B=10で`2.2432E-02`
（tolerance `1E-08`）のため、いずれもfixed-H前で停止した。B>6のPにはphysical
canonical aliasが含まれるため、P全点を数える現outer-shell/total normは重複像を含み、
buffer収束指標として再定義が必要である。gateを緩和して先へ進めてはいけない。

少なくとも以下をrank-local値とcollective aggregateの両方で確認する。

- occupied Wが16/fragment、128/globalであり、stable ownerが一意である。
- core値と隣接fragmentへ通信したWannier tailから得るW normが有限かつ妥当である。
- routed overlapとdirect overlapのtotal/W/P captureが一致する。
- assembled Sと実空間 `A^dagger A` のWW/WP/PP分解が一致する。
- `captured_norm`、source/projected/normalized density、chargeが物理的な許容範囲に入る。
- occupied rank、source Gram、normalization-density identityが合格する。
- fixed-Hへ入る場合、H0/interface分離、frozen value/shape/fingerprint、transport IDが不変である。
- failure時に全rankが同じ場所で停止し、部分的context/checkpointをpublishしない。

## 全RHS収束条件の扱い

metric solveの全RHS `1e-10`条件を今すぐ緩和・変更しない。これまでのblock Jacobi診断では
condition estimateは約`1.06e5`から`2.31e3`へ改善したが、256反復でworst RHSは
`O(1e-7)`だった。その後に見つかった巨大charge/densityは、PCGだけでなく旧W表現と
metric/operator contractの不一致が原因だった。

occupied-W routeで物理計算をfixed-H到達点まで進め、その結果を得てから、全RHS条件、
iteration cap、overlap-1 Additive Schwarzの必要性を再検討する。physical gate前に
cutoffを緩めたり、反復数だけを増やして問題を隠さない。

## レビューポイント

レビューは次の順序で行う。

1. source provenance: descriptorの再構築やgaugeの二重化がないか。
2. overlap routing: 各coreがnonzeroな全support W/P rowを評価し、stable source IDとbasis-row IDでcanonical ownerへrouteしているか。
3. metric/operator identity: WW/WP/PPが同じtail-carrying fieldから組み立てられているか。
4. MPI collective: rank依存分岐の外でcollective順序が一致し、missing/duplicate/wrong-imageをcollective failureにするか。
5. transaction/rollback: allocation・validation失敗時に既存contextとcheckpoint stateを変えないか。
6. frozen invariant: 値、shape、stable ID、H0/interface、projector、face項を実データとfingerprintの双方で検査するか。
7. publication boundary: diagnostic continuationからcheckpoint/manifest/RTへ到達できないか。

## 停止条件と後続順序

以下のいずれかで停止して診断・レビューする。

- W count/owner/tail coverage不一致
- assembled/real-space metric identity不一致
- 非物理的なcapture、charge、density、rank
- provenance曖昧、frozen invariant破壊、collective不一致
- nonpublishable diagnostic結果しか得られない状態

全physical gateが通った場合のみ、既存計画に従ってzero-interface fixed-H、interface continuation、
field-off、20-step laser smokeの順で進む。Si64長時間productionやHHG解釈には進まない。

## 2026-07-22 規格化後の非停止診断とfixed-H到達結果

ユーザー方針に従い、規格化済みoccupied-Wの近似品質を実計算で評価するため、有限で
構造的に妥当なouter-tail超過とsource-to-DC density差をwarningへ分離した。非有限値、
非正norm、射影・規格化・電荷・rank不良は引き続きcollective fatalであり、fixed-H、
frozen-state、checkpoint/RTの各gateは変更していない。実装は`06e0b0c`と`6622842`、
設計・計画は`da152bc`と`1dae824`。二回目のコードレビューでCritical/Importantなし。

fresh B=6 runは
`stage2d_wpw_runs/20260722_task13_dc_density_warning_b6/run.log`。DC-SCFは87反復、
`diff=9.8498E-10`、電子数256を維持した。128 Wのtail ratio `7.2686E-03`はwarningで
通過した。続くdensity seedは`captured_norm=1`、projection residual
`7.7960E-11`、normalization residual `1.3946E-15`、projected charge error
`5.9380E-11`で構造gateを通過し、source-to-DC residual `6.0445E-02`をwarningとして
初めてzero-interface fixed-H algebraへ到達した。

fixed-Hの最初の`solve_window`は160 inner反復で`info=40`となった。これは反復上限で
`residual < 1E-8`を満たさなかったことを表す。orthogonalityは約`1.36E-14`を維持したが、
最大generalized residualは初回`3.3809E-01`から反復150で`1.8466E-03`、反復160で
`1.6650E-03`までしか低下せず停滞した。したがって現在の次の境界はWannier tailや
source-density gateではなく、fixed-H window eigensolverの全160 retained stateに対する
最大残差である。runはWPW checkpoint/manifest/RTをpublishせず、full LCFOへ安全に
fallbackしてSALMON自体はexit 0で終了した。

次は許容値や反復上限を直ちに変更せず、occupied 128 stateとextra 32 stateの残差を分離し、
最大残差がoccupied側かextra側か、またpreconditioner後にどのstateが停滞しているかを
診断する。その結果により、occupied-projector収束を判定対象とするか、extra-state seed/
preconditionerを改善するかを選ぶ。

### 2026-07-22 fixed-H状態別残差診断

`920381e`で収束判定を変えずにoccupied/extra blockの最大残差、raw-worst state、
preconditioned norm/ratioを追加した。fresh B=6 run
`stage2d_wpw_runs/20260722_task14_window_state_residual_b6/run.log`は同じdensity seedを
再現し、fixed-H反復160で再び`info=40`となった。

- inner 1: occupied `6.2988E-02` (state 94)、extra `3.3809E-01` (state 160)
- inner 32: occupied `2.9049E-03` (state 128)、extra `5.3288E-03` (state 150)
- inner 80: occupied `1.1683E-03` (state 64)、extra `3.5465E-03` (state 144)
- inner 160: occupied `4.8659E-04` (state 66)、extra `1.6650E-03` (state 146)

extra blockが全160状態の最大残差を決めるが、occupied blockだけでも`1E-8`より約4桁
大きく、extra statesを判定から除くだけではfixed-Hをqualifyできない。さらにinner 160の
raw-worst stateで、現diagonal preconditionerはoccupied残差normを`2.2726E-02`
（ratio `46.704`）、extraを`5.6867E-01`（ratio `341.54`）へ増幅する。selected history
全体でもratioはoccupied約33--563、extra約8--910であり、次の主対象は許容値や
Wannier tailではなく`H_ii-epsilon S_ii` diagonal preconditionerの小分母・符号・scaleである。
checkpoint/manifest/RTはpublishされていない。

### 2026-07-22 無preconditioner比較

`a59c5f9`でdefault=`y`の明示比較controlを追加し、`n`のときだけfixed-H/continuationの
optional preconditioner callbackを省略する。fresh B=6
`stage2d_wpw_runs/20260722_task15_no_precondition_b6/run.log`は同じnormalized seedを
使い、反復160の`info=40`まで実行した。

無preconditionerではinner 96のoccupied/extra残差が`3.6079E-04`/`1.1093E-03`
（有りでは`1.0673E-03`/`3.4053E-03`）、inner 160が
`2.5607E-04`/`2.0289E-04`（有りでは`4.8659E-04`/`1.6650E-03`）となった。
従って現diagonal preconditionerは明確な悪化要因であり、特にextra statesをinner 160で
約8.2倍悪化させる。しかし無しでも`1E-8`より約4桁大きく、inner 112以降は非単調なので、
原因をpreconditionerだけに限定できない。残る候補はW/P metricの悪条件性とLOBPCG
search-space更新である。checkpoint/manifest/RTはpublishされていない。

### 2026-07-22 search metric-mode診断

`e7ca9ca`でselected iterationの`Z^dagger S Z` spectrumと
`V^dagger(Z^dagger R)` weightを診断した。fresh無preconditioner B=6 run
`stage2d_wpw_runs/20260722_task16_search_metric_no_precondition_b6/run.log`は同じ
residual historyと`info=40`を再現した。

480-column search metricのeffective rankはinner 32/96/160で437/292/213まで低下し、
search spaceの冗長化は明白である。しかしinner 160のdiscarded residual fractionは
occupied `1.0231E-03`、extra `1.6596E-03`、cutoff直上のlowest retained decadeも
occupied `9.3282E-03`、extra `3.6089E-02`に留まる。従って残差weightの99%弱
（occupied）および96%以上（extra）はcutoff近傍外のretained方向にあり、metric modeの
discardだけでは`O(1E-4)` plateauを説明できない。

現在の証拠は、physical occupied-W spanの不足よりLOBPCG recurrence/restartの手続き不良を
優先して示す。次はsearch履歴を毎回3-blockで累積する現更新と、明示restart
（例えば`Z=[Q,R]`または一定間隔で`search=0`）を同一基底上で比較する。

### 2026-07-22 完全search restart比較

`54f334b`でdefault=`y`の`yn_dg_wpw_search_history`を追加し、明示`n`では各reduced
update後にsearch blockを直接zero化して、次反復を`Z=[Q,R,0]`とした。無効時は捨てる
history `matmul`自体を実行しない。preconditioner有無、fixed-H、continuation、通常の
production algebra routeへ同じcontrolを伝播した。関連テスト、MPI fixture、全体
MPI/EigenExa buildが通り、コードレビューはCritical/Importantなしだった。

fresh B=6、無preconditioner、完全restart runは
`stage2d_wpw_runs/20260722_task17_search_restart_no_precondition_b6/run.log`。同じseedを再現し、
inner 160で同じ`info=40`となった。履歴ありTask 16に対し、inner 96のoccupied/extraは
`1.1676E-03`/`1.8784E-03`（履歴あり`3.6079E-04`/`1.1093E-03`）、inner 160は
`8.0773E-04`/`1.2379E-03`（履歴あり`2.5607E-04`/`2.0289E-04`）である。完全restartは
最終occupiedを約3.15倍、extraを約6.10倍悪化させた。

完全restartのeffective rankはinner 32/96/160で320/315/305、履歴ありは437/292/213
なので、冗長性を抑えても残差は改善しない。従って累積search historyはplateau原因ではなく、
むしろ収束に必要である。現時点ではphysical occupied-W空間が悪いという証拠より、
retained Rayleigh--Ritz/operator recurrenceまたはcorrection direction生成の手続きに問題が
残るという証拠が強い。ただし単純な完全restartは対策にならない。次は履歴を保持したまま、
reduced solve後の明示残差最小化・Ritz pairの再計算整合性を診断する。cutoff、tolerance、
basis、publication gateは変更しない。

### 2026-07-23 Ritz update整合性診断

`f9fe8b4`で、post update 31/95/159の`Q=ZC`、`HQ=HZC`、`SQ=SZC`から得る
physical Ritz residualを保持し、次のdirect H/S適用32/96/160と同一stateの残差vectorを
比較した。同時に`Hred*C-Sred*C*lambda`のstate別reduced residualと
`Q^dagger S Q-I`を測定した。全rankでpending transactionとfinite/Gram検証を揃え、
solver state、収束判定、operator call回数は変更していない。2-rank runtime fixtureは
3組全てを実行し、関連テスト、全体build、二回のレビューはblockerなし。

fresh B=6、無preconditioner、履歴ありrunは
`stage2d_wpw_runs/20260723_task18_ritz_consistency_b6/run.log`。31→32、95→96、159→160の
reduced occupied/extra residualはそれぞれ`3.21E-08/7.37E-08`、
`3.72E-09/9.71E-09`、`1.99E-09/1.48E-09`であるのに対し、physical residualは
`1.97E-03/4.04E-03`、`3.57E-04/1.03E-03`、`2.56E-04/2.09E-04`だった。
一方、post/direct残差vectorの相対defectは`3.29E-12/9.71E-13`、
`1.29E-11/2.21E-11`、`3.90E-11/1.02E-10`、post metric orthogonalityは
`1.62E-11`、`2.34E-12`、`1.98E-13`である。

従ってstale H/S image、Ritz update、次反復direct operatorとのrecurrence不整合はplateau
原因ではない。reduced problemはphysical残差より4--5桁以上高精度に解けており、残差は
現在のRayleigh--Ritz search spanの外にある。現時点の証拠はoccupied-W物理表現の不足より、
retained residual correction directionを十分速く取り込めない手続き、特に適切なmetric-aware
preconditionerの欠如を優先して示す。既存diagonal preconditionerは悪化が実測済みなので、
それを戻さず、S-metricを考慮したcorrection equation/preconditionerを次に設計する。

### 2026-07-23 fragment-local metric block correction比較

`803f6df`でdefault-offの`yn_dg_wpw_metric_preconditioner`とcollective production adapterを
追加した。既存diagonal controlとの同時`y`は拒否し、fixed-H/continuationだけを
metric/diagonal/noneに分岐する。通常production algebraはcallback-freeのままである。
adapterはbounded-operator communicatorで固有値/RHS metadataをcollective検証してから
既存fragment-local block inverseへ委譲し、因子のepoch/fingerprint拘束を保持する。
2-rank実fixtureはexact block solve、interface lambda変更後の再利用、operator rebuild後の
stale拒否を確認した。関連回帰、MPI/EigenExa build、レビューはすべてPASSし、
Critical/Importantはなかった。

fresh B=6 runは
`stage2d_wpw_runs/20260723_task19_metric_block_b6/run.log`。Task 16と同一seed、
diagonal=`n`、metric=`y`、history=`y`で8 rank実行した。fragment blockはdimension 138、
condition `4.5069E+05`。inner 32/96/160のoccupied/extraはそれぞれ
`7.6809E-05/6.2909E-03`、`2.4120E-06/1.2688E-03`、
`1.9816E-06/9.5958E-04`だった。Task 16に対してoccupiedは約25.6/149.6/129.2倍改善したが、
extraは約1.56/1.14/4.73倍悪化した。effective rankは337/210/175（Task 16は
437/292/213）。inner 160のpreconditioned/residual比はoccupied`1.5280E+02`、
extra`8.4597E+02`で、局所逆によるextra過増幅が明確である。

Ritz post/direct defectは最大`1.81E-09`、metric orthogonalityは
`1.97E-12`から`1.08E-13`で、既存更新の
整合性は保たれた。fixed-Hは`info=40`、publicationなし。これは局所
`S_block^{-1}r` block-Jacobi近似の失敗であり、metric-aware correction一般やoccupied-W
基底空間の失敗とは解釈しない。次はextra stateの過増幅を抑えるregularization/state別適用、
または`H-epsilon S`を近似するmetric-consistent correctionを同一条件で設計・比較する。

### 2026-07-23 S-orthogonal PW complement比較

`0087f37`でdistributed `A=S_WW^{-1}S_WP`と
`T=[[I,-A],[0,I]]`を実装し、`db67872`でdefault-offかつfixed-H限定の
`yn_dg_wpw_s_orthogonal_pw` routeを統合した。H/S actionは`T^dagger H T`/
`T^dagger S T`で包み、seed・density境界の係数mapを明示した。metric-only fingerprintは
H/interface更新で維持し、S/ID変更ではcollectiveにstale拒否する。関連fixture、source
contract、MPI/EigenExa build、コードレビューはPASSした。

Task 16 restartからfresh B=6を8 rankで実行し、diagonal/metric preconditioner=`n`、
history=`y`、complement=`y`とした。complement solve residualとcross-metric defectは
ともに`6.4725E-13`、raw PW 1024列を削除せずnumerical rankは976、係数round-trip defectは
`4.2352E-22`だった。

inner 32/96/160のoccupied/extra residualは
`1.8708E-03/3.4680E-03`、`2.9849E-04/8.7240E-04`、
`2.0653E-04/1.5936E-04`、effective rankは`426/284/202`である。Task 16比の改善は
約1.05--1.27倍に留まり、`info=40`、publicationなしだった。Ritz post/direct defectは
最大`1.29E-10`、metric orthogonalityは最大`6.80E-12`で、座標mapと更新整合性に異常はない。

matching untransformed runとのselected generalized eigenvalueの絶対差は
`1.13E-08/8.33E-08/2.66E-08/1.61E-06`で、非収束runの`O(1E-4)` residualと整合する。
従ってS-orthogonal complementはmetric overlapを除去してもplateauを質的に改善せず、
span-preservingな座標変換だけでは不十分である。normal SCF/checkpoint/RTへは昇格させず、
次はstate-partitioned/regularized correctionまたは`H-epsilon S`型の
metric-consistent correctionを設計する。

## 新セッション開始文（コピー用）

```text
/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md
を読み、executing-plansを使ってoccupied-W bootstrap計画のTask 6から開始してください。
dirty worktreeと既存のstage2d_wpw_runsを保護し、最初にworktree/branch/HEAD/status/diff-checkとfocused verificationを確認してください。
commit・push・PRは新たに明示指示するまで行わないでください。
```
