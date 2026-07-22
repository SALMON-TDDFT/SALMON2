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

## 新セッション開始文（コピー用）

```text
/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md
を読み、executing-plansを使ってoccupied-W bootstrap計画のTask 6から開始してください。
dirty worktreeと既存のstage2d_wpw_runsを保護し、最初にworktree/branch/HEAD/status/diff-checkとfocused verificationを確認してください。
commit・push・PRは新たに明示指示するまで行わないでください。
```
