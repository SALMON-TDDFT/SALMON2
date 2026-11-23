# Niの短時間SOC検査

[English](NICKEL_SOC.en.md)

NiをSOC用の小規模CI試験に使います。fcc基本セル1原子、SALMON本体の標準2×2×2（8点）を使用します。Pythonで独自のk点を生成しません。

## 擬ポテンシャルと条件

- PseudoDojo PBE v0.4の完全相対論NC-UPF `Ni_std.upf`、SOC・NLCCあり。
- カタログID：`pseudodojo-pbe-0.4-fr-standard-ni`。固定コミットとSHA256で照合します。
- **3s・3pセミコアを価電子に含む18電子**。原子参照配置は3s² 3p⁶ 3d⁸ 4s²です。固体中の軌道占有数をこの整数配置に固定する意味ではありません。
- 36スピノル状態、16³格子、Libxc-PBE、電子温度300 K。
- 初期波動関数は乱数シード0のランダムスピノル。磁化方向や磁気モーメントを指定した初期状態ではありません。
- 混合は `simple_dm`、係数0.2。SCFは3回だけ実行します。

格子・k点・温度・格子定数は短い接続検査用です。収束したNiの磁性、磁気異方性、光学物性を保証しません。

## 明示的に実行する

Libxc付きでビルド・禁止固有値ソルバー監査を済ませた `build.json` が必要です。ソース版では次を実行します。

```sh
python tools/test_nickel_soc.py BUILD_JSON NEW_OUTPUT --potential NI_UPF
```

`BUILD_JSON`、`NEW_OUTPUT`、`NI_UPF` は実際のパスに置き換えます。入力準備だけで終了する処理とは異なり、**このテストコマンドは短い計算を開始します**。管理ノードでは実行せず、割当済み資源かローカルの検証環境を使います。既存の保存先は上書きしません。

固定UPFの取得も明示的に許可する場合だけ、`--potential NI_UPF` の代わりに `--allow-download` を使えます。計算ノードから外部通信できない場合は、管理ノード等で[UPFを取得](PSEUDOPOTENTIALS.md)して持ち込みます。

wheelにはソース資材の `salmonpy/_runtime_source/tools/` と `tests/runtime/` に同じ試験を含めます。通常の `import salmonpy` や最初のLD/SBE計算でこのビルド・計算は始まりません。

## 判定すること

Fortranで密度行列・XC行列・混合履歴の有限性、密度行列のエルミート性、電子数18、非共線配列の独立した複製と正確な復元を検査します。1／2スレッドで実行します。

SCF3回での未収束は想定した状態として記録します。未収束GSからRTへ自動で進まず、波動関数・密度・GS再開データをディスクへ出力しません。`result.json` と `execution/test-record.json` に普通のログ・残差・判定を保存します。

**検査対象は非共線スピンの作業配列です。SCF全体やTDDFTの直接保存／再開を検証したという意味ではありません。** これらのSOC対応は別の作業として残っています。

## CI

既存の手動ワークフロー `NLCC native consistency` のLibxcビルドを再利用します。Ni専用の重複ビルド・自動物性計算・新規ワークフローは追加しません。2026-11-23にApple Silicon Macでローカル実行を確認しました。GitHub Actionsの遠隔実行はまだ行っていません。

## SCF状態の独立性とメモリ上の保存・復元

`tools/test_nickel_scf_mpi.py` は、セルと混合条件の異なる2状態の交互実行を個別実行と比較します。
18電子・36状態、300 K、SALMON標準の2×2×2 k点で、密度・占有数・固有値・軌道・残差履歴を比較します。

```sh
mpiexec -n 3 python tools/test_nickel_scf_mpi.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
mpiexec -n 6 python tools/test_nickel_scf_mpi.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

2026-11-23のMac検証では内部MPI 1・2ランクとも配列が完全一致し、復元時のSCF実行回数は0でした。
射影子構築の一時配列は呼び出しごと、非共線の可変配列はSCF状態ごとに保持します。
波動関数・密度・GS再開ファイルは保存しません。短い未収束SCFの接続試験であり、
SOCのディスク再開・RT・力・高速HSEの検証完了を意味しません。

## プロセス終了後のSCF再開

```sh
python tools/test_nickel_scf_disk.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 2
```

この手動試験は電子状態・混合履歴と比較用配列をディスクに保存します。保存先と保存の許可を先に確認します。
SCFを2回進めて保存し、別のMPIプロセスで復元してさらに2回進め、連続4回と比較します。
2026-11-23のMac検証では内部MPI 1・2ランクとも密度・占有数・固有値・軌道・残差履歴が完全一致し、
復元時のSCF実行回数は0でした。非共線配列を含むストリーム形式2を使い、旧形式へは戻しません。
同じ実行ファイル・同じMPI分割での短い未収束SCF試験です。SOCのランク数変更、RT、力、高速HSEは未検証です。

## SOCのA→電流と直接保存／復元

```sh
mpiexec -n 3 python tools/test_nickel_tddft_mpi.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
mpiexec -n 6 python tools/test_nickel_tddft_mpi.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

300 K・18電子・36状態・2×2×2 k点のNiを使う手動の接続試験です。
GSは4回だけ進めた未収束状態をメモリ内で使い、物理結果として利用しません。
2つの独立状態をx/y/z成分を持つ異なるAで8ステップ駆動し、個別実行と比較します。
ステップ4の電子状態をディスク保存し、新しい状態へ直接復元してステップ5～8を比較します。
Mac上の内部MPI 1・2ランクで電流は完全一致しました。SiのSOCなし電流も変更前と完全一致しました。

この試験は電子状態を保存します。保存先と許可を確認し、保存しない比較には `--without-checkpoint` を使います。
対象は固定イオン・全k点・CPU・PBEのNiで、非共線配列を状態ごとに保持します。
結晶の駆動はA/c、返す値は電荷電流です。SOCなしの時間発展式・更新順は変更していません。
TDDFT保存形式は3で、旧形式を暗黙変換しません。
同じプロセス内の復元に加え、以下の別プロセス再開も確認しました。MPI分割変更・
深部異常からの復帰・高速HSE・物理的収束は未検証です。WSe₂実材料の検証は別途行います。

## プロセス終了後のTDDFT再開

```sh
python tools/test_nickel_tddft_disk.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 1
python tools/test_nickel_tddft_disk.py BUILD_JSON NEW_OUTPUT_MPI --potential NI_UPF --size 2
```

この手動試験はGSとTDDFTの電子状態を保存するため、保存先と許可を先に確認します。
MPI分割ごとに、最初のプロセスでSCFを4回だけ進めた検証専用GSを一度保存し、連続8ステップの電流・電子数を記録します。
次のプロセスでは同じGSを読み込み、TDDFTをステップ4まで進めて保存・終了します。
最後のプロセスでGSとTDDFTを読み込み、ステップ5～8だけを進めます。

Mac上の内部MPI 1・2ランクで、保存・再開の各プロセスのSCF実行回数は0回、
TDDFTの読み込み時はステップ4、追加実行は4ステップでした。
電流と実測電子数は連続実行と完全一致しました。同じ実行ファイル・同じMPI分割での検証です。
初期化時の作業配列・通信・ポテンシャルの準備は行いますが、SCF反復や保存時点までの時間発展は再実行しません。
未収束GSを使う短いソフトウェア試験であり、Niの物理的な収束を保証するものではありません。


## 出力も保持する巻き戻し・再試行

```sh
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_OUTPUT --potential NI_UPF --size 2
```

未収束GSはソフトウェア比較専用としてメモリ内だけで扱い、GSや電子チェックポイントを保存しません。
毎ステップのスピン出力を指定し、ステップ後のホスト例外を発生させます。
内部MPI 1・2ランクで、電流・エネルギー・スピンの3ファイルが巻き戻し前とバイト一致し、
再試行後の電流・磁化・スピン流が連続実行と一致しました。別状態の出力は保持されます。
これは処理可能なホスト例外の検証であり、深部数値異常・MPI停止・プロセスクラッシュからの復旧ではありません。
[読み込みと保存の例](OUTPUTS.md)も参照してください。

## SOCの監査付きモードとパフォーマンスモード

```python
# model = gs.tddft_model(dt=..., max_steps=...)
audited = model.new_state('audit', execution_mode='checked')
fast = model.new_state('fast', execution_mode='performance')
```

| モード | ステップ内の処理 | 異常時 |
|---|---|---|
| `checked`（既定・監査付き） | 入力・MPI整合、波動関数・密度・Hartree・XC/SOC配列・観測値の有限性検査と状態コピー | 対応する数値異常は直前の状態へ自動復帰。電流・エネルギー・スピン出力は変更せず、`NumericalError` を返す。原因を直して同じステップを再試行できる |
| `performance` | ステップ内の診断検査・異常合意通信・保護用の状態コピーを省く | 自動復帰は保証しない。必要な操作だけ明示的な `guarded_step` やチェックポイントで保護する |

監査付きでも物理的収束を保証するものではありません。ここでの監査はステップの数値検査・状態保護であり、
禁止固有値ソルバーのソース・オブジェクト監査は両モードともビルド時に行います。
変更していない計算ライブラリの再利用時は記録されたSHA256を照合します。ジョブ開始時や各ステップにソース監査を繰り返しません。

SOCの非共線配列の検査は既存のXC段階の異常合意へ、磁化・スピン流の検査は既存の観測値段階の異常合意へまとめています。
磁化・スピン流は出力前に計算・検査し、その値を出力に使うため、同じ観測量を二度計算しません。
性能モードではこれらの追加検査と保護用コピーを呼びません。

本番用のライブラリで、次の明示的な試験を実行できます。

```sh
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_CHECKED --potential NI_UPF --size 2 --execution-mode checked
mpiexec -n 2 python tools/test_nickel_tddft_retry.py BUILD_JSON NEW_FAST --potential NI_UPF --size 2 --execution-mode performance
```

故障注入は `SALMONPY_RECOVERY_TESTING` を定義して別途ビルドした試験専用ライブラリだけを使い、
上の試験に `--deep-faults` を付けます。本番用ライブラリにはこの入口を含めません。
Mac上のNi Libxc-PBE・内部MPI 1/2で、実オーバーフローと波動関数・密度・Hartree・XC・観測値の故障を検証しました。
失敗時の状態・3出力ファイルの保持、再試行後の電流・磁化・スピン流の一致、別状態の独立性を確認済みです。
性能モードは正常時に監査付きモードと一致し、明示的に保護したホスト例外からの復元・再試行も一致しました。
これは固定イオン・全k点・CPU・k点MPIでの検証で、MPI停止・プロセスクラッシュからの復旧ではありません。
GSは短い未収束の接続試験用であり、物理的収束の検証には使いません。

補助磁場を段階的に弱めるSCF入口と1/2ランクの保存試験は、[SOCのSCFを補助磁場で整える](NATIVE_GS.md#socのscfを補助磁場で整える)を参照してください。長い固定入力で整列の挙動は確認しましたが、最終ゼロ磁場SCFの収束は未確認です。

## 1000回単位のSCF比較（手動実行）

`tools/test_nickel_scf_mpi.py --alignment-study` は、通常の短い状態独立性試験とは別に、
磁場なしと補助磁場の段階的解除を比較します。これは**計算を実行するコマンド**であり、
管理ノードで行う入力検査ではありません。計算可能な環境／ジョブ内で明示的に実行します。

```sh
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 python tools/test_nickel_scf_mpi.py \
  BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF --size 2 \
  --alignment-study --stage-limit 1000
```

18電子・36スピノル・300 K・16³格子・SALMONの2×2×2 k点を使い、
`simple_dm`の混合率0.3と任意の自動低下を明示します。乱数種0の初期条件を揃え、
初期配列の差を許容誤差とともに記録します。多スレッドでは丸め差があり得ます。
補助場はBz=0.002→0.0005→0.0001→0（原子単位）、各段階最大1000反復です。磁場なしの比較も最大1000反復です。
磁場なしと最終ゼロ磁場のSCF残差の基準は1e-9、磁気モーメントの変化は
直前の観測点から最大成分差1e-5 μBを診断の目安にします。中間の非ゼロ磁場は
整列準備として残差1e-5・磁気モーメント変化1e-3 μBを目安に次へ進みます。
この二つは`--preliminary-residual`と`--preliminary-moment-change`で指定できます。
中間段階の目安に達しても「最終SCFが収束した」とは記録しません。観測間隔も記録します。
残差が1e-6より小さくなれば1反復ずつ観測し、それ以外は25反復間隔です。
各段階の目安を満たさなければ、次の小さい磁場へ自動で進みません。
最終ゼロ磁場の基準は緩めません。同じ条件の対照は`--control-record PATH_TO_ZERO_CONTROL_STAGES_JSON`で
明示的に再利用でき、元の入力・擬ポテンシャル・実行条件を照合して由来を記録します。

`history.csv`、`progress.json`、各段階のJSON、残差CSV／任意のPNG、
比較図`scf-alignment.png`と`result.json`を保存します。
波動関数・密度・再開データはメモリ内だけに保持し、終了後の再開用には保存しません。
粗い固定入力でのSCF安定性の調査であり、k点・格子やNiの物性の収束保証ではありません。
この長い計算はCIへ追加しません。

### 2026-11-23の固定入力での観測

MPI2・各4スレッドで長いSCF比較を実行しました。磁場なしは1001回時点で残差約3.0e-2、
Bz=0.002の固定補助場では1001回時点で約4.0e-6でした。
中間段階の目安を分けると、Bz=0.002→0.0005→0.0001→0まで継続できました。
ただし、ゼロへ戻して1000回進めた2801回時点の残差は約5.9e-6で、最終基準1e-9には未達です。
磁気モーメントは約(0.2844,0.2877,0.2893) μB、直前25回での最大成分差は約1.5e-4 μBでした。
整列・段階的解除の挙動は確認できましたが、最終ゼロ磁場のSCF収束は確認できていません。
自動混合率が磁場なしでは0.01、補助場から戻した場合は0.15となるため、
これは補助場と自動調整を含む手順全体の比較であり、磁場だけの効果を分離した比較ではありません。
反復数を増やすだけで完了扱いにせず、残差停滞の切り分けが残っています。
