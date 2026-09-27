# DC-HSE・MLWF・ACE：実装と測定結果

更新：2026-09-27。対象ブランチ：`dc-hse-mlwf-ace`。最新変更：HSE/LCFOの制御をnamelistへ統一。前版：`3c4ff577`（半径診断のループ整理）。半径namelist：`e6c54a87`。分散seedのメモリ測定：`133e8b37`。

**時間発展はTaylor4。局所交換と2段階MPIは実装済みですが、直接WF伝播の疎な局所化は未完了です。** このノートを開発状況の入口とし、詳細な時系列記録・図・数値データを下記にまとめています。

- [最新状況・検証範囲・富岳パッチ](#latest-status)
- [Si 3D弱スケーリング入力・半径namelist](docs/reports/si-3d-weak-scaling/README.md)（4³/6³/8³/10³、未実行）
- [Diamond 3D弱スケーリング入力](docs/reports/diamond-3d-weak-scaling/README.md)（4³/6³/8³/10³、未実行・大規模メモリ制約あり）
- [富岳：通常のCMakeビルド](#fugaku-build)
- [現在の実装と制約](#implementation)
- [最新：seed QRのMPI分散](#distributed-seed)
- [root係数の二重保持を除去](#streamed-seed)
- [Gamma seedの2次元化](#seed-memory)
- [初期MLWFのリンクメモリ削減と本番制約](#initial-memory)
- [交換FFTの計画最適化](#fft-measure)
- [Gram検査の演算・通信削減](#packed-gram)
- [交換ソースの支持領域処理](#source-support)
- [前回：Diamond 32・64・128原子の弱スケーリング](#weak-scaling)
- [局所積分半径と精度](#radius)
- [ACE・U輸送と高速化](#reuse)
- [検証と残る課題](#remaining)
- [入力・再実行・詳細記録](#records)

<a id="latest-status"></a>

## 最新状況：2026-09-27

### 富岳DC初期化のTofuマッピングSTOP

4×4×4/MPI64/OMP12で報告された`logical error: get # of process/node`は、
フラグメント内通信グループのノード内rank数を全体Tofu形状から求めた値と
比較して停止する処理に対応する。未対応Tofu次元のfallback判定もこの比較の後だった。
[修正パッチ](tools/patches/fugaku-dc-tofu-fallback.patch)は、1D/2D等の未対応次元と
DCフラグメントを既存の通常マッピングへ戻す。MPI64/OMP12を変更する対策ではない。
全体系のDC初期化（内部`yn_dc='t'`）と通常3Dの最適化は維持する。

Fujitsu APIを模擬して実際のマッピング手続きをコンパイルするテストで、旧版の同一STOPを
再現。修正後はDCの1D/3D、全体系の1D、通常3Dの4ケース成功。GNU通常ビルド成功。
富岳実機でのビルド・再実行は未確認。適用は
[before](tools/patches/fugaku-dc-tofu-fallback-before.sha256)照合→dry-run→適用→
[after](tools/patches/fugaku-dc-tofu-fallback-after.sha256)照合→再ビルド。
この1ファイルの修正はnamelist移行と独立して適用できる。


**HSE/LCFOの追加制御はすべて`&functional`へ統一しました。**
`SALMON_LCFO_RT_MLWF=1`などを設定する必要はありません。環境変数による
アルゴリズム制御は削除し、`yn_hse_wannier='y'`をGS・LCFO RT共通の有効化にしました。
LCFO RTでは併せて`yn_hse_lcfo_rt='y'`を指定します。
ACE/U間隔、直接WF伝播、FFT、分散seed、診断出力もnamelistです。
[設定一覧](docs/inputs/lcfo-rt-development.md)。OMP/MPI/BLASの標準実行環境設定は別です。
旧測定記録中の環境変数名は当時の条件であり、現在の実行手順ではありません。


| 項目 | 最新の状態 | 確認範囲 |
|---|---|---|
| 制御方法 | HSE/LCFOの追加設定をnamelistに統一。旧カスタム環境変数は読まない | 小規模MPIで旧バイナリと半径3/ACE4/U2の電流・エネルギー差0。旧環境変数混入の影響なし |
| 半径指定 | `&functional hse_lcfo_wf_radius`、常にbohr。0=全範囲（既定）、正値=固定半径、負値は拒否。環境変数fallbackなし | 小規模MPIで旧指定と電流・密度・エネルギー差0 |
| 99.9% Warning | 初期MLWFの各球内ノルム/全ノルムが0.999未満なら通知。自動半径変更・再規格化なし | 3 bohrの試験で最小0.9440599349375678、Warning発生。全範囲1、Warningなし |
| 診断出力 | `lcfo_mlwf_radius.dat`に各WFの保持率とprotectedフラグ。protected WFは従来どおり切らない | 初期時点のみ。RT全時刻や電流精度を保証する指標ではない |
| 不要処理の削減 | 全範囲時のWF再走査・追加MPI集計を削除。半径判定/二乗をループ外へ移し、非rootの保持率計算も省略 | OMP1/2/4、ビルド、Taylor4・フラグメント×軌道MPI・ACE/U再利用の回帰試験成功。速度/RSSの追加実測なし |
| Si 3D入力 | 4³/6³/8³/10³のGS/RT、座標、Si擬ポテンシャル、設定・適用パッチ | 静的検証・アーカイブ展開後の再生成に成功。大規模ジョブ未投入 |
| 富岳 | 累積seedメモリ版のコンパイル・リンク成功を利用者ログで確認 | `-Nalloc_assign`修正後、半径namelist、ループ整理、namelist統一の実機ビルド/数値は未確認 |

Si入力はa=10.26 bohr、core16³、buffer各方向8点、RT dt=0.02、16step、
半径9 bohr、ACE1/U1。原子数512/1728/4096/8000、MPI64/216/512/1000。
9 bohrで初期Si3D-WFの99.9%を満たすか、SCFが収束するかはまだ測定していません。
短時間の性能・動作確認用で、収束した誘電関数の計算ではありません。

**Si 8³のメモリは未実測。** 占有8192軌道で複素密行列1枚は1 GiB。
初期局在化rootのリンク6枚＋U＋勾配＋勾配作業配列＋ローカル全占有WFの
既知容量小計は約9.5 GiB。基底、他の軌道配列、Hamiltonian、一時配列、MPI領域は
この外側にあるため、全体ピークの上限ではありません。GS/RTの別段階でさらに増える
可能性があり、同一ノードの他rank分も加算されます。交換半径の縮小や今回の診断整理では、
この密な初期局在化配列は小さくなりません。

### 富岳のGit管理されていないソースへの適用

| ソースの状態 | 適用する差分 |
|---|---|
| 初期リンク削減前（利用者の4ファイルSHAと照合済み） | [累積メモリパッチ](tools/patches/gamma-memory-fugaku-verified.patch)。[before](tools/patches/gamma-memory-fugaku-before.sha256)→[after](tools/patches/gamma-memory-fugaku-after.sha256)を確認 |
| allocatable代入の自動確保が無効 | [富岳設定](tools/patches/fugaku-alloc-assign.patch)。正確なオプションは`-Nalloc_assign`（アスタリスクなし） |
| 累積メモリ版まで適用済み | [半径namelist](tools/patches/lcfo-radius-namelist.patch)。[before](tools/patches/lcfo-radius-before.sha256)→[after](tools/patches/lcfo-radius-after.sha256)を確認 |
| 半径namelist適用済み（e6c54a87/6e07618d相当） | [ループ整理](tools/patches/lcfo-radius-loop-cleanup.patch) |
| ループ整理適用済み（3c4ff577/eb23bfea相当） | [namelist統一](tools/patches/lcfo-namelist-only.patch)。[before](tools/patches/lcfo-namelist-only-before.sha256)→[after](tools/patches/lcfo-namelist-only-after.sha256)を確認 |

各差分はソース直下で`patch --batch --forward --fuzz=0 --dry-run -p1`に成功した場合だけ
本適用します。失敗したら後続を実行しません。半径パッチのafter照合はループ整理を適用する
**前**に行います（同じ2ファイルが次段で変わるため）。namelist統一はその後です。適用後はCMake再設定・再ビルド。
古い`gamma-memory-from-pre-root.patch`や3枚の手順は、SHAで特定した富岳旧版には使用しません。

[Si入力・再実行手順](docs/reports/si-3d-weak-scaling/README.md) ／
[半径機能の検証データ](docs/reports/si-3d-weak-scaling/verification.json) ／
[全測定と変更の時系列](docs/reports/diamond64-mlwf-support/report.md)

## 半径namelistと初期保持率のWarning

`e6c54a87`で`&functional hse_lcfo_wf_radius=9d0`を追加（常にbohr）。
0は全範囲（現在の既定値）。旧環境変数による指定は廃止しました。
初期MLWFで球内ノルム/全ノルムが99.9%未満のWFがあればWarningを出し、
`lcfo_mlwf_radius.dat`へ各WFの保持率を保存。半径は自動変更しません。
小規模MPIでnamelist/旧指定の電流・密度・エネルギーが一致し、Warning分岐を検証。
[Si入力一式・適用手順](docs/reports/si-3d-weak-scaling/README.md)／
[検証結果](docs/reports/si-3d-weak-scaling/verification.json)。富岳でのこの追加機能は未検証。

<a id="fugaku-build"></a>

## 富岳：通常のCMakeビルド

[公式マニュアル](https://salmon-tddft.jp/webmanual/current/html/install_and_run.html#build-and-install)と同じ入口を使えます。ソース直下から：

```sh
mkdir build
cd build
python3 ../configure.py --arch=fujitsu-a64fx-ea --enable-scalapack
make -j 8
```

実行ファイルは `build/salmon`。インストール先を指定する場合はconfigure.pyに `--prefix=/absolute/path/to/install` を追加し、続けて `make install`。今回追加したHSE・MLWF・ACEも同じ手順で組み込まれます。FFTW/Libxcは利用可能なものをリンク検査し、なければ対象コンパイラで自動ビルドします（初回ダウンロードにはネットワークが必要）。直接CMakeを呼ぶ場合の富岳自動選択も残しています。

ローカルでは設定選択の自動試験と、追加のライブラリ指定なしのMPI/HSE・HSE無効ビルドを確認。**富岳では累積seedメモリ削減版までコンパイル・リンク完了を利用者ログで確認。GS/RT・3次元の数値実行は未確認**です。ログでallocatable代入の自動確保が無効と判明したため、実行前に[Fortran設定の追加パッチ](tools/patches/fugaku-alloc-assign.patch)を適用し再ビルドします。この追加設定の実機確認はまだです。[手順・設定の優先順位・実機検証範囲](docs/hse-platforms.md#fugaku)

<a id="implementation"></a>

## 現在の実装と制約

| 項目 | 実装内容・対応範囲 |
|---|---|
| DC-SCF HSE | フラグメント＋バッファー＋周期境界で遮蔽交換を計算し、LCFO基底・初期状態を構築。複素数経路に対応 |
| Native LCFO RT | 固定・直交LCFO部分空間にHamiltonianを射影。密度・Hartree・半局所XC・擬ポテンシャル・電流は既存ルーチンを使用 |
| 時間積分 | 既存のTaylor4と予測・修正。今回の直接WF開発ではPT-CNを導入していない |
| MLWF | 初期局在化後、polar分解で占有空間のUを輸送。毎回のMV最適化を避ける |
| 局所交換 | 初期WF中心からの3次元最小像距離でソースWFを球状マスク。密度・Hartreeは切らず、WF再規格化なし |
| ACE | 係数空間で構築・作用。物理ステップ間の再利用を実装。impulse初回は必ず再構築 |
| MPI | LCFO RTはフラグメント×軌道の2段階MPI。1フラグメントあたり空間1ランク。DC-HSE GSの軌道MPI制約とは別 |
| 分散処理 | 交換・ACEの局所基底行、必要なWF列のhalo通信。U輸送にはなお全占有軌道依存が残る |
| 厳密な仕事量削減 | 必要なWF列・球内格子・厳密に非ゼロの基底ブロックだけ再構築。厳密ゼロの交換ペアFFTを省略 |
| 直接WF係数Taylor4 | `&functional yn_hse_lcfo_direct_wf='y'`で試験可能。係数を保持してTaylor多項式を蓄積し、受理時にWFへ回転。伝播範囲の切断はまだない |
| 係数の直接受け渡し | ACE入力の再射影・出力の再構築・Taylor側の再射影を省略。予測・修正込み1ステップで24回の格子/基底行列積を削減 |

実験的LCFO RTは**Gamma・非スピン分極・直交セル・完全占有の固定占有数**に限定。一般のk点RT、restart/checkpoint、フラグメント内の追加空間分割はこの経路で未対応です。通常SALMONやDC-SCFの対応範囲と混同しないでください。

有限半径のソースマスクは非変分的近似です。エルミート性と全電子数保存だけでは局所連続性・交換電流との整合性まで保証しません。現時点で収束した誘電関数を得たとは主張していません。

実装詳細：[LCFO RT開発仕様](docs/inputs/lcfo-rt-development.md)、[HSE入力](docs/inputs/hse.md)、[ビルド](docs/hse-build.md)。

<a id="distributed-seed"></a>

## 最新：seed QRのMPI分散

`133e8b37`で分散seed QR（現在は`yn_hse_lcfo_seed_distributed='y'`）を追加。ScaLAPACKのピボット付きQRを使い、係数を所有rankからQR所有rankへ直接転送します。rootに全体QR配列を置かず、snapshotも列ごとに保存。**既定値は'n'（従来経路）**で、MPI＋ScaLAPACKが必要です。

seed単独、合成複素係数・基底数=16×占有数、OMP/BLAS各1。各条件は独立プロセスで1回ずつ測定。rootのピークを分散する変更であり、全rank合計の削減率ではありません。

| 占有数・MPI | root peak RSS MiB（旧→新） | 非root peak RSS MiB（旧→新） | seed秒（旧→新）／速度比 | U最大差 |
|---|---:|---:|---:|---:|
| 256・2 | 49.17→38.42 | 23.44→32.70 | 0.261→0.221／1.182 | 0 |
| 512・2 | 141.27→102.55（27.4%減） | 50.78→83.81 | 2.120→1.848／1.147 | 0 |
| 512・4 | 126.31→71.84（43.1%減） | 34.91–35.58→52.36–52.92 | 2.215→1.614／1.372 | 0 |

| C128/MPI16、同一GS・R6・ACE1/U1・dt0.02・16step | 結果 |
|---|---|
| 初期snapshot／リンクsnapshot | byte単位で一致。U・中心・norm・広がり・勾配差0 |
| 最大電流差／E∞（基準ピーク規格化） | 1.97e-17 a.u.／1.87e-10% |
| 終端電流差DT（基準ピーク規格化） | −1.01e-10% |
| 最大密度差／出力エネルギー差 | 5.00e-15／0 |
| RT秒・速度比（旧/新） | 53.377→55.984、0.953 |
| 背景CPU平均（%） | 242.3→228.7 |

RTは各1回で負荷差があり、速度改善の主張はしません。QRは初期化時だけです。一般には同率pivotの選択順が変わる可能性があるため、今回の一致を全入力へ一般化せず選択式を維持します。

MPI1/2/3/4、空のroot・不均等行分割・端数ブロック・特異行列、非対応ビルドの拒否を検証。分散seedを有効にしたnative Taylor4、フラグメント×軌道MPI、局所半径、ACE/U再利用も通過。**富岳での新経路のビルドと数値計算は未検証**です。密SVD/U/ACEとrootの6リンクは残り、8³・10³本番の保留は継続します。

[詳細](docs/reports/diamond64-mlwf-support/report.md) ／ [測定データ](docs/reports/diamond64-mlwf-support/distributed-seed-memory-results.json) ／ [適用用差分](tools/patches/distributed-gamma-seed.patch)

富岳側の実ファイル4個のSHA-256を照合し、初期リンク削減前の版と特定。[実ファイル照合済み累積パッチ](tools/patches/gamma-memory-fugaku-verified.patch)と[適用前](tools/patches/gamma-memory-fugaku-before.sha256)・[適用後](tools/patches/gamma-memory-fugaku-after.sha256)の検証表を使用します。以前の3枚とpre-root用差分はこの版には適用しません。経緯と手順は詳細ノート末尾に追記。

<a id="streamed-seed"></a>

## 前段：root係数の二重保持を除去

`64d147d8`で全係数をQR配置へ直接集約し、元配置の`full_coeff`を廃止。snapshotは列ごとに書き出し、QR後は選択行だけを各rankから回収してSVDへ渡します。QR行列は行回収前に解放します。

MPI2・占有512/基底8192のseed単独peak RSSは **root204.25→141.27 MiB（30.8%減）**、非root50.19→50.64 MiB。Uと係数snapshotは一致。前段の1プロセス測定とは条件が異なり、削減率は合算しません。

| C128/MPI16、同一GS・16step | 結果 |
|---|---|
| 初期snapshot／リンクsnapshot | byte単位で一致 |
| 最大電流差／密度差／出力エネルギー差 | 1.46e-17 a.u.／4.00e-15／0 |
| RT秒・速度比（旧/新） | 52.042→56.467、0.922（各1回・背景負荷差あり） |

この標準経路ではrootにQR行列1枚は残り、QR自体は未分散。密SVD/U/ACEも残るため、RT全体や8³のピーク削減率・本番可否を保証する結果ではありません。今回の変更は富岳では未検証です。

[詳細](docs/reports/diamond64-mlwf-support/report.md) ／ [測定データ](docs/reports/diamond64-mlwf-support/streamed-seed-memory-results.json) ／ [適用用差分](tools/patches/streamed-gamma-seed.patch)

<a id="seed-memory"></a>

## 前段：Gamma seedの2次元化

`841bafe6`で全係数のreshapeコピーを除去し、QRの大配列を解放してからSVD行列を確保。seed単独のpeak RSSは占有512・基底8192で **220.78→156.55 MiB（29.1%減）**。Uは全要素一致しました。これは合成係数の単独測定で、RT全体や富岳8³のピークではありません。

| C128/MPI16、同一GS・16step | 結果 |
|---|---|
| 初期snapshot／リンクsnapshot | byte単位で一致 |
| 最大電流差／密度差／出力エネルギー差 | 2.12e-17 a.u.／3.00e-15／0 |
| RT秒・速度比（旧/新） | 53.122→55.642、0.955（各1回・負荷変動あり） |

初期化のメモリ削減であり、RT速度改善の主張はしません。rootの全係数＋転置QR配列、全体ピボット選択、密U/ACEは残るため、8³・10³の本番は引き続きピーク検証待ちです。

[詳細](docs/reports/diamond64-mlwf-support/report.md) ／ [測定データ](docs/reports/diamond64-mlwf-support/gamma-seed-memory-results.json) ／ [適用用差分](tools/patches/gamma-seed-memory.patch)

<a id="initial-memory"></a>

## 初期MLWFのリンクメモリ削減と本番制約

`bee87f63`ではrootのGamma局在化を更新済みリンクの直接評価に変更し、リンク関連の保持量を18 N²→6 N²複素数へ削減。係数のroot集約バッファは64列に制限し、全係数もリンク生成前に解放します。占有1024のGamma初期評価単独でpeak RSS **364→156 MiB（57%減）**。これはRT全体のピークではありません。

| C128/MPI16、同一GS・16step | 変更前→後 | 検証 |
|---|---:|---|
| RT秒 | 52.222→54.895（速度比0.951） | 各1回・負荷変動あり、速度改善の主張なし |
| 最大電流差 | 1.32e-17 a.u. | E∞=1.25e-10%（基準ピーク規格化） |
| 密度／出力エネルギー差 | 3.00e-15／0 | 初期U・WF中心も丸め差の範囲 |

**大規模本番はまだ保留。** 8³のリンク容量は18→6 GiB、10³は68.7→22.9 GiB。ただしseedの全係数集約・QR/SVDとRTの密U/ACEが残るため、root全体のピークを保証しません。前回の非rootリンク削減と今回のroot削減を分けて記録しています。

[詳細](docs/reports/diamond64-mlwf-support/report.md) ／ [root測定データ](docs/reports/diamond64-mlwf-support/root-gamma-memory-results.json) ／ [前回のリンク構築測定](docs/reports/diamond64-mlwf-support/initial-links-memory-results.json)

<a id="fft-measure"></a>

## 交換FFTの計画最適化

`2c996869`では `yn_hse_lcfo_fft_measure='y'` でFFTの実測計画を選択可能。初回に専用scratchで計画を作り、反復中は再利用。物理近似は追加していません。既定0（ESTIMATE）を維持します。

| 系 | 通常→実測計画 RT秒 | 速度比 | 最大電流差 a.u. |
|---|---:|---:|---:|
| C64/MPI8 | 33.140→33.096 | 1.001 | 3.48e-17 |
| C128/MPI16 | 52.831→51.472 | 1.026 | 1.91e-17 |

同一GS・同一バイナリ、R6/ACE1/U1/16steps、各1回・逐次実行。密度差最大5.00e-15、出力エネルギー差0。C64基準の弱効率62.73→64.30%は参考値。負荷変動があり、大幅・確定的な改善は確認できません。計画方式の選択だけでは残る全占有軌道処理を解決しません。

[初期費用を含む詳細記録](docs/reports/diamond64-mlwf-support/report.md) ／ [測定データ](docs/reports/diamond64-mlwf-support/fft-measure-results.json)

<a id="packed-gram"></a>

## Gram検査の演算・通信削減

`9508f8e9`では、直接WFの直交性検査をHermitian上三角に限定。C128で通信要素65536→32896、単独検査1.792→0.936 ms（約1.9倍）。検査頻度・閾値は同じで新しい近似はない。漸近次数と作業メモリ全体はほぼ変わらない。

同じGS・R6・ACE1/U1・16steps・MPI16のRT比較は47.113→50.944秒（各1回）、電流差最大6.60e-18 a.u.、密度差5.00e-15、出力エネルギー差0。**全体の高速化は未確認。** Gram削減は17検査で約0.015秒の規模で、背景負荷も変動した。主要な弱スケーリング課題は全占有列の伝播・U輸送・交換FFTに残る。

[詳細記録](docs/reports/diamond64-mlwf-support/report.md) ／ [単独検査・C128比較データ](docs/reports/diamond64-mlwf-support/packed-gram-results.json)

<a id="source-support"></a>

## 交換ソースの支持領域処理

`8fa826da`では、厳密に非ゼロのソース格子点だけでペア密度の積と交換結果の加算を行います。広い支持領域は従来演算を維持。積の格子点数45.6%、加算34.6%減。FFT格子・実行ペア7424・R6・Taylor4は変更せず、新しい切断近似は加えていません。

| 系 | 旧→新RT秒 | 旧/新比 | 旧→新弱効率 % |
|---|---:|---:|---:|
| C32/MPI4 | 21.388→20.506 | 1.043 | 100→100 |
| C64/MPI8 | 40.171→32.216 | 1.247 | 53.2→63.7 |
| C128/MPI16 | 52.892→52.032 | 1.017 | 40.4→39.4 |

同一GS、直接係数Taylor4、R6、16steps、ACE1/U1、OMP1/BLAS1。各1回、各バージョンのC32時間で効率を規格化。背景負荷の変動があり速度比は参考値。**C64では改善した一方、C128の弱効率改善は確認できません。** 電流差最大4.36e-17 a.u.、密度差最大4.00e-15、出力エネルギー差0。全占有軌道の伝播・密なU・FFT本体の課題は残ります。

[図](docs/reports/diamond64-mlwf-support/source-support-scaling.png) ／ [条件・測定データ](docs/reports/diamond64-mlwf-support/source-support-results.json) ／ [詳細記録](docs/reports/diamond64-mlwf-support/report.md)

<a id="weak-scaling"></a>

## 前回：Diamondの弱スケーリング

`88385397`、直接WF係数Taylor4。各ランク8原子、コア16³、バッファー8×0×0、R=6 bohr、dt=0.02 a.u.、16ステップ、ACE1/U1、OMP1/BLAS1、FFT batch1。4×1×1→8×1×1→16×1×1の**1次元拡張**です。MPI4→8→16、軌道MPI分割なし。C32は新規DC-HSE GS（1147反復、残差9.9556585e-8）、C64/C128は収束済みGSを再利用しました。

弱効率は **100×T32/TN**、理想100%。GS作成・初期MLWF局在化を除くRT反復時間を比較します。

| 原子数 / MPI | RT秒 | T/T32 | 弱効率 % | Hamiltonian秒 | WF source秒 | 交換構築秒 | ACE構築秒 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 32 / 4 | 22.856 | 1.000 | 100.0 | 1.379 | 0.558 | 19.650 | 0.026 |
| 64 / 8 | 39.193 | 1.715 | 58.3 | 3.588 | 0.808 | 31.744 | 0.155 |
| 128 / 16 | 53.758 | 2.352 | 42.5 | 9.498 | 1.546 | 34.320 | 0.642 |

WF source/交換/ACEは初回を除いた34構築のrank最大時間の合計。Hamiltonianとは計測範囲・集計方法が異なるため、列を単純加算してRT時間とはしません。

![Diamondの弱スケーリング](docs/reports/diamond64-mlwf-support/weak-32-64-128.png)

各フラグメントの局所基底64行、halo基底192行、source WF60本、FFT実行ペア7424/11520は全サイズで一定。一方、各空間ランクの全占有軌道列は64→128→256本に増えます。Uのpolar分解は128軌道以上で分散経路へ切り替わります。局所仕事量が一定でも、全占有列の伝播・通信・同一ノード内の資源競合が残ります。

**各サイズ1回、同一ホスト、CPU固定なしの参考測定です。** 平均背景CPUは66.1/87.2/84.1%。異なる性能クラスのCPUを含むホスト上の測定であり、均一ノードを増やすクラスタの弱スケーリングとは区別します。3次元で各辺L倍に伸ばせば原子数はL³倍ですが、現実装の速度改善がL³倍になることを意味しません。

[詳細・機械可読データ](docs/reports/diamond64-mlwf-support/weak-32-64-128-results.json) ／ [図PDF](docs/reports/diamond64-mlwf-support/weak-32-64-128.pdf)

<a id="radius"></a>

## 局所積分半径と精度

Diamond C64の半径試験は、初期実装`b1773c78`、MPI8/OMP2/BLAS1、dt=.02、64ステップ、ACE1。上の最新16ステップ性能測定とは別です。各半径の励起・無励起電流の差を ΔJ とし、全範囲を基準に比較しました。

- E∞ = 100 maxₜ|ΔJ_R−ΔJ_full| / maxₜ|ΔJ_full|
- D_T = 100[ΔJ_R(T)−ΔJ_full(T)] / maxₜ|ΔJ_full|

| 半径 bohr | 破棄WFノルム最大 % | E∞ % | D_T % | impulse全時間 秒 | 全時間比 |
|---|---:|---:|---:|---:|---:|
| 全範囲 | 0 | 0 | 0 | 317.32 | 1.000 |
| 8 | 0.02783 | 0.01830 | −0.01830 | 196.50 | 1.615 |
| 6 | 0.10799 | 0.07578 | −0.07578 | 170.44 | 1.862 |
| 4 | 1.1102 | 0.64699 | −0.64699 | 150.33 | 2.111 |
| 3 | 4.2712 | 1.65013 | −1.65013 | 161.15 | 1.969 |
| 2 | 13.3871 | 1.95983 | −1.95983 | 134.76 | 2.355 |

この短時間試験では6 bohrが精査の候補、8 bohrが高精度側の比較点です。4 bohr以下では無励起電流も大きくなります。全範囲のdt半減差はE∞=0.000403%。総時間は約0.031 fsで、長時間誘電応答の許容半径を確定する試験ではありません。

Si128での初期9/8/7 bohr試験は**x方向だけの半幅**であり、現在の球状半径とは異なります。Siの「9 bohr・ACE2〜4」をDiamondの球状マスクへそのまま移植しません。[Si時系列記録](docs/reports/si128-mlwf-support/report.md) ／ [Diamond半径試験と時系列](docs/reports/diamond64-mlwf-support/report.md)

<a id="reuse"></a>

## ACE・U輸送と高速化

impulseでは最初の予測伝播前にACEを必ず再構築。滑らかなレーザーでは初期GS軌道から構築したACEを再利用できます（保存済みACE因子をファイルから読む機能ではない）。エネルギー診断の頻度は既存namelistの`out_rt_energy_step`で設定します。

U輸送間隔はACEとは独立です。間隔kでは初期と物理ステップ1,1+k,…で更新し、それ以外は保存Uを使用します。輸送基準は「最後に輸送したWF」を保存し、保持中の位相変形が基準に累積しないようにしています。

Diamond R6/16steps/ACE1、U1を基準とした比較（各1回）：

| 系 | U間隔 | RT秒 | U1/RT比 | source秒 | E∞ % |
|---|---:|---:|---:|---:|---:|
| C64 | 1 | 39.299 | 1.000 | 0.8276 | 0 |
| C64 | 2 | 38.835 | 1.012 | 0.7332 | 0.00974 |
| C64 | 4 | 39.511 | 0.995 | 0.6889 | 0.02863 |
| C128 | 1 | 60.551 | 1.000 | 1.6949 | 0 |
| C128 | 2 | 57.708 | 1.049 | 1.3483 | 0.00974 |
| C128 | 4 | 54.030 | 1.121 | 1.0505 | 0.02863 |

有限半径ではゲージ保持による追加誤差があるため、U間隔の既定値は1。[比較データ](docs/reports/diamond64-mlwf-support/u-cadence-results.json)

主な実装段階：

| コミット | 変更 | 根拠・結果 |
|---|---|---|
| `35d1d38e` | 必要WF列のみhalo通信 | C128で受信係数49152→11520/rank。自己・ノード内転送を含む数であり、ネットワーク実測バイト数ではない |
| `beccb272` | 厳密ゼロの交換ペアを省略 | FFT実行7424/11520。数値しきい値による近似なし |
| `19b9333c` | 非ゼロブロック投影とFFT batch | batch既定1。4/8で一貫した高速化なし |
| `eb13fdc6` | 球内格子・非ゼロ基底だけ再構築 | 積数94,371,840→13,238,272/fragment。source時間約29%減の単回比較 |
| `073830ef` | ACE独立のU輸送間隔 | 基準WF固定、predictor巻き戻し・101ステップ位相回帰検証 |
| `78519d44` | 直接WF係数Taylor4の密な基準 | 通常経路と数値一致、当初は約8〜14%遅い。既定は通常経路 |
| `88385397` | 係数Hamiltonian受け渡し | 重複射影・再構築を削減。下表 |

直接WF経路の変更前→後を同じ入力で比較（各1回、R6/16steps）：

| 系 | 旧→新RT秒 | 旧/新比 | 旧→新Hamiltonian秒 | 最大電流差 a.u. |
|---|---:|---:|---:|---:|
| C64/MPI8 | 42.083→40.527 | 1.038 | 5.6648→3.5785 | 1.13e-17 |
| C128/MPI16 | 64.000→52.698 | 1.214 | 15.154→9.2638 | 1.04e-17 |

最大密度差5.00e-15、出力エネルギー差0、ACE構築回数は同じ。背景負荷が異なるため、全体速度比をコードだけの効果とは断定しません。この表は上の3サイズ再測定とは別の測定系列です。[比較データ](docs/reports/diamond64-mlwf-support/coefficient-handoff-results.json)

過去のC128 285.550→57.364秒という見かけの4.98倍は、無変更部分も大きく変動したため**性能比較として無効**です。時系列記録には残しますが、成果として合算・引用しません。

<a id="remaining"></a>

## 検証と残る課題

検証済み：複素係数の密参照一致、不均等MPI行分割、MPIフラグメント×軌道、有限半径、ACE/U再利用、impulse/laserの初回処理、半dt、密度・電流・エネルギー・Gram、球内再構築の微小成分/ゼロ支持。最新の係数ACEは因子ランク2/4/6で検証。GNU Fortran15/AArch64で`-fno-tree-loop-vectorize`を使用し、HSE有効/無効のMPIビルドを確認しました。

| 残る仕事 | 現状 |
|---|---|
| 直接WFの疎な格納・Hamiltonian/ACE作用 | 全占有列・密なゲージ補正が残る。現在は正しさの基準実装 |
| 伝播半径R_prop | 未実装。交換半径R_intと独立に収束を検証する必要がある |
| ベクトルポテンシャル位相 | exp(i A·r)の補正は未実装。固定LCFO部分空間への射影は一般にunitaryでなく、基底外漏れ・周期境界を検証する |
| 光学応答 | 長時間、励起/無励起、dt、k点、buffer、局所連続性と電流整合性の収束が未完 |
| 3次元・複数ノード | ここでの性能は細長い1次元系列・単一ホスト。ノード間測定は保留 |
| その他 | 最新の直接halo fallback（ACE構築失敗時）の強制試験、GPU構成は未検証 |

<a id="records"></a>

## 入力・再実行・詳細記録

- [32/64/128入力・原子座標・再実行手順](docs/reports/weak-scaling-inputs/README.md)
- [Diamond詳細ノート・図・JSON](docs/reports/diamond64-mlwf-support/report.md)
- [Si128詳細ノート・ACE間隔と局所範囲](docs/reports/si128-mlwf-support/report.md)
- [LCFO RT設定と実装の詳細](docs/inputs/lcfo-rt-development.md)
- [単体・MPIテスト](testsuites/unit_lcfo_rt)
- [実装計画の履歴](docs/plans)

本ノートのJSONは既存測定をGitHubで読めるように整理したものです。ローカル絶対パスは`<local>/…`に置換し、バックグラウンドのプロセス一覧は除外しました。数値結果、集約負荷、入力・実行ファイルのハッシュは保持。実行バイナリ、大きなGS/WFデータは含めません。
