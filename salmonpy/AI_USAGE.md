# SALMONpy: AI向け利用・実装ガイド

LD単独利用のPython試作APIを提供しています。実際の入口は `LorentzDrude` の `respond`、`new_state`、`permittivity` です。構築方法・単位・時間位置は README.md と examples/ld_response.py を参照してください。FDTDには外部実行ファイルを使うバッチ入口があります。WORKFLOWS.mdを参照してください。結晶は共通A駆動インターフェースを持ち、内蔵SBE速度ゲージと初期状態ローダーを提供します。TDDFTはPZ／TB-mBJ／Libxc-PZ／Libxc-PBE・中点伝播のプロセス分離型接続を提供します。存在しないAPIを例示しないでください。

## 物質応答の境界

結晶の公開インターフェースはA駆動のみです。TDDFTとSBEは、内部状態と新旧Aを受けて時間発展しJを返す交換可能な物質応答として設計します。結晶にE駆動の選択肢を追加しないでください。SBEは既存の演算順を保った共有Fortranカーネルに接続しています。

E,Hを扱うFDTDとA駆動の結晶の間には結合処理を置きます。単位、ゲージ、時間位置、符号、座標回転、Jの体積規格化は現行マルチスケール実装に準拠して確定し、アダプターで明示します。Aの積分履歴は結合側、波動関数・密度行列は物質側が保持します。

ドルーデ＝ローレンツはE駆動の構成則として独立させます。結晶用A駆動APIに合わせるための数式変更はしません。

## 今回利用できるFortran境界

`src/maxwell/lorentz_drude.f90` の `s_ld_state` がP、J、旧J、旧E、係数、次の電場更新用の補助量Qを所有します。公開関数は `ld_allocate`、`ld_coefficients`、`ld_store_old`、`ld_update`、`ld_finalize` です。具体的な配列形状・呼び方は `testsuites/unit_ld/test_ld.f90` と `src/maxwell/fdtd_eh.f90` を参照してください。

`ld_coefficients` の累積先は呼び出し前にゼロ初期化します。`ld_store_old` と `ld_update` は既存演算順を保つため分割E成分を受け取ります。`rj*_fdtd_ld` は物理電流の単純和ではなく補助量Qです。

`fdtd_eh_kernel.f90` は差分演算 `eh_fd` を提供します。`fdtd_eh` は引き続き時間ループ、境界、通信、観測、再開ファイルを調整します。Python向け単一ステップAPIはまだありません。

## 変更してはいけない更新順

旧J/E保存 → E更新（LD係数を含む）→ 光源 → 旧Q注入 → E通信 → LDのJ/P/Q更新 → 線形応答 → 旧H保存 → H更新 → H通信 → TTM → 観測・保存。

LDは電流を返すだけでは分離できません。E更新係数への寄与も必要です。係数を固定したまま時間刻みを変更しないでください。独立計算には別々の状態インスタンスを使います。

## 検証と禁止事項

ビルド前に `python tools/check_forbidden_eigensolvers.py --source-root src` を実行します。構成・増分ビルド時の検査を外しません。新規実行ファイルは全アプリケーションオブジェクトを検査し、再利用時は記録されたSHA256と照合します。禁止された固有値ドライバーへのフォールバックはありません。PZHEEVDの検証には固有値、残差、直交性、ワークスペース容量を含めます。

IEEE判定は各要素へのスカラー呼び出しに限定します。GSデータを収集しません。計算結果は利用者が指定した書き込み可能な出力先へ保存し、コードツリーに格納しません。

変更後は独立LDの旧演算との比較とFDTD全体の比較を行います。単一スレッドのビット一致と並列の数値許容差検証を区別します。再開互換性や未検証の計算機での動作を推定だけで保証しないでください。

## 入力の寿命

INPUT_CONTRACT.md に従います。材料係数はモデル作成時、dtと初期条件は状態作成時、Eは各stepで指定します。`material.config` と `state.config` は確定値の読み取り専用レコードです。設定変更は新しいモデル・状態の作成で行い、private属性を書き換えません。`respond` は毎回新規状態です。設定レコードはチェックポイントではありません。

## 配布

wheelの構築と修復・OS別検証は PACKAGING.md に従います。検証済み環境とCI設定のみの環境を混同しません。wheel利用者にはFortranコンパイラを要求せず、ソース構築時のみ要求します。

構造・波形・解析の単独入口は WORKFLOWS.md に記載しています。FDTDの一ステップAPI、周期構造、斜め入射、未対応のTDDFT計算方式を対応済みとして例示しないでください。共通Acos2パルスとFDTD側のCEP規約の差はアダプターで0.25周期補正します。


## SBEの接続

`salmonpy.sbe.SBEInitialState.from_salmon` は原子単位・単一スピンの
`_k.data`, `_eigen.data`, `_tm.data` を読みます。`bands` はファイルの全バンド数、
`active_bands` は任意の伝播バンド数です。占有数は現行SBEと同じ閉殻規約です。
`new_state(dt=..., initial_A=...)` ごとに独立した密度行列を生成します。
結晶はA/c駆動のみ。中点Aで四次Taylor伝播、新しい端点AでJ評価の順を保持します。
補正次数と非局所補正は初期モデルで固定し、dtは状態作成時に固定します。
SBEの長さゲージ・バイナリ再開は未対応です。四次Taylor法なので厳密な
ユニタリー性や自動誤差制御を保証しません。APIは刻み幅を自動変更しません。


## TDDFTの接続と寿命

`salmonpy.tddft.TDDFTInitialState.from_salmon` は自己完結した入力ディレクトリを
複製し、ハッシュ付きの初期スナップショットを所有します。初期GSのバイナリ解釈は
SALMONの既存ローダーが行います。各 `new_state` は専用プロセスと作業コピーを所有し、
波動関数、Vh履歴、予測子修正子、XCや通信の共有変数を他の計算と共有しません。
`with` で解放してください。共通 `respond` は成功・失敗とも自動的に解放します。
PZ／TB-mBJ／Libxc-PZ／Libxc-PBE、周期系、固定イオン、非スピン分極、中点伝播、CPU一MPIランクが対象です。
実行ファイルは別途構築・監査し、SHA256を指定します。wheelだけでTDDFTが動くとは
説明しないでください。dt、最大ステップ数、初期Aは状態作成時に固定します。
共通JはSALMONの電子物質電流です。Maxwell側で電荷電流が必要なら符号は−Jです。
異常終了時は受理済みA履歴を再実行して復元し、復元失敗時は状態を閉じます。
実GSによる検証と合成波動関数による接続回帰を区別してください。

Libxc-PZ/PBEには `USE_LIBXC=ON` の外部SALMONが必要です。Libxc全汎関数の対応ではありません。GSとRTで汎関数・擬ポテンシャルを一致させてください。現ローダーはGS生成条件との整合を自動判定しません。PBEの検証はヘッダーでPBE・NC・スカラー相対論・NLCCありを確認したSi UPFを使用しました。詳細は docs/libxc-validation-2026-11-13.md を参照してください。

## 後処理用アーカイブ

OUTPUTS.md に従います。`crystal.respond(..., output_directory=...)`、
`TDDFTState.finish(directory)`、`outputs.OutputArchive` と `outputs.save_arrays` が
実在する入口です。TDDFTのfinishは状態を閉じ、初期ファイルハッシュと選択した
本体出力を残します。大きい分布は初期入力のanalysisで選択・間引きます。
波動関数チェックポイントとprojection/unfoldingは未対応です。NPZアーカイブを
restartファイルとして渡さないでください。別ノード通信やコア固定は未実装・未検証です。

## バンド図

BANDS.mdを参照してください。`bands.BandStructure.from_salmon` と `.plot` が入口です。
固有値の単位、実際のCartesian k、セルの選択、余剰点、収束対象バンド、エネルギー基準を
明示してください。経路上の価電子帯上端を全域のVBMと呼ばないでください。

バンドの計算入口は `BandInitialState.from_salmon` → `run(BandPath, ...)`。
初期束の外の新規保存先、対応した監査済み実行ファイルを使い、未収束を成功扱いしません。

ASE連携は任意の独立モジュール `salmonpy.ase`。ASE.mdを参照。
`run_bands` は構造を照合し、ASEのk点列をそのまま計算します。Calculator・力は独立したASE_CALCULATOR.mdの対応範囲を参照してください。

GS作成は `salmonpy.dft.DFT`、ASE変換は `salmonpy.ase.dft_from_atoms`。DFT.mdに従い、電子数・擬ポテンシャル・格子を明示します。
収束済み `GroundStateResult` からバンド／TDDFT初期状態を作れます。PZ/KYのCalculatorと本体の力はASE_CALCULATOR.mdを参照してください。

`salmonpy.ase_calculator.SalmonCalculator` は任意依存の独立モジュール。energy/forcesのみ、PZ/KY・単一ランクが対象です。
粗い格子のエネルギー有限差分はegg-box効果を含むため正解と仮定せず、剛体並進と格子収束を調べます。

## 初心者への最初の案内

manual.jsonの `beginner_intake` を必ず適用します。新しい計算を選ぶ前に、未回答の
「使ってよい総メモリ」「使ってよいCPU並列数」「やりたいこと」を短くまとめて聞きます。
既に回答された条件は再質問せず、「不明」も受け付けます。装置の全メモリ・全コア数を
自動的に利用許可と解釈しません。目的に応じて構造規模・精度・所要時間の希望を補います。

最初に、今できること／追加準備が必要なこと／未対応と、予算内で始められる小規模例を示します。
計算規模の上限は実測・暫定見積もり・不明を区別します。未測定の原子数上限などを断言しません。
必要なら予算に収まる小規模試行でピークメモリ・時間・保存量を測定し、条件を記録して更新します。
見積もりはメモリ強制制限ではありません。自動容量推定器・ハードメモリ制限は現状未実装です。

プロセスとスレッドを合計したCPU使用量、全子プロセスを含むメモリ、保存容量を扱います。
メモリの余裕は許可された上限の内側に確保します。現在のPython起動が単一MPIランクである
機能を、任意のMPI並列数で実行可能と案内しません。収まる計算規模と物理的な収束精度は別です。

実行可能な用途別Pythonサンプルと必要な入力は [examples/README.md](examples/README.md) を参照してください。AI検索の各taskの `sample_files` からも辿れます。

## 初見の環境診断と検索

`python -m salmonpy.doctor --json` でライブラリの読み込み可否、追加依存の所在、
任意指定した外部実行ファイルのハッシュ照合を確認できます。計算・コンパイル・導入は行いません。
外部SALMONの `hash_matches` を機能対応や実行成功と読み替えないでください。
`manual.search` は常に `getting_started` と `task_catalog` を返します。一致しない検索語では
目的別一覧を示して目的を確認し、対応不能とは即断しません。
