# 公開Python APIと互換性

[English](API.en.md) · Topic: `public_api` · Paired revision: `2026-11-14.3`

版は **0.1.0a5** です。`salmonpy.__version__` と配布メタデータに同じ版を示します。
プレビュー版であり、安定した1.0 APIを宣言するものではありません。物質C ABIのプレビュー1は[C_ABI.md](C_ABI.md)を参照してください。

## 公開範囲

`salmonpy.api.PUBLIC_API` がモジュール別の公開名一覧です。指定モジュールからimportします。
偶然importされている名前や、先頭がアンダースコアの実装名は公開APIではありません。
戻り値の状態・結果オブジェクトは各ガイドに記載したメソッド・プロパティを公開します。
内部の場状態クラスを直接構築しません。`StateSnapshot` は型を示す名前で、利用者が構築するものではありません。

| モジュール | 主な用途 |
|---|---|
| `ld`, `sbe`, `crystal`, `materials` | 物質・初期状態・A駆動応答 |
| `maxwell`, `timing` | 場の生成・結合・時刻宣言 |
| `fields`, `geometry`, `analysis`, `outputs` | パルス・構造・解析・出力 |
| `dft`, `bands`, `tddft`, `fdtd`, `runtime` | 外部ソルバーの明示的実行・準備 |
| `ase`, `ase_calculator`, `salmon_vr` | 任意の連携機能 |
| `distributed`, `material_response`, `embedded_tddft`, `coupled` | 既存の寿命制限がある実験的MPI・埋め込みAPI |
| `errors`, `snapshots`, `api` | 例外・不変の状態トークン・公開一覧 |

トップレベルの `__all__` は版、`LorentzDrude`、例外型に限定します。
任意依存を一括importしません。import時のコンパイル・取得・ネイティブ読込み・MPI初期化はありません。
snapshot／restoreは標準出力へ書きません。必須Python依存は引き続きNumPyです。
LD/SBEカーネルや外部ランタイムの明示的準備、任意依存の要件は従来どおりです。

## 互換方針

配布するAPI・挙動の変更には新しいパッケージ版を付けます。同じリリース内で単位・電流符号・
時刻配置・配列名を黙って変更しません。プレビュー版間では明記した破壊的変更を行うことがあります。
ホスト側は検証した版を固定し、移行説明を確認してください。将来の安定版では廃止予告期間を
設けますが、現在の実験的MPI API全体に安定性を保証するものではありません。
公開一覧への掲載と、機能の成熟度は別です。

アーカイブ・時刻・ランタイム・snapshotのスキーマ版は独立です。パッケージ版だけから
保存形式の互換性を推測しません。メモリ上のsnapshotは同じパッケージ版・同じ元モデルを要求します。

a5への移行で更新式・電流規約は変更しません。結晶は **A/c → 電荷電流密度j**、
内部の電子物質電流は J = -j です。過去のa4開発中に規約が変更されたため、古いデータは
記録されたメタデータで解釈し、一律に符号反転しません。基本はtransverseで、Longitudinalは明示的選択です。

## 例外の契約

| 型（`salmonpy.errors`） | 意味 | 従来の捕捉型 |
|---|---|---|
| `SalmonpyError` | ライブラリが明示的に送出する例外の親 | `Exception` |
| `InputError` | 値・互換性・状態遷移の不正 | `ValueError` |
| `InputTypeError` | 型の不正。`InputError`の子でもある | `TypeError` |
| `SnapshotError` | 不正・不適合な復元。状態は保持 | `ValueError` |
| `NumericalError` | 検出した非有限値などの数値異常 | `FloatingPointError` |
| `RuntimeFailure` | ランタイムの欠如・失敗・結果の未確認 | `RuntimeError` |

既存の `RuntimeUnavailable` と `OutputSaveError` は名前を維持し、`RuntimeFailure`を継承します。
従来の組込み例外型での捕捉も利用できます。OS・ファイルの例外（`FileExistsError`を含む）、
プロセスの時間切れ、NumPy変換、任意依存の不足、利用者コールバックの例外は元の型を保つ場合があります。
外部例外を一律に包みません。収束報告がないことはランタイムの失敗であり、数値発散の証明ではありません。

例外の型だけから計算全体の巻き戻しを推測しません。内蔵LD・SBE・Pythonの場のstep／restoreは
失敗時に受理済み状態を保持します。MPI応答・結合の失敗はセッションを無効化し得ます。
FortranのSTOP／abortを回復可能なPython例外にする変更ではありません。

プロセス分離`TDDFTState`は`save_checkpoint(directory)`と
`restore_checkpoint(directory, max_bytes=None)`も公開します。
同じ内容の初期モデルを明示的に用意し、Zarr履歴から別プロセスで再開できます。
初期GS本体は保存しません。詳細は[状態保存](SNAPSHOTS.md)を参照してください。

直列モデルの共通契約、適合チェッカー、MPI点分散への接続は[MATERIALS.md](MATERIALS.md)を参照してください。

`fast_hse` の `FastHSEMaterialProvider`／`FastHSEResponseBackend` は明示ビルドした高速PR版serialワーカーを外側MPIの物質応答へ接続するプレビュー入口です。COUPLED.mdを参照してください。HSE内部MPIを初期化しません。

`fast_hse_mpi.NativeFastHSEBackend` は高速PR版HSEの実験的なホストMPI入口です。k点のみ1／2ランクの物質群と、同じ分割での直接保存／復元に対応します。メインスレッド・1プロセス1回のnative寿命に限定します。COUPLED.mdを参照してください。

`checkpointing.CheckpointCapability` と `check_checkpoint_compatibility` は、自作バックエンドにも使える保存／復元前の共通宣言・検査です。[CHECKPOINTING.md](CHECKPOINTING.md)を参照してください。

`checkpointing.CheckpointTarget` と `execute_checkpoint` は検査後にバックエンド独自の保存／復元関数を呼び、不透明な戻り値をプロセス内で返します。

`CheckpointState`／`restore_checkpoint` は、独立候補を全ランクで準備した後にまとめて公開する保護付き復元です。`AtomicRestoreError`で元のホルダー維持と候補の破棄エラーを区別します。

メモリ助言の公開入口はsalmonpy.memoryのMemoryEstimate／MemoryReport／memory_advice／process_memoryです。予算超過は警告のみで、実行可否や上限の保証ではありません。詳細は[MEMORY.md](MEMORY.md)。

`stepping.StepTarget`・`guarded_step`・`StepRecoveryError` は明示した完全退避契約によるステップ保護のプレビュー入口です。STEPPING日英ガイドを参照してください。

`native_steps.native_step_targets` は、監査済みnative状態と通常出力を共通の保護付きステップへつなぐプレビューAPIです。NATIVE_STEPPING日英ガイドを参照してください。

### ソルバー初期化時の原子入力

標準ソルバーと高速PR版HSEは、原子数・元素数・座標（Cartesian／分数）・元素IDを、内部の
値オブジェクト `s_ion_initial_config` に取り込みます。`initialize_ions` は選択した状態の
格子から座標を解決し、独立した所有配列を作ります。形状不一致・非有限座標・範囲外元素ID・
下限が1でない入力配列は、対象状態を変えずに拒否します。明示的な処理は共有入力を読み書きしません。
既存互換経路では下流処理のための分数座標→Cartesian共有入力正規化を維持し、正常な計算結果を保ちます。

これは内部の所有境界であり、新しいPython APIや安定したC ABIではありません。セル・k点・
対称性・コミュニケータなどの依存は残り、異なる全入力のGS／TDDFTを同一プロセスへ自由に
配置できることは意味しません。メモリ確保失敗や初期化中のSTOPからの復帰も未対応です。

### セル・実空間格子の明示的な解決

標準ソルバーにも高速PR版HSEと同じ値型・処理 `s_cell_config`、`s_cell_geometry`、
`s_rgrid_config`、`resolve_cell`、`init_grid_whole_config` を揃えました。両方とも
入力の取り込みだけを行う `cell_config_from_input` を使い、明示的な解決処理は共有入力を
読みません。正常入力のセル指定の優先順位・演算順・周期／孤立系の添字・epsilonを含む丸め・
Maxwell用整合長を維持します。拒否時は対象の形状・格子添字配列を変えません。
これはセル値と全体格子添字の解決で、全座標・並列分割・対称性・FFT資源の初期化ではありません。

既存互換経路には共有入力の取り込みと、下流処理用の `al` 正規化が残ります。異なるセルの
GS／TDDFTを同一プロセスで自由に初期化できる段階ではなく、新しいPython／C公開入口もありません。
実験用RT設定版にも同じ値型を持たせ、準備時に更新した置換前後のハッシュを検査します。

### 初期軌道構成の明示的な選択

標準版・高速PR版の初期化で、計算モード・周期性・スピン／SOC・状態数／電子数・温度・
k点判定を内部の `s_orbital_initial_config` に取り込みます。`resolve_orbital_layout` は
共有入力を読まずにスピン数・軌道数・実数／複素数格納を選びます。セルの直交性を明示的に渡し、
高速PR版は所有済みの占有設定・実際のk格子と、従来と同じ取り込み時のhybrid判定を使います。
標準版のファイル指定k点も、この判定には従来の `-1` 格子印を維持します。

RTの `temperature < -1e-12` 分岐など、正常入力の既存GS／RT選択規則は維持します。
不明なモード・非有限温度・使えない構成は、結果オブジェクトを変えずに拒否します。
これは初期配列構成の選択であり、占有容量・物理的妥当性・SCF入力全体の所有や、native初期化からの
異常復帰を保証するものではありません。新しいPython／C公開入口もなく、k点・対称性・
コミュニケータなどの共有依存は引き続き残ります。

HSEは高速PR版のみです。`NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)` で
対称性を明示的に許可できます。内部kランク1・周期・非スピン分極・CPU・full-q・factor-off
Taylorに限定し、共有の対称化設定は不変、1プロセス1寿命です。既定値Falseと既存C入口は
対称性を拒否します。学習／PT-CNの対称性と内部k MPI2は未対応です。

新C入口 `sp_fast_hse_init_symmetry` は既存initと同じ引数で、
`sp_fast_hse_symmetry_protocol()==1` を宣言します。全ランクの有限なstatus34とこの宣言を
確認した場合だけ、Pythonは状態を変えない `InputError` を返し、Aを直した再試行を許します。
未知・旧ABIの失敗は継続しません。FDTD結合では全体の再試行に `step_protected` が必要です。
直接／候補復元は共通checkpoint経路を使います。

[対称性ガイド](FAST_HSE_SYMMETRY.md)には、混在サンプル、対応範囲と既知の極小勾配XC感度を
記載しています。接続・保存復元のテストを物理収束の証明と解釈しないでください。
