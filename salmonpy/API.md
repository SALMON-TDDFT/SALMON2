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
