# nativeを含む混在ステップの保護

[English](NATIVE_STEPPING.en.md) · Topic: `native_protected_step` · Paired revision: `2026-11-17.2`

`salmonpy.native_steps.native_step_targets` は対応したnative部品を
[共通のステップ保護](STEPPING.md)へ接続します。明示的に新しくビルドした対応ホストMPI
ライブラリが必要です。個別の数値復帰・checkpoint契約から全体保護を推定しません。
旧ライブラリは時間発展前に拒否し、自動fallbackやimport時のビルドは行いません。

```python
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepRecoveryError

targets = native_step_targets(coupled)
try:
    frame = coupled.step_protected(targets)
except StepRecoveryError as error:
    if not error.original_preserved or error.cleanup_errors:
        raise  # この生存期間では継続しません。
    # 受理済み状態を維持しています。原因を確認・修正してから同じステップを再試行します。
    # 拒否された試行でループの時計を進めないでください。
```

全ランクで参加します。対応は既知の組み込み型に限定します。NativeEHField（PMLと場内の
全LDを含む）、1ランクのMultiscaleField、周期・非スピン分極・固定イオン・中点伝播の
標準native TDDFT、高速PR版native HSE、通常サンプルの組み込みSBE端点ラッパーです。
各native群が数値退避を保持し、FFTプラン・コミュニケータ・所有資源のハンドルは複製しません。
波動関数、密度／ポテンシャル、XC状態／係数、ACE／履歴／キャッシュ、場／PML／LD／電流と
受理済み時計を直接戻し、過去の伝播を再実行しません。PythonのA／E／電流／時計／フラグの
参照も戻します。nativeの任意の複数コンテキストや安定した公開C ABIを追加するものではありません。

自作バックエンドには `custom_targets=(StepTarget(...), ...)` を明示し、ownerを実際の
local状態に一致させます。未知のラッパー・派生クラスをメソッド名だけで認証しません。
コールバックは状態と副作用全体を扱い、全群ランクで通信を終えて戻る必要があります。
checkpointの反映で状態オブジェクトが置き換わった場合、対象を作り直します。古いownerは
退避前に拒否します。native checkpoint候補とステップ退避の同時使用はPython／nativeの
両方で拒否します。

## 出力と失敗時の扱い

対象出力は、各owner専用作業フォルダ直下の通常RT電流／エネルギー時系列と、E/Hの通常の
点観測です。退避時に、受理済みの未flush行をflushしてから基準位置を記録します。
巻戻しではnativeのRT出力を閉じ、失敗試行の追記を基準位置まで戻し、同じユニットを
受理済み末尾へ開き直します。正常終了時は受理した行をflushします。通常の `step()` は
従来の経路と演算順を維持します。

監査するのは専用の通常ファイル・単一リンクの `.data` だけです。別の出力先／サブフォルダ、
密度・微視的電流・ELF・静電場・projection・TTM・ASE／ART・周波数観測・restart出力などの
未監査設定は、該当する経路で更新前に拒否します。設定を黙って無効化しません。
ファイルの置換・削除や想定外の新規時系列は、勝手に削除せず失敗扱いにします。
実行中のフォルダ移動や同時の外部読み書きはしないでください。巻戻しまでの間、仮の追記が
外から見える場合があり、永続化の原子的公開やMPI／プロセス停止からの復帰ではありません。
ディスク障害、メモリ枯渇、初期化／外部ライブラリのSTOP、nativeから戻らない異常は対象外です。
標準native TDDFTの未知エラーも、E/H・HSEと同様にfailedを保持し、後の正常らしい戻り値で
生存期間を復活させません。

StepRecoveryErrorのoriginal_preservedがTrueかつcleanup_errorsが空の場合だけ、原因修正後に
再試行します。自動再試行・波形の減衰・刻み変更はしません。全体退避と部品内部の復帰用退避が
同時に存在するため、メモリ助言では未較正の追加費用と実測ピークを確認します。予算は警告のみです。

## 通常サンプルで使う

`salmonpy.examples.mpi_mixed_native`、または推奨の1次元経路
`salmonpy.examples.mpi_multiscale_tddft` に `--protect-steps` を追加します。
通常CLIは未処理の例外では停止します。独自ホストは安全に拒否されたステップを捕捉して扱えます。

同梱の `examples/mixed_tddft_hse_layers.json` はSBE・TDDFT・HSEの3電子点と、場内2種類のLDを
配置します。場所は `importlib.resources.files('salmonpy').joinpath('examples/mixed_tddft_hse_layers.json')`
で取得できます。HSEの初期状態が位置引数BUNDLEと異なる場合、`--hse-bundle` で別に指定します。
そのファイル群はHSE固有の識別情報とcheckpoint互換性へ入り、HSE側だけを変更してもTDDFTの
物質ID情報は変わりません。DFT実行・ダウンロード・形式変換は自動では行いません。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 7 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD_JSON SERIAL_RUNTIME PBE_BUNDLE NEW_OUTPUT \
  --layers LAYOUT_JSON --field-grid 1 2 1 --tddft-ranks 2 \
  --fast-hse-build SERIAL_HSE_BUILD_JSON --fast-hse-mpi-build HOST_HSE_BUILD_JSON \
  --hse-bundle HSE_BUNDLE --hse-ranks 2 --protect-steps \
  --steps 32 --rank-budget 7 --memory-budget-gib 48
```

この40³の3次元例は接続回帰用で、収束した実材料／多層膜計算ではありません。実用TDDFTの
マルチスケールは引き続き1次元から始め、大きな3次元の電子応答にはSBE／自作コードを推奨します。
並列数と予算は実際の割当に合わせ、大きなMacの値をそのまま流用しません。
ローカル検証は標準PZ・TB-mBJ自動／固定係数・Libxc-PZ/PBEと、高速PR版HSE Taylor4／学習PT-CNの
内部MPI2、ならびに1次元の独立2 TDDFT群です。測定範囲は検証記録を参照し、今回の小規模系だけで
新たな他OS対応・物理収束が保証されたとは判断しないでください。


## 大きめの格子で測った保護の負担

64GB・18コアのApple Siliconで、96³のMaxwell＋局所LD＋SBE＋Libxc-PBE2＋高速PR版HSE
PT-CN2を64ステップ進め、保護なし／ありを各2回比較しました。場2／全7ランクの平均は
12.00／14.38秒、場4／全9ランクでは7.08／7.66秒でした。各ランクのピークRSSの合計の最大値は
それぞれ2.012／2.692GiB、2.949／4.417GiBです。全配列と通常出力が完全一致しました。
これはローカルの測定例で、メモリ上限や推奨並列数ではありません。拡大したのはMaxwell格子で、
電子状態の精度は変えていません。未知のメモリ量を確定値に置き換えず、実際の計算環境で比較してください。
