# 計算の一ステップが失敗したときに元の状態を保つ

[English](NATIVE_STEPPING.en.md)

FDTD・LD・電子状態を一緒に更新するとき、一部だけが先へ進むと再試行できません。この機能は更新前の状態を退避し、失敗時に電磁場・物質状態・時刻をまとめて戻します。

`salmonpy.native_steps.native_step_targets` で対応する部品を `StepTarget` にまとめ、`step_protected` で一ステップ進めます。呼び出し元と同じMPIを使ってビルドした共有ライブラリが必要です。対応しないライブラリは時間発展を始める前にエラーを返します。

```python
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepRecoveryError

targets = native_step_targets(coupled)
try:
    frame = coupled.step_protected(targets)
except StepRecoveryError as error:
    if not error.original_preserved or error.cleanup_errors:
        raise  # この計算状態は継続して使いません。
    # 受理済み状態を維持しています。原因を確認・修正してから同じステップを再試行します。
    # 拒否された試行でループの時計を進めないでください。
```

全ランクで参加します。対応は既知の組み込み型に限定します。NativeEHField（PMLと場内の
全LDを含む）、1ランクのMultiscaleField、周期・非スピン分極・固定イオン・中点伝播の
標準Pythonから直接呼ぶTDDFT、高速PR版Fortran製HSE、通常サンプルの組み込みSBE端点ラッパーです。
各Fortran/Cで計算するMPI群が数値退避を保持し、FFTプラン・コミュニケータ・所有資源のハンドルは複製しません。
波動関数、密度／ポテンシャル、XC状態／係数、ACE／履歴／キャッシュ、場／PML／LD／電流と
最後に成功した計算時刻を直接戻し、過去の伝播を再実行しません。PythonのA／E／電流／時計／フラグの
参照も戻します。同じプロセス内の任意の複数計算状態や安定した公開C ABIを追加するものではありません。

自作の計算コードには `custom_targets=(StepTarget(...), ...)` を明示し、ownerを実際の
local状態に一致させます。未知のラッパー・派生クラスをメソッド名だけで認証しません。
コールバックは状態と副作用全体を扱い、全群ランクで通信を終えて戻る必要があります。
checkpointの反映で状態オブジェクトが置き換わった場合、対象を作り直します。古いownerは
退避前に拒否します。計算処理 checkpoint候補とステップ退避の同時使用はPython／計算処理の
両方で拒否します。

## 出力と失敗時の扱い

対象出力は、各owner専用作業フォルダ直下の通常RT電流／エネルギー時系列と、E/Hの通常の
点観測です。退避時に、受理済みの未flush行をflushしてから基準位置を記録します。
巻戻しでは計算処理のRT出力を閉じ、失敗試行の追記を基準位置まで戻し、同じユニットを
受理済み末尾へ開き直します。正常終了時は受理した行をflushします。通常の `step()` は
従来の経路と演算順を維持します。

監査するのは専用の通常ファイル・単一リンクの `.data` だけです。別の出力先／サブフォルダ、
密度・微視的電流・ELF・静電場・projection・ASE／ART・周波数観測・restart出力などの
未監査設定は、該当する経路で更新前に拒否します。設定を黙って無効化しません。
ファイルの置換・削除や想定外の新規時系列は、勝手に削除せず失敗扱いにします。
実行中のフォルダ移動や同時の外部読み書きはしないでください。巻戻しまでの間、仮の追記が
外から見える場合があり、永続化の原子的公開やMPI／プロセス停止からの復帰ではありません。
ディスク障害、メモリ枯渇、初期化／外部ライブラリのSTOP、計算処理から戻らない異常は対象外です。
標準Pythonから直接呼ぶTDDFTの未知エラーも、E/H・HSEと同様にfailedを保持し、後の正常らしい戻り値で
初期化から終了までの期間を復活させません。

StepRecoveryErrorのoriginal_preservedがTrueかつcleanup_errorsが空の場合だけ、原因修正後に
再試行します。自動再試行・波形の減衰・刻み変更はしません。全体退避と部品内部の復帰用退避が
同時に存在するため、メモリ助言では未較正の追加費用と実測ピークを確認します。予算は警告のみです。

## 通常サンプルで使う

`salmonpy.examples.mpi_mixed_native`、または推奨の1次元経路
`salmonpy.examples.mpi_multiscale_tddft` に `--protect-steps` を追加します。
通常CLIは未処理の例外では停止します。独自呼び出し元プログラムは安全に拒否されたステップを捕捉して扱えます。

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

## 同じ物質を複数の点へ配置する

同梱の `examples/mixed_multipoint_layers.json` は、SBE4点・Libxc-PBE TDDFT2点・
高速PR版HSE2点と、Drude／Lorentz各1,152セルを配置する接続例です。電子物質の各セルに
独立した状態とMPI群を作り、同じ物質IDを繰り返しても時間発展中の状態は共有しません。
各点の論理IDと担当情報は保存・再開にも引き継がれます。Maxwell格子は既定の40³ではなく64³です。

上の7ランク例の `LAYOUT_JSON` をこのファイルに変更し、`mpiexec -n 16` と
`--field-grid 2 2 1 --field-shape 64 64 64 --steps 64 --rank-budget 16`
を指定します。場4ランク、SBE4ランク、TDDFT2ランク×2群、HSE2ランク×2群です。
LD専用群は不要です。格子・形状・並列数は対応しているため、層を厚くするだけでなく
必要並列数も再確認します。予算は実際の割当に合わせてください。

直接再開は、最初に `--direct-save-at 31` で保存し、同じ入力・バイナリ・MPI分割の
別実行から `--direct-resume PREVIOUS_OUTPUT/direct-checkpoint` を指定します。
出力先は新しくし、両実行とも `--steps 64` を維持します。残り33ステップだけを進め、
受理済み伝播は再実行しません。新しい実行の通常RTファイルはその区間分です。応答配列の
全履歴は混在checkpointに保持されますが、過去のRTファイルが結合されるとは解釈しません。

ローカル16ランクで連続実行・保護付き実行・直接再開と、1物質群の失敗／再試行を検証しました。
これは接続・資源検証であり、収束した3次元TDDFTの物理結果ではありません。実用TDDFTは
引き続き1次元から始め、利用者の資源に合わせます。自作物質は明示的なbindingと仕様で接続します。
