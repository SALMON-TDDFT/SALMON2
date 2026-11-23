# Si入門例：基底状態から電流へ

[English](SI_INTRO.en.md)

固定Si入力でSCFを計算し、A/cで駆動した21時刻の電荷電流を期待値と比べます。
2026-11-23に、現行ソースから作った電子数付き応答の計算プログラムで再検証しました。
入力・擬ポテンシャル・期待値・許容誤差は以前と同じです。

## 必要な配布物

- OS・CPU・Python版に合う、更新済みSALMONpy wheel。Python依存はNumPyです。
- `provenance.json` の現在の `runtime_artifacts` に登録された外側ZIP。

現時点で新しく検証したのはApple Silicon Macで作成したローカル配布物です。公開Releaseにはまだ置いていません。配布者からwheelと対応ZIPを受け取ってください。保存先URLを推測しないでください。新ZIPがない環境では、まずwheelだけの [LD/SBEテスト](README.md)を使います。

旧 `runtime-2026.11.14` Releaseは現在の電子数付き応答と合いません。過去のハッシュと出典は `historical_runtime_artifacts` 等に保持し、現在の入口では受け付けません。Intel Mac・Linux用の新しい配布物については再検証が必要です。WindowsのDFT/TDDFT配布を保証しません。

## 導入・準備・実行

```sh
python -m pip install --only-binary=:all: ./SALMONPY_WHEEL.whl
python -m salmonpy.quickstart
python -m salmonpy.examples.si_intro prepare new-si --runtime-artifact RUNTIME_ARTIFACT.zip
python -m salmonpy.examples.si_intro run new-si --threads 1 --memory-budget-gib 4
```

仮名のwheelとZIPを受け取った実ファイル名に置き換えます。外側ZIPは自動展開せず、そのまま渡します。`prepare` はZIP・実行ファイル・依存ライブラリ・入力を検査して配置するだけで、計算を始めません。管理ノードでここまで実行し、`run` は割当済み計算資源で実行します。暗黙の通信・ビルドはありません。既存フォルダを上書きしません。

`--threads` は1〜4。メモリ4 GiBはジョブ全体の申告予算で、OSによるハード上限ではありません。計測不能を計算失敗とは扱わず、未計測値で予算内とも判定しません。

## 確認する結果

`new-si/result/report.json` の `calculation` と `validation` が共に `passed` なら成功です。資源計測は `resources` で別に報告します。

- GSと再開データ：`result/gs/`、ログ：`result/gs/native/stdout.log`
- 電流：`result/current.csv`、`result/current.npz`、応答アーカイブ：`result/rt/`
- CSV列：`time_au,Ac_x_au,Ac_y_au,Ac_z_au,j_x_au,j_y_au,j_z_au`

原子単位の**電荷電流密度j**です。ソルバーの電子物質電流Jとは `J = -j` の関係で、返値を再度反転しません。

期待GSエネルギーは約−31.7465390575 Hartree。エネルギー絶対誤差1e−7 Hartree、電流各成分は `atol=1e-9, rtol=1e-5`、時刻は丸め誤差だけを許します。反復数は診断値で完全一致を要求しません。2026-11-23のソース実行では54反復、全21時刻の電流・エネルギーの期待値との差はゼロでした。

## この例が検証する範囲

立方Si 8原子、セル長10.26 bohr、32電子・32状態、16³格子、2³ k点、PZ、SCF基準1e−10です。dt=0.02 a.u.で20ステップ（0.4 a.u.）だけ進めます。これは導入・接続の確認で、実際のレーザーパルスや収束した光学物性ではありません。

次の実計算では波形を数周期、観測を目的に必要な時間へ伸ばし、[収束診断](CONVERGENCE.md)で格子・k点・dtを調べます。格子はdtを収束させた電流の相対RMS差5%をスクリーニング基準とし、dtは電子数減少率1e−5を目安にします。誘電率や高次高調波の精度保証とは区別します。

擬ポテンシャルは既にwheel内に同梱されています。個別取得や他の参照データの準備は[擬ポテンシャル取得](PSEUDOPOTENTIALS.md)を参照してください。
