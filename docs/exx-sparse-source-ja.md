# 交換用MLWFの支持点保存

2026-09-29。native Γ点の有限支持交換で、球状マスク後のMLWFを格子点番号・複素値・列オフセットで保持する。実時間発展の波動関数と、次のゲージ輸送に使う`previous`は従来の全メッシュ表現を維持する。

## 表現と処理

`exx_sparse_orbitals`の`s_sparse_orbitals`は、担当空間の格子点番号、非ゼロ値、各軌道の開始位置を持つ。既存マスクの非ゼロ成分は微小値もすべて保存し、新たな閾値や半径近似は加えない。格子点番号は列内で昇順。空軌道・空の支持領域を許し、不正な添字や非有限値は検査する。

- 既存の99.9%規則／明示半径でマスクを決めた後に圧縮し、密な交換用sourceを解放する。明示半径の優先順位・warningは維持する。
- 軌道MPI間では支持点の添字と値だけを配信する。FFTに渡す際は1軌道分だけ担当メッシュへ展開する。
- source ACEの構築対象は32軌道ずつ展開し、交換作用を求める。各ブロックの仕事量・対数を合算する。
- ACE計量行列の左側MLWFは支持点だけで内積を計算する。ScaLAPACK分散とLAPACKの両経路で、全訓練MLWFへの復元を避ける。
- source ACEが拒否された場合は圧縮sourceを維持して従来のoccupied-vector ACEへ戻る。受理後は不要なsourceを解放する。
- full-supportの場合は密経路を維持する。途中で支持条件が変われば更新時に古い圧縮データを破棄し、現在のsourceから作り直す。

全MPIランクへの軌道・メッシュ全体の集約は導入しない。局所化と球の判定の間には、密sourceを一時生成する処理が残る。最終交換作用`w`、輸送用previous、実時間波動関数なども残り、全計算のメモリが線形化したわけではない。

## 検証

GNU MPI/ScaLAPACK版と逐次版のビルド成功。HSE06/PBE0/PBEhのsource ACE 7試験、対省略SCF/RT 6試験に合格。疎な訓練軌道を使う計量行列も、空ランク・ゼロ作用・悪条件行列を含め、MPI 1/2/4で密訓練と比較。ScaLAPACKなしのLAPACK経路でも比較した。読み取りによる独立レビューで重大な指摘なし。富岳・NVHPC/GPUは今回の変更を未検証。

```sh
python3 experiments/unit_exx_sparse_orbitals/run.py
python3 developer_tests/651_hybrid_exchange/metric/run.py --build /path/to/build
python3 developer_tests/651_hybrid_exchange/metric/run.py --build /path/to/build --no-scalapack
```

## 実測

同一の保存GS・入力を使う128 H₂、MPI 8・OMP 1・BLAS 1、PBEh(40)、0.999支持、source ACE、impulse 16ステップ。前回の採用済み構成の保存結果を再利用し、今回の計算を1回行った。各1回の時間差は実行環境のばらつきを含む。無外場計算や電流差引きは行わない。


| 指標 | 採用済み構成 | 圧縮source |
|---|---:|---:|
| 最大ランクRSS | 2,917.06 MiB | 3,012.72 MiB |
| RT時間/step | 12.050秒 | 11.659秒 |
| RT 16ステップ | 192.80秒 | 186.54秒 |

全体ピークRSSは3.28%増、時間は3.25%減だった。時間差は各1回の観測であり、有意な高速化とは断定しない。電流の最大差は1.02e-20、エネルギーの最大差は9.95e-13 Ha。測定値と入力ハッシュは[集計JSON](reports/exx-sparse-source/summary.json)に保存した。

ランク0の交換用sourceは8,388,608複素成分（128 MiB）から64,192支持点（値・整数添字・オフセットで約1.225 MiB）になった。この値はランク0の単一配列の比較であり、全ランク最大RSSではない。圧縮保存自体は機能したが、プロセス全体のピーク削減にはつながらなかった。密sourceの一時生成などは残るが、ピーク増加の原因はこの計測だけでは特定していない。

## 採用状態と再現

全体のメモリ削減が確認できなかったため、native経路の`use_sparse_source=.false.`を既定とし、今回の経路は実験用として保持する。新しいnamelist入力は追加していない。既定バイナリは従来の密sourceを使う。

以下は共有オブジェクトを再利用して、密source、圧縮source、1軌道ブロックの境界試験用バイナリを別ディレクトリに作る。既定ビルドは変更しない。出力先には未使用のディレクトリを指定する。

```sh
python3 benchmarks/benchmark_matrix_memory/build_sparse_variants.py \
  --build /path/to/build --output /path/to/new-comparison-directory
```

`before/salmon`が既定、`after/salmon`が圧縮source、`boundary/salmon`が圧縮sourceを1軌道ずつ処理する版。各ディレクトリのmanifestにコマンド・ソースと実行ファイルのSHA256・共有オブジェクトのSHA256を記録する。通常の32軌道版と1軌道版のsource ACE回帰をそれぞれ実行し、ブロック境界も検証する。
