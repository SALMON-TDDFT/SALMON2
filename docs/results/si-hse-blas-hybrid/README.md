# Taylor4＋ACE：BLAS内部OpenMP

時間発展は `hse_taylor4`（ACE使用）。従来のOpenBLAS内部1スレッド固定を廃止し、交換作用・ACE構築・ACE適用のBLAS呼び出しをアプリケーションOpenMP領域の外に置く。位相乗算、転置用データ整形、カーネル乗算はSALMON側でOpenMP並列化する。MPIとFFTW呼び出しは並列領域外。軌道配列をスレッド数分複製しない。BLAS自身の作業領域はライブラリ依存。

富岳を含むベンダーBLASに固有APIを要求しない。OpenMP並列版BLASをリンクし、実行環境側で内部スレッド数を設定する。`OMP_NUM_THREADS=12` は逐次版BLASを並列化する指定ではなく、行列サイズによってBLASが使用するスレッド数を減らすこともある。富岳でのビルド・性能は未検証。

ローカル測定はMac、GNU Fortran/OpenMPとOpenBLAS 0.3.33 OpenMP版。`OMP_NUM_THREADS` と `OPENBLAS_NUM_THREADS` を同じ値に設定、`OMP_DYNAMIC=FALSE`、MPIのbindingなし。停止済みPython参照計算は再開していない。PythonはFortran実行の時間・RSS採取と結果比較にのみ使用。

Si、12³実空間格子、4³ k点、dt=0.16 a.u.、5ステップ、各条件2回。測定は同一MPI数でスレッド数を変えた比較。RSSは0.1秒間隔の観測値であり、厳密なピークや共有ページを除いた専有メモリではない。短時間計測から長時間の励起子スペクトル精度は評価しない。

数値検証：1728格子×16軌道×8 k点のACE作用を独立MATMULと比較し、OMP 1/2/4/12で相対誤差1e-11以下。ネイティブ交換の4テスト合格。実時間計算はOMP1の最終軌道と相対差1e-10未満を要求する。

生データ：`runs.json`。再現用計測器：`samples/hse_mlwf_reference/benchmark_blas_hybrid.py`。

入力の `directory_read_data` は既存のHSE基底状態 `calculations/si_hse_native/reference_restart` を使用する。保存した入力は測定時の絶対擬ポテンシャルパスを含むため、別環境ではパスを変更する。BLAS設定照会の `blas-thread-config.txt` は「GNUスレッド番号、OpenBLAS設定スレッド数、GNUチーム数」。設定12を確認したもので、各ZGEMM内部の実際の稼働人数の計測ではない。

## 測定結果（MPI 1固定）

| OMP | 平均 秒/step | OMP1比 | RSS MiB |
|---|---:|---:|---:|
|1|5.046|1.00×|597.0|
|2|3.648|1.38×|597.6|
|4|3.193|1.58×|597.9|
|12|4.687|1.08×|599.1|

この小さいSi計算ではOMP12よりOMP4が速い。全体時間の変化にはSALMON側OpenMPも含まれ、BLAS単独の加速率ではない。FFTWは直列のまま。富岳実機での最適MPI×OMPは未測定。

追加検証：実MPIの独立交換作用テスト合格。MPI×OMP=(1,1),(3,1),(8,1),(3,2),(3,4),(1,12)、不均等分割・空行担当・任意の試行軌道・rank局所NaN/shape異常を検証。最終Taylor軌道の最大相対差は4.18e-12。
