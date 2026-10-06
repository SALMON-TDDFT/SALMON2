# 性能測定用スクリプト

旧testsuites内のbenchmark_*.pyをここへ分離した。標準CTestの合否判定には使用しない。実行はリポジトリルートから `python3 benchmarks/unit_lcfo_rt/benchmark_gram.py` 等。元のtestsuitesにある共通probe/stubを参照するため、それらは移動していない。

MPI・BLAS/ScaLAPACK・コンパイラ設定は各スクリプトを確認する。既存Mac固有のライブラリパスは今回変更していない。重いジョブと重ねず一件ずつ測定し、速度/RSSは物理精度やCI合格を意味しない。
