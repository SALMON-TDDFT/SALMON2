# 通常回帰試験とローカル開発検証

2026-10-06。`testsuites`にはSALMON入力、参照データ、CMake登録、準備・合否判定だけを置く。420番台の代表GS→RTと129/130のDC-GS回帰を維持する。

内部演算のドライバ・Python単体検証165ファイルはGitHub配布から除外し、ブランチ外のローカル開発記録へ保存した。`BUILD_DEVELOPER_TESTS`とCMake登録も削除。過去の報告・計画にあるパスと検証結果は当時の履歴であり、現行配布物のコマンドではない。

性能測定は`benchmarks`に維持する。性能測定に必要なドライバのみ各ベンチマークへ移し、RSS取得のC実装は`benchmarks/peak_rss.c`へ置いた。計算ソースと通常testsuitesの合否条件・許容誤差は変更しない。

標準のpreparationはsh、verificationはPythonによる出力検査、CMakeはcreate_test/create_mpi_test。追加420番台は明示的な失敗表示とsys.exit(-1)、成功時sys.exit(0)を使用する。
