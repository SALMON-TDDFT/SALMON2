# 回帰試験と開発検証の分離

`testsuites`には通常のSALMON入力、参照データ、CMake登録、準備・合否判定を置く。内部演算用FortranドライバとPython単体検証はすべて`developer_tests`へ移動した。

開発検証は交換/ACE・DC-LCFO・汎関数の3機能に集約。`BUILD_DEVELOPER_TESTS`は既定OFF。ON時だけLCFO・学習ACE数学試験をビルド・CTest登録する。手動機能検証の手順はdeveloper_tests/README.md。

420番台の代表GS→RTと129/130のDC-GS回帰は維持。性能測定はbenchmarks、GPU試作検証はexperiments、機種設定はtools/validation。過去の計算記録中の当時のパスは保存する。

検証: 開発検証OFFでは4内部数値試験が未登録、ONではLCFO・学習履歴2件・LAPACKの4試験合格。Python構文・固定パス・420番台fixture依存とtestsuites内Fortran/Cドライバゼロを確認。代表GS→RT全計算とGPU実機検証は今回未実施。

既存111/112の形式に準拠: preparationはsh、verificationはPythonによる出力検査、CMakeはcreate_test/create_mpi_test。追加420番台のassert判定を明示的な失敗表示とsys.exit(-1)へ変更し、成功時sys.exit(0)を明記。判定条件・許容誤差・物理入力は変更していない。Python最適化実行でも合否判定が消えない。基準値のない正常終了・有限値検査を物理精度の基準値回帰と混同しない。
