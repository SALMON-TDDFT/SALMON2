# 開発検証の分離

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

承認済み方針: testsuitesには入力・参照データ・標準準備/合否判定を残す。3機能のFortran/Python内部検証はdeveloper_testsへ同じ階層深さで移動。CMake登録はBUILD_DEVELOPER_TESTS既定OFFで分離。参照、構文、ON/OFF登録、ビルド、数値検証を確認する。既存出力と固定計算バイナリは変更しない。
