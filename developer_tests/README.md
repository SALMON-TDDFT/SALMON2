# 開発用数値検証

Fortranドライバ、単体試験、ビルド・実行スクリプトを保管する。通常の入力ベース回帰試験はtestsuitesに置く。

CMake単体試験は既定OFF。`-DBUILD_DEVELOPER_TESTS=ON`でCPU用LCFO時間発展・学習ACE履歴を登録する。GPU実機試験を登録するものではない。

機能単位の手動検証:

```sh
python3 developer_tests/651_hybrid_exchange/check_function.py --build /path/to/build
python3 developer_tests/652_dc_lcfo/check_function.py --build /path/to/build
python3 developer_tests/653_functional/check_function.py --build /path/to/build
```

ハイブリッド機能窓口は学習履歴CTestも使うため開発検証ONのビルドを指定する。`--detailed`で並列条件を拡大する。通常testsuitesはSALMON本体、入力、参照データ、準備・合否判定のみで実行する。
