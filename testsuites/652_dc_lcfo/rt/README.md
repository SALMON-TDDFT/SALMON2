# 652_dc_lcfo/rt

機能単位の数値検証。実行窓口は一つ:

```sh
python3 testsuites/652_dc_lcfo/rt/check_function.py --build /path/to/build
```

内部の複数ケースは同じ機能の数値一致・境界条件・並列条件を検証する。未選択の旧開発スクリプトは標準窓口では実行しない。GS/RTの統合試験は420番台を使用する。GPU試験は含まない。
