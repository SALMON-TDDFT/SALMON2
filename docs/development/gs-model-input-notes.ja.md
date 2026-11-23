Implementation notes recorded 2026-11-24. Current supported Python interfaces: [retained GS](../../salmonpy/NATIVE_GS.md).

## 開発用：モデル所有の標準GS入力

標準計算プログラムに内部Fortran型 `salmonpy_gs_model::s_gs_model` を追加しました。
`load(path,status)` がモデル固有のパーサーと入力値を保持し、`run(status)` がその入力を
有効化して既存のSCFを順次実行します。`clear()` は入力と作業領域を解放します。
未対応の入力へ置き換えようとしても、既存のモデルは保持します。

対象は周期・全k点・非偏極・固定イオンのPZ／TB-mBJ／Libxc-PZ/PBEです。
呼び出し元プログラムがMPIを初期化し、`comm_attach`、数学定数・並列設定・タイマーを準備した後に、
同じ群の全ランクで同じモデルを同じ順序で呼びます。入力は原子座標の読み込み後に保存し、
実行前後で呼出側の入力ビュー・単位表・calc_mode・対称性の実スイッチ・最適化表・XC作業領域を
切り替えます。モデル用の呼出しは通信キャッシュをローカル配列より先に終了し、
各回の子通信資源を解放します。既定の run(status) は状態を保持しない順次GS入口です。run(status,state) による保持と
Pythonハンドルは[保持GS](NATIVE_GS.md)で利用できます。同時GS／TDDFT、
スレッド並行、HSE／SOC／+U／DCのGS独立化は別作業です。

statusは0が正常復帰、1が空モデル／使用中の子通信スコープ、2がこの内部入口の未対応条件です。
SCF未収束は別扱いで、従来どおり警告して最後の反復結果を出力します。
パーサーの致命的エラー・STOP・メモリ不足などの一般復帰は保証しません。
相対パスの基準は呼出側の作業ディレクトリです。通常のPython `DFT.run` は従来どおり
別プロセスを使用し、この内部入口へ自動的には切り替わりません。

開発者用の再現試験は `tools/test_standard_gs_models.py NEW_BUILD OLD_BUILD NEW_OUTPUT
--reference-sha256 SHA --size 2 --split spatial` です。明示的に監査済みの計算プログラムを指定し、
新しいプローブのソース・全リンク対象オブジェクトを検査してから実行します。
8³／12³、異なるセル・原子数・電子数・XCのA→B→Aを旧版の個別実行と比較し、
エネルギー・残差・反復数・固有値・k点出力、対称性復元、反復後の通信資源数を検査します。
`write_gs_restart_data='no'` を使い、再開GSバンドルは生成・収集しません。
4反復の診断用入力であり、収束した物性の検証ではありません。

