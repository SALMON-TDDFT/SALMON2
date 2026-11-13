# SALMON VR公開版とのファイル連携

[English](SALMON_VR.en.md)

独立モジュール `salmonpy.salmon_vr` は密度差とA/cのCUBE時系列を公開Mac版向けに出力します。
SALMON本体・ASE・描画ライブラリは不要です。録画操作やリアルタイム通信は対象外です。

```python
from salmonpy.salmon_vr import export_cubes
result = export_cubes(
    ['drho_step0.cube', 'drho_step8.cube'], 'new-view',
    steps=[0, 8], dt_au=0.02,
    missing_ac='zero',  # A未出力のときだけ、表示用ゼロを明示的に選択
    # ac_files=['Ac_step0.cube', 'Ac_step8.cube'],  # 実際のA/cがある場合
    # application='/path/to/SALMON VR.app',  # 指定したアプリを新規保存先へコピー
)
# result.launch()  # application指定時のみ、Macで明示起動
```

`dt_au` は出力間隔ではなく計算の1ステップ幅です。stepsは0～999999の厳密昇順整数。
間引いた時刻は保持し、表示側が線形補間します。入力順序で対応付けるため、ファイル名から
時刻や物理量を推定しません。入力は原子単位、密度は初期状態からの差、AはA/cであることを
利用者が確認してください。格子、原点、原子配置は全フレーム・両物理量で照合します。
全密度を密度差に変換したり、EからAを積分したりする処理はありません。

出力は `config.json`、`Ac_Je_drho_CUBE/it000008_drho.cube` と対応する `it000008_Ac_x.cube`、
由来・時刻・ハッシュを記録する `salmonpy-vr.json` です。元ファイルは変更しません。
既存保存先は拒否し、失敗時は今回の新規出力だけを削除します。

Aなしでは通常エラーにします。`missing_ac='zero'` の場合は画面用説明と記録に
「A unavailable: zero placeholder」を残します。これは物理的なA=0の計算結果ではありません。
公開Mac版でゼロ補完の読み込み、密度等値面と時刻更新を確認しました。任意のデータでの表示を保証するものではありません。

公開Mac版はアプリと同階層の設定・データを読みます。`application` を省略した場合は、
公開アプリを生成フォルダーへ配置してFinderから起動してください。コピー容量も確保します。
起動時にOSの確認が出た場合は利用者が判断します。このモジュールは隔離属性やOS設定を変更しません。
公開マニュアルの対象はApple Silicon・macOS 14.5以上です。書き出し自体は他OSでも利用できますが、
Windows版の配布物・起動は今回検証していません。

既存 `OutputArchive` では、索引のCUBEファイルを選び、そのパスを渡します。
物理量・時刻は実行条件から明示してください。CUBE保存を選んでいない計算では密度動画を復元できません。
大きいCUBEは読み込みメモリと保存容量に注意してください。公開アプリは全フレームを読み込みます。

サンプル: `python -m salmonpy.examples.salmon_vr_export --help`

仕様: [公開Macマニュアル](https://salmon-tddft.jp/download/macos_manual.pdf)、
[公開アプリ](https://github.com/SALMON-TDDFT/SALMON-VR)、
[Otobe et al., Virtual Reality 30, 52 (2026)](https://doi.org/10.1007/s10055-026-01319-0)。

初期対応は正方向の軸に沿った直交格子です。斜交・回転格子は拒否します。単一フレームや先頭ステップが0以外の表示は未検証です。公開tarの実行権限が不足する場合は公式Macマニュアルの起動手順を参照してください。

`launch()` はInfo.plistの実行ファイル指定、ファイルの存在、実行権限を起動前に確認します。
不足時は起動せず、対象パスを引用した `chmod u+x` の手順をエラーで案内します。
操作は利用者がアプリを確認して行います。アプリコピーや起動で権限を自動変更しません。
