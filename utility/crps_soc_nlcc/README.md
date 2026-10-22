# CRPS SOC+NLCC 擬ポテンシャル生成器

2026-10-22に標準選択したW・Seを生成するソースと入力です。

- `crps_soc_nlcc.f`：元のcrpsAu.fから拡張した多元素・j分解SOC生成器。
- `crps_output.f90`：`KY_SOC_NLCC_V1`出力。NLCCの原点近傍一次微分の符号修正を含みます。
- `inputs/W.nml`：W radius2、s1.59/p1.01/d1.38/f1.59 Å、価電子12、NLCC芯半径0.65 Å。
- `inputs/Se.nml`：Se s1.38/p1.43/d1.22 Å、価電子6、NLCC芯半径0.50 Å。

## ビルドと生成

GNU Fortranが必要です。以下はこのディレクトリから実行する例です。生成先は新しい空ディレクトリを使ってください。出力済みファイルがある場合は上書きせず停止します。

```sh
gfortran -std=legacy -O0 -fcheck=all -fbacktrace -Wno-align-commons \
  crps_output.f90 crps_soc_nlcc.f -o crps_soc_nlcc
mkdir W-run Se-run
(cd W-run && ../crps_soc_nlcc ../inputs/W.nml > stdout.txt 2> stderr.txt)
(cd Se-run && ../crps_soc_nlcc ../inputs/Se.nml > stdout.txt 2> stderr.txt)
```

Wは`W_soc_nlcc_rps.dat`、Seは`Se_soc_nlcc_rps.dat`が出力されます。共通配布先の`testsuites/pseudo/W_rps.dat`、`Se_rps.dat`とはファイル名だけが異なります。標準ファイルへの自動コピーや置換は行いません。

c2探索範囲は入力の`c2lo/c2hi`（現在−10.5〜5.0）です。各チャネルで一意に符号変化する区間のみ採用し、解がない／複数のときは停止します。範囲を広げて分岐を自動選択する処理はありません。

生成器は元の相対論的LDA交換・PZ型相関を継承しています。PBE用生成器ではありません。原子の価電子構成やカットオフ半径を変更した場合は別の候補として検証してください。動径ファイルは0.01 Å刻みなので指定半径は丸められます。

## 使用と検証範囲

生成ファイルの使用には`KY_SOC_NLCC_V1`対応readerとSOC設定が必要です。通常の数値ヘッダーだけの旧KY形式ではありません。本ブランチのSALMONには対応readerを統合しています。ファイルのヘッダーで旧KY／拡張形式を判別し、拡張形式だけはSOC有効時に読みます。SOC無効での暗黙の平均化は行わず拒否します。新readerを使うにはSALMON本体の再ビルドが必要です。稼働中の既存実行ファイルは変更しません。

標準W・Seとの再生成照合結果と全ファイルのSHA256は`provenance.json`を参照してください。固体・メッシュ・光学精度の検証状況は`testsuites/pseudo/crps_soc_nlcc/README.md`に記録しています。元のcrpsAu.fと過去の失敗出力はwork配下に保存しており、変更していません。

Coding Rules対応として出力writerをprivate既定のmoduleに置き、追加した探索範囲・jチャネル配列はCOMMONではなく明示引数とintentで渡しています。450はparser・生成器、451は本体serial／MPI読込の回帰試験です。通常CTestは有限値・構造・数値許容差で比較し、同一ツールチェーンでのバイト一致は来歴に別記します。
