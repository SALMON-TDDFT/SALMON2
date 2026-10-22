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


## 入力制御の自動化（2026-10-22）

SALMONでは拡張KYのヘッダーを検出すると、省略した`spin`と`yn_spinorbit`を
`noncollinear`と`y`へ解決します。`spin='noncollinear'`または`yn_spinorbit='y'`の
片方だけの指定でも他方を補完します。旧KYで両方省略した場合は従来どおり
unpolarized/SOCなしです。明示的なSOC無効化や不整合なスピン指定は拒否し、
自動的なスカラー平均はしません。解決済みの値をvariables.logへ記録します。
他形式のSOC情報の自動検出はこの変更の対象外です。
拡張KYの`lloc_ps`指定も不要です。ファイル内に格納された局所ポテンシャルを使います。
旧KYの局所参照チャネルは物理的な選択なので、この省略規則は適用しません。

生成器の`&crps`では、以下の制御指定を省略できます。

- `Lmax`: 先頭から連続する3/5/7チャネルの半径リストから1/2/3を導出。
  部分・疎なリストを上書きする場合はLmaxを明示してください。
- `NA`: 有効なPtoAの最大参照軌道番号。占有数合計と核電荷の検査は維持。
- `lloc`: 解決後のLmax。W標準のlloc=0は物理的選択なので明示を維持。
- `output_file`: 入力namelistの拡張子を除いたbasenameに`_soc_nlcc_rps.dat`を追加。
- `c2lo`: 既存の探索下限-10.5。標準W/Seの上限c2hi=5は引き続き明示。

明示した値は上書きせず、引数なしの従来Au初期値も維持します。
元素、占有数、参照軌道、半径、NLCC芯半径、必要な局所チャネルは物理入力として残します。
W/Seの省略入力と明示入力は、同じ生成器でbyte一致し、配布済み標準表とも数値照合します。
