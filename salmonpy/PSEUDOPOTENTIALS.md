# 擬ポテンシャルを明示的に取得する

[English](PSEUDOPOTENTIALS.en.md)

計算中には取得しません。インターネットにつながる管理ノード等で、利用者の許可を得て実行し、入力と一緒に計算ノードへ転送します。importだけでは通信しません。

同梱Si例と同じ、固定コミットのApache-2.0データを取得する例です。

```sh
python -m salmonpy.pseudopotentials NEW_DIRECTORY --silicon-example --allow-download
```

Pythonからは次のように呼びます。保存先の親フォルダは先に用意してください。

```python
from salmonpy.pseudopotentials import download, silicon_reference
path = download(silicon_reference(), 'NEW_DIRECTORY', allow_download=True)
```

ファイルと `download.json` に出典・ライセンス・URL・SHA256を保存します。既存の保存先は上書きしません。HTTPS（転送先を含む）、64 MiB以下、参照と同じSHA256を確認します。不一致なら計算用ファイルを公開せず、参照値も更新しません。

他のデータは `PseudopotentialReference(filename, atomic_number, url, sha256, license, source)` を指定します。同じフィールドを持つローカルJSONは `--reference FILE.json` でも使えます。SHA256は信頼できる出典から取得し、ダウンロードした未知のファイルから自動採用しないでください。ライセンス欄は権利者の宣言を記録するもので、任意のURLの再配布権を保証するものではありません。利用者は利用条件を確認します。

取得はファイルの同一性だけを確認します。形式・XC・価電子数・局所チャネル・格子収束の適合を保証しません。計算に渡すときは `Pseudopotential.from_file(path, ...)` で条件を明示します。GS/2は引き続き参照のみを保存し、取得を自動実行しません。

## PseudoDojoの選択とSOC

現在のカタログはPBE v0.4の固定版です。スカラー相対論版はH・C・N・O・Al・Si・Ti・Se・W、完全相対論版はWSe₂用のW・Seと、CI用のNiを登録しています。Niは3s・3pセミコアを含む18電子です。[Niの短時間SOC検査](NICKEL_SOC.md)。カタログはファイルの同一性とヘッダーの記録であり、全元素の実計算検証済みという意味ではありません。

```sh
python -m salmonpy.pseudopotentials --list --element W --xc PBE --relativity fully-relativistic
python -m salmonpy.pseudopotentials NEW_W_DIRECTORY --catalog-id pseudodojo-pbe-0.4-fr-standard-w --allow-download
```

```python
from salmonpy.pseudopotentials import catalog_reference, download
from salmonpy.dft import Pseudopotential
reference = catalog_reference('pseudodojo-pbe-0.4-fr-standard-w')
path = download(reference, 'NEW_W_DIRECTORY', allow_download=True)
potential = Pseudopotential.from_file(path, atomic_number=74, local_channel=3,
                                     reference=reference)
```

版・XC・相対論条件・価電子数・NLCC・SOCフラグ・出典コミットを記録します。取得時とPythonへの読込み時にヘッダー・ハッシュを照合します。原子番号や元素名だけで別の版へ置き換えません。物理的な収束や局所チャネルの適合は別途確認します。

PseudoDojoのデータはCC-BY-4.0です。利用時はPseudoDojoとONCVPSPに帰属を示し、関連文献を引用します：[公式プロジェクト](https://github.com/PseudoDojo/pseudodojo)、[PseudoDojo論文](https://doi.org/10.1016/j.cpc.2018.01.012)。固定UPFをwheelへ一括同梱せず、参照カタログだけを保存します。

**SOCの全経路対応は整備中です。** 完全相対論ファイルの取得・読込みとSOC用DFT入力の生成に対応しています。[Ni検証](NICKEL_SOC.md)ではSCFの状態独立性・別プロセス再開と、TDDFTの状態独立性・直接復元・同じMPI分割での別プロセス再開を確認しました。MPI分割変更などは未検証です。SOCが必要な指定をスカラー相対論計算へ切り替えることはありません。最初の検証対象は単層WSe₂です。
