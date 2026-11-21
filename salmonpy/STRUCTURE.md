# 構造と最大格子間隔の指定

[English](STRUCTURE.en.md)

`salmonpy.structure` に `StructureInput`、`GridSelection`、`grid_from_spacing` を用意しました。
ASEは構造読込・コピー時だけ必要な任意依存です。通信・計算起動、繰り返し、セル標準化、
座標の折り返しを追加しません。スーパーセルは別途作成したファイルで指定します。

```python
from salmonpy.structure import StructureInput, grid_from_spacing
structure = StructureInput.from_file("supercell.xyz", format="extxyz", index=0)
atoms = structure.atoms  # Å単位の独立コピー。変更しても元の構造を変えない
selection = grid_from_spacing(atoms.cell.array, spacing=0.4, unit="angstrom")
print(selection.as_dict())
print(structure.record())
```

数値は指定方法の例であり、特定物質の収束済み推奨値ではありません。
メモリ上の入口は `StructureInput.from_atoms(atoms)` で、こちらもコピーします。
ファイルは `cif` または `extxyz` を明示し、非負整数のフレーム番号を指定します。
既定0は先頭で、ASE標準の最終フレームとは違います。同じファイル内容からhashと構造を取得し、
ファイル名の`@`をフレーム指定として解釈しません。ASEのMPI配信を無効にして読みます。
CIFの対称性展開は読込そのものに含まれるため、非対称単位のサイト数と原子数は一致しない場合があります。
追加のスーパーセル化はしません。原子配置を明示列挙した構造にはextxyzを使えます。
ASE/書式やOSの例外は伝播し、別形式への自動切替・単位の推測はしません。

`record()`はセル（Å）・周期性・原子数と、ファイルの場合は名前・形式・フレーム番号・
全ファイルのSHA256を返す独立辞書です。完全な構造アーカイブではなく、元ファイルを保持します。
`atoms`は読込後のセル・原子位置・種・ASE設定のコピーです。ローダーはDFT適合を認定せず、
電荷・磁性・非周期系に関する既存DFTの制限は維持します。

## 格子の仕様

希望する最大幅をスカラーまたは3成分で、`angstrom`／`bohr`の単位付きで指定します。
セルの入力単位は常にÅです。Bohr換算はSALMONの0.52917721067 Åを使います。
各軸でceil(L/h)から、最低点数（既定8）と明示した倍数（既定1）を満たすように切り上げます。
実際の幅L/Nは指定幅を超えず、入力セルは変更しません。戻り値`GridSelection`は
ヘルパーからだけ作成する不変tupleの集合で、`as_dict()`に希望幅・実際の幅（Å）、点数、
最低点数、倍数、総点数を返します。5×6×7 Å、最大幅0.4 Åなら(13,15,18)です。

初版は正の軸方向直交セルに限定し、非直交・回転セルの幅解決は拒否します。
セルだけを受ける単独ヘルパーでは呼出側が周期性を確認し、DFTのASE入口では3軸周期を検査します。
各軸2**20以下という入力上限は、メモリ予算の保証ではありません。倍数は利用者/backendが
明示する制約で、ランク数から推測しません。特別な整列が必要なbackendでは明示的に解決し、
`selection.counts`を`grid`へ渡してください。

## 既存DFT／TDDFTへの接続

`salmonpy.ase.dft_from_atoms` と `tddft_from_atoms` で、`grid=(...)`の代わりに
`grid_spacing=...` と `grid_spacing_unit=...`を指定できます。点数と幅の両方、両方なし、
幅なしの単位指定は拒否します。従来の点数指定の挙動・幾何の対応範囲は変えません。
`states`、`kmesh`、擬ポテンシャル、XCは別途指定し、収束や計算プログラム対応を推測しません。
DFTモデルは採用点数を保存するので、詳細レポートが必要なら作成前に格子ヘルパーを使います。
`SalmonCalculator`へ幅指定キーワードを追加したわけではなく、そちらにはヘルパーの点数を渡します。

## 物質IDを付けてオフライン準備する例

```sh
python -m salmonpy.examples.prepare_structure supercell.xyz new-inputs --format extxyz --material-id film --spacing 0.4 --unit angstrom --kmesh 2 2 2
python -m salmonpy.preparation new-inputs
```

共通物質registryへ構造・格子の検証関数を登録し、指定ファイルをコピーします。
物質ID・元ファイルhash・フレーム番号・セル・採用格子・明示kmeshをresolved.jsonへ保存します。
Pythonでは同じregistry・準備APIを使い、物質ごとに条件を変えられます。
サンプルは1回1ファイルの指定で、DFT backend登録や擬ポテンシャル選択、GS計算、
科学入力の完全性認定ではありません。ASEは事前にオフライン環境へ導入します。計算出力は入力一式の外へ置きます。

混合占有・部分占有のCIF/ASE構造は、代表元素や完全占有へ黙って変換しないよう読込時とDFT入口で拒否します。明示した秩序構造を用意してください。ASE infoの入れ子も独立コピーします。
