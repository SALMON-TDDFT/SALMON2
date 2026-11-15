# 組み込みLD・SBEの物質設定

[English](BUILTIN_MATERIALS.en.md) · Topic: `builtin_materials` · Paired revision: `2026-11-15.1`

`builtin_material_registry()`で、`salmonpy.ld`／`ld`と`salmonpy.sbe`／`sbe`の
検証関数を持つMaterialRegistryを明示作成します。アダプター版は`1`です。
登録・設定解決・モデル作成ではnative状態の初期化・取得・コンパイルを行いません。
`new_state`では既存の導入済みLD/SBEライブラリが必要です。

```python
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
registry = builtin_material_registry()
catalog = load_materials("inputs/materials.json", registry=registry)
material = material_from_catalog(catalog, "film", base_directory="inputs")
with material.new_state(dt=0.02) as state:
    print(state.material_id, state.definition_digest)
```

モデルに応じて状態引数を選びます。LDは**E駆動**で`initial_field`、SBEは**A/c駆動**で
`initial_A`と`time0`を使います。数値計算・電流・snapshot/restore・contextでの解放は
既存状態へ委譲し、契約や更新順を変えません。LDをA駆動に変換する機能ではなく、
その場合は既存の明示的なLD結合アダプターを使います。

## LDパラメータ

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"backend": "salmonpy.ld", "model": "ld", "parameters": {
      "units": "atomic_gaussian", "plasma_frequency": 0.5,
      "strength": [1.0], "resonance": [0.3], "damping": [0.1]
    }}
  }
}
```

5項目すべて必須です。単位はGaussian原子単位・角周波数を明示し、受動的な等方極と
背景誘電率1という既存契約です。LorentzDrudeの検証を使い、未知キー・他の単位は拒否します。
値は書式例であり実測物性値ではありません。eV/fsからの自動換算は行いません。

## SBEパラメータ

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "silicon": {"backend": "salmonpy.sbe", "model": "sbe", "parameters": {
      "units": "atomic", "directory": "sbe", "sysname": "si",
      "bands": 2, "kpoints": 1, "electrons": 2,
      "lattice": [[4,0,0],[0,4,0],[0,0,4]]
    }}
  }
}
```

これは書式例で、実際のSiバンドデータではありません。表示した項目は必須です。
任意指定は`active_bands`（既定は全バンド）、`correction_order`（既定0）、
`nonlocal_correction`（既定false）。格子は**列がベクトル・bohr単位**です。
ASEの行ベクトル・Åを無変換で渡さないでください。directoryは入力基準からの相対パス
（`.`可）、sysnameは単一ファイル名の接頭辞です。対応する3ファイルは
`si_k.data`、`si_eigen.data`、`si_tm.data`です。単一スピン・原子単位のテキスト、
閉殻占有（出力ファイルの占有数列は使いません）・行列形状等は既存[SBE.md](SBE.md)の制約を維持します。
バイナリcheckpointや汎用NPZローダーを追加したわけではありません。

解決時に配列を読込・検証し、読込前後で3ファイルのhashを照合します。
採用parametersには`file_sha256`を記録し、モデル作成時も再読込してカタログと照合します。
変更されたファイルを同じ物質として黙って使用せず拒否します。毎stepでなく初期化時の処理です。
ファイル容量・メモリ予算の自動推定はなく、カタログ作成と後のモデル作成で再度読み込みます。
この2回の読込ではnativeハンドルは作りません。

## IDと独立状態

`material_from_catalog`の戻り値`ConfiguredMaterial`は、`material_id`、
`definition_digest`、`model`（従来のLD/SBEモデル）、ID情報を返す`record()`を持ちます。
`new_state(...)`は同じ不変ID情報を持つ`ConfiguredState`を返し、既存状態の操作を使えます。
別の呼出しで独立状態を作り、それぞれ解放します。他の既存モデル公開属性も委譲します。
具体的な型を要求する利用先には`.model`を渡せますが、その場合は呼出側がID対応を保持します。

snapshotは従来の同一モデル・固定dtのメモリ保存で、移設可能なcheckpointではありません。
初期状態はLDのゼロ履歴、SBEの閉殻状態です。カタログの任意`initial_states`のsource/prepareは
このfactoryでは**選択・消費しません**。励起状態を暗黙に読み込むこともありません。

自作resolverをregistryに追加できますが、このfactoryは既知版の組み込みモデルだけを作成します。
backend名から任意importを行いません。任意の保存済みTDDFTアダプター、一般的な能力スキーマ、
点/run ID、MPI配置、出力/restartへの自動ID継承は別項目です。
出力関数はwrapperを自動認識しないので、結果と一緒に`material.record()`を明示保存します。

## オフライン準備と実行例

同じregistryを`prepare_inputs`へ渡し、SBEの3ファイルを入力内の相対パスで列挙します。
検証関数はコピー先で動くため、参照不足なら保存先確定前に拒否します。
転送後は一式を照合し、このregistryでmaterials.jsonを読み、転送先をbase_directoryにして
モデルを作ります。カタログが物理収束や実物質データの正しさを保証するわけではありません。

```sh
python -m salmonpy.examples.configured_ld new-ld-run
```

合成の二物質例です。入力一式を準備して独立状態を確認し、ID・電流を入力一式の**外**の
response.jsonに保存します。完全版SALMONランタイムや通信は不要です。

DFT設定とGS→TDDFTへのID継承も[DFT_CONFIG.md](DFT_CONFIG.md)で利用できます。LD/SBEのE/A契約は変わりません。

保存済みGS/1に限定した明示的なID選択は、別入口の
`tddft_from_catalog`で利用できます。[INITIAL_STATES.md](INITIAL_STATES.md)参照。
上記の汎用material構築自体は初期状態を選びません。

直列の点/run IDと結果対応表は、別層の[PointRun](POINT_RUNS.md)で保持・保存できます。
上記の未対応範囲は従来API単体の自動継承、MPI、restartに関するものです。
