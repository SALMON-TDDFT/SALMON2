# 物質ごとの設定カタログ

[English](CONFIG.en.md) · Topic: `material_config` · Paired revision: `2026-11-15.1`

`salmonpy.config` は設定だけを扱うプレビューです。Python辞書とUTF-8 JSONを
`salmonpy.materials/1` で読みます。設計中の全計算スキーマとは別です。
ソルバー生成、構造/GS読込、ダウンロード、MPI初期化、能力の認定は行いません。
組み込みバックエンドの自動登録もありません。
`python -m salmonpy.examples.material_config` はソルバー不要の実行例です。

## 入口

- `MaterialRegistry()`：明示的なローカル登録簿。
- `registry.register(name, version=..., models=[...], resolve_parameters=...)`：
  信頼するPython検証関数を登録。名前やモデルの重複は拒否します。
- `resolve_materials(data, registry=..., base_directory=...)`：PythonのJSON値を
  検証します。基準ディレクトリは必須。tuple、NumPy配列、Path、自作オブジェクトは
  JSON値ではないため、呼び出し前に明示変換してください。
- `load_materials(path, registry=...)`：同じ検証経路。相対参照は入力ファイルの
  ディレクトリを基準にします。JSONから任意importを行いません。
- 戻り値 `MaterialCatalog` の `materials`、`initial_states`、`definition_digests`
  は再帰的に不変です。`as_dict()` は独立したJSON保存用コピーを返します。
  その形式は `salmonpy.materials-resolved/1`。入力レシピやrestartではなく、
  `load_materials` へ再投入できません。

## 入力契約

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"model": "lab:response", "backend": "local",
             "parameters": {"coefficient": 2.0}, "label": "First layer"}
  },
  "initial_states": {
    "film-ground": {"material_id": "film", "prepare": {"kind": "equilibrium"}}
  }
}
```

`initial_states` と `label` は省略可能です。初期状態の各項目には、実在する
`material_id` と、辞書である `source` または `prepare` の片方だけを指定します。
この辞書はレシピの記録だけです。内部キー、外部ファイルの存在・ハッシュ、物理的整合性は
まだ検証しません。準備済み状態として使用しないでください。配置やMPI割当など未知の
最上位キーは拒否します。IDの大小文字・Unicodeはそのまま保持し、空ID、前後空白、
制御文字、重複JSONキー、非有限数値は拒否します。IDはファイルパスではありません。

## 自作の検証関数

署名は `resolve_parameters(model, parameters, *, base_directory)`、戻り値は通常の
JSON辞書です。ホストが明示的に信頼したPythonコードであり、サンドボックスではありません。
未知パラメータ・単位・範囲・組合せを検証し、単位変換と既定値補完、由来の記録をして
採用値を返します。ここでソルバー起動、状態確保、共有設定の変更を行わないでください。
入力と戻り値はコピーされ、元データと採用設定は可変領域を共有しません。
スキーマ・ID・参照の不正は検証関数を呼ぶ前に検出します。後続関数の失敗時に先行関数の
副作用は巻き戻せません。関数のValueError/TypeErrorは物質の設定パスと元例外を伴う
InputErrorにし、それ以外の例外やOSエラーはそのまま伝えます。

SchemaからAIがパラメータを生成する機能は未実装です。登録しただけでMPI、保存復元、
スレッド、複数状態対応にはなりません。現在の計算モデルの契約は
[MATERIALS.md](MATERIALS.md)を参照。組み込みLD/SBEへの接続は末尾のガイドを参照してください。

## 定義のハッシュ

backend名・版、model、採用parametersのSHA256を計算し、物質IDとlabelは除きます。
正規化規則v1はキー整列・空白なしJSON・UTF-8・有限数値で、`1`と`1.0`は区別します。
解決規則を変えたらbackend版を更新してください。ファイル内容や関数コードは自動hashしません。
内容識別が必要なら、検証関数でデータhashと移設可能な参照を含め、機種固有の絶対パスを
識別に使わないでください。設定のハッシュであり、物理精度や全機種の結果一致の保証ではありません。

## 独立した構造・格子アダプター

ASE AtomsとローカルCIF/extxyzをSTRUCTURE.mdの独立アダプターへ渡せます。
スーパーセル・複雑な構造は別途ファイル指定とし、そのセルと原子配置をそのまま使用します。
暗黙の繰り返し・回転・標準化はしません。出典・形式・必要なフレーム番号・内容hash・単位を記録します。
MateriAppsは検索の入口で、データは選択した提供元（COD/Materials Project等）の利用条件に従い
取得し、明示的アダプターで読み込みます。構造だけからXC、擬ポテンシャル、電荷、磁性、
収束条件を自動決定しません。

登録した検証関数からgrid_from_spacingを使い、物質ごとに単位付きの最大格子間隔を指定できます。対応する直交周期セルでは
N_i = ceil(L_i/h_i) を出発点に、選択バックエンドの最低点数・必要な倍数へ切り上げます。
希望幅・採用点数・実際のL_i/N_iを表示・保存します。MPI資源変更で確定済み格子を黙って変えません。
点数の直接指定も残しますが、幅との同時指定はASE DFT/TDDFT入口で拒否します。
非直交・孤立系は専用の幾何検証が必要です。FDTD格子とは分離し、切り上げを精度保証としません。
直接DFTは点数指定です。ASEのDFT/TDDFT入口には、下記の単位付き幅指定も追加されました。

## オフライン環境での構造指定

計算ノードからのデータベース検索・取得を前提にせず、接続可能な準備環境で
構造・スーパーセルファイルと全参照データを揃え、出典・版・SHA256を記録して転送します。
計算入力ではローカル参照を指定し、転送後にハッシュと相対パスを確認します。
MateriApps等での検索とジョブ実行は分離します。欠けたファイルを実行時に自動取得する
設計にはしません。依存wheel・ランタイムも含む準備は
[RUNTIME.md](RUNTIME.md)のオフライン計算ノードの項を参照してください。
構造ファイルはSTRUCTURE.mdのアダプターで扱い、カタログ自体は自動収集・検証しません。

入力一式の明示的なローカル収集・転送後照合は[PREPARATION.md](PREPARATION.md)で利用できます。構造解釈・依存自動探索は別です。

構造・幅解決の独立ヘルパーと明示登録例が追加されました：[STRUCTURE.md](STRUCTURE.md)。本カタログ自身の自動読込ではありません。

LD/SBEは[BUILTIN_MATERIALS.md](BUILTIN_MATERIALS.md)の明示登録とモデル作成を利用できます。

保存済みGS/1に限定した明示的なID選択は、別入口の
`tddft_from_catalog`で利用できます。[INITIAL_STATES.md](INITIAL_STATES.md)参照。
上記の汎用material構築自体は初期状態を選びません。
