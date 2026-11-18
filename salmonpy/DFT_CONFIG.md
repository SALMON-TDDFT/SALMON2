# 物質設定からDFT・TDDFTへ接続する

[English](DFT_CONFIG.en.md) · Topic: `dft_material_config` · Paired revision: `2026-11-15.1`

`builtin_material_registry()`にbackend `salmonpy.dft`、model `dft`、アダプター版`1`を
明示登録します。ローカル構造・擬ポテンシャルを既存DFT APIへ渡すもので、設定解決時に
ランタイム起動・取得・コンパイルはしません。このresolverには任意依存のASEが必要ですが、
LD/SBEの登録はASEなしでも利用できます。

## 入力レシピ

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"backend": "salmonpy.dft", "model": "dft", "parameters": {
      "structure": {"path": "structures/cell.xyz", "format": "extxyz", "index": 0},
      "pseudopotentials": [{"path": "pseudo/Si_rps.dat", "atomic_number": 14, "local_channel": 2}],
      "xc": "pz", "states": 32, "kmesh": [2,2,2],
      "grid": {"spacing": {"value": 0.4, "unit": "angstrom"}}
    }}
  }
}
```

数値は書式例で、収束済み条件ではありません。構造、擬ポテンシャル、XC、states、kmesh、gridは必須。
indexの既定は0です。任意electronsは選択した擬ポテンシャルの価電子数から求める中性電子数と照合します。
gridは`{"counts":[16,16,16]}`か`{"spacing":{"value":0.4,"unit":"angstrom"}}`の片方だけです。
幅は3成分やbohrも指定できます。幅解決はSTRUCTURE.mdの正の軸方向直交・周期セルに限定し、
点数指定では既存DFTの幾何の契約を維持します。ランク数から格子を変えず、未知キーは拒否します。

参照はbase_directory内のローカル相対パスです。構造は明示CIF/extxyzで、CIFの対称性展開を除き
指定スーパーセルを保持します。混合/部分占有、電荷/磁性、非周期ASE条件は既存入口の範囲で拒否します。
擬ポテンシャルは既存KY/NC-UPF読込、原子種との一対一対応、中性電子数の検証を使います。
追加でNC-UPFヘッダーの元素名とatomic_numberを照合します。KYの元素メタデータや
擬ポテンシャルの汎関数/XC由来まで認定するものではなく、検証済みの適合した組を指定してください。
XCは既存のpz/libxc_pz/libxc_pbeで、Libxc全般やHSE対応を追加するものではありません。

採用parametersに元参照・正規化条件・電子数・structure_record・grid_resolution・file_sha256を保持します。
読込前後のhashと、モデル作成時の再読込を照合し、変更時は明示的な再解決を要求します。
点数モードは点数を、幅モードは希望幅と実際の幅も保存します。大きい波動関数をJSONへ格納しません。

## 実行とID継承

```python
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
from salmonpy.dft import SCFSettings
from salmonpy.runtime import load_runtime

registry = builtin_material_registry()
catalog = load_materials("inputs/materials.json", registry=registry)
material = material_from_catalog(catalog, "film", base_directory="inputs")
engine = load_runtime("runtime", required=["dft", "tddft", "pz", "single-rank"])
gs = material.run(engine, "new-gs", scf=SCFSettings(max_iterations=300, tolerance=1e-9), threads=1)
with gs.tddft_initial_state(executable=engine, max_steps=100, threads=1) as initial:
    with initial.new_state(dt=0.02) as state:
        j = state.step([0.001, 0.0, 0.0])  # 原子単位のA/c。戻り値は電荷電流j
        identity = state.material_id, state.definition_digest
```

ランタイム、SCF停止条件、threads、timeout、warm startは実行引数で、固定の物質定義から分離します。
選んだXCに必要なランタイム機能（必要ならLibxc）を指定してください。既存DFTの検査、収束報告、
GSからTDDFTへの検証済みローダーを使い、機能・物理条件の制限を解除しません。
結晶の基本はtransverse・A/c駆動です。2回のnew_stateで既存の独立子プロセス状態を作れます。

material.runはConfiguredGroundStateを返し、resultが従来のGroundStateResultです。他の属性は委譲します。
initial_guessにはラッパー結果も従来結果も渡せ、展開後に既存のwarm start互換性検査を通します。
tddft_initial_stateは従来TDDFTInitialStateを持つConfiguredMaterialを返し、with終了時に所有初期状態を
閉じます。closeを持たないモデルのwith終了は何もしません。時間発展状態も通常どおりwithで解放します。

物質IDと定義hashは、実行中のDFT→GS→TDDFT→stateで継承します。material/GSのrecord()を結果へ明示保存できます。
このhashは固定物質入力を識別し、SCF実行・GS実体・波形・dtの識別子ではありません。
GS/TDDFT自身のmanifestとデータhashは従来どおり別に保持します。汎用のJSON TDDFT backendや
catalog.initial_statesのレシピ選択ではありません。ディスクへのID自動保存、point/run ID、MPI配置、
任意の保存済みGSの選択は後続です。model/resultを他のAPIへ渡す場合、呼出側がID対応を保持します。

## オフラインSiサンプル

```sh
python -m salmonpy.examples.configured_si RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY
```

同梱Si入門データと明示したローカルランタイムを使い、構造・擬ポテンシャルを一式にまとめます。
短いGS/TDDFT計算で逆符号の波形を持つ2状態の交互実行と完全一致の再生を確認し、
既存期待値との許容誤差付き比較、物質ID・電流の記録を入力一式の外へ保存します。
各プロセス1スレッド、TDDFTワーカーは最大2つです。メモリ強制制限はありません。
接続確認であり、格子/k/時間収束済み光学計算ではありません。通信やランタイムビルドは行いません。

保存済みGS/1に限定した明示的なID選択は、別入口の
`tddft_from_catalog`で利用できます。[INITIAL_STATES.md](INITIAL_STATES.md)参照。
上記の汎用material構築自体は初期状態を選びません。

直列の点/run IDと結果対応表は、別層の[PointRun](POINT_RUNS.md)で保持・保存できます。
上記の未対応範囲は従来API単体の自動継承、MPI、restartに関するものです。

## 開発用：モデル所有の標準GS入力

標準ランタイムに内部Fortran型 `salmonpy_gs_model::s_gs_model` を追加しました。
`load(path,status)` がモデル固有のパーサーと入力値を保持し、`run(status)` がその入力を
有効化して既存のSCFを順次実行します。`clear()` は入力と作業領域を解放します。
未対応の入力へ置き換えようとしても、既存のモデルは保持します。

対象は周期・全k点・非偏極・固定イオンのPZ／TB-mBJ／Libxc-PZ/PBEです。
ホストがMPIを初期化し、`comm_attach`、数学定数・並列設定・タイマーを準備した後に、
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
--reference-sha256 SHA --size 2 --split spatial` です。明示的に監査済みのランタイムを指定し、
新しいプローブのソース・全リンク対象オブジェクトを検査してから実行します。
8³／12³、異なるセル・原子数・電子数・XCのA→B→Aを旧版の個別実行と比較し、
エネルギー・残差・反復数・固有値・k点出力、対称性復元、反復後の通信資源数を検査します。
`write_gs_restart_data='no'` を使い、再開GSバンドルは生成・収集しません。
4反復の診断用入力であり、収束した物性の検証ではありません。
