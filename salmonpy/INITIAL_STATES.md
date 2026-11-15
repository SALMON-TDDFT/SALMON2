# 保存済み基底状態をIDで選択

`tddft_from_catalog`は、明示設定したDFT物質に対してSALMONpyで生成したローカルの
`salmonpy.gs/1`基底状態を選びます。SCF・ダウンロード・コンパイル・native workerの起動は
行いません。`new_state(dt=...)`で独立workerを開始します。任意依存のASE/f90nmlと、
検証済みのローカルDFT/TDDFTランタイムを先に用意してください。

`salmonpy.materials/1`の`materials`と並べて次を記述します。

```json
"initial_states": {
  "Si-ground": {
    "material_id": "Si-starter",
    "source": {
      "format": "salmonpy.gs/1",
      "path": "saved",
      "manifest_sha256": "GROUND_STATE_JSONのSHA256に置換"
    }
  }
}
```

pathは明示した基準ディレクトリ内の相対パスです。物質IDやURLではありません。
`ground-state.json`とその`files`に列挙された全ファイルを保存してください。物質設定の
構造・擬ポテンシャル・ライセンス・出典も保持します。`prepare_inputs`に転送先名→元ファイルの
一覧を明示して収集し、転送後は別途記録したmanifestのhashで`verify_inputs`を実行します。
[オフライン準備](PREPARATION.md)を参照。計算ノードで通信する必要はありません。
native restartの機種間互換性は保証しません。この入口はGS作成時と同じ実行ファイルのSHA256を
要求し、実行先自体もランタイムに対応する必要があります。

```python
from salmonpy.builtin_materials import builtin_material_registry
from salmonpy.config import load_materials
from salmonpy.initial_states import tddft_from_catalog
from salmonpy.runtime import load_runtime

engine = load_runtime(runtime_directory, required=('dft', 'tddft', 'single-rank'))
catalog = load_materials(bundle / 'materials.json', registry=builtin_material_registry())
with tddft_from_catalog(catalog, 'Si-ground', base_directory=bundle,
                       executable=engine, max_steps=20, threads=1) as model:
    with model.new_state(dt=0.02) as first, model.new_state(dt=0.02) as second:
        j = first.step([0., 0., 1.e-5])  # 原子単位のA/c。secondは変わらない
    metadata = model.record()
```

manifestのhash、nativeファイル、SCF完了、原子の順序を含む構造、セル、電子数、状態数、
格子、k点、XC、PP実体、実際のGS科学入力、コピー後のsnapshotを検査します。
原子順やセルを暗黙に変換せず、不整合や不正なmetadataを拒否します。
hashはバイト列を固定するもので、物理的収束や出所の信頼性を単独で保証しません。
短いSi入門例は接続検査で、光学物性の収束例ではありません。

返すモデルと状態に`material_id`/`definition_digest`、`initial_state_id`/`initial_state_digest`
を保持します。後者のdigestはGS manifestのSHAで、時間発展中のcheckpointではありません。
結果と一緒に`record()`を明示保存します。point/run ID、ディスクへの自動ID保存、MPI割当は別作業です。
汎用material作成は初期状態を選びません。`prepare`レシピ、任意のSALMON restart、GS/2は未対応です。

実行サンプルには、原子単位のscalar `dt`と有限値の`A`（形状`(steps+1, 3)`、`A[0]=0`）を
NPZで渡します。

```sh
python -m salmonpy.examples.catalog_tddft inputs/materials.json Si-ground \
  /path/to/local/runtime waveform.npz /new/output/directory
```

電流配列とID metadataを保存します。電流は電荷電流密度j = -J_SALMONです。
符号・時刻は[FIELDS.md](FIELDS.md)参照。状態ごとにプロセスとsnapshot保存領域を使い、
共有メモリTDDFTやメモリ予算の保証ではありません。スレッド数・同時worker数を資源に合わせてください。

直列の点/run IDと結果対応表は、別層の[PointRun](POINT_RUNS.md)で保持・保存できます。
上記の未対応範囲は従来API単体の自動継承、MPI、restartに関するものです。
