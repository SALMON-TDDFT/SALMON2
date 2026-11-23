# IDを持つA駆動モデルをMPIへ接続

`bind_recorded_material(model, key=...)`は既存の端点Aアダプターを再利用し、
`model.record()`から識別情報を作ります。設定済みSBE、およびbackend/new_state/timingの
仕様を持つ自作Pythonモデルが対象です。適合試験やMPI初期化は暗黙に行いません。
準備済みモデルは借用するため、応答セッションを全て閉じてからモデルを閉じてください。

```python
from salmonpy.materials import bind_recorded_material, response_identity
from salmonpy.material_response import ResponseSession

binding = bind_recorded_material(model, key='prepared-film')
bindings = {10: binding, 90: binding}
point_map = {10: 'right', 90: 'left'}
record = response_identity(bindings, point_map, run_id='experiment-1')
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=point_map)
```

全ランクで同じrun IDと対応表を使い、[MPI_IDENTITY.md](MPI_IDENTITY.md)の集団呼出し・
解放順に従います。整数通信番号順に並べ、利用者のIDはそのまま保持します。
同じキャッシュkeyを共有する点には同じbindingオブジェクトを渡しますが、可変状態は点ごとに作ります。
物質群は各1 MPIランクに限定します。直列計算ルーチンを内部MPI並列化する機能ではありません。

`record()`は[PointRunの仕様](POINT_RUNS.md)に従い、物質ID・定義SHA256、任意で
初期状態IDとdigestの両方を返します。bindingは記録をコピーして不変化します。
MaterialBindingを直接作る場合の任意material_recordも検査・コピー・不変化します。
prepareと状態生成時にモデルの記録を再確認し、ResponseSessionもコールバック前に
点の記録とbindingの記録の不一致を拒否します。ただしhashは宣言値であり、記録を変えずに
内部状態を変更する自作モデルまでは検出できません。設定済みSBEの既存ローダーによる
ファイル検査と配列所有は維持し、ファイルシステムの常時監視は追加しません。

`response_identity`は記録付きbindingを要求し、整数の対応表と点間の物質・初期状態IDを
検査して独立したpoint-run/1記録を返します。worker起動・MPI通信はしません。
既存bind_materialと記録なしbindingは従来APIで引き続き使えます。

LDはE駆動のままです。この入口はA端点からEへ暗黙変換しません。
その変換が必要な場合は明示的なLDMaterialProvider/LDResponseBackend経路を使います。
TDDFTはConfiguredMaterialで包んだ場合もここでは拒否し、検証済み管理計算プログラムと
TDDFTMaterialProvider/TDDFTResponseBackendを使います。設定済みTDDFTをこの専用providerへ
つなぐfactoryは未実装です。集団Pythonから直接呼ぶTDDFTとrestart IDはこの入口の対象外です。

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_recorded_material /new/output
mpiexec -n 3 python -m salmonpy.examples.mpi_recorded_material /new/sbe-output \
  --catalog /local/inputs/materials.json --material my-sbe
```

前者は合成の自作モデル、後者は設定済みSBEで、後者には計算処理 SBE入りwheelが必要です。
MPI電流を直接の直列応答と比較し、保存したIDを読み戻します。通信によるデータ取得やコンパイルは行いません。
完全一致は回帰確認の観測であり、物理的収束や異なる機種でのビット一致の保証ではありません。

上記の専用factoryは別入口として追加しました。[CONFIGURED_NATIVE.md](CONFIGURED_NATIVE.md)参照。
汎用bind_recorded_materialからのTDDFT/LD自動選択は行いません。
