# 設定済みLD・TDDFTの専用binding

数値・計算プログラムの仕様を守るため、入口を分けます。

```python
from salmonpy.materials import (
    bind_configured_ld_endpoints, bind_configured_tddft, response_identity,
)

# 対応するConfiguredMaterialに対して選択する。
ld_binding = bind_configured_ld_endpoints(ld_material, key='ld-cache')
tddft_binding = bind_configured_tddft(tddft_initial, key='gs-cache', runtime=runtime_directory)
```

LDは既存LDResponseBackendの`E_new = -(A_new-A_old)/dt`で端点AをEへ変換します。
これは計算処理の半陰的E/H結合ではなく、分散FDTDランク内LD経路を置き換えるものでもありません。
初期の場・履歴はゼロです。

TDDFTはConfiguredMaterial内のTDDFTInitialStateを要求します。設定済みGSの橋渡し、または
`tddft_from_catalog`から取得してください。明示した管理計算プログラムが直列TDDFTに対応し、
実行ファイルのhashがsnapshotと一致することを、binding作成・prepare・既存backendのworker生成で
確認します。bindingとprepareはsnapshotを借用し、追加コピー、SCF、プロセス起動、通信による取得、
コンパイルはしません。状態生成では既存の独立子workerを使い、計算処理更新順、A/c入力、
電荷電流j=-J_SALMONを維持します。全セッションを閉じるまで準備モデルを開いたままにし、
その後ConfiguredMaterialを明示closeしてください。このproviderは借用モデルを閉じません。

両bindingは物質・任意初期状態の不変記録を保持します。providerの識別情報には、さらにLDの数値
パラメータ、またはTDDFTのsnapshot manifest/入力・実行ファイルhash・max_steps・threads・timeoutを
含めます。絶対計算プログラムパスやオブジェクトアドレスは同一性に使いません。
点/run記録のdigestは従来の物質・初期状態のものです。backend設定はprovider identity側なので、
完全な出典記録が必要なら設定も別途保存してください。

```python
bindings = {10: tddft_binding, 90: tddft_binding}
point_map = {10: 'positive', 90: 'negative'}
record = response_identity(bindings, point_map, run_id='experiment-1')
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=point_map)
```

集団呼出しと出力は[MPI_IDENTITY.md](MPI_IDENTITY.md)に従います。物質群は各1 MPIランクです。
群内並列Pythonから直接呼ぶTDDFTや同一プロセスの複数TDDFTを有効にする機能ではありません。
資源予算・物理収束・restartの機種間互換性は認定しません。LD直接E、SBE/自作の端点A経路は変更しません。

3ランクの例。NPZには正のscalar dtと、有限のA（形状(steps+1,3)、A[0]=0、原子単位）を保存します。

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_configured_native \
  inputs/materials.json metal-A waveform.npz /new/ld-output --mode ld-endpoints
mpiexec -n 3 python -m salmonpy.examples.mpi_configured_native \
  inputs/materials.json Si-ground waveform.npz /new/td-output \
  --mode tddft --runtime /local/runtime
```

TDDFTにはinitial_state_id、LDにはmaterial_idを渡します。ローカルのcatalog/GS、対応計算プログラム、
計算処理 LD入りwheelを事前準備してください。2点に逆符号の波形を与え、ID/電流を保存し、MPIセッションを
閉じてから独立した直接状態と比較します。TDDFT workerは各1スレッド、同時最大2個です。
メモリ上限の強制はありません。Si入門例は短い接続検査で、収束した光学計算ではありません。
