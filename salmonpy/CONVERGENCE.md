# 格子幅・時間刻みの診断

`salmonpy.convergence` は計画と実行を分けます。計画だけならNumPyで作成でき、
計算・コンパイル・通信は始まりません。実行時の有限値検査、電子数の保存誤差、
応答の比較はwheelに同梱するFortranカーネルが担当します。Pythonは試験の
呼び出しと結果の記録を担当します。古いwheelに必要ABIがなければ、明示的な更新を案内します。

まず少数のk点を固定し、物理的な擬ポテンシャル半径の1/3から格子幅を小さくします。
既定の調査下限は0.1 Åですが、さらに細かい指定も可能です。半径はspeciesごとに
値・定義・PPのSHA256を明示します。半径の自動抽出はまだ行いません。
動径表の末端や刻みをcutoff半径として代用しないでください。
既存の幅→ceiling格子変換を使い、同じ実効格子の重複は除きます。

```python
from salmonpy.convergence import plan_grid_scan, plan_dt_scan
spatial = plan_grid_scan(
    [[5.43,0,0],[0,5.43,0],[0,0,5.43]], material_id="silicon",
    radii={"Si": {"radius_angstrom": 1.2,
                  "definition": "利用者が選択したPPチャネルのcutoff半径",
                  "pp_sha256": pp.sha256}},
    metadata={"xc": "pz", "kpoints": [1,1,1], "ranks": 4})
temporal = plan_dt_scan(material_id="silicon", duration_fs=3,
                        impulse=[1e-4,0,0], expected_electrons=8)
```

`run_grid_scan(plan, evaluate, atol=..., rtol=...)` のコールバックは
`{"value": energy_per_atom, "scf_converged": bool, "metadata": {...}}` を返します。
単位と観測量、セル・原子・PP・XC・バンド・SCF条件は全試験で揃えます。
SCF未収束は警告して継続し、判定は暫定と明記します。隣接する細かい格子との
比較で目安を満たさなければ、選択値はNoneです。k点収束の証明にはなりません。

格子とGSを固定した後、dtを0.2 auから半減します。物理時間は最初のdtの
整数ステップへ切り上げ、全試験で同じ終了時刻を使います。結晶では初期A/c=0から
一定のA/cへ跳ぶ弱いimpulseを使い、積分したEは `-impulse` です。
x・y・z・任意の偏光ベクトルを指定できます。

```python
from salmonpy.convergence import run_dt_scan
report = run_dt_scan(
    temporal, lambda dt: initial.new_state(dt=dt),
    electron_number=backend.measured_electron_number,
    rtol=1e-5, atol=0.0)
```

上は接続契約の例です。`measured_electron_number` は利用者バックエンドの関数名で、
SALMONpyの既存TDDFT状態に自動追加されたメソッドではありません。
コールバックはFortran/C側で軌道の占有・k点重み・格子体積を含めて測った電子数を
返す必要があります。設定のnelecを返すだけでは誤差を検出できません。
初期状態のcurrentと各stepの戻り値は長さ3の電荷電流とし、各試験を独立に作成・closeします。
既存native TDDFTへの電子数測定口の接続と実材料での検証は次段階です。

主指標は初期電子数からの最大減少率で、基本の相対目安は **1/100000** です。
増加率と絶対ずれも併記し、許容値内判定には両側の最大ずれを使います。
`expected_electrons` を指定すると、それを基準にして初期値のずれも検出します。
目安を超えても試験は止めません。NaN／Inf・非正のノルム・数値例外だけ、その試験を
直ちに終え、細かいdtへ進みます。入力・実装の契約違反を数値破綻と偽装しません。
外部処理・I/O失敗はruntime_failure、close失敗はcleanup_errorとして分けます。

- `nonbreaking_dt_au`: 試験時間を破綻せず完走した最大dt。
- `norm_acceptable_dt_au`: 電子数の保存誤差が目安内だった最大dt。該当なしはNone。
- `current_converged_dt_au`: オプションの電流比較で隣接する細かいdtと一致した最大dt。

ノルム保存だけで位相・電流・光学応答の精度は保証されません。
`current_atol` と `current_rtol` を両方指定すれば、共通の物理時刻で電流も比較します。
最小dtを自動的に収束済みとは扱いません。最後にk点収束と本計算条件を別途確認します。
MD、HSE等で適用可能な診断かは各バックエンドの測定定義も確認してください。
