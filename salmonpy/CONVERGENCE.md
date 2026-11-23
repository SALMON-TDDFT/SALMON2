# 格子幅・時間刻みの診断

`salmonpy.convergence` は計画と実行を分けます。計画だけならNumPyで作成でき、
計算・コンパイル・通信は始まりません。実行時の有限値検査、電子数の保存誤差、
応答の比較はwheelに同梱するFortran計算ルーチンが担当します。Pythonは試験の
呼び出しと結果の記録を担当します。古いwheelに必要ABIがなければ、明示的な更新を案内します。

まず少数のk点を固定し、物理的な擬ポテンシャル半径の1/3から格子幅を小さくします。
既定の調査下限は0.1 Åですが、そこへ到達する必要はありません。SCF収束済みの
隣接2格子で指定した誤差内なら、その時点で終了します。さらに細かい指定も可能です。半径はspeciesごとに
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

格子の主判定は `run_current_grid_scan(plan, evaluate)` を使います。
各格子で、まず電子数保存1e-5を満たすdtを `run_dt_scan` で求め、そのdtの電流を
隣接する細かい格子の電流と比較します。許容する相対誤差は既定で **5%以下** です。

\[
\epsilon_J=\frac{\sqrt{\sum_{t,\alpha}(j^{\rm coarse}_{t,\alpha}-j^{\rm fine}_{t,\alpha})^2}}
{\sqrt{\sum_{t,\alpha}(j^{\rm fine}_{t,\alpha})^2}}\le0.05
\]

全時刻・全3成分を使う相対RMS（相対L2）誤差です。各時刻の電流で割らないため、
ゼロ交差で相対誤差が発散しません。細かい格子を基準にし、同じ物理時間のサンプルで
比較します。dtが異なる場合は入れ子の時刻列を間引いて揃え、補間はしません。
電流がゼロの基準波形では相対判定不能です。絶対許容値atolを明示した場合だけ
絶対RMS誤差で判定できます。

```python
from salmonpy.convergence import run_current_grid_scan
# evaluate_gridは各格子のGS/RTをFortran/C計算処理で計算する関数です。
# {"scf_converged": bool, "dt_report": run_dt_scanの結果} を返します。
report = run_current_grid_scan(spatial, evaluate_grid, rtol=0.05, atol=0.0)
```

セル・原子・PP・XC・k点・励起・物理時間を揃えます。dtが未解決の格子は
比較対象に使いません。SCF未収束は警告して継続し、判定は暫定と明記します。
両格子のSCFが収束し電流誤差5%以下なら、0.1Åより粗い段階でも終了します。
`stop_when_converged=False` なら予定した全格子を調べます。
一致する格子対がなければ選択値はNoneです。k点収束の証明にはなりません。
`run_grid_scan` によるエネルギー等のスカラー比較は補助診断として残しています。

格子とGSを固定した後、dtを0.2 auから半減します。物理時間は最初のdtの
整数ステップへ切り上げ、全試験で同じ終了時刻を使います。結晶では初期A/c=0から
一定のA/cへ跳ぶ弱いimpulseを使い、積分したEは `-impulse` です。
x・y・z・任意の偏光ベクトルを指定できます。

```python
from salmonpy.convergence import run_dt_scan
report = run_dt_scan(
    temporal, lambda dt: initial.new_state(dt=dt),
    rtol=1e-5, atol=0.0)
```

新しく明示ビルドした計算プログラムでは、プロセス分離 `TDDFTState` と通常XCの
計算処理／共有モデル／異種モデル登録の状態が `electron_number` を持ちます。
Fortranの密度積分（占有・k点重み・スピンを含む電子密度×格子体積）を返します。
入力のnelecで代用しません。`run_dt_scan` は省略時にこのプロパティを読みます。
自作Fortran/C計算処理は `electron_number=callback` でも接続できます。
計算処理の取得は全内部MPIランクが同じ順序で呼ぶ集団操作です。状態ハンドルごとに
取得し、時刻・電流・軌道は更新しません。保存／再開後も復元した密度から測定します。
対応範囲は結晶A駆動のPZ／TB-mBJ／Libxc-PZ/PBEです。HSEと孤立系E駆動の
測定口は未接続です。既存計算プログラムや旧プロトコルは通常計算に使えますが、
実測値が取れなければ診断は具体的なエラーを返します。新wheelだけでは旧本体は更新されません。

```sh
mpiexec -n 2 python -m salmonpy.examples.tddft_convergence BUILD_JSON NEW_OUTPUT --grid 16
```

この例はGSをメモリ内で保持し、GS再開データの保存を明示的に無効にします。
各dtは同じGSから独立した物質モデルを作ります。SCF未収束なら警告して続け、
結果はprovisionalです。小さな低k格子での実行例は、光学物性の収束証明ではありません。
既存の電子数が10倍以上ずれた際の本体停止判定は、回復対応APIでは状態を戻して
その試験の数値エラーとして返します。1e-5の目安超過とは別の既存判定です。

主指標は初期電子数からの最大減少率で、基本の相対目安は **1/100000** です。
増加率と絶対ずれも併記し、許容値内判定には両側の最大ずれを使います。
`expected_electrons` を指定すると、それを基準にして初期値のずれも検出します。
目安を超えても試験は止めません。NaN／Inf・非正のノルム・数値例外だけ、その試験を
直ちに終え、細かいdtへ進みます。入力・実装の仕様違反を数値破綻と偽装しません。
外部処理・I/O失敗はruntime_failure、close失敗はcleanup_errorとして分けます。

- `nonbreaking_dt_au`: 試験時間を破綻せず完走した最大dt。
- `norm_acceptable_dt_au`: 電子数の保存誤差が目安内だった最大dt。該当なしはNone。
- `current_converged_dt_au`: オプションの電流比較で隣接する細かいdtと一致した最大dt。

ノルム保存だけで位相・電流・光学応答の精度は保証されません。
`current_atol` と `current_rtol` を両方指定すれば、共通の物理時刻で電流も比較します。
最小dtを自動的に収束済みとは扱いません。最後にk点収束と本計算条件を別途確認します。
MD、HSE等で適用可能な診断かは各計算処理の測定定義も確認してください。

電流の相対RMS5%という既定値は格子のスクリーニング基準であり、誘電関数や高次高調波スペクトルの精度保証ではありません。
用途に応じて`rtol`／`atol`を厳しくし、目的のスペクトル自体の収束も確認してください。採用した許容値は結果に記録します。既定の数値基準は変更しません。
