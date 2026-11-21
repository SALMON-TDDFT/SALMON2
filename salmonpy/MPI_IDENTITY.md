# MPI応答セッションで利用者のIDを継承

直列のpoint-run記録を既存の二段階MPI応答サービスへつなぐ、明示選択式の機能です。
整数の通信番号、ランク割当、A/c→電荷電流の更新・加算、場のソルバーは変更しません。
計算処理 backendの能力・制約は従来どおりです。import時にMPIを初期化しません。

全contextランクで同じ`point_record`と`point_map`を
`ResponseSession(..., dt=..., point_record=record, point_map=mapping)`または
`CoupledStepper.from_multiscale(..., point_record=record, point_map=mapping)`へ渡します。
記録は`salmonpy.point-run/1`（[POINT_RUNS.md](POINT_RUNS.md)参照）です。
plain dataとして作るか、直列状態を作らずに`PointRun.record()`から得られます。
run_idは呼び出し元プログラム側で決め、セッション作成前に共有します。各ランクで別々にUUIDを生成しないでください。

```python
record = {'schema': 'salmonpy.point-run/1', 'run_id': 'experiment-1', 'points': [
    {'point_id': 'right', **right_material.record()},
    {'point_id': 'left', **left_material.record()},
]}
# 既存plan.point_ids == (10, 90)。リストと対応表の順序を明示して一致させる。
mapping = {10: 'right', 90: 'left'}
session = ResponseSession(context, bindings, dt=0.02,
                          point_record=record, point_map=mapping)
```

整数の全点を一度ずつ対応付け、記録の順序は`plan.point_ids`をこの対応表で変換した順と一致させます。
利用者のIDをソート・正規化したりランク名・ファイル名に変換したりしません。
物質定義・任意の初期状態hashの内部整合性と、全ランク間の記述子・記録の一致を
prepare/createより前に検査します。ただし宣言されたIDであり、providerが実際にそのデータを
読み込んだ証明ではありません。アダプター側の実体hash・互換性検査も必要です。
MaterialBindingのキャッシュkeyはmaterial_idとは別です。同じbinding/keyで異なる物質・初期状態を
参照することは拒否します。

provider/backendへ渡る`SolverContext.point_records`は、自群内の整数点番号→不変な記録です。
各記録にrun_id、point_id、material_id/definition_digest、任意のinitial_state_id/digestが入ります。
既存コールバックのpoint_id引数は整数のままです。metadataを省略するとpoint_recordsと
`session.identity_record()`はNoneで、従来の利用を維持します。指定時のidentity_recordは
point_runとmpi_point_map（整数点、利用者点、群名、worldランク）の独立コピーを返します。
集団Fortran/Cで保持する計算状態を新たな状態ラッパーで包むことはしません。

## 集団での出力

全contextランクで`session.save_arrays`を呼びます。引数を渡すのは場のleaderだけで、他は引数なしです。
事前検査・書込エラーを全ランクへ通知し、場leaderだけにOutputArchiveを返します。

```python
if context.group.role == 'field' and context.is_leader:
    archive = session.save_arrays(new_directory,
        {'current': current_history, 'time': time}, point_axes={'current': 1},
        metadata={'units': 'atomic', 'current': 'charge j'})
else:
    session.save_arrays()
```

配列はplan.point_ids順です。時刻は推定しないため明示保存してください。
point_run、point_axes、mpi_point_mapは予約キーです。`open_point_arrays`で配列hash・点軸・
利用者ID・通信点と群所有の整合性を検査できます。MPI状態の復元や計算再開は行いません。
全ランクで同じ順序の呼出しが必要です。プロセス死亡や計算処理 MPI abortは回復不能なので、
ランチャー・ジョブに制限時間を設けます。通常の出力失敗では数値状態を進めたり無効にしたりしません。
新しい出力先で再試行できます。既存の出力先は上書きしません。
metadataの真正性には、別経路で信頼できるmanifestのhashが必要です。

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_point_identity /new/output
```

4ランク（片方の物質群が2ランク）でも動き、直接式・従来セッションと比較して保存IDを読み戻します。
合成Python応答であり、TDDFTの物理検証や新しいFDTDソルバーではありません。
各ランク1スレッドとし、MPI/mpi4pyは事前導入してください。暗黙の取得・コンパイルは行いません。
