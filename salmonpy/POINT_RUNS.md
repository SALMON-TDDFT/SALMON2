# 直列計算の点IDと結果配列

`PointRun`に順序付きの点ID→準備済みモデルの対応を渡します。点ごとに一つの状態の寿命を管理し、
指定したrun ID（省略時はUUID）と物質・任意の初期状態ID/hashを保持します。
IDは大文字小文字を区別する不透明なUnicode文字列で、ファイル名やランク番号に変換しません。

```python
from salmonpy.point_runs import PointRun, open_point_arrays

# modelはConfiguredMaterial。tddft_from_catalogの戻り値も使える。
with PointRun({'right': model, 'left': model}, run_id='experiment-1') as run:
    right = run.new_state('right', dt=0.02)
    left = run.new_state('left', dt=0.02)
    right.step([0., 0., 1.e-5])
    # rightだけ進む。LDの入力はE、SBE/TDDFTはA/c。
    identity = right.record()
    # 配列の点順はrun.point_ids == ('right', 'left')に合わせる。
    run.save_arrays(new_output, {'current': [[right.current, left.current]]},
                    point_axes={'current': 1},
                    metadata={'units': 'atomic units', 'note': 'points at different times'})
archive = open_point_arrays(new_output)
point_ids = [p['point_id'] for p in archive.metadata['point_run']['points']]
```

更新順の制御、時刻の一致、駆動の変換、単位推定、run全体の巻き戻しは行いません。
物理解析には時刻配列も明示保存してください。一つのrunで各点を作れるのは一度で、close後も
再作成できません。new_state自体が失敗した場合は再試行できます。分岐・反復計算には別runを使います。
run.closeは全状態のcloseを試み、最初の失敗を再送出します。失敗したcloseを自動再試行しません。
準備済みモデルは呼出側所有で、全状態を解放するまで維持してください。
new_stateは初期化済み状態を返す契約で、状態のcontext-entry hookは呼びません。
資源を持つ自作状態はclose()を実装してください。contextメソッドだけで資源管理する状態は未対応です。
資源を持たない軽量な状態ではclose()を省略できます。
状態のcontextはその点だけを解放し、run.closeは二度目以降何もしません。
この層には同時呼出し、MPI割当、checkpoint、restart対応表はありません。
委譲されるsnapshot/restoreも既存backendの制約を維持し、native checkpointへIDは追加しません。

## 自作モデル

継承は不要です。`new_state(**kwargs)`と`record()`を実装します。recordは厳密に
`material_id`と`definition_digest`（小文字SHA256）、任意で`initial_state_id`と
`initial_state_digest`の両方を返します。digestには検証済み採用設定と必要な外部データhashを使い、
オブジェクトのメモリアドレスは使いません。状態生成前にもidentityを照合します。
他の出典情報は出力metadataへ保存します。同じ物質IDの定義digest、同じ初期状態IDの物質と内容は
一致する必要があります。複数点に同一の状態オブジェクトを返す場合は拒否しますが、隠れた共有配列や
nativeグローバル変数までは検出できません。独立性は自作backend側でも保証してください。
identityから物理的妥当性・実行能力は認定しません。JSON設定の解決時は引き続きregistryの検証が必要です。
このPythonのbindingは任意JSONを読み込む仕組みではありません。

## 保存と読込

`run.save_arrays`は既存`salmonpy.outputs/1`のmetadataに`point_run`として順序付きの
`salmonpy.point-run/1`記録を追加します。`point_axes`で配列名→非負の点軸番号を明示します。
少なくとも一つの配列で、点軸の長さが`len(run.point_ids)`と一致する必要があります。
他の配列は解釈しません。点順は呼出側で合わせてください。長さの検査だけでは点の入れ替えを検出できません。
二つのmetadataキーは予約済みです。`open_point_arrays`は配列hash、形状と対応表、IDと定義の
一貫性を検査します。自動並べ替えや計算再開はしません。metadataの意味を検査しても真正性は保証しないため、
信頼境界を越える場合はmanifestのhashを別途固定してください。通常のOutputArchiveでも読み取れますが、
追加の点schema検証は行いません。

```sh
python -m salmonpy.examples.point_ld /new/output/directory
```

導入済みLD wheelが必要です。2点の独立性と直接呼出しとの一致を確認し、電流/E/時刻を保存して
対応表を読み戻します。合成材料であり、実験で検証したパラメータではありません。MPI・通信は不要です。

MPI応答サービスでの明示的なID継承・出力は[MPI_IDENTITY.md](MPI_IDENTITY.md)を参照。
PointRun自体の直列専用という制約やrestart未対応は維持します。

管理TDDFTの点ID付きディスク復元を別途追加しました。[POINT_CHECKPOINTS.md](POINT_CHECKPOINTS.md)参照。
上記のrestart未対応はMPI全体・汎用backendの制約です。

## 点ID付きのメモリsnapshot

`PointState.snapshot()` は内部バックエンドのsnapshotを `salmonpy.point-snapshot/1` の不変な点情報で包みます。`state.restore(saved)` は実験・点・物質定義・初期状態のIDとパッケージ版を照合してからバックエンドを呼びます。別の点や実験、裸のバックエンドtokenは拒否されます。明示した同一run/point IDであっても、バックエンドのモデル・時間刻み等の整合性検査は引き続き必要です。

これはプロセス内の巻き戻し用であり、分岐・複製やディスク保存、MPI全体のトランザクションではありません。snapshot未対応の自作バックエンドでは `SnapshotError` となります。バックエンドの復元失敗時の保証はその実装に従います。点ID不一致の拒否では内部のrestoreを呼びません。裸の状態に対する従来のsnapshot APIは変わりません。

`python -m salmonpy.examples.point_ld /new/output/directory` は別の点への復元拒否と、同一点での復元後の電流一致も確認します。
