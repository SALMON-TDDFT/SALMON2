# 点IDを照合してTDDFTを再開

PointState.save_checkpoint / restore_checkpointは既存の管理TDDFT再開履歴にID検査を重ねます。
形式は`salmonpy.point-checkpoint/1`で、point.jsonと、既存TDDFT履歴archiveを入れたstate/です。
GS波動関数は含みません。呼出側で互換なローカル初期状態とruntimeを用意します。
archive内のパスからコードをimport・実行することはありません。

```python
from salmonpy.point_runs import PointRun

with PointRun({'p': configured_initial}, run_id='experiment-1') as run:
    state = run.new_state('p', dt=0.02)
    state.step([0., 0., 1.e-5])
    state.save_checkpoint(new_checkpoint_directory)

# 同じID・定義・初期データ・数値条件で明示的に再開する。
with PointRun({'p': rebuilt_configured_initial}, run_id='experiment-1') as run:
    state = run.new_state('p', dt=0.02)
    state.restore_checkpoint(new_checkpoint_directory)
```

外側のローダーはrun_id、point_id、material_id/definition_digest、任意initial_state_id/digestの
完全一致を要求します。別の点・runはbackendの復元前に拒否します。同じrun IDを使うのはそのrunの
再開を意味し、分岐ではありません。モデルオブジェクトは新しく作れますが、計算処理側の内容・runtime・
dt・容量・threadsの互換性検査は維持します。モデル生成や任意restart変換は行いません。

外側で内部manifestのhashを照合し、既存ローダーが配列を検証して、受理済みA/電流履歴を再計算します。
復元失敗時は既存TDDFTの原子的復元仕様に従います。履歴に比例して時間がかかり、復元用の別状態workerと
対象workerが共存する場合があるため、両方のメモリを確保してください。ID不一致では再計算を始めません。
保存先は上書きせず、失敗時は今回新しく作ったディレクトリだけを削除します。

symlink、余分な最上位項目、JSON重複キー、schema/version不一致を拒否します。
ID JSONは64 KiB、内部manifestは16 MiBまでです。max_bytes指定時は正の整数で、
内部checkpointのファイル/復号payloadの予算と、manifest事前hash読込の上限に使います。
総RSS上限ではありません。hashは破損検出であり真正性の保証ではないため、必要なら信頼できる別経路で
point.jsonのhashを固定してください。読込と同時にarchiveを変更しないでください。

既知の管理TDDFTState（ConfiguredStateで包んだ場合を含む）だけが対象です。
LD/SBEディスク形式、自作serializer、別runへの分岐、MPI全体の再開は未対応です。
既存のメモリ上snapshot/restoreと、直接TDDFTのcheckpoint APIは変更しません。
PointStateのディスクメソッドはこの外側形式を使い、point.jsonのない直接TDDFT checkpointは受け付けません。

```sh
python -m salmonpy.examples.point_tddft_restart \
  inputs/materials.json Si-ground /local/runtime waveform.npz /new/output
```

NPZは実数有限のA（steps+1,3、A[0]=0）と正のscalar dtです。途中保存、別点への復元拒否と非変更を
確認し、初期モデルを作り直して同じ点/runを再開、継続電流と完全一致を比較します。
各worker 1スレッド、同時最大2個です。短いSi入門では接続検査であり、光学収束・機種間互換性の証明ではありません。
