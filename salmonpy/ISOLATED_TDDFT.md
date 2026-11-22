# 孤立系TDDFT：E入力と双極子

PZ・固定イオン・非偏極・スカラー擬ポテンシャル・中点伝播を、明示的な別プロセスで使います。結晶のTDDFT/SBEは従来どおりA駆動です。返す量は電荷双極子であり、電流密度ではありません。格子・時間刻み・計算箱の収束は利用者が確認します。

## 計算プログラムと初期状態

```python
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState

runtime = load_isolated_runtime('salmon', 'validation/report.json')
with IsolatedTDDFTInitialState.from_salmon(
        'existing-bundle', 'rt.inp', runtime=runtime,
        max_steps=100, threads=1, timeout=60) as initial:
    with initial.new_state(dt=.02) as state:
        saved = state.savepoint()
        frame = state.step([0., 0., 1.e-5])
        print(frame.time, frame.dipole, frame.electron_number)
        state.restore_savepoint(saved)  # 同じワーカーへ直接コピー、時間発展なし
        state.release_savepoint(saved)
```

E・時刻は原子単位、双極子はe·bohrです。`dipole` はイオンと電子を合わせた全電荷双極子、`electronic_change` は電子双極子変化です。配列は読み取り専用です。初期E・時刻はゼロ。`new_state(position=...)` は場を評価する位置（既定(0,0,0) bohr）を指定し、分子座標の移動・回転は行いません。

読込みは64 KiB以下のローカル検証記録と、実行ファイルのSHA256を照合します。現行の記録形式は `salmonpy.isolated-worker-validation/4`。数値応答はE3の9成分のままですが、同一プロセス・時間発展追加0回の直接復元と、古い退避の拒否の実証を必須にしています。旧revision3の記録は受け付けません。読込み時にダウンロード・ビルド・検証計算・MPI初期化・ワーカー起動は行いません。記録はローカル検証の証拠で、電子署名や禁止ソルバー監査の代用ではありません。

既存の参照入力を使う `tools/test_isolated_worker.py EXECUTABLE EXISTING_REFERENCE NEW_VALIDATION` で、通常応答と直接復元を検証した記録を作成できます。これはRT計算を実行する開発用ツールです。新しいSCFは行いません。対応環境への配布が完了したという意味ではなく、手元の検証済み計算プログラムを明示して使います。

初期bundleには孤立系RT入力、`data_for_restart/info.bin`・`occupation.bin`・`wfn.bin`、相対パスで参照する擬ポテンシャル等を用意します。ローダーは不変の初期コピーを保持し、各状態が独立した作業コピーと子プロセスを持ちます。初期ローダーを閉じても既存状態は使えます。周期系GSの変換や孤立系GSの自動作成はしません。

## 同じ計算プロセス内での直接退避

`savepoint()` はFortran側の波動関数両バッファ、密度・ポテンシャル履歴、Poisson作業状態、双極子・外場履歴、エネルギー、更新番号を保持し、`IsolatedSavepoint` を返します。一つの状態につき一つだけ保持できます。`restore_savepoint(token)` は出力を元の長さへ戻して同じプロセスで復元し、時間発展を呼びません。トークンは復元後も残り、繰り返し復元できます。不要になったら `release_savepoint(token)` で明示的に解放します。解放だけなら現在の計算結果は変えません。

別状態のトークン・解放済みのトークン・再起動前のトークンは、復元前に拒否します。トークンはディスクへ保存する形式ではありません。退避を保持したまま履歴再生の `snapshot()/restore()` や `save_checkpoint()` は使えません。計算状態と退避の寿命は呼び出し元が管理し、状態を閉じれば退避も終了します。

`worker_statistics()` は明示的な問い合わせで、PID、実際の時間発展回数 `updates`、受理番号 `step`、直接復元回数 `direct_restores` を返します。通常の更新に全体通信を追加しません。復元だけでは `updates` は増えません。

直接退避中の失敗で、別ワーカーの起動や履歴再生は行いません。出力の巻き戻しに一時失敗した場合は `usable=False` として更新を止め、同じトークンで復元を完了できます。ワーカー喪失や復元結果を確定できない障害は使用不可として返します。`usable` を確認せずに更新を再開しません。任意のOS/MPI障害からの復帰を保証するものではありません。

退避の対象は現在検証している主要二出力です。未対応の空間出力を要求した場合は退避作成だけを拒否し、通常のE駆動計算は引き続き使えます。退避には電子状態一式を保持する分のメモリが必要です。

## 共通の場入力と観測

```python
from salmonpy.field_input import FieldPacket, FieldValue
frame = state.step_fields(FieldPacket(E=FieldValue(
    [0., 0., 1.e-5], time=state.time+state.dt,
    position=(0.,0.,0.), unit='E:au')))
```

この例は上記の生きている状態内で実行します。Eの時刻・位置・単位を照合し、未提供とゼロを区別します。A/Hは消費せず、AからEを暗黙に再構成しません。単独の `step(E)` は既存の定数区間Eを受ける積分順を維持します。

`frame.electron_number` はFortranが積分した実測電子数です。減少率1e-5は診断の目安であり、それだけで計算を停止しません。`frame.internal_energy` は既存の計算頻度で新しく評価された端点だけ値を返し、それ以外は `None` です。計算頻度を勝手に増やしたりゼロ外場ドリフトを自動で差し引いたりしません。

## 履歴の保存と再実行（直接復元とは別）

`snapshot()/restore()` は受理済みEと観測量の履歴を保存・再実行する入口です。同じ初期モデル、dt・位置・容量・スレッド数・バイナリが必要です。候補プロセスの観測量を照合してから交換し、復元失敗時は元の状態を保持します。費用は履歴長に依存し、検証中は新旧二ワーカー分の資源を使います。直接退避を使わない通常の `step()` の異常回復もこの履歴再実行方式です。

```python
state.save_checkpoint('new-history-checkpoint')
resumed.restore_checkpoint('new-history-checkpoint', max_bytes=16*1024*1024)
```

`resumed` は同じ初期データ・設定で呼び出し元が作成した状態です。ディスク形式は `salmonpy.isolated-tddft-checkpoint/3`。E0、初期／現在の9観測量、E履歴と観測履歴をハッシュ付きで保存するZarr v3互換形式です。GSは別途用意します。`max_bytes` は追加のファイル読み込み予算であり全RSS上限ではありません。これは波動関数の直接チェックポイントではありません。直接の耐久ディスク保存は残件です。

## 出力とFDTD

`python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT` は、E[steps,3]とスカラーdtを原子単位で持つNPZから双極子アーカイブを作ります。`--checkpoint` は明示した履歴再実行用の保存です。応答・単位・ソースの出典を明示し、結果はソースツリー外へ保存します。

双方向FDTD結合は [分子とFDTD](DIPOLE_RADIATION.md) の `MolecularFeedback` と、完全例 `examples/molecular_feedback.py` を使います。取消に `savepoint/restore_savepoint/release_savepoint` を使い、電子履歴を再計算しません。真空自己場の差し引きは真空の放射反作用も除きます。

既存C2H2 GS・PZ、64場ステップの6ランク検証では、正常時は電子更新63回、1区間の取消・再試行時は64回で、各実行内のPIDは変わりません。復元による追加更新は0回、受理した全保存配列・六面記録は以前の結果と完全一致しました。これは短時間の接続検証で、光学物性の収束保証ではありません。
