# 明示したバックエンド契約によるステップ保護

[English](STEPPING.en.md) · Topic: `protected_step` · Paired revision: `2026-11-17.1`

`salmonpy.stepping` は、ホストが明示した複数の状態を一つの操作として保護します。
`StepTarget(key, owner, capture, rollback, release=None)` でバックエンド側の完全な
退避・復元契約を指定します。トークンは各ランク内に保持し、SALMONpyはコピー・直列化・
MPI転送・内容解釈をしません。`capture()` は状態を変更しません。
`rollback(token)` は、数値状態だけでなく時計・キャッシュ・失敗フラグ・再試行に必要な
副作用も戻します。`release(token)` は退避資源だけを解放します。
一部ランクで退避に失敗しても、その対象へ入った全ランクがreleaseを呼び、失敗した
ランクはNoneを渡します。releaseはNoneを受け、途中の準備も安全に後片付けできる必要が
あります。正常なトークンとしてNoneを返すことも可能です。releaseを省略する場合、退避失敗時の
後片付けはcapture自身が行い、トークンの明示的な解放が不要であることが条件です。

```python
from salmonpy.stepping import StepTarget, StepRecoveryError, guarded_step

# 自作バックエンドが完全な巻戻しを保証するメソッドを明示します。
target = StepTarget("material", state, state.capture, state.rollback, state.release)
try:
    result = guarded_step(lambda: state.step(A), [target])
except StepRecoveryError as error:
    if not error.original_preserved or error.cleanup_errors:
        raise  # 安全を確認できない生存期間では再試行しません。
    # 受理済み状態に戻っています。原因を修正してから明示的に再試行します。
```

MPIでは全参加ランクが `context=parallel_context` を渡します。同じ物質群内では
対象のキー順とreleaseの有無を一致させます。群どうしの対象数は異なって構いません。
コールバックは自群の通信だけを使い、ホストのleaders通信へ入りません。操作・
コールバックが全ランクで通信を終えて戻ることが前提です。デッドロック、MPI／
プロセス停止、メモリ枯渇、native STOPからの復帰ではありません。宣言だけで
隠れた可変状態や外部I/Oの正しさが認証されることもありません。

事前検査後、退避を一対象ずつ全群で合意します。操作の途中で失敗したら、成功した
部品を含め全退避対象を逆順で戻し、トークンを解放します。トークンや計算結果は
合意通信に載せません。復元失敗ではoriginal_preservedがFalseになり、残りの復元と
後片付けも続けます。正常操作後のrelease失敗はphaseが`release`、original_preservedが
Falseです。この場合は既にステップが受理された可能性があり、安全な拒否と区別します。

## 場・物質をまとめて保護する

`CoupledStepper.step_protected(targets)` には、各ランクの生存中の場または各物質状態に
一つずつ対象を渡します。target.ownerは実際のオブジェクトにします。対象の不足・
重複・余分な所有者は、退避や場更新前に全ランクで拒否します。結合側の電流・時計・
応答A／失敗フラグ・端点アダプター・計測メタデータは結合器が追加で保護します。

部品の物理状態や副作用の保証は登録したコールバックに依存します。復元・後片付けが
不完全なら応答をfailedにし、場のあるランクでは場も無効化します。通常の`step()`は
従来動作のままです。自動fallback、刻み変更、自動再試行はしません。

導入済みwheelの例を6ランク、OMP／BLAS1で実行できます。

```sh
mpiexec -n 6 python -m salmonpy.examples.mpi_protected_step NEW_OUTPUT
```

Multiscale1D 1ランク、LD 1ランク、SBE 1ランク、自作モデル3ランク群（独立した2点）を
結合します。退避対象はA／電流だけでなく、ラッパーとLD／SBEの内部ハンドル状態全体です。
最初・途中・最後の失敗、場更新後の失敗、LDの失敗で元を保持し、再試行後の80受理
ステップが通常実行と完全一致します。過去の伝播を再計算しません。事前拒否・一部退避
失敗・復元失敗・解放失敗も検査します。合成模型による実装確認で、実材料の予測ではありません。

native E/H／TDDFT／高速PR版HSEの監査済み全体保護は、別の版付き契約を持つ
`native_steps.native_step_targets` で接続します。[NATIVE_STEPPING.md](NATIVE_STEPPING.md)
を参照してください。この軽量／自作例自体はnative全体の出力保護の検証ではなく、個別の
数値復帰・checkpoint能力から全体の原子性は推定しません。退避費用はバックエンドに依存し、
予算は警告のみです。import時のMPI初期化・nativeロード・ビルド・ダウンロードはありません。
