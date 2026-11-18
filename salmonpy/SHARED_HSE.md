# Pythonから共有する高速PR版HSEモデル

[English](SHARED_HSE.en.md) · Topic: `shared_fast_hse` · Revision: `2026-11-18.1`

`NativeFastHSEBackend.open_model(initial, dt, context)` が `NativeFastHSEModel` を返します。
不変の全k点入力を一度だけ検証・コピーし、ホストのMPIコミュニケータに一度接続します。
実際のnative物質初期化は最初の `new_state` です。以後は元のGSから独立した状態を作り、
他状態の励起履歴を引き継ぎません。隠れた電子状態テンプレートや、状態ごとのGSディレクトリコピーはありません。

```python
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend

initial = FastHSEMaterialProvider(gs_directory, max_steps=100).prepare(context)
model = NativeFastHSEBackend(build_json, work_root).open_model(initial, .001, context)
states = []
try:
    with model:
        a = model.new_state('point-a')
        states.append(a)
        b = model.new_state('point-b', dt=.002, capacity=100)
        states.append(b)
        ja = a.step([0., 0., 1e-4])
        jb = b.step([1e-4, 0., 0.])
        # MPIランクごとに異なる共有ファイル名を使う。
        a.checkpoint(rank_checkpoint_path, False)
        restored = model.new_state('point-a')
        states.append(restored)
        restored.checkpoint(rank_checkpoint_path, True)
        j_next = restored.step([0., 0., 2e-4])
finally:
    for state in reversed(states):
        state.close()
    model.close()
```

`context` は `solver_comm` と `point_ids` を持ちます。呼び出し元のParallelContextも利用できます。
全ランクが同じ順序で、メインスレッドから呼び出してください。全ての状態・モデルのleaseを閉じるまで
コミュニケータを維持します。入力・ライブラリの検査はnative寿命消費前です。
importだけではMPI・nativeライブラリの読み込み、コンパイル、ダウンロードを行いません。

`new_state(point_id=None, *, dt=None, capacity=None)` は、既定でモデルのdtと準備済み最大ステップ数を使います。
point_idを省くと連番整数、指定した文字列・整数IDは `state.point_id` に保持します。
nativeハンドルはこれと別です。各状態はA・電流・時計・容量・HSE履歴・出力basenameを所有し、
`state.model_identity` が共通物質スナップショットを識別します。dtやcapacityの変更は兄弟状態や
モデルの既定値を変更しません。元の入力は再利用でき、検証済み作業コピーは変更しないでください。

状態の `step(A)`、`advance(Aold,Anew)`、`checkpoint(path,False)`、
`checkpoint(path,True)`、`close()` は従来と同じ意味です。
直接読込は、dt・容量・バイナリ・GS／入力・MPI分割が対応する新規状態へ行います。
過去ステップの再計算なしに内部履歴を復元します。保存ファイルと同一性の照合はホスト側で行います。
ログは新しい状態のラベルで開始し、数値復元だけでは過去ログを統合しません。
通常の保護付きステップは、その状態の2本の出力だけを退避するため、共通ディレクトリでも
兄弟状態のログを巻き戻しません。

`model.close()` は新規生成を止め、既存状態は継続できます。最後の状態を閉じた時にMPI接続を解放します。
モデルを開いたまま全状態を閉じても、隠れた完全電子状態を保持せずに再生成できます。
モデルのコンテキストマネージャーは生存する状態を閉じません。各状態を `finally` で閉じてください。
native初期化を一度も試していないモデルの終了は、未使用のPython予約も解放します。
PythonのGCからMPI終了処理は呼びません。

対象は、プロセスごとに一つの不変な**対称性なし・全k点**モデル、既存backendのk点MPI1／2、
CPU・メインスレッドの範囲です。従来の `create_state` は1寿命のままで、同じプロセスでは
こちらか `open_model` のどちらかを選びます。handle protocol1とmodel protocol1を備えた
ランタイムを明示的にビルドしてください。縮約対称性の複数状態、異なる物質入力の同一プロセス実行、
複数点の自動ランク再配置、スレッド実行は別作業です。子コミュニケータは最後のdetachまで接続単位で保持し、
無制限の生成・終了を保証しません。資源予算と物理収束の助言は警告扱いで、入力を自動変更しません。
AはCartesian A/cの原子単位、jは電荷電流密度（SALMON物質電流の負）、dtは原子時間、energyはHartreeです。
x／y／z・円偏光を指定できます。

対応MPIで [examples/shared_fast_hse.py](examples/shared_fast_hse.py) を明示実行できます。
全k点・内部パルスなしの準備済みGSを指定します。

```sh
mpiexec -n 2 python shared_fast_hse.py BUILD_JSON GS NEW_OUTPUT
```

2種類の偏光を進め、新規状態への直接復元を確認する例です。実材料の物理的収束を保証する例ではありません。
下位C ABIと残る共有資源については[nativeハンドル](NATIVE_HANDLES.md)を参照してください。
