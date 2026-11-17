# メモリ見積もりと助言予算

メモリ予算は警告用です。超過しても計算・保存／復元を停止せず、OSやジョブスケジューラの制限も変更しません。
見積もり以下でも実行可能という保証ではありません。4GiBは4×2³⁰ bytesです。

## Pythonから使う

```python
from salmonpy.memory import MemoryEstimate, process_memory, memory_advice

live, basis = process_memory()
# candidate_bytesは自分のバックエンドの追加バッファ量。512MiBは記述例です。
estimate = MemoryEstimate(live_bytes=live, candidate_bytes=512 * 2**20,
                          reserve_bytes=128 * 2**20, basis=basis)
report = memory_advice(estimate, budget_bytes=4 * 2**30)
print(report.to_dict())
```

自作バックエンドも同じ入口へ、実測や配列サイズから求めた量とbasis（根拠）を渡せます。
取得できない値はNone、未計上の部分はunknownで示します。ゼロへ置き換えません。
known_peak_bytesは分かる部分の合計、estimated_peak_bytesは不明部分があればNoneです。
exceeds_budgetのFalseは「与えた見積もりが予算以下」であり、安全認証ではありません。
Trueなら超過量excess_bytesも返します。不明で判定できなければNoneです。

MPIでは既存のcommをmemory_advice(..., comm=comm)へ渡し、全ランクで呼びます。
予算はそのcomm全体の合計で、ランク毎ではありません。ノード別に見る場合は呼出側でノード別commを用意します。
予算がランク間で違う場合は不明として報告します。警告はcommのrank0だけです。
バックエンドの見積もりが失敗したら、各ランクで例外を処理してestimate=None等を渡し、集合呼出を省略しないでください。

process_memoryはLinux/macOSでプロセスの最大RSSを読むだけで、import時の計測・MPI初期化・コンパイルはしません。
最大RSSは現在使用量ではなく、子プロセスを含みません。共有ページをランク毎に数えるため合計は重複し得ます。
未対応環境や計測失敗はNoneと理由を返します。Windows等では利用者が見積もりを与えられますが、RSS取得を推測しません。
警告が例外化される設定でも警告表示失敗をwarning_errorへ分け、予算超過を計算失敗にしません。
不正なAPI引数はInputErrorです。

## 通常の混在サンプル

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 \
  python -m salmonpy.examples.mpi_mixed_native BUILD_JSON RUNTIME GS_BUNDLE NEW_OUTPUT \
  --memory-budget-gib 4 --rank-budget 4
```

4ランク×1スレッドで4コアを使います。予算を指定するとNEW_OUTPUT/memory-advice.jsonに次を残します。

- 初期化前：格子／LD極数、SBE行列、入力波動関数サイズからの概算。不明な初期化ワーク領域も明示。
- 初期化後と計算区間の終了時：各プロセスの最大RSSとランク毎128MiBの余裕分。
- 保護付き復元前：最大RSS＋native保存データの4倍＋JSONの32倍＋snapshot候補の追加領域＋余裕分。

復元候補と読込用コピーを同時に保持する分を含みます。倍率は内部配列の全てを証明した上限ではなく、
同じ状態を複数回数える場合もある保守的な概算です。HSE等の子ワーカー使用量は不明として残します。
各ランクの内訳・根拠と不明点を読んでください。予算は物理入力ではなく、再開時に変更できます。
資源計測・警告表示・助言ファイル出力の失敗は診断へ分け、計算のfailedフラグを立てません。
実際の入力破損・数値エラー・MPI障害とは区別します。未指定ならサンプルの助言処理は無効です。
この自動見積もりは通常混在サンプルのCPU経路向けで、全DFT/ASE/外部バックエンドを自動で推定するものではありません。
