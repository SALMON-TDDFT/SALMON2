# 二階層MPI応答の基盤（実験的API）

[English](DISTRIBUTED.en.md) · Topic: `distributed` · Paired revision: `2026-11-13.4`

Pythonが起動済みMPIジョブを計算群に分けます。各物質群は一つ以上のマクロ点IDを担当し、各点が独立した状態を持ちます。
場を担当する群がA/cを配り、点ID順に物理的な電荷電流密度Jを回収します。単位は原子単位、時刻は固定dtの端点、初期Aは0です。
このサービス単独では所定のAへの応答を計算します。**内部MPIのTDDFTとE/H結合**には追加の[ネイティブ結合](COUPLED.md)を使います。

同じMPI環境に任意依存の `salmonpy[mpi]` を導入します。通常のimportではMPIの読み込み、コンパイル、ランク起動を行いません。
mpiexecまたはスケジューラでPythonを一度起動し、全ランクが同じ集団呼び出しを実行します。

```sh
mpiexec -n 5 python -m salmonpy.examples.mpi_response NEW_JSON
```

例は、場役1ランクと、各2ランクの物質群二つです。各物質群に一つのマクロ点（ID 10と20）があり、異なる係数の材料を使います。
利用者が書いた履歴付きの人工的な応答則を使い、実材料やFDTDを計算する例ではありません。

## 利用者コードの交換口

`MaterialBinding(key, provider, backend)` で物質準備と応答を分離します。同一材料キーの点では同じbindingインスタンスを再利用します。
providerの `identity` は版と**設定内容**を識別し、全ランクで一致させます。`prepare(context)` は読み取り専用の初期データを返します。
独自形式を読み込め、SALMONの内部型に合わせる必要はありません。backendのidentityにも同じ条件があります。

- `execution`：全群メンバーが実行する `collective`、または現在1ランク群のみで使う `leader_only`。
- `reduction`：各ランクの寄与を一度だけ加算する `sum`、または群リーダーが完成した電流を持つ `leader`。後者でも他ランクは有限な3成分を返します。
- `create_state(initial, point_id, dt, context)`：毎回独立した状態を作成。その状態は `advance(A_old,A_new)` と `close()` を持ちます。
  advanceの戻り値は有限な3成分の**物理的な電荷電流密度**です。電子数電流や全電流と混同しないでください。

コールバックに渡す `SolverContext` は群のコミュニケータ、群名、担当点IDだけを含みます。全メンバーが同じ順序で点を処理します。
群内計算にMPI_COMM_WORLDを使ったり、MPIを初期化・終了したり、共有初期データを変更したり、一つの可変状態を複数点へ返したりしないでください。
必要なら状態closeも群内で集団実行します。prepare/createが失敗するときは、コールバック自身が構築途中の資源を群内で解放してから例外を返します。
サービス側はそれ以前に作成を完了した状態を解放します。

`ParallelPlan` はランク・点の重複を禁止し、場役を一群に限定します。`ParallelContext.create(world,plan)` は分割前に全ランクの計画一致を確認します。
`ResponseSession(context,bindings,dt=...)` が状態を準備します。全ランクからadvanceを呼び、場群のリーダーだけが `(点数,3)` のA配列を渡し、他はNoneを渡します。
場群の全ランクには `plan.point_ids` 順のJ、物質群にはNoneが返ります。使用中のprovider/backendの設定やidentityは変更しません。
コンテキスト・コミュニケータ・状態のスレッド共有は未対応です。

不正な入力Aは伝播前に拒否します。応答計算中の例外では全群のセッションを使用不可にし、全体の巻き戻しは保証しません。
closeして既知の初期状態から再開してください。ネイティブプロセスの故障や利用者コード内部のMPI停止はPython例外だけでは回復できず、
起動側のタイムアウト・ジョブ停止が必要です。セッションとFortran接続を閉じてからParallelContextを閉じます。MPI自体の終了はPythonジョブ側が担当します。

`LDMaterialProvider(model)` と `LDResponseBackend()` は既存LD状態を使い、E_new=−(A_new−A_old)/dtをLDState.stepへ渡して物理Jを返します。
方式名は `salmonpy.ld.a-endpoints/1`。**既存FDTDの半陰的なLD係数・Q結合ではなく、その更新処理には置き換えられません**。
TDDFT/SBEの電子数電流も各アダプターで符号変換が必要です。このサービスは利用者のJの意味を推測して符号を変えません。

## Fortranへの群の受け渡し

mpi4pyと同じMPIを使い、新しく準備した外部ソースから明示的にビルドします。現在はGNU互換mpifort、macOS/Linuxを対象とします。

```sh
python tools/prepare_salmon_runtime.py UPSTREAM_CHECKOUT NEW_SOURCE
python tools/build_parallel_bridge.py NEW_SOURCE NEW_BRIDGE --fc mpifort
```

wheelにも `salmonpy/_runtime_source/tools` 配下にツールを含めます。場所はimportlib.resourcesで確認できます。
ソース・オブジェクト監査、build.json、ハッシュ、ログを残します。MPIブリッジは実行環境用に作り、直列ランタイムへの同梱や自動コンパイルは行いません。
Python/FortranのMPI ABI互換性が必要です。実装名と初期化状態の照合は追加検査であり、あらゆるABI互換性の保証ではありません。

生成されたライブラリとbuild.jsonのSHA256を `NativeParallelContext(context,library,sha256)` に渡します。
例の `--bridge PATH --sha256 HASH` でも試せます。mpi4pyのプロセス内Fortranハンドルを受け、整数表現を確認し、SALMON通信層が群を複製します。
既定の通信・バリアはその群を使い、切り離しは複製だけを解放してMPIを終了しません。info/sum/barrierは接続検証用で、ソルバー入口ではありません。
現在は各プロセス・ライブラリに一つの接続です。ネイティブ試験では接続中の既存comm_init/comm_finalizeも確認します。

一群一点の組み込みTDDFTと限定したE/H↔A/J結合は[COUPLED.md](COUPLED.md)に記載しています。残る作業は、組み込みDFT、同一プロセス内の複数電子状態、層からマクロ格子への自動割り当て、再開と性能検証です。
従来のバッチFDTDとプロセス分離型TDDFTは、これまでの対応範囲で利用できます。

collective型のcreate_stateは、群の全ランクで成功するか、全ランクで途中の資源を解放して例外を返す契約です。一部ランクだけの成功は対象外です。コンテキストのcloseは全群で集団実行し、どの群にも使用中の状態・接続がないことを確認してから共通のリーダー通信路を解放します。


## プロセス分離型TDDFTの接続

`TDDFTMaterialProvider(bundle,input_file,runtime=managed_runtime,max_steps=...,threads=1)`
は管理済みランタイムのtddft・single-rank機能を確認します。入力・資産・バイナリのハッシュとワーカー設定を材料の識別情報に含めます。
各物質群で必要になった時点で初期状態を複製し、識別情報の作成後にファイルが変わっていれば拒否します。
使用後はproviderを閉じます。すでに作った点の状態は専用コピーを持つため、その後も計算できます。
波動関数、Hartree履歴、出力先は点ごとの既存TDDFTワーカーが所有し、伝播アルゴリズムは変更しません。

`TDDFTResponseBackend()` はleader_onlyで、**物質群はそれぞれMPI 1ランクに限定**します。
一群で複数点を順に計算できますが、各点には別々の直列子プロセスが常駐し、メモリとファイルを持ちます。
全点分のメモリと稼働ワーカーのスレッド数を見積もってください。MPIランク数とOSプロセス数は同じではありません。
この直列ワーカー入口にはTDDFT内部MPIはありません。群を直接渡す計算とE/H結合には[別のネイティブ入口](COUPLED.md)を使います。
電流は `TDDFTState` がすでに電荷電流 j = −J として返すので、ここでは符号を変えません。汎関数・伝播法は既存TDDFTワーカーと同じ対応範囲で、
準備済みGSの読み込みを行い、DFTを自動実行して作ることはありません。

```sh
mpiexec -n 3 python -m salmonpy.examples.mpi_tddft RUNTIME GS_BUNDLE RT_INPUT NEW_JSON
```

例は二つの1ランク物質群で計2状態を保持し、短い人工的なAランプを与えます。実際のレーザー実験を再現する条件ではありません。
実計算ではGS、空間格子、k点、バンド数、時間刻みを収束させてください。例はエラー時にMPIジョブ全体を停止し、他ランクの待ち残りを防ぎます。

独自ローダーから既存のTDDFTInitialStateを渡す場合は `TDDFTPreparedState(initial, managed_runtime)` で包みます。バックエンド側でも直列ランタイムの機能とバイナリハッシュを再確認し、未検証・不一致の初期状態は起動前に拒否します。

当面の電子計算は、一つのMPI計算群に一物質・一マクロ点を割り当てます。将来用の独立複数状態機能は残しますが、標準TDDFT例は一群一点です。LD中心のMaxwell計算では、既存FDTDの空間MPI分割の各ランクが多数のLD状態を管理します。FDTD_MATERIALS.mdを参照してください。両者は別の実行経路です。
