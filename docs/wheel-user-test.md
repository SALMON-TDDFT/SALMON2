# wheelから始める利用テスト

[English](wheel-user-test.en.md)

目的は、配布wheelの導入からLD・合成SBE・波形解析・形状作成・結果保存の5例まで、
SALMONpy本体を編集せずに実行できるかを確認することです。今回は実行・入出力の
スモークテストです。先のNotebookにあった解析解・保存則の追加検証は含みません。
DFT/TDDFT/FDTD伝播/MPI、実材料の予測、物理収束、Jupyter UIも対象外です。

## テスト参加者に渡すもの

- 成功したLD wheels CIから取得した、OS・CPU・Python版に対応する `.whl`。
  Actionsの外側ZIPは展開します。配布元URL・CI run・対象commitを記録してください。
- この手順書と `tools/test_wheel_user_workflow.py`。
  スクリプトはソース配布側のテスト用です。wheel内のモジュールではありません。
- 作業用Python（対応版）、書き込み可能な新規保存先、ネットワーク接続。
  Fortran、CMake、外部SALMON、Matplotlibは不要です。
  オフラインでは、対応するNumPy wheelもwheelhouseに用意します。

AIへの依頼例：

> このwheelを新しい環境へ導入し、指定の5例を動かしてください。
> 最初に未回答のメモリ予算とCPU予算を確認してください。数値処理は1スレッド・逐次実行です。
> 導入、計算、数値検証、資源計測を区別して記録してください。ソースビルドや本体修正へ
> 自動で切り替えないでください。失敗は残し、再試行には新しい出力先を使ってください。

## 実行

wheelに合うPythonで、次を実行します。大文字の引数は実際のパス・予算に置き換えます。
保存先はソースツリー外の未作成ディレクトリです。Windowsでも同じPythonコマンドを使えます。

```sh
python tools/test_wheel_user_workflow.py WHEEL_PATH NEW_OUTPUT --memory-budget "2 GiB"
```

予算が不明なら `--memory-budget unknown`。これは記録でありメモリ制限ではありません。
オフラインの場合は `--wheelhouse DEPENDENCY_WHEEL_DIRECTORY` を追加します。

スクリプトは新しいvenvを作り、バイナリ配布物だけをpipで導入します。
ソースからのコンパイルは許可しません。通常はNumPyをパッケージ索引から取得します。
導入後は `-I` と新規作業ディレクトリで、読み込んだSALMONpyがvenv内にあることを検査します。
その後、次の公式例を各1回、別プロセスで実行します。

1. `ld_response`
2. `sbe_response`
3. `pulse_spectrum`
4. `geometry_preview`
5. `outputs_roundtrip`

独立した例は前の例が失敗しても続行します。導入・読込元確認に失敗した場合は計算しません。
自動再試行はありません。実行時間の上限・メモリの強制制限はないため、異常に長い場合は
中断して記録を保存してください。

## 記録と判定

`report.json` にwheelのSHA256、環境、予算、各工程の成否を記録します。
各工程の `.stdout.txt`、`.stderr.txt`、`.json` と、SBE・格子・保存例の出力も残します。
`installed_identity.stdout.txt` に実際のパッケージ版と読込元を記録します。

計算の終了コードと計測の状態は別です。時計の取得に失敗しても、成功した計算を失敗にしません。
計算自体の非ゼロ終了は失敗のままです。ピークメモリは全工程で `not_measured` / null とし、
0や推定値で埋めません。全体成功もメモリ予算内だったことの証明にはなりません。
終了コード0はこのスモークテストの完了であり、独立した科学的検証の完了ではありません。

参加者は最初のAI依頼、AI名・版、手動介入、迷った点、準備も含む所要時間を別途残してください。
再試行は初回と分離し、変更内容を記録します。共有前にログのユーザー名・パス等を確認し、
venv全体を共有せず、ログ・report・計算結果を共有してください。
