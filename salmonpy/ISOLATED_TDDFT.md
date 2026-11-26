# 孤立系TDDFT：E入力と双極子

PZ・固定イオン・非偏極・スカラー擬ポテンシャル・中点伝播を使います。分子を別プロセスで動かす方法と、利用者が用意したMPI通信グループで動かす方法があります。結晶のTDDFT/SBEは従来どおりA駆動です。返す量は電荷双極子であり、電流密度ではありません。格子・時間刻み・計算箱の収束は利用者が確認します。

## 用途から選ぶ

- **FDTDの近接場で分子を励起する**：[一方向の励起](#one-way-grid)から始めます。現在推奨する入口です。一度の起動で両方を進める方法と、保存した場を読む方法があります。分子電流はFDTDへ返しません。
- **一様な電場で分子を励起する**：下の「計算プログラムと初期状態」を使います。
- 双方向の自己場差し引きを扱う後半の説明は、検証用の別経路です。一方向計算には不要です。

<a id="one-way-grid"></a>
## 空間依存FDTD E/Hによる一方向の励起

`python -m salmonpy.examples.molecular_grid_response WORKER VALIDATION BUNDLE input.inp FIELDS.npz NEW_OUTPUT` で、保存したFDTDの場から分子だけを時間発展させます。FDTDの起動、分子電流の返却、自己真空参照の差し引き、放射計算は行いません。分子が照明場やナノ構造に及ぼす影響は無視します。

既存GSの入力一式には、同じ格子形状の `salmonpy-near-eh.dat` を準備します。計算プログラムの検証記録には `grid_eh` と `grid_stream` が必要です。FIELDS.npzには `origin`, `spacing`, `times`, `E`, `H`, `molecule_position`, `center` と任意の `time_origin` を保存します。E/Hの形は `[フレーム,z,y,x,3]`。Yee格子の各成分の原点は保持します。HはFDTDの成分取得機能で整数ステップに時刻を揃えた値を使い、未処理の半ステップ値を渡しません。単位は原子単位、分子位置の媒質は真空です。

一定dtでフレーム0..N+1を用意します。分子のステップn-1→nには場のフレームn,n+1を渡します。`history.window(step_index=n,dt=dt)` がその2枚を取り出します。フィット、空間平均、時間補間は加えず、空間補間と経路積分による結合はFortranで処理します。`time_origin` は元データの時刻の記録であり、それ以前の分子の励起は再現しません。

`--stop-after 2 --direct-checkpoint NEW_CHECKPOINT` で途中保存し、同じコマンドの出力先を変えて `--resume-direct NEW_CHECKPOINT` で続けます。全履歴、準備済み形状、dt、計算容量、計算プログラムを揃えます。入力場も保存状態に保持され、再計算せず次の分子ステップから再開します。出力の `restore_statistics` と `final_statistics` で更新回数を確認できます。ステップ番号、双極子、電子分極の変化、電子数を保存します。内部エネルギーは全フレームで取得できる場合に保存し、取得可否も記録します。追加の一様電場ファイルは不要です。

### 分子内部もMPI並列にする

外側のPythonで場と分子にランクを分け、分子用のグループをFortranへ渡せます。分子の中では軌道と実空間の格子をMPIで分割できます。計算用共有ライブラリに分子MPIの入口が必要です。以前の実行ファイルと検証記録を渡す方式とは別で、対応するライブラリの`build.json`と既存のC2H2入力一式を使います。共有ライブラリは配布物を使うか、[計算プログラムの準備](RUNTIME.md)に従って明示的にビルドしてください。import時には準備しません。

次のコマンドは計算を開始します。管理ノードでの検査コマンドではなく、計算資源の割当内で実行します。

```sh
OMP_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-mpi --field-ranks 2 --molecule-ranks 2 --steps 96 --ld
```

MPIの総数は`field-ranks + molecule-ranks`です。上では場2＋分子2の計4プロセスを使い、子プロセスや内側の`mpiexec`は起動しません。両方のグループを一度の起動で進め、分子へは空間分布を保持したE/Hを渡します。分子電流の返却や自己真空参照はありません。現在の検証は、固定形状のLD＋C2H2例、内部1／2ランクの軌道並列と、x・y・z各方向の空間2分割、x/z方向の空間4分割、1ランク1スレッドです。軌道2×空間2／4の混合並列も検証しています。多ノード、これらより大きなランク数、他の計算条件は別途検証が必要です。

`--molecule-grid NX NY NZ`で分子の空間分割数を指定します。省略時は`1 1 1`で、これまでどおり軌道を分割します。空間分割数の積は`--molecule-ranks`を割り切る必要があります。残りの倍率が軌道分割数です。これはTDDFTの分割であり、FDTDの分割や計算格子の総数を変更する指定ではありません。

```sh
# 場2ランク＋分子のz方向2分割。計4ランク。
OMP_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-space --field-ranks 2 --molecule-ranks 2 --molecule-grid 1 1 2 --steps 96 --ld
# 場2ランク＋分子のx/z方向4分割。計6ランク。
OMP_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-space4 --field-ranks 2 --molecule-ranks 4 --molecule-grid 2 1 2 --steps 96 --ld
```

空間分割で保存／再開する時は、上と同じ分割指定を両方のコマンドに付け、最初の起動に`--stop-after 48`、次の起動に新しい出力先と`--resume SAVED_OUTPUT/checkpoint`を指定します。最大ステップ数は両方とも96のままです。検証では、空間分割と独立した1ランク実行の差は観測量で4×10⁻¹³以内でした。z方向2分割とx/z方向4分割の直接再開は、連続実行の全配列と通常出力4ファイルに完全一致し、復元時の時間発展は0回でした。これは接続と再開の検証であり、物理的収束の保証ではありません。

場と分子がアクセスできる共有作業フォルダを出力先に使います。場の代表ランクがE/H入力を書き、分子の代表ランクが一度だけ読み、その分子グループ内へ配ります。毎ステップの通信にアプリケーション全体の集約は使いません。mpi4pyと計算ライブラリには同じMPI実装が必要です。分子ランクは計算中専用で、同じプロセス内で他の計算やファイル操作を並行させません。

```sh
OMP_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-mpi-first --field-ranks 2 --molecule-ranks 2 --steps 96 --ld --stop-after 48
OMP_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-mpi-resumed --field-ranks 2 --molecule-ranks 2 --steps 96 --ld --resume run-mpi-first/checkpoint
```

分子48／場49で、場・LD・分子をまとめて直接保存します。再開時も最大ステップ数96、入力、バイナリ、MPI分割、スレッド設定を揃えます。分子の保存データは内部ランク別で、入力・バイナリ・MPI分割とファイルの整合性を確認して復元します。波動関数からの再計算や、経過したステップの再実行はしません。異なるMPI分割への移行はこの入口では扱いません。

結果は`response/arrays.npz`、実行情報は`response/manifest.json`の`metadata`です。分子の配列は`molecular_step`、E/Hは`field_step`で対応付けます。分子のステップn−1→nには場のn,n+1を使うため、場のフレームが分子より1枚多くなります。再開後の履歴には保存前も含まれ、前ジョブの配列を再結合する必要はありません。`restore_statistics`の更新回数0と、`worker_statistics`の新規更新回数を確認してください。保存／再開を含む同じ設定で全履歴と通常出力の一致を確認していますが、光学物性の収束を保証する試験ではありません。

軌道と空間を同時に分割する場合も、同じ入口を使います。例えば、分子4ランクに`--molecule-grid 1 1 2`を指定すると、空間2分割×軌道2分割です。分子8ランクに`--molecule-grid 2 1 2`なら空間4分割×軌道2分割になります。場に2ランクを使う場合、MPIの総数はそれぞれ6／10です。

```sh
# 場2＋分子（軌道2×空間2）の計6ランク。
OMP_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-mixed --field-ranks 2 --molecule-ranks 4 --molecule-grid 1 1 2 --steps 96 --ld --stop-after 48
OMP_NUM_THREADS=1 mpiexec -n 6 python -m salmonpy.examples.molecular_fdtd_mpi_response BUILD_RECORD EXISTING_C2H2_BUNDLE run-mixed-resumed --field-ranks 2 --molecule-ranks 4 --molecule-grid 1 1 2 --steps 96 --ld --resume run-mixed/checkpoint
```

同じ固定C2H2例で、両方の混合分割による96ステップ計算を独立した1ランク実行と比較しました。双極子・電子数の差は4×10⁻¹³以内で、FDTDのE/Hは完全一致しました。両方の混合分割について、分子48／場49からの直接再開も、全履歴と通常出力4ファイルが連続実行に完全一致し、復元時の時間発展は0回です。これらは1台のMac・PZ・固定イオン・スピンなし・1ランク1スレッドの検証です。空間分割数と軌道分割数を再開時に変える機能ではありません。

部品として使う場合は、分子グループ全員で`run_molecular_mpi(..., decomposition=(1, N_ob, NX, NY, NZ))`を呼び、場側の代表ランクが`MolecularMPIController`で入力します。前者は分子の寿命の間Fortran内で待機・計算する呼び出しで、後者の`close()`で終了します。`N_ob`は軌道分割数で、分子グループのランク数は`N_ob*NX*NY*NZ`です。`MolecularMPIController`にも同じ`decomposition`を渡します。通信グループと作業先は呼び出し元が所有し、終了後もMPIは利用できます。この入口は1プロセスにつき1寿命です。`controller.save_direct_checkpoint()`／`restore_direct_checkpoint()`は同じ配置での直接保存／復元です。プロセス終了・MPI異常・深部数値異常からの一般的な復旧や、復元途中の全グループの一括巻き戻しは保証しません。異常時はジョブを終了し、保存状態から新しい起動で再開します。密度・軌道磁化の読み出しは従来の1ランク側の機能で、分子MPI側にはまだ接続していません。

### 一度の起動でFDTDと分子を進める

FDTDと分子TDDFTを一つのPythonプログラムから順に呼ぶ例もあります。場をファイルへ書いてから別の計算を起動する必要はありません。一方向とは、分子の電流をFDTDへ返さないという意味です。

```sh
OMP_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.molecular_fdtd_response FIELD_BUILD WORKER VALIDATION EXISTING_C2H2_BUNDLE run-live --field-ranks 2 --steps 96 --ld
```

`FIELD_BUILD`はFDTD共有ライブラリのbuild.json、`WORKER`と`VALIDATION`は分子実行ファイルとその検証記録です。`EXISTING_C2H2_BUNDLE`には既存GSの再開データ、C/H擬ポテンシャル、`input.inp`を置きます。サンプルが実際のパッチ形状に合わせた初期入力を作るので、FIELDS.npzや`.dat`の事前生成は不要です。GS計算・ビルド・ダウンロードは行いません。

これは固定形状の接続例です。40³・幅2bohrのFDTD格子、任意で合成LD箱、格子点からずれたC2H2位置、24³の場パッチ、dt=0.003auを使います。任意のナノ構造をこのコマンドだけで指定する入口ではありません。FDTDは1または2ランク、ほかに調整用1ランクを使います。上の3ランクの例では分子ワーカーが別プロセスで1つ起動するため、計4プロセス分のCPU資源を見込みます。`--field-ranks 1`ならMPIは2ランクです。分子へのE/Hは成分ごとの空間分布を保持し、Fortranの補間・経路積分で入力します。自己真空参照は起動しません。

```sh
OMP_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.molecular_fdtd_response FIELD_BUILD WORKER VALIDATION EXISTING_C2H2_BUNDLE run-first --field-ranks 2 --steps 96 --ld --stop-after 48
OMP_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.molecular_fdtd_response FIELD_BUILD WORKER VALIDATION EXISTING_C2H2_BUNDLE run-resumed --field-ranks 2 --steps 96 --ld --resume run-first/checkpoint
```

分子48／FDTD49で場・LD・分子を保存して終了し、同じランク配置・入力・最大ステップ数で再開します。`--save-step 48`なら保存後も継続します。この例の保存先は出力先直下の`checkpoint`で、一回の起動で1つだけ作ります。

`summary.json`に実際の分子・FDTDステップ番号、更新回数、返却と自己真空参照の有無を記録します。`history.npz`の`observations[k]`は分子ステップkの9成分（電子分極変化3、全電荷双極子3、電子数、内部エネルギー、エネルギー取得フラグ）です。`total[j]`は`field_step[j]`で示すFDTDステップj+1のE/H。0ステップ目は`field-step-0.npz`の`EH`にあります。再開した実行の履歴には保存前の履歴も含まれるため、前の出力とさらに結合する必要はありません。分子の通常出力は`molecule-streams`に保存します。大きな計算向けのI/O量削減や任意のMPI障害からの復旧を保証する例ではありません。

以下は、FDTDを既に計算済みで、保存した場から分子だけを計算する場合の手順です。

### 入力を同じ場データから準備する

`molecule-input` は、既存GSの再開データ・擬ポテンシャル・孤立系RT入力 `rt.inp` を揃えた準備先です。ここでGS計算は行いません。`salmon` は使用する検証済み実行ファイル、`validation.json` はその検証記録に置き換えます。まず次の形式でFDTDの場を `FIELDS.npz` に保存してください。

| 項目 | 形と意味 |
|---|---|
| `E`, `H` | `[N+2,nz,ny,nx,3]`。電場・磁場とも原子単位。全フレームの格子形状は固定 |
| `times` | `[N+2]`、原子時間。0から一定dtで増える。Nは分子を進めるステップ数で1以上 |
| `origin` | 共通格子なら`[3]`、Yee格子なら`[2,3,3]`＝`[E/H,成分,座標軸]`、bohr |
| `spacing` | `[3]`、各軸の格子幅、bohr |
| `molecule_position` | `[3]`、FDTD座標での分子の積分中心、bohr。格子点に一致する必要はない |
| `center` | `[3]`、同じ積分中心のTDDFT座標、bohr |
| `time_origin` | 任意のスカラー。元のFDTDデータの開始時刻、原子時間。再開ステップ番号とは別 |

FDTDの`field.component_fields(indices)`を保存する場合、戻り値は`[2,点数,3]`です。indicesをz→y→x（xが最速）の順に並べて矩形パッチを読み、`[2,nz,ny,nx,3]`へ整形します。各ステップの先頭軸を重ねてE/Hを取り出します。`field.component_field_geometry()`の成分別originにパッチ開始インデックス×spacingを足した原点を記録します。この取得機能のHは整数ステップに揃え済みです。配置が異なるE/Hを同じ原点とみなさないでください。サンプリングは確定したFDTDステップで行います。

次の準備コードを一度実行すると、**計算に使うFIELDS.npzそのもの**から初期化用ファイルも作れます。計算は起動しません。既存の`salmonpy-near-eh.dat`は上書きせず、新しい準備先を使います。

```python
from pathlib import Path
import numpy as np
from salmonpy.near_field import GridFieldSeries

with np.load("FIELDS.npz", allow_pickle=False) as data:
    names = ("origin", "spacing", "times", "E", "H", "molecule_position", "center")
    history = GridFieldSeries(
        **{name: data[name] for name in names},
        time_origin=float(data["time_origin"]) if "time_origin" in data else 0.)
history.write(Path("molecule-input") / "salmonpy-near-eh.dat")
```

`prepare_grid_field.py`は点座標付きサンプルから`.dat`を作る別の準備例です。その入力NPZは上表の格子形式とは異なるため、そのまま今回の計算例には渡しません。ここでは上のコードで2つの入力を揃えます。

### 計算・途中保存・再開

分子を2ステップ以上進める場データを用意した場合の例です。出力先と保存先は未使用の名前にします。2回目も同じFIELDS.npzと入力一式を使います。

```sh
python -m salmonpy.examples.molecular_grid_response salmon validation.json molecule-input rt.inp FIELDS.npz result-1 --stop-after 2 --direct-checkpoint checkpoint-2
python -m salmonpy.examples.molecular_grid_response salmon validation.json molecule-input rt.inp FIELDS.npz result-2 --resume-direct checkpoint-2
```

途中停止しない連続計算では、1行目から`--stop-after`と`--direct-checkpoint`を省き、別の出力先にします。ここで起動するのは分子計算のみで、MPIでFDTDを起動するコマンドではありません。この保存済み場を読む経路はPZ・固定イオン・非偏極スカラー擬ポテンシャル・CPU1ランクです。

### 出力を読む・再開前後を結合する

各出力先の`arrays.npz`に数値配列、`manifest.json`の`metadata`に実行条件と更新回数があります。再開後の最初のフレームは保存点そのものです。同じ保存点が一致することを確認し、後半の先頭を1枚だけ除いて、最後に一度結合します。

```python
import json
from pathlib import Path
import numpy as np

metadata = json.loads(Path("result-2/manifest.json").read_text())["metadata"]
print(metadata["restore_statistics"], metadata["final_statistics"])
with np.load("result-1/arrays.npz", allow_pickle=False) as first, \
        np.load("result-2/arrays.npz", allow_pickle=False) as second:
    for name in ("molecular_step", "time", "dipole", "electronic_change", "electron_number"):
        np.testing.assert_array_equal(first[name][-1], second[name][0])
    combined = {
        name: np.concatenate((first[name], second[name][1:]), axis=0)
        for name in ("molecular_step", "time", "dipole", "electronic_change", "electron_number")
    }
```

`restore_statistics["updates"]`は復元直後の新規時間発展回数（通常0）、`final_statistics["updates"]`はそのプロセスで進めた回数です。`internal_energy`は全出力フレームで取得できる区間にだけ存在します。結合する両区間にある場合だけ同様に扱い、存在しない値を0で補いません。`metadata["internal_energy_available"]`で取得可否を確認します。

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

読込みは64 KiB以下のローカル検証記録と、実行ファイルのSHA256を照合します。現行の記録形式は `salmonpy.isolated-worker-validation/5`。数値応答はE3の9成分のままですが、同一プロセス・時間発展追加0回の直接復元と、古い退避の拒否の実証を必須にしています。旧revision4以前の記録は受け付けません。読込み時にダウンロード・ビルド・検証計算・MPI初期化・ワーカー起動は行いません。記録はローカル検証の証拠で、電子署名や禁止ソルバー監査の代用ではありません。

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

## プロセス終了後の直接ディスク再開

```python
with model.new_state(dt=dt) as state:
    for E in first_fields:
        state.step(E)
    state.save_direct_checkpoint('new-electronic-checkpoint')
# 元のワーカーは終了。同じGSと設定から新しい状態を作る。
with model.new_state(dt=dt) as resumed:
    resumed.restore_direct_checkpoint('new-electronic-checkpoint')
    assert resumed.worker_statistics()['updates'] == 0
    for E in remaining_fields:
        frame = resumed.step(E)
```

`salmonpy.isolated-direct/1`はFortranの数値配列・履歴・Poisson・XC作業状態、
受理ステップ、観測値、出力ファイルを保存します。各ファイルのSHA256と、
モデル・実行ファイルの一致を読み込み前に確認します。元のGSは別途指定し、
初期化時に読み込みます。SCFや経過済み時間ステップの再計算はしません。
別ワーカーで読み込みと出力再接続を完了してから切り替えるため、
候補の復元失敗時も元のワーカーを保持します。復元した状態の異常時に、
履歴再生へ暗黙に切り替えません。復元中は電子状態が追加で必要です。

実行ファイル・GS／入力・dt・最大ステップ数・位置・スレッド数は同一にします。
異なるコンパイラや配列配置をまたぐ可搬形式ではありません。
保存点を解放してからディスク操作を行います。max_bytesは任意の読み込み量制限で、
メモリ使用量の制限ではありません。既存ディレクトリは上書きしません。
対象は既存の確認済み出力モードです。ワーカー検証記録は実際のプロセス終了後の
再開を含むrevision5を要求し、古い記録は受け付けません。

応答サンプルで`--stop-after 7 --direct-checkpoint NEW_CHECKPOINT`を指定すると、
7ステップで終了して保存します。同じE NPZと設定で別途
`--resume-direct NEW_CHECKPOINT`を指定すると、残りだけを進めます。
これは分子単体の保存です。FDTD結合全体の再開には二つの場と結合時刻も
協調して保存する必要があります。

## 履歴の保存と再実行（直接復元とは別）

`snapshot()/restore()` は受理済みEと観測量の履歴を保存・再実行する入口です。同じ初期モデル、dt・位置・容量・スレッド数・バイナリが必要です。候補プロセスの観測量を照合してから交換し、復元失敗時は元の状態を保持します。費用は履歴長に依存し、検証中は新旧二ワーカー分の資源を使います。直接退避を使わない通常の `step()` の異常回復もこの履歴再実行方式です。

```python
state.save_checkpoint('new-history-checkpoint')
resumed.restore_checkpoint('new-history-checkpoint', max_bytes=16*1024*1024)
```

`resumed` は同じ初期データ・設定で呼び出し元が作成した状態です。ディスク形式は `salmonpy.isolated-tddft-checkpoint/3`。E0、初期／現在の9観測量、E履歴と観測履歴をハッシュ付きで保存するZarr v3互換形式です。GSは別途用意します。`max_bytes` は追加のファイル読み込み予算であり全RSS上限ではありません。これは波動関数の直接チェックポイントではありません。直接の数値保存には上記のsave_direct_checkpointを使います。

## 出力とFDTD

`python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT` は、E[steps,3]とスカラーdtを原子単位で持つNPZから双極子アーカイブを作ります。`--checkpoint` は明示した履歴再実行用の保存です。応答・単位・ソースの出典を明示し、結果はソースツリー外へ保存します。

双方向FDTD結合は [分子とFDTD](DIPOLE_RADIATION.md) の `MolecularFeedback` と、完全例 `examples/molecular_feedback.py` を使います。取消に `savepoint/restore_savepoint/release_savepoint` を使い、電子履歴を再計算しません。真空自己場の差し引きは真空の放射反作用も除きます。

既存C2H2 GS・PZ、64場ステップの6ランク検証では、正常時は電子更新63回、1区間の取消・再試行時は64回で、各実行内のPIDは変わりません。復元による追加更新は0回、受理した全保存配列・六面記録は以前の結果と完全一致しました。これは短時間の接続検証で、光学物性の収束保証ではありません。

## 入力で指定する双極子源の近接場

SALMON入力の `num_dipole_source`、`cood_dipole_source`、`vec_dipole_source`、
`rad_dipole_source` を使う入口は、`from_salmon(..., native_near_field=True)` で明示します。
既定は無効です。この入口には、使用する実行ファイルについて近接場の比較結果を記録した
検証記録が必要です。従来の一様電場だけの検証記録では有効になりません。
C₂H₂の短時間検証で本体との双極子比較、独立状態、同じプロセス内の復元と、プロセス終了後の直接ディスク再開を確認しました。物理量の収束検証とは区別します。

対応範囲はPZ・固定イオン・スカラー擬ポテンシャル・1ランクです。
双極子源は最大2個、半径は正の有限値、波形は `Acos2` または `Ecos2` です。
座標・半径はSALMON入力の単位系に従い、Pythonの `step(E)` は引き続き原子単位です。
入力の双極子源の空間分布と時間波形を保持し、`step(E)` は**追加の一様電場**を渡します。
源の数が0なら一様電場だけです。源の入力や有効化設定が異なる状態間では復元できません。

これは入力で定めた双極子源の機能です。一般の空間依存E/Hを受け取る入口とは区別します。

## 空間依存E/Hを使う近接場モード

`from_salmon(..., near_field=True)` は、分子周囲のE/Hを経路積分して使います。
現在の対象は **PZ・固定イオン・非偏極・スカラー擬ポテンシャル・CPUの1ランク**です。
磁気結合は軌道に対するものです。スピン・SOCの磁気結合、分子からFDTDへの放射の戻しは含みません。
分子が置かれる領域は真空とし、Gaussian原子単位で `B=H` を使います。

### FDTD格子から直接渡す（展開しない方法）

`GridFieldSeries` は、保存したE/Hの格子データを直接補間し、Iwasaらの電場の経路積分へ渡します。
多極子の次数や減衰率を選ぶ必要はありません。現在は等間隔の3次元直交格子で、セル内は3重線形補間です。
経路が格子セルを跨ぐ場所で積分を分割します。磁場による軌道結合も同じ補間場から計算します。

```python
from salmonpy.near_field import GridFieldSeries

# E_samples、H_left、H_rightは、同じセル中心で保存したSampledField。
# H_leftはn-1/2、H_rightはn+1/2。相対時刻0の元のFDTD時刻も保持する。
grid = GridFieldSeries.from_yee_samples(
    E_samples, H_left, H_right,
    counts=(16, 17, 17), origin=(7., -15., -15.), spacing=(2., 2., 2.),
    molecule_position=(21.2, .7, .9), center=(0., 0., 0.), time_origin=.303)
grid.write('NEW_INPUT_BUNDLE/salmonpy-near-eh.dat')  # 保存GSを再利用する入力準備のみ

with IsolatedTDDFTInitialState.from_salmon(
        'NEW_INPUT_BUNDLE', 'rt.inp', runtime=runtime,
        max_steps=20, near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(7):
            state.step([0., 0., 0.])  # 追加の一様E。格子の場は保存時系列から取得。
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

位置・間隔はbohr、時刻は原子単位です。`counts`は各軸のサンプル点数で、FDTD全体の格子数とは区別します。
`origin`は最初のサンプル位置です。`molecule_position`は積分中心のFDTD座標、`center`は同じ点のTDDFT座標で、
Fortranが `FDTD位置 = molecule_position + (TDDFT位置-center)` として座標を対応させます。
どちらの位置も格子点と一致する必要はありません。座標軸は共通とし、回転は自動では行いません。原子だけでなく、TDDFT計算箱・差分の外側格子点と
非局所ポテンシャルの経路まで覆うパッチを用意します。格子外の推測や暗黙の周期折り返しはしません。
分子領域は真空としてH=Bを使います。

HをEの時刻へ揃えて保存済みなら `from_samples(E_samples, H_samples, ...)` を使います。
構造化配列を渡す場合は `GridFieldSeries(origin, spacing, times, E, H, ...)` で、E/Hは
`[時刻,z,y,x,成分]` です。成分ごとの格子位置が異なる場合は、`origin`を `[E/H,成分,座標軸]` の
配列で指定します。各成分は同じ点数・間隔を持つ共通のパッチに切り出します。
半ステップの構造化H配列には `from_yee_arrays(..., H_left_times=..., H_right_times=...)` を使えます。
いずれも時刻合わせ・空間補間・積分はFortranで行い、Hや端点を省略して推測することはありません。

入力サイズは時刻数×格子点数に比例するため、必要な局所パッチと区間を切り出します。
格子の時系列は相対時刻0から一定dtで、入力フレーム0〜max_steps+1を用意します。入力のdtとTDDFTのdtを揃えます。
ステップnではフレームnとn+1を使い、`step(E_extra)`で進めます。`step_near_field`は展開係数用なので使いません。
元のFDTD時刻は`time_origin`として入力ファイルにも記録します。この情報で分子の過去の時間発展を
再現したことにはなりません。初期状態と照射開始の条件を明示してください。
保存／再開では格子データ・座標対応と電子状態を保持し、同じ入力モデル・dt・計算プログラムで
`restore_direct_checkpoint`後に残りのステップを進めます。

準備用サンプルは `python -m salmonpy.examples.prepare_grid_field --help` です。準備だけで計算は始まりません。
直接格子入力の検証記録を持つ実行ファイルが必要です。磁気結合は元のIwasaらの電気相互作用とは
別に検証しており、スピン磁気結合や分子からFDTDへの放射の戻しは対象外です。
補間してもFDTDで解像できない場は復元できないので、格子と積分・分子応答の収束を確認します。

直接補間の解析場試験では、成分ごとにずれた格子で線形場・経路積分の絶対誤差は最大約2.1×10⁻¹⁷でした。
滑らかな場を格子幅0.5→0.25→0.125 bohrで調べると、電場の相対RMS誤差は6.60→1.46→0.328%、
結合ポテンシャルは3.42→0.646→0.0901%へ減少しました。実際のLD/FDTD保存場をC2H2へ渡す試験では、
分子位置をFDTD格子からずらし、ステップ7の保存から20までの全9観測量と出力が連続実行と完全一致しました。
これは接続と再開の検証であり、実ナノ構造・分子応答の物理的収束を保証するものではありません。

### フィット時に中心を決める

```python
from salmonpy.near_field import NearFieldSeries

# E_samples/H_samples: 同じ位置・時刻のSampledField。単位はE:au/H:au、位置はbohr。
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=1)
print(near.center, near.relative_residual)  # 残差の列はE、H
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # 準備だけ。計算は始めない。
```

中心は分子中心などを明示指定することを勧めます。省略時はサンプル点の重心、半径の省略時は
そこから最も遠い点までの距離を、最初のフィットでFortranが決めます。
決まった中心・半径を全フレームと保存／再開で共有し、実行途中に選び直しません。
半径は基底の規格化尺度です。FDTDの分解能を細かくする操作ではありません。

E/Hの各成分を、半径方向の項を含む実Ylmで展開します。`lmax` は総次数dの上限です。
係数は `d=0..lmax`、各dで `l=d,d-2,...`、各lで `m=-l..l` の順に並び、
`(r/radius)^d Y_lm` を掛けます。d=l+2kなので半径因子r^(2k)も含みます。
項数は `P=(lmax+1)(lmax+2)(lmax+3)/6`、配列は `[時刻,P,xyz]` です。
電気ポテンシャルだけの入口とは係数の単位が異なります。
`lmax=0` は空間一様場、`lmax=1` は勾配を含む場です。1次には電気四重極までが含まれ、
一様Hは軌道磁気双極子に結合します。より高い次数も使えますが、有限次数で任意の場を厳密に表せるわけではありません。
フィットの相対残差が既定の5%を超えると警告し、計算可能な結果は返します。

FDTDのHが半ステップずれる場合は `NearFieldSeries.from_yee_samples(E_samples, H_left, H_right, ...)`
を使います。H_leftはn−1/2、H_rightはn+1/2、Eはnです。
位置・時刻の対応を確認し、Hの平均をFortranで計算します。初期の前後のHも明示的に渡してください。
欠けた端点の推測や外挿はしません。

準備用の例は `python -m salmonpy.examples.fit_near_field --help`。
実行例 `examples/isolated_response.py` には `--near-field` があります。
フィット・保存・importはTDDFTを起動しません。


分子位置と展開中心は実座標で指定し、FDTD格子点との一致は要求しません。
周囲の既知のE/Hサンプルから、その中心での展開を作ります。格子点に分子を移動したり、
最近傍の場だけを採用したりしません。現在は保存した場をTDDFT入力として使う一方向の処理です。
分子応答をFDTDへ戻す補間・電流配置・エネルギー整合性は別の設計課題です。
中心から分子の電子密度が広がる領域をサンプルが十分に覆うようにし、金属界面をまたぐ
単一の滑らかな展開や、サンプル領域外への外挿を精度保証と扱わないでください。

TipやBow Tieに必要な次数は未確定です。次数だけでなく、先端半径・ギャップ・分子位置・
フィット領域を変え、フィットに使わない点で絶対・相対残差を確認します。
現在の基底は指定した総次数までの完全多項式を表せます。ただしFDTDで解像できない情報は
復元せず、Maxwell方程式の整合性も自動では保証しません。独立した位置の残差とFDTD分解能を確認します。


従来の診断では角度成分だけの表現と完全多項式を比較しました。TDDFT入力も現在は半径項を含みます。
格子点での残差が小さくても、格子間の場や分子応答の精度は別途確認が必要です。
総次数8では各成分165係数になり、少なくとも165点の独立なサンプルが必要です。
点数・独立性不足は収束として扱いません。入力ファイルV3・直接保存V7で表現を識別し、
以前のファイルは作り直します。

半径項付きの接続検証では、C2H2の既存GSを使い、非ゼロの2次半径項を含むE/Hで20更新しました。
更新7で直接保存し、プロセスを終了してから再開した13更新の全9観測量は連続実行と完全一致し、
復元時の追加時間発展は0回でした。独立した解析場では総次数2・4・6・8の格子間点のE/H誤差は
最大約1.3×10⁻¹⁵ auです。これは表現・接続の検証であり、Tip/Bow Tieの物理的収束保証ではありません。


### 指数減衰を含む場

既知の方向に急減衰する場には、同じ半径項付き基底へ指数因子を掛けられます。

```python
near = NearFieldSeries.from_samples(
    E_samples, H_samples, center=(0., 0., 0.), radius=8., lmax=4,
    decay_vector=(0., 0., .25))  # 単位bohr^-1。正のz方向に減衰。
near.write('GS_BUNDLE/salmonpy-near-eh.dat')  # 入力準備だけ。
```

各E/H成分の基底は `exp(-decay_vector·(r-center)) (|r-center|/radius)^d Y_lm` です。
半径は多項式部分の規格化尺度で、減衰ベクトルの単位は変えません。
既定の `(0,0,0)` は通常の多項式展開です。E/Hで共通の減衰ベクトルを使い、中心・半径・次数と同様に、
全フレームと保存／再開で固定します。`NearFieldFrame.from_samples` にも同じベクトルを指定します。
減衰率の自動選択は行いません。別の点で選んだ率を使う場合も、最終検証点は分けて残差を確認します。
選ばれた率が物理的な減衰長を一意に示すとは限りません。単一因子は一般のTip/Bow Tieの全近接場を
表すものではなく、逆方向には増大するため、分子を覆う検証済み領域で評価してください。

場の近似と分子への結合積分はFortranで分離しています。指数因子を含む場合も、電場は経路上のEを、
磁場は経路パラメータを掛けたHを積分して結合ポテンシャルを作ります。指数積分では小さい引数の
桁落ちを避け、磁気の非局所ポテンシャル位相は指数依存に対応した数値積分を使います。
別の場表現を後から渡すための、経路上のサンプルから積分する共通処理もあります。
FDTDからの直接補間や双方向接続まで実装したという意味ではありません。


自作の場表現には `salmonpy.near_field.couple_path_samples(displacement, parameters, weights, E, H)`
からも共通積分を呼べます。中心からの変位はbohr、経路パラメータは0〜1、重みはその積分に対する重み、
E/Hは経路上の同じ位置にある `[点数,3]` 配列です。戻り値は電子のスカラーポテンシャルエネルギーとA/cです。
経路上の場と積分点・重みは呼び出し元が用意し、積分点数の収束も確認します。この関数だけではTDDFTは進みません。

### 係数を更新ごとに渡す

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        near_field=True) as initial:
    with initial.new_state(dt=.003) as state:
        saved = state.savepoint()
        frame = state.step_near_field(near.frame(1), near.frame(2))
        state.restore_savepoint(saved)
        state.release_savepoint(saved)
        state.step_near_field(near.frame(1), near.frame(2))
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

更新nでは、SALMONの既存順序に合わせ、場のステップnとn+1を渡します。
各フレームには中心・半径・次数・時刻・ステップ番号・E/H係数が入ります。
`NearFieldFrame.from_samples(..., step_index=n, center=..., radius=...)` で1回分の空間E/Hから
係数を作ることも、`NearFieldFrame` に自作コードの係数を直接渡すこともできます。
隣接更新で共有する場のフレームを修正する場合は、その場を使う前のsavepointへ戻して更新し直します。

将来の場がまだない場合は、最初の `NearFieldFrame` の `write_live_input(path,max_steps=...,dt=...)`
で初期値を保持する入力を準備し、各更新で `step_near_field` を呼びます。
この入力は将来の場の予測ではありません。
通常の `step(E_extra)` は入力済みのE/H系列を使い、引数を追加の一様Eとして扱います。

保存は `savepoint` または `save_direct_checkpoint`、再開は対応する直接復元を使います。
ライブ係数はEだけの再実行履歴では保存できないため、近接場モードのreplay snapshotは提供しません。
係数系列・中心・半径・磁気用配列・非局所位相まで直接保存します。復元時に時間発展を再計算しません。

### 結合と確認範囲

Fortran側では、展開中心からの変位rについて
`V_e = r·integral_0^1 E(s r) ds`、`A/c = -r×integral_0^1 s H(s r) ds / c`
を使います。`V_e` は電子のポテンシャルエネルギーです。
対称化した運動エネルギーの磁気項とA²項を保持し、非局所擬ポテンシャルには
原子中心からの線積分位相を使います。これは[動径方向のゲージによるポテンシャル](https://arxiv.org/abs/physics/0204034)
に基づく、指定された古典場との結合です。

C₂H₂で電気四重極・磁気の効果、入力済み系列とライブ係数の一致、初期状態への復元、
プロセス終了後の直接再開を確認しています。これは短時間のソフトウェア検証です。
実材料の応答、展開次数、FDTD/電子の格子、dt、場の近似残差は目的に応じて確認します。E/Hを独立にフィットしただけでMaxwell方程式との整合性が保証されるわけではありません。
場のないGSから開始し、先頭から非ゼロのE/Hを与えた場合は、その初期状態からの時間発展です。
過去の照射履歴や場ありのGSを自動で作ることはありません。

## 保存した電場の準静的ポテンシャル近似

`salmonpy.near_field.ElectricPotentialSeries` は、同じ位置・時刻で取得した電場を、
実球面調和関数の正則な固体調和関数で近似します。数値処理はFortranです。
**一般の近接場のE/H入力ではなく、スカラーポテンシャルで表せる電気成分の近似です。**
Hの結合、磁気多極子、分子からFDTDへの放射の戻しはこの入口に含みません。
FDTDの格子で失われた空間情報も復元しません。

```python
from salmonpy.near_field import ElectricPotentialSeries

# electric_samples: SampledField('E', times, positions, values, 'E:au')
# positionsはbohr。Eの3成分の位置・時刻を揃えてから渡す。
potential = ElectricPotentialSeries.from_samples(
    electric_samples, center=(0., 0., 0.), radius=8., lmax=2)
print(potential.relative_residual)
potential.write('GS_BUNDLE/salmonpy-ylm.dat')  # 準備だけ。計算は開始しない
```

同じ入射波・位置・時刻の真空参照があれば、`reference=vacuum_samples` で差の電場を近似できます。
差の計算もFortranで行います。相対残差が既定の5%を超えた場合は警告し、値を返します。
電場がほぼゼロの時の相対値だけで判断せず、`electric_field(...)` と元の電場も比較します。
次数を上げても磁場や回転する電場成分の欠落を補えるとは限りません。

実関数は複素Ylmの[NISTの規格化](https://dlmf.nist.gov/14.30#E1)から、正のmは
`sqrt(2)*(-1)^m*Re(Ylm)`、負のmは対応する正のmの虚部を同じ係数で使います。
`l=1..lmax`、各lの中で`m=-l..l`の順です。`lmax=1..8`に対応します。
電子のポテンシャルは `V_e=sum(c_lm*(r/radius)^l*Y_lm)`、電場は `E=grad(V_e)`。
係数はHartreeで、定数ポテンシャルl=0は省きます。領域内に外部電荷がない近似を前提とします。

実行は明示的に別の操作で行います。対応する実行ファイルのYlm検証記録が必要です。

```python
with IsolatedTDDFTInitialState.from_salmon(
        'GS_BUNDLE', 'rt.inp', runtime=runtime, max_steps=80,
        electric_potential=True) as initial:
    with initial.new_state(dt=.003) as state:
        for _ in range(40):
            state.step([0., 0., 0.])  # このEは追加の一様電場
        state.save_direct_checkpoint('NEW_CHECKPOINT')
```

係数表の時刻は0から一定dtで、`0..max_steps+1`の行が必要です。
更新番号nではSALMONの既存処理に従いnとn+1の係数を参照します。
dt・入力・実行ファイルが同じ状態へ、通常と同じ`restore_direct_checkpoint`で復元できます。
FDTDの途中の時間窓を使う場合、元の時間ステップ・時刻と分子側の時刻原点の対応を保存してください。
場のないGSからその窓で開始する計算は、前の照射履歴を持つ分子の継続とは異なります。

3次元LD/FDTDの27点から分子へ渡し、80更新のうち40更新でプロセスを終了して直接再開した試験で、
連続実行と完全一致しました。一様場のYlm表現は本体と約4e-12 e*bohr以内で一致しています。
ただし実FDTD試験では電場の近似残差が大きく、近接場の物理精度を検証した結果ではありません。
空間依存E/Hを使う場合は、上記の近接場モードを使って電気・磁気を一緒に渡します。

## 空間依存電場に対応する電荷応答の分配

`GridFieldSeries.charge_polarization` と `charge_current` は、電場の経路積分と対になる分配をFortranで行います。
分子内の各点の電荷密度から、E成分の格子上に分極または区間平均電流を返します。
入力の `charge_density` は **符号付きの物理的電荷密度**（e/bohr³）です。正の電子数密度ではありません。
電子のみなら電子数密度に負号を付け、`volume_elements` は各点の積分体積（bohr³）を渡します。
前後の原子配置・点・体積が同じ場合に限り、固定イオンの電荷は電流の差分から消えます。

```python
J = field_series.charge_current(
    molecular_positions, old_charge_density, new_charge_density,
    volume_elements=volumes, dt=dt)
P = field_series.charge_polarization(
    molecular_positions, new_charge_density, volume_elements=volumes)
```

分子座標の点を展開中心からFDTD上の分子位置へ平行移動し、各経路を格子セル境界で区切って分配します。
返り値は `[z,y,x,component]`、Pはe/bohr²、Jは原子単位です。成分ごとに入力Eと同じ格子を使います。
空間の仕事量は `dt Σ_grid E·J ΔV = Σ_mol Δrho_charge V_e Δv` となります。
ここで `V_e(r)=(r-R)·∫E(R+s(r-R))ds` は電子ハミルトニアンに入る経路積分です。
この等式は指定した同じ電場と電荷変化に対する空間の対応であり、時間発展の精度や全領域の収支を保証しません。

この入口は配列を受け取る計算部品です。自作コード等が用意した密度に利用できますが、密度の明示的な取得と成分格子への直接電流注入は次節の入口を使えます。
逐次入力と自己真空参照を組み合わせた電気応答の双方向試験は、下記の検証用ドライバーで確認しています。汎用の本番用統括APIは未提供です。
このJを双極子用の `DipoleMap` やマクロ点の電流入口へそのまま渡すと、重複した補間になります。
軌道磁化電流の別入口は後述します。全結合の離散的な電荷連続則の検証は残っています。
`python -m salmonpy.examples.grid_electric_response --help` で、保存配列からの準備例を確認できます。

## TDDFTから密度を取得して電流を渡す

`state.charge_density()` は、現在の時間ステップの電子電荷密度を明示的に取得します。
時間発展・再初期化・再計算を行いません。返す `ChargeDensitySamples` の `step_index`、`time`、
`positions`（bohr）、`charge_density`（負の電子数密度、e/bohr³）、`volume_elements`（bohr³）は独立した読み取り専用データです。
イオン電荷は含みません。固定イオンでは前後の差分に影響しません。

```python
old = state.charge_density()
state.step(E)  # 指定した駆動で1ステップ進める
new = state.charge_density()
J = series.charge_current(new.positions, old.charge_density, new.charge_density,
                          volume_elements=new.volume_elements, dt=state.dt)
```

この操作は対応検証記録を持つ分子PZ・固定イオン・単一ランクの計算プログラムで利用します。
公開ローダーが、実行ファイルのSHA256と `charge_density` の検証記録を確認します。
毎ステップの自動出力ではありません。要求時だけ、ワーカーの私有作業フォルダーでバイナリ転送用ファイルを更新します。
大規模計算ではこの明示的なI/Oの量を考慮してください。恒久的な密度保存は利用者が選択します。
直接保存・復元の電子状態には元々密度が含まれ、復元後も同じ時間ステップの密度を取得できます。

FDTDの `field.advance_edges(indices, current)` は、E成分に既に配置した電流密度を直接使います。
`indices[N,3]` は全体格子のゼロ始まりの添字、`current[N,3]` は原子単位で、体積や負号を追加しません。
格子ごとの所有ランクだけが値を消費し、ランクごとにリスト長を変えられます。重複行は加算されます。
現在この入口は監査付き（`checked`）モードで、既存の異常復帰・E→LD→H更新を使います。
成分格子と異なる配置のJを渡したり、マクロ点の補間を二重に適用したりしないでください。
`electric_edge_geometry()` が全体格子数、E成分ごとの原点（bohr）、格子間隔を返します。
`electric_edges(indices)` は指定したE成分を読み、FDTD群内だけで所有データを合算します。全体MPIへの通信はありません。
電流はこれから更新する区間の平均値です。結合方法によっては試行・巻き戻し・自己真空参照が必要で、これらの入口だけでは自己無撞着な双方向計算にはなりません。

C₂H₂の18,432点で、密度取得の有無による観測量の差は0、ステップ2で保存・終了した電子状態の再開後も密度と観測量が完全一致しました。
その密度変化をFortranで電流へ変換し、1／2ランクのFDTDで放射を確認しました。
FDTDへの直接注入は従来の分配と一致し、NaN拒否・巻き戻し再試行・別プロセスの直接ディスク再開でも配列が一致しています。
これらは接続試験で、物理的収束や完全な空間依存の双方向結合を検証したものではありません。

## 逐次グリッド入力と時間ステップ番号

対応検証記録を持つ実行ファイルでは、`state.step_grid_field(window, step_index=n)` でE/Hを逐次渡せます。
`window` は2フレームの `GridFieldSeries` です。相対時刻は `[0,dt]`、`time_origin=n*dt` とし、
E/HはFDTDのステップ **nとn+1** を表します。分子はこの呼び出しで **n−1からn** へ進みます。
格子の原点・間隔・点数、分子位置、展開中心は準備した入力と同じにします。

```python
window = GridFieldSeries(origin, spacing, [0., dt],
                         [E_n, E_next], [H_n, H_next],
                         molecule_position=position, center=center,
                         time_origin=n*dt)
state.step_grid_field(window, step_index=n)
```

Fortranが対象ステップ、固定した幾何条件とdt、両フレームの全値を検査してから反映します。
不正な幾何条件、NaN/Inf、不完全な転送は、元の数値状態を保って拒否します。
Eのみなら電子格子上の積分経路を覆います。途中からHを有効にする場合は、差分用の周辺格子と
非局所擬ポテンシャルの積分経路も必要です。不足する場合はHを有効にする前に拒否し、Eのみの計算は継続できます。
Hの入力には軌道磁気結合が含まれますが、スピン/SOCの磁気応答を追加するものではありません。

初期化には `salmonpy-near-eh.dat` が引き続き必要です。未来の場が不明なら、同じ幾何条件で
`max_steps+2` フレームのゼロ場などを明示的に準備し、各ステップで必要な2フレームを上書きします。
現在の保持配列は `max_steps×格子点数` に比例します。全領域でなく分子を覆うパッチを用い、
長時間・大領域の計算ではメモリと転送ファイルのI/Oを見積もってください。

既存サンプル `isolated_response` の `--grid-series FIELDS.npz` でも使えます。
FIELDS.npzには `origin, spacing, times, E, H, molecule_position, center` を入れ、
フレーム0からN+1を用意します。従来のE_NPZはN個の追加一様電場とdtです。
再開には `--direct-checkpoint` / `--resume-direct` を使い、Eだけを再計算する保存方式は使いません。
ステップ番号はこの実行での番号です。原データの開始時刻を記録しても、GSが過去の照射を受けた状態にはなりません。

## 空間依存電場による双方向結合の検証範囲

`tools/test_grid_feedback_mpi.py` は、公開の密度取得・逐次入力・直接電流注入を組み合わせた検証用ドライバーです。
本計算と、対象分子だけが放射する自己真空参照のFDTDを別の群で進め、両者の差をFortranで計算します。
分子の密度変化から分極電流を作り、各FDTD群の所有ランクへ必要な電流だけを配ります。
全グリッドを全MPIランクへブロードキャストしません。

最初にFDTDだけをゼロ電流で1ステップ進めます。その後は次の対応です。

| 分子の更新 | 渡すFDTDのE | 受理後のFDTDステップ |
|---|---|---|
| 0→1 | 1、2 | 2 |
| 1→2 | 2、3 | 3 |
| 95→96 | 96、97 | 97 |

試行した場から分子応答を求め、場を戻して電流を注入した後の選択電場が一致することを確認します。
この検証では電流の分配先を真空に置き、Hを分子へ渡さない **電気応答のみ** の模型を使います。
磁気応答をFDTDへ返す完全なE/H結合ではありません。

C₂H₂を格子点からずれた位置に置き、外側4／6ランクで確認しました。LDを含む192分子ステップでは、
環境から戻る場の差は約4.0e-7 au、予測場と電流注入後の選択場の差は約1.4e-17 auでした。
途中拒否後の再試行と、分子96／FDTD97で保存して別MPIプロセスへ直接再開した結果は、
連続実行の配列・通常の出力ファイルと完全一致し、復元時の追加時間発展は0回でした。
これは接続の検証であり、物理的収束や任意の材料配置を保証するものではありません。
このドライバーは検証条件外の予期しないエラーではMPI試験を終了します。汎用的な失敗回復を備えた本番APIとは区別してください。

## スピンなし分子の軌道磁気応答を返す

電気相互作用は長さゲージです。一様電場なら電子のスカラーポテンシャルは `r·E` です。磁場も扱う場合は、磁場を表すベクトルポテンシャルを残し、運動項に `(p + A_mag/c)^2/2` として入れます。ここで `A_mag = -r × ∫₀¹ s B(R+s r) ds` です。これは電気項を速度ゲージへ戻す処理ではありません。一様磁場では `A_mag = B × r/2`、逆向きには `B = curl A_mag` です。後者を入力場との一致条件として使うには、磁場の発散ゼロと補間精度も考慮します。

対応する検証記録のある実行ファイルでは、`state.orbital_magnetic_response()` で、現在の波動関数から軌道磁化密度を読み出せます。時間発展は進めません。分子のスピンなし・固定イオン・直交セル・直接E/H格子入力・電子計算1ランクが対象です。スピン磁化は含みません。

Fortranで、実際に使っているハミルトニアンの磁場微分 `M = -δ〈H〉/δB` を計算します。局所項の常磁性・反磁性成分と、非局所擬ポテンシャルの経路位相の微分を含みます。入力のHは分子を置いた真空のBと同じ値です。原子単位・Gaussian表記で `j_M = c curl M` です。

```python
previous_magnetic = state.orbital_magnetic_response()  # 分子ステップ n-1
# 対応する場を渡して分子をステップ n へ進める
state.step_grid_field(window, step_index=n)
magnetic = state.orbital_magnetic_response()
# electric_current は密度差から得た電気分極電流だけ
current = magnetic.total_current(electric_current, previous=previous_magnetic)
# current を対応するE成分格子の所有ランクへ送り、advance_edges へ一度だけ渡す
```

`magnetization[z,y,x,3]` は入力Hの成分格子上です。`current_density()` は瞬時の磁化電流だけ、`current_density(previous=...)` は隣接ステップの磁化をFortranで平均した区間電流です。`total_current(...,previous=...)` はその区間電流を電気分極電流へ一度だけ加えます。既に磁気成分を含む全電流を渡すと二重計上になります。完全な軌道輸送電流をさらに加えません。

電流への変換はYee配置を要求します。E成分の半格子ずれと、H成分の残り2方向の半格子ずれを指定してください。磁化の支持領域に加え循環電流が入る余白を取り、パッチの正方向端面では磁化をゼロにします。切断したパッチの境界電流を失う変換は拒否します。分子位置は格子点上でなくても構いません。空間微分・経路積分の逆演算・区間平均・電流の加算はFortranです。

磁化は波動関数と適用済みの場から再取得するため、追加の物理的な履歴を保存しません。直接チェックポイントから復元した直後に前ステップの磁化を取得し、次の分子更新後の磁化と組み合わせれば、再計算なしで区間電流を作れます。C₂H₂で電子ステップ2から4の直接再開後に磁化・区間電流・通常観測量の完全一致を確認しました。

`tools/test_orbital_magnetic_response.py` は実ハミルトニアンの磁場微分と保存／再開を検証します。`tools/test_electric_edges_mpi.py response --response-history ...` は生成された電気・磁気・合計電流をFDTDへ渡します。これらは接続検証です。`test_grid_feedback_mpi.py` の既定は電気応答のみで、磁気応答を含む自己真空参照は後述の `--orbital-magnetic` で選びます。一般的な異常復帰保証と物理的収束の検証は別に残ります。

## 磁気応答を含む自己真空参照の双方向結合

`tools/test_grid_feedback_mpi.py` の `--orbital-magnetic` は、直接E/H格子入力と軌道磁気電流の返却を、本計算・自己真空参照の両方へ接続します。スピンなし・固定イオン・PZのC₂H₂、電子計算1ランク、各FDTD群1または2ランクが現在の検証対象です。自己真空参照にも同じ分子電流を注入し、本計算からその参照のE/Hを差し引きます。真空中の自己場を除き、環境で散乱して戻った場を残す方針です。真空の放射反作用まで含めた計算ではありません。

`field.component_field_geometry()` はE/Hそれぞれの成分格子の原点を返し、`field.component_fields(indices)` は配列 `[2,N,3]` にEとHを返します。Hは前後の半ステップの平均で、Eと同じ整数ステップ番号です。成分ごとの空間配置は保ちます。Fortranは観測出力用配列を流用せず、前の生のHを独立に保持します。その情報は巻き戻し・ディスク保存の対象です。E/Hは担当FDTD群内の一度の加算通信で取得し、アプリケーション全体へ格子を複製しません。

更新順は次のとおりです。

1. 最初に両方のFDTDを電流ゼロで1回進めます。
2. 分子ステップn−1からnへの更新には、FDTDのnとn+1のE/Hを使います。まず電流ゼロの試行でn+1を予測します。
3. 分子の前後の状態から電気分極電流と区間磁化電流をFortranで求め、一度だけ加算します。
4. 両方のFDTDを試行前へ戻し、同じ合計電流を、それぞれ担当するE成分格子へ一度ずつ注入します。差し引き後のE/Hが予測と一致することを確認して更新を受け入れます。

192分子ステップの試験では、分子96／場97で両FDTD・LD・分子を直接保存します。別のMPI起動で復元し、保存状態から前ステップの磁化を読み出して継続するため、復元のための時間発展は不要です。E/H・分子観測量と通常出力6ファイルは連続実行と完全一致しました。64ステップ目の意図的な失敗・再試行、および各FDTD群1／2ランクの比較も一致しました。磁気なしの既存試験は、同じ入口で `--orbital-magnetic` を省略します。

```sh
mpiexec -n 6 python tools/test_grid_feedback_mpi.py feedback \
  FIELD_BUILD E_WORKER VALIDATION_JSON EXISTING_GS NEW_OUTPUT \
  --orbital-magnetic --ld --field-ranks 2 --steps 192 --save-step 96
```

再開は同じ入力条件で、新しい出力先と `--resume SAVED_OUTPUT/checkpoint` を指定します。`VALIDATION_JSON` は使用する分子実行ファイルについて、`grid_eh`・`grid_stream`・`charge_density`・`orbital_magnetic` の検証を含む必要があります。インストール済みwheelでは、付属ツールは `salmonpy` の `_runtime_source/tools/` にあります。検証用の固定形状・光源を使う入口であり、任意の構造に対する一般的な実行・異常復帰の保証ではありません。予期しないMPI異常はその試験ジョブを終了します。物理的収束、スピン磁化、真空放射反作用の扱いは別途確認が必要です。
