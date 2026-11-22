# 孤立系TDDFT：E入力と名前付き双極子

salmonpy.isolated_tddftは、PZ・固定イオン・非偏極・スカラー擬ポテンシャル・中点伝播を
明示的な別プロセスで使う入口です。結晶はsalmonpy.tddftのA駆動を維持します。
返すのは電荷双極子で、電流密度ではありません。格子・時間刻み・計算箱の収束は利用者が確認します。

## 計算プログラムと初期状態

```python
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
runtime = load_isolated_runtime('/path/to/salmon', '/path/to/validation/report.json')
with IsolatedTDDFTInitialState.from_salmon(
        '/path/to/bundle', 'rt.inp', runtime=runtime, max_steps=100,
        threads=1, timeout=60) as initial:
    with initial.new_state(dt=.02) as state:
        saved = state.snapshot()
        frame = state.step([0., 0., 1.e-5])
        print(frame.time, frame.dipole, frame.electronic_change, frame.electron_number)
        state.restore(saved)
```

E・時刻は原子単位、双極子はe*bohrです。dipoleはイオンと電子を合わせた全電荷双極子、
electronic_changeは計算処理の電子双極子変化です。配列は読取専用で、物質電流とは区別します。
初期E・時刻はゼロです。new_state(position=...)は場を評価する位置（既定(0,0,0) bohr）を
指定しますが、分子座標の移動・回転は行いません。

計算プログラム読込みは64 KiB以下のローカル検証記録について、E3プロトコル・対応範囲・
比較／独立性／履歴再実行の成功記録と実行ファイルのSHA256を確認します。
ダウンロード・コンパイル・MPI初期化・worker起動はしません。
記録はローカル検証の証拠で、電子署名による認証やソース監査の代わりではありません。
ソースを監査して明示的にビルドしたstandard runtimeについて、ソース内の
 tools/test_isolated_worker.pyを小さな参照fixtureとともに実行して記録を作ります。
形式はsalmonpy.isolated-worker-validation/1です。対象バイナリ・検証範囲に限定され、
任意の物理入力の精度を保証しません。通常のA専用計算プログラムのmanifestでは代用できません。
今回の変更でE対応のビルド済みcatalog配布を公開したわけではありません。

bundleには孤立系RT入力、data_for_restart/info.bin・occupation.bin・wfn.bin、
相対パスで参照する擬ポテンシャル等を揃えます。実行ファイルと互換な計算処理記録を使ってください。
読込みで不変の初期コピーを作り、各状態が独立した作業コピーと子プロセスを持ちます。
元のファイルを変更しても初期コピーは変わりません。未対応の入力条件は計算状態の初期化時に
伝播前に拒否します。周期系GSを変換したり、孤立系GSを自動生成したりする機能ではありません。

## 共通の場入力と復元

```python
from salmonpy.field_input import FieldPacket, FieldValue
# Inside the live state context above.
frame = state.step_fields(FieldPacket(E=FieldValue(
    [0., 0., 1.e-5], time=state.time+state.dt,
    position=(0.,0.,0.), unit='E:au')))
```

Eの時刻・位置・単位を照合し、未提供とゼロを区別します。このbackendはA/Hを消費しません。
step(E)は原子単位を明示的に使う簡単な入口です。Eは計算処理中点伝播のstep indexで使われ、
現行の更新順を保ちます。Aから暗黙に再構成しません。

snapshot()は受理済みEと双極子履歴のメモリ上のトークンです。restore()には同じ初期モデル、
dt・位置・容量・スレッド数・バイナリが必要です。候補プロセスで履歴を再実行し、観測量が一致
（rtol=1e-12、atol=1e-14）してから状態を交換します。復元失敗時は元の状態を保持します。
復元中は新旧2つのworkerが存在するため、2つ分の資源を確保してください。
初期ローダーを閉じても既存状態が参照する初期コピーは残り、復元できます。新規状態の作成は止まります。

step異常時は直前の最後に成功した計算状態の再構成を試み、RuntimeFailureを返します。
回復も失敗すると状態を閉じます。入力検査だけの失敗は状態を進めません。
外側のFieldAdapterはbackend例外で無効化される従来仕様なので、必要なら作り直します。
計算処理波動関数snapshot、MPI全体の復元、同一Fortranプロセス内の
複数contextを提供するものではありません。復元コストは履歴の長さに依存します。

## 出力サンプルとFDTD

```sh
python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE RT_INPUT E_NPZ NEW_OUTPUT
```

E_NPZにはdt,2dt,...のE[steps,3]と正のスカラーdtを原子単位で保存します。
出力archiveには時刻、初期ゼロを含むE、名前付き双極子と単位を保存します。
一方向放射にはframe.dipoleをDipoleMapへ渡します（DIPOLE_RADIATION.md）。
周期系電流の符号反転や、TDDFTの真空セル体積での除算は行いません。
自己場の帰還モデルは別途選択するもので、自動的には有効にしません。

## ディスクへの保存と再開

```python
# 初期bundleの外にある新規ディレクトリを指定します。
state.save_checkpoint('new-isolated-checkpoint')
# 同じ初期データから別途作った状態へ復元できます。
resumed.restore_checkpoint('new-isolated-checkpoint', max_bytes=16*1024*1024)
```

形式は `salmonpy.isolated-tddft-checkpoint/3`、非圧縮Zarr v3の限定形式です。
E0・observables0・observables・E_history・observable_historyを保存します。観測配列は9成分で、
双極子変化3・全電荷双極子3（単位e*bohr）・実測電子数・内部エネルギー（Hartree）・鮮度フラグです。
波動関数や後処理archiveではなく、受理済み履歴の再実行用記録です。
GSと検証済みruntimeは別途必要です。初期ファイル・実行ファイルのハッシュ、dt、容量、
スレッド数、位置、パッケージ版を照合します。初期モデルのオブジェクトは作り直せます。
周期系Aのcheckpointは拒否します。ハッシュ・配置・サイズ検査後に候補workerで
再実行し、失敗時は元の状態を保持します。max_bytesはファイル・復号データの上限であり、
RSS上限ではありません。省略時は設定容量から上限を計算します。保存先の上書きはしません。

応答サンプルの `--checkpoint NEW_DIRECTORY` で最終状態も保存できます。
容量いっぱいまで計算した記録はその終了点までの復元専用です。続きを計算する用途では、
自分のコードで必要な総容量を最初に確保し、同じ容量で復元してください。


### 実測電子数

E2では `frame.electron_number` と `state.electron_number` が占有電子密度の積分値を
返します。入力nelecからの推定ではありません。Fortranの既存双極子観測が計算した値を
再利用するので、観測のための全体通信は増やしません。応答は双極子変化3成分、全電荷
双極子3成分、電子数・内部エネルギー・鮮度の各1成分で9値です。読取専用の配列と有限なスカラーとして返します。
E2の検証記録は `salmonpy.isolated-worker-validation/3`、チェックポイントは
`salmonpy.isolated-tddft-checkpoint/3`。観測配列の先頭6列は双極子、7列目は電子数、8列目はエネルギー、9列目は鮮度です。
旧E1・旧チェックポイントは受け付けません。電子数も履歴再生・復元の照合対象です。
分子FDTD開発試験はFortranで電子数減少率を評価し、目安1e-5を診断結果へ記録します。
目安の超過だけでは計算を止めません。値が非有限ならその試行を拒否します。


### 内部エネルギーと仕事量

`frame.internal_energy` と `state.internal_energy` は外場の相互作用エネルギーを含まない
内部エネルギー（Hartree）です。既存の固有値計算が更新された時刻だけ値を返し、その他は
`None`です。計算を増やすために `out_rt_energy_step` を暗黙に変更しません。

開発用分子FDTD試験では、Fortranで実際に与えた終点Eによる区間一定の仕事
Σ E[n]·(p[n]−p[n−1]) と、中点Eによる仕事を別々に積算します。内部エネルギーとの差は
更新済みの時刻だけ比較します。無電場のエネルギードリフトも別に測定し、主結果から
自動的に差し引きません。仕事の残差は診断であり、有限な計算を閾値で止めません。
NPZの `internal_energy_raw` は鮮度なしでは使用せず、`energy_fresh` または
`energy_observation_steps`・`energy_observations` を使って更新済み値を選んでください。
全電磁エネルギー収支の証明とは区別します。旧E1/E2・旧保存形式は受け付けません。

分子FDTD開発試験の既定は中点E結合です。`step(E)`に渡す区間一定Eを結合側でFortran平均し、単独TDDFTの更新を変えません。終点方式は明示的な比較用です。DIPOLE_RADIATION.mdの時間配置・仕事量・自己場の範囲を参照してください。
