# 構造から基底状態を計算する

GS計算本体の入口はASEに依存しない `salmonpy.dft` です。
`salmonpy.ase.dft_from_atoms` はASE構造をこのモデルへコピーする薄い変換です。
ASEを使わない場合は `DFT` に原子番号、行ベクトルのセル（bohr）、分数座標を直接渡せます。

## SOCと金属の明示的な入力

`DFT(..., spin_orbit=True)` は完全相対論NC-UPFを使い、非共線スピノルの入力を生成します。状態数はスピノル本数で、省略時は電子数の2倍です。スカラー擬ポテンシャルとの暗黙の置換はしません。正の電子温度は `electronic_temperature_k=300`（K）または `electronic_temperature`（Hartree）の一方でモデルごとに指定します。

`SCFSettings` では、反復数・許容差に加え `mixing_method`、`mixing_rate`、`initial_wavefunction`、`random_seed` を明示できます。混合率の既定値は0.3です。方法は省略時にBroydenを使い、明示した係数を優先します。`simple_dm` は非共線密度行列の混合です。

[NiのSOC検証](NICKEL_SOC.md)には、セミコア付きUPF・2×2×2・300 Kの短いSCF試験、SCFの独立性・別プロセス再開、TDDFTの独立性・直接復元・同じMPI分割での別プロセス再開を含みます。収束した物性、SOCの力・HSE・TDDFTのMPI分割変更は未検証です。

## ASEからの使用例

```python
from ase.build import bulk
from salmonpy.ase import dft_from_atoms
from salmonpy.dft import Pseudopotential, SCFSettings
from salmonpy.fdtd import SalmonExecutable

atoms = bulk('Si', 'diamond', a=5.43, cubic=True)
pseudo = Pseudopotential.from_file('Si_rps.dat', atomic_number=14, local_channel=2)
model = dft_from_atoms(
    atoms, pseudopotentials=[pseudo], electrons=32, states=32,
    grid=(16, 16, 16), kmesh=(2, 2, 2), xc='pz',
)
engine = SalmonExecutable(executable_path, validated_sha256)
gs = model.run(engine, 'new-gs',
               scf=SCFSettings(max_iterations=300, tolerance=1e-9), threads=2)
print(gs.energy_hartree, gs.iterations)
```

この格子・k点は接続用の例であり、任意の物性に対する収束を保証しません。
ASE入口では電子数を省略すると擬ポテンシャルの価電子数から中性系の値を決定し、明示値との一致を確認します。直接DFTを作る場合は電子数を明示します。原子番号の総和からは推定しません。
対応形式はKYの `_rps.dat` とスカラー・ノルム保存UPFです。UPFのPAW、ultrasoft、SOCは拒否します。
擬ポテンシャルの元素・価電子数・生成汎関数と計算条件の整合は利用者が確認してください。
擬ポテンシャルをダウンロードする処理はありません。

ASEの初期電荷・磁気モーメント指定はこの入口では受け付けません。
閉殻・非スピン分極として偶数の電子数と十分な軌道数を指定してください。
セルの変換や原子の並べ替えは行わず、セル・座標をASEのÅから原子単位へ変換します。

## 設定の寿命

- モデル作成時に固定：構造、擬ポテンシャルのバイト列、電子数、軌道数、格子、k点、汎関数。
  元のAtomsや擬ポテンシャルファイルを後から変更しても、モデルは変わりません。
- 計算ごとに指定：SCFの累積反復上限、`rho_dne` 収束許容値、実行ファイル、新規保存先、スレッド数、タイムアウト。
- `model.run` は既定で新規GS計算。`initial_guess=completed_gs` を明示すると、検証済みの波動関数を初期値としてSCFを再実行します。既存の保存先を上書きしません。

収束判定は本体の既存フラグを使います。固有状態の更新・混合・計算順は変更しません。
正常終了だけでは成功とみなさず、収束報告と再開データを確認します。
`iterations` は本体のMiter（更新回数）で、収束を検出した次のループ番号とは異なる場合があります。
失敗時も `ground-state.json` の失敗状態と `native/stdout.log` を残します。

## 保存・再利用

```python
from salmonpy.dft import GroundStateResult
from salmonpy.bands import BandPath

gs = GroundStateResult.open('new-gs')
with gs.band_initial_state(executable=engine) as initial:
    result = initial.run(
        BandPath([[0, 0, 0], [.5, 0, 0]], [10], ['Γ', 'X']),
        'new-bands', bands=24,
    )

with gs.tddft_initial_state(executable=engine, max_steps=1000) as initial:
    with initial.new_state(dt=0.02) as state:
        j = state.step([1e-4, 0., 0.])  # 新しい端点のA/c、原子単位
```

`ground-state.json` に設定・入力・実行ファイルの識別情報・収束結果・計算処理ファイルの
ハッシュを記録します。`native/` に入力・擬ポテンシャル・出力・再開データを保持します。
ディレクトリ全体を移動して `open` できます。計算処理内の変更・追加・欠落ファイルや
シンボリックリンクを拒否します。自分のメモや後処理結果は計算処理の外へ置いてください。
ハッシュは改変検出用であり、第三者による改ざんに対する電子署名ではありません。

初期状態への変換も独立コピーを作り、保存したGSを変更しません。
TDDFT用入力は元の物理条件を保持して生成します。dtは `new_state`、A/cは各stepで指定します。
現在のTDDFT初期化は中点伝播・予測子修正子なしです。

## 環境と対応範囲

任意のASE依存の導入方法はASE.mdを参照してください。外部SALMON実行ファイルは別途必要で、
`SALMONPY_GS_V1` 収束報告と `SALMONPY_FORCES_V1` 力報告を出す本体を監査して指定します。
バンド／TDDFTへ接続する実行ファイルにも各入口の対応機能が必要です。
バイナリ再開形式の互換性は任意の異なるビルド間で保証しません。

固定イオン・3次元周期・非スピン分極、PZ／Libxc-PZ／Libxc-PBE、CPU単一MPIランクが対象。
Libxcには対応ビルドが必要です。今回の実GS一貫テストはmacOS arm64・PZ Siです。
本体の力は `gs.forces_hartree_per_bohr` で取得できます。力の格子収束は別途確認します。
ASE Calculatorと固定セル最適化への接続は独立したASE_CALCULATOR.mdに従ってください。
TB-mBJの新規GSは直列の `DFT(xc="tbmbj", cval=None).run(...)` で作成できます。`cval=None` は自動係数、正の固定係数（1e-10より大きい値）も指定できます。分子・孤立系、応力、MPI起動はこの入口には含みません。計算処理並列SCFの対応範囲は変更していません。

## 波動関数を初期値として再利用する

`model.run(engine, NEW_DIRECTORY, initial_guess=previous_gs)` は、セル・格子・k点・元素順・電子数・軌道数・汎関数・擬ポテンシャルの内容と局所チャンネル・実行ファイルのSHA256が一致する場合に利用できます。原子位置のみ変更可能です。不一致やファイル破損は本体を起動する前に拒否します。新しい保存先は元GSの外に置きます。

`info.bin`、`occupation.bin`、`wfn.bin` を専用領域へコピーし、`read_gs_restart_data='wfn'` と反復数リセットを使います。旧原子位置や混合履歴は読み込みません。SCFは新しい原子位置で密度を再計算して始め、元のGSは変更しません。出力の `ground-state.json/initial_guess` に元GSとハッシュを記録します。warm startで必ず高速化するとは限らず、別の局所解への収束や物理量の収束確認は別途必要です。

## HSEの計算処理 GS入口

DFT(..., xc="hse06")はHSEのNativeGSBackend session.run用に使用できます。
HSE06は全k点・k点方向のMPI分割が対象で、計算機能を宣言した新しいビルドが必要です。
詳細はNATIVE_GS.mdを参照してください。バッチDFT.runのHSE出力は未接続です。
HSE入力ではPBE予備SCFを無効にし、指定反復をHSEで実行します。未収束時は確認待ちです。

## 再開用配列を保存せず、エネルギー・力を得る

バッチDFT.runは既定で従来どおり再開データを保存します。エネルギーと力だけ必要なら
実行時にsave_restart=Falseを指定します（モデル設定は変わりません）。

```python
gs = calculation.run(engine, "new-output", scf=SCFSettings(max_iterations=4),
                     save_restart=False)
print(gs.energy_hartree, gs.forces_hartree_per_bohr, gs.restart_saved)
```

対象はバッチ入口のPZ／Libxc-PZ／PBEです。4反復は接続確認用で、物理的な収束を保証しません。
未収束時は結果を保持して確認待ちにします。restart_savedは実際に再開用3ファイルが揃ったかを
示し、省略した場合はFalseです。エネルギー・力の読み込みは可能ですが、バンド初期化、
TDDFT初期化、warm start、GSアーカイブ出力には再開用配列を保存した結果が必要です。
save_restartを記録していない従来の結果は保存ありとして扱います。
計算プログラムが保存省略を無視しても計算失敗にはせず、警告して実際の保存状態を記録します。

高速版ソースもgs_status.dat（SALMONPY_GS_V1、Hartree）と
gs_forces.dat（SALMONPY_FORCES_V1、Hartree/Bohr）を出します。SCFや力の計算順は変えません。
HSEの計算処理入口とCLIもこのレポートを出しますが、DFT.runのHSE選択は
HSEの来歴確認が未接続のため、引き続きNativeGSBackendを使います。

Siサンプルは--no-restartと--max-iterationsに対応します。
`python -m salmonpy.examples.dft_si EXECUTABLE SHA256 PSEUDO NEW_OUTPUT --no-restart --max-iterations 4`

### TB-mBJの直列GS入口

`DFT(..., xc="tbmbj", cval=None)` または `dft_from_atoms(..., xc="tbmbj", cval=1.2)` を使います。JSONの物質レシピでも `xc` と任意の `cval` を指定できます。係数はGS記録、warm startの整合性検査、TDDFT入力へ引き継ぎます。TB-mBJはモデルポテンシャルのため、全エネルギーと力は未検証の診断値として警告し、計算は続行します。ASE Calculatorによる最適化の対象にはしません。サンプルは `python -m salmonpy.examples.dft_tbmbj --help` で確認できます。今回の確認は合成密度の計算処理XC評価とインターフェース検査であり、実材料の物理精度保証ではありません。


管理ノードで入力保存だけを行う場合はpreparation.prepare_dftを使います。実行とは独立した入口で、自動的に計算へ進みません。initial_guessを明示すれば、同じwarm start検査と入力生成を使って3ファイルの初期推定も準備できます。日英のPREPARATION.md／PREPARATION.en.mdを参照してください。


SCF残差が振動して下がらない場合は、混合率を下げて試します。Broydenでは`alpha_mb`、Pulayでは`beta_p`が対象です。Pythonの指定例と保持状態からの継続は[NATIVE_GS.mdの混合率変更手順](NATIVE_GS.md#残差が振動するときは混合率を下げる)を参照してください。反復数を増やすだけで押し切らず、残差の変化を確認します。

## SCF未収束時の確認待ち（2026-11-21）

`DFT.run()` と `NativeGSSession.run()` は、反復上限で未収束なら警告し、
`SCFReviewRequired` を送出して後続計算への自動継続を止めます。
失敗した実行としてデータを消すことはしません。バッチの保存済み結果／計算処理の生存状態は
例外の `result` に保持し、`report_directory` に残差CSV・任意のPNG・report.jsonを保存します。
MatplotlibがなくてもCSVは利用できます。図の生成失敗と計算失敗を混同しません。
レポートの保存に失敗しても確認待ちは維持し、元の出力を確認するよう警告します。
通常のGS保存設定を維持し、`save_restart=False` なら新たなGS再開データは保存しません。

```python
from salmonpy.errors import SCFReviewRequired
try:
    gs = model.run(engine, 'new-gs-output')
except SCFReviewRequired as review:
    gs = review.result
    print('Review SCF residuals:', review.report_directory)
    # ここで終了し、残差の推移・閾値・入力条件を利用者が確認する。
    # このexcept節からTDDFTを自動で開始しない。
```

残差が下がっているか、停滞／振動しているかを確認し、反復上限・混合条件などを見直します。
確認後に未収束状態を使うと利用者が明示した場合だけ `allow_unconverged=True` を指定できます。
AIが自動的にこの指定で再実行してはいけません。既存結果の `GroundStateResult.open()` や
残差の読取・保存・warm startは引き続き可能です。既知の未収束GSからの
`band_initial_state()`／`tddft_initial_state()` にも同じ確認待ちを適用します。

計算処理の `run(..., steps=K)` と `continue_scf(steps=M)` で意図的に反復境界へ
止めた `paused=True` の状態は正常な途中状態です。全反復上限に達した未収束終了だけが
確認待ちになります。確認待ち後もFortran/Cで保持する計算状態は閉じず、保存や診断に利用できます。
ランク間で `allow_unconverged` の指定を揃えます。計算の更新順・収束閾値は変更しません。
この方針は以前の「未収束でも警告のみで自動継続」を置き換えます。

## 混合率の既定値

SALMONpyの混合率は**0.3**を既定とします。Pythonの`SCFSettings()`、
Fortran入力パーサー、状態ごとの設定で揃えています。
Broydenは`alpha_mb=0.3`、Pulayは`beta_p=0.3`、simple系は`mixrate=0.3`です。
明示的な`mixing_rate`や対応するSALMON入力係数があれば、その指定を使います。
以前のBroyden／Pulayの0.75、simple系の0.5から下げたため、既定値を使うSCFの
反復経路・反復数は変わり得ます。収束許容差は変更しません。
振動が残る場合は0.1などを試し、残差の推移を確認します。任意の自動低下は既定では無効です。

### 磁性体のSCF閾値

磁性体では `SCFSettings.for_magnetic(max_iterations=4000)` を明示的に選ぶと、残差の閾値は `1e-5` になります。`tolerance=` で変更できます。一般の `SCFSettings()` は `1e-9` のままです。SOCの有無だけでは磁性体と判断しません。

残差だけでなくエネルギーと磁気モーメントの変化も確認してください。振動が続く場合は混合率を段階的に下げます。Niの検証では最終無磁場段階の混合率0.15、観測間のエネルギー変化 `1e-5 Ha`、磁気モーメント成分の最大変化 `1e-4 μB` を目安とします。これは入力・観測間隔に依存する診断で、物理精度の保証ではありません。電荷残差が閾値に達するとSCFは停止するため、観測量が安定していなければ要確認として報告します。
