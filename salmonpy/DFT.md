# 構造から基底状態を計算する

GS計算本体の入口はASEに依存しない `salmonpy.dft` です。
`salmonpy.ase.dft_from_atoms` はASE構造をこのモデルへコピーする薄い変換です。
ASEを使わない場合は `DFT` に原子番号、行ベクトルのセル（bohr）、分数座標を直接渡せます。

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
電子数は擬ポテンシャルの価電子数に合わせて明示します。原子番号の総和からは推定しません。
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
- `model.run` は毎回新規GS計算。既存の保存先を上書きしません。再開・warm startのAPIではありません。

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
        J = state.step([1e-4, 0., 0.])  # 新しい端点のA/c、原子単位
```

`ground-state.json` に設定・入力・実行ファイルの識別情報・収束結果・nativeファイルの
ハッシュを記録します。`native/` に入力・擬ポテンシャル・出力・再開データを保持します。
ディレクトリ全体を移動して `open` できます。native内の変更・追加・欠落ファイルや
シンボリックリンクを拒否します。自分のメモや後処理結果はnativeの外へ置いてください。
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
TB-mBJの新規GS作成、分子・孤立系、応力、MPI起動はこの入口には含みません。既存のTB-mBJ GSを使うTDDFT入口は引き続き利用できます。
