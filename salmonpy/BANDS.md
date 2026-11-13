# バンド図

`salmonpy.bands.BandStructure.from_salmon('band.dat')` はSALMONの
`theory='dft_band'` のブロック形式を読みます。`k_cartesian` は実際に計算した
Cartesian k (bohr⁻¹)、`energies` は (k, band, spin) のHartree値です。
`plot(reference_ev=..., ticks=[(index,label), ...])` はMatplotlibのAxesを返します。
描画には追加依存のMatplotlibが必要です。対応する配布wheelをインストールした環境で
`python -m pip install matplotlib` を実行します。PyPI公開は未実施です。
ソースからの構築では `python -m pip install ".[plot]"` も使用できますが、
Fortranなどのビルド環境が別途必要です。import時のコンパイルはありません。

```python
from salmonpy.bands import BandStructure
b = BandStructure.from_salmon('band.dat')
ax = b.plot(reference_ev=0.0)  # 明示的に絶対エネルギーのまま表示
ax.figure.savefig('bands.png', dpi=160)
```

## 計算

現在は別途構築したSALMON実行ファイルを使用します。既存のGS入力と擬ポテンシャルを
独立した作業ディレクトリへコピーし、収束済み `data_for_restart` を `seed/` にコピーします。
GSと同じ格子、セル、状態数、汎関数、擬ポテンシャル、k点数を維持し、次を変更します。

```
&calculation
 theory='dft_band'
/
&control
 sysname='Si', yn_restart='y', directory_read_data='./seed/',
 yn_reset_step_restart='y'
/
&band
 lattice='non', nref_band=24, num_of_segments=4,
 ndiv_segment(1:4)=4,4,4,4,
 kpt(:,1)=0,0,0, kpt(:,2)=.5,0,0,
 kpt(:,3)=.5,.5,0, kpt(:,4)=.5,.5,.5, kpt(:,5)=0,0,0,
 tol_esp_diff=1d-5
/
```

この例は8原子の立方Siセル用。kptはそのセルの逆格子基底の分数座標です。
原始FCCセルの経路とは異なり、折り畳まれたバンドになります。
通常のSCFと異なり、バンド反復では読み込んだ密度からのポテンシャルを固定します。
終了ログだけでなく、各k点ブロックの必要バンドの収束を確認してください。
`nref_band=24` なら最初の24バンドまでが収束判定対象です。

今回の本体修正では、最後の頂点の読み落としと逆格子変換の重複を修正し、
最後の終点を必ず含めました。計算単位に合わせた余剰点は終点の複製です。
4区間×4分割なら17点を描画し、残りの余剰点を除きます。
旧実行ファイルにはこの修正がないため、同じ入力から同じ経路を得られません。
ローダーは旧形式も読めますが、誤った経路で計算済みの固有値を修復することはできません。

`examples/plot_bands.py band.dat bands.png` はこのSi例を描きます。
32電子・非スピン分極の占有16バンドの経路上の最高エネルギーを0 eVにします。
これは全Brillouin域のVBMやFermi準位の自動推定ではありません。
今回の12³格子・粗い経路は接続確認用であり、定量的なギャップの収束計算ではありません。
この入口自体は対称点探索・unfolding・金属の占有判定を行いません。
GS作成はDFT.md、ASEによる経路生成はASE.mdを参照してください。

## 初期状態から計算する

`BandInitialState` は読み込み時の自己完結したGSディレクトリをコピーして保持します。
入力はそのGSを生成した元のDFT入力です。再開データは `data_for_restart/`、
擬ポテンシャルなどは同じ束の中に相対パスで置いてください。
バイナリの生成条件と入力の物理的一致は利用者が確認します。ハッシュ検査だけでは
異なる計算由来のGSを判別できません。

```python
from salmonpy.bands import BandInitialState, BandPath
from salmonpy.fdtd import SalmonExecutable

# この環境の監査済み実行ファイルと記録済みSHA256を利用者が指定する。
engine = SalmonExecutable(executable_path, validated_sha256)
path = BandPath(
    vertices=[[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [.5, .5, .5], [0, 0, 0]],
    intervals=[4, 4, 4, 4], labels=['Γ', 'X', 'M', 'R', 'Γ'],
)
with BandInitialState.from_salmon('gs-bundle', 'gs.inp', executable=engine) as gs:
    result = gs.run(path, 'new-band-run', bands=24,
                    tolerance=1e-5, max_iterations=300, threads=2)
# 任意の明示的なエネルギー基準。0.0なら絶対エネルギー表示。
ax = result.plot(reference_ev=0.0)
ax.figure.savefig('bands.png', dpi=160)
```

- `vertices` は当該セルの逆格子の分数座標、`intervals` は各区間の分割数。
  最大10区間。始点と最後の終点を含め、描画用の余剰点は自動的に取り除きます。
- `bands` は必要な収束バンド数。`tolerance` は固有値変化の許容値（Hartree）。
  これは離散化・k点数・擬ポテンシャル由来の誤差限界ではありません。
- 保存先は初期束の外にある新規ディレクトリ。`native/` にコピーした初期ファイル、
  実際の入力、ログ、未加工の出力を、`result/` に配列と単位・条件を保存します。
  追加の図は `result.plot` で生成します。図の失敗で完了済みの計算を失いません。
- 同じ `gs` から異なる `path` を何回でも計算できます。前の計算の状態は再利用しません。
  `with` の終了で保持した初期コピーを解放します。計算結果の保存先は残ります。
- 正常終了に加え、全k点ブロックの `bands` 本が収束し、実際のk座標が指定と一致した
  場合だけ結果を返します。未収束・中断・旧実行ファイルで判定情報がない場合はエラー。
  失敗時も `band-run.json` と `native/stdout.log` を残します。

この入口の対応範囲は固定イオン・周期系・非スピン分極のPZ／Libxc-PZ／Libxc-PBEです。
Libxcには対応する別途構築のSALMONが必要です。TB-mBJ、SOC、DC、構造最適化、
MPI起動はこの入口では未対応です。元入力の解析出力と並列配置は引き継がず、
単一MPIランクのバンド計算として実行します。ソルバーの反復設定は引き継ぎ、
`nscf` は `max_iterations` に置き換えます。この回数は本体の規約どおり
全k点ブロックにわたる累積反復数の上限であり、ブロックごとの回数ではありません。

本体には `band_status.dat` の `SALMONPY_BANDS_V1` 出力を含む今回の修正が必要です。
旧wheelの対応OS一覧は、この外部実行ファイルの動作保証ではありません。
実GS結合テストはmacOS arm64で実施しました。斜交セルの検証は本体の経路生成・
k設定処理の単体テストであり、斜交セルのGS全体を検証したものではありません。

## 任意のk点列

`ExplicitKPoints(kpoints)` はASEに依存しない、分数座標 `(nk,3)` の読み取り専用コピーです。
`BandInitialState.run` に `BandPath` の代わりに渡せます。本体には
`&band/file_band_kpoints` を使い、点数に続く3列の分数座標を読み込ませます。
再補間は行わず、末尾の計算用余剰点だけを追加・除去します。
この列自体は経路の切れ目やラベルを持たないので、`BandRun.plot` では描画せず、
ASEなどの経路情報を持つ側で描画します。ASE接続は独立した [ASE.md](ASE.md) を参照。
