# ASEとの構造・バンド・光応答連携

`salmonpy.ase` は任意の連携モジュールです。ASEとの構造・経路の照合と結果変換を
担当し、計算の実行・初期状態の保護はASEに依存しない `salmonpy.bands` が担当します。
LD・SBE・TDDFT・FDTDの通常利用にはASEもf90nmlも不要です。
`salmonpy.ase` のimportだけでもこれらの依存ライブラリは読み込みません。

対応するSALMONpy wheelをインストールした環境で、追加依存を導入します。

```
python -m pip install "ase>=3.26,<4" "f90nml>=1.4,<2"
```

ソースから構築する場合は `python -m pip install ".[ase]"`。
この操作には別途ビルド環境が必要です。PyPI公開済みとはしていません。
本体の `&band/file_band_kpoints` 対応と `SALMONPY_BANDS_V1` 報告を含む、
監査済みの外部SALMON実行ファイルを用意してください。

## 使い方

```python
from salmonpy.bands import BandInitialState
from salmonpy.fdtd import SalmonExecutable
from salmonpy.ase import atoms_from_initial, run_bands

# 利用環境で監査した実行ファイルとSHA256を指定。
engine = SalmonExecutable(executable_path, validated_sha256)
with BandInitialState.from_salmon('gs-bundle', 'gs.inp', executable=engine) as gs:
    atoms = atoms_from_initial(gs)  # 元のGS入力から新しいASE Atomsを作る
    # 8原子の立方Siセルの例。別のセルではそのセルに適した経路を選ぶ。
    path = atoms.cell.bandpath('GX,MR', npoints=17)
    result = run_bands(gs, atoms, path, 'new-ase-bands',
                       bands=24, reference_ev=0.0, threads=2)
result.plot(filename='bands.png', show=False)
# ASE標準の形式を直接利用可能。
ase_bands = result.band_structure
```

GUIのない環境では描画前にMatplotlibの `Agg` backendを選びます。
例の0 eVは絶対エネルギー基準です。`reference_ev` は明示的に指定し、
Fermi準位や価電子帯上端を推定して設定することはしません。

利用者がASEで用意した `atoms` も渡せます。その場合、元のGS入力と同じセル、
原子順・元素、原子位置、3次元周期境界であることを確認します。位置は整数倍の
格子並進を同一視します。セル成分と原子位置の照合許容差は1e-8 Åで、計算後のk点照合にも同じセル精度を反映します。
セルの変換、原始胞化、原子の並べ替えは自動では行いません。
セル・位置を変えた系は、それに対応するGSを用意してください。

## 経路と結果

- ASEの `path.kpts` の順序・値をそのまま本体へ渡します。再補間しません。
- `GX,MR` のコンマによる区間の切れ目、対称点名、セルをASE結果に保持します。
  計算単位に合わせた末尾の余剰点だけを確認して除きます。
- `ase-bands.json` はASEの `BandStructure` JSON、`ase-input.json` は構造・経路・
  エネルギー基準・換算係数の記録です。`native/` と `result/` も保存します。
- `to_ase_band_structure(bands, path, reference_ev=...)` は既存のSALMONpy結果を
  変換するだけの入口です。実際のCartesian kとASE側のセル・k点が一致しなければ拒否します。
- ASEの可変な構造・経路・エネルギー配列をコピーし、呼び出し側のオブジェクトを変更しません。

ASEではセルは行ベクトル（Å）、k点は逆格子の分数座標、エネルギー配列は
`(spin, k, band)` のeVです。SALMONpyの `(k, band, spin)` Hartreeから変換します。
逆格子の2πと転置は連携側で明示的に扱います。長さはSALMON入力の換算値
1 bohr = 0.52917721067 Å、エネルギーは既存SALMONpyの
1 Hartree = 27.211386245988 eVを使い、記録します。ASE既定のCODATA値に暗黙置換しません。

## 現在の範囲

- 元入力からの構造読み込み：`al`／`al_vec1..3`／`dl*num_rgrid`、原子単位または
  `A_eV_fs`、inlineまたは束内ファイルのCartesian／分数原子座標。明示した原子配置が必要です。
- 固定イオン・周期・非スピン分極。PZ／Libxc-PZ／Libxc-PBEは既存バンド入口に準拠。
- バイナリGS自体の生成条件を解析して入力と照合する機能ではありません。
  GSと元入力・擬ポテンシャルの一致は利用者が確認します。
- PZ／Libxc-PZ／Libxc-PBE・KY／NC-UPFのASE Calculator・力・固定セル最適化は独立した [ASE_CALCULATOR.md](ASE_CALCULATOR.md) に対応。応力・MPI起動は未実装。
  構造からのGS計算と全エネルギー取得は別途 [DFT.md](DFT.md) に対応。
- ASE 3.29.0とmacOS arm64上で実GSの結合を検証。斜交セルは変換と本体経路処理の
  単体検証であり、斜交セルの実GS全体の検証ではありません。

仕様の参照先：
https://docs.ase-lib.org/ase/dft/kpoints.html
https://docs.ase-lib.org/_modules/ase/spectrum/band_structure.html

## 構造からGSを用意する

`dft_from_atoms` でASE非依存のDFTモデルを作れます。詳細は [DFT.md](DFT.md)。

## 中性系の電子数と光応答

`dft_from_atoms` と `SalmonCalculator` の `electrons` は省略できます。
KYの先頭行の第4項、またはNC-UPFの `PP_HEADER/z_valence` を読み、原子数を掛けて
中性系の電子数を決めます。`Pseudopotential.valence_electrons` でも参照できます。
明示した電子数がこの合計と異なる場合、計算開始前に拒否します。読めないヘッダー、
非整数の価電子数、奇数電子の閉殻系も拒否し、原子番号からの推測や丸めは行いません。
KYには元素名の自己確認情報がないため、`atomic_number` と実際のファイルの対応、
生成汎関数は引き続き利用者が確認します。

`states`・`kmesh` は必須です。格子は`grid`または`grid_spacing`＋`grid_spacing_unit`で指定します（[構造・格子の準備](STRUCTURE.md)）。資源予算と目的から候補を選び、格子・k点・
非占有状態数・時間刻みを変えて対象物理量の収束を調べてください。中性電子数が正しい
ことと、物性が収束していることは別です。以下の16³格子はSi接続試験の条件です。

```python
from salmonpy.ase import tddft_from_atoms
from salmonpy.crystal import respond
from salmonpy.dft import SCFSettings

# atoms、pseudo、監査済みengine、A/c波形driveを先に用意する。
with tddft_from_atoms(
    atoms, 'new-gs', executable=engine, pseudopotentials=[pseudo],
    states=32, grid=(16, 16, 16), kmesh=(2, 2, 2),
    scf=SCFSettings(max_iterations=500, tolerance=1e-10),
    max_steps=20, threads=1,
) as initial:
    response = respond(initial, drive, output_directory='new-current')
    # 別の計算は独立状態から始まる。stepも同じA→j契約。
    with initial.new_state(dt=0.02) as state:
        j = state.step([0.0, 0.0, 1e-4])
```

`tddft_from_atoms` はGS計算と既存の初期状態ローダーを呼び、
`TDDFTInitialState` を返します。`with` で初期スナップショットを解放しても、
`new-gs` の基底状態とログは残ります。再利用には `GroundStateResult.open` を使います。
実行ファイルにはDFTとTDDFT workerの両方が必要です。実行時の自動取得・コンパイルはありません。
周期・中性・非スピン分極・固定イオン・一MPIランク、PZ/Libxc-PZ/Libxc-PBEが対象です。
構造変換 `dft_from_atoms(..., xc="tbmbj", cval=None)` はTB-mBJの直列GS生成に対応します。Calculatorの力・構造最適化にはTB-mBJを追加していません。孤立系、non-collinear、イオン運動、応力は今回追加しません。Calculatorの明示的なwarm startはASE_CALCULATOR.mdを参照してください。
Aは原子単位のA/c、戻り値は**電荷電流密度j**（SALMON物質電流Jの負号）です。

実行可能な例は [examples/ase_tddft.py](examples/ase_tddft.py)。
[Si入門例](SI_INTRO.md)で取得・展開した計算プログラムを再利用できます。

```sh
python -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY --verify-independence
```

8原子Si・16³格子・2³ k点・20ステップを1スレッドで計算します。
`gs/`、`rt/`、`current.npz`（time, A, current）、`report.json` を保存し、
許容誤差付きの既存Si期待値と比較します。独立性検査は逆波形と元波形を順番に再計算します。
短い接続試験であり、収束した光学スペクトルの例ではありません。
資源計測は別扱いで、サンプル自体はピークメモリを測定・強制制限しません。

結晶TDDFTの標準は **transverse**（`trans_longi='tr'`）です。Longitudinalは明示的に場の部品を選ぶ場合だけ使います。孤立系・non-collinear・Ehrenfestは初版の目標に含まれますが、現時点のこの入口では未対応です。

### SOC入力の準備

`dft_from_atoms(..., spin_orbit=True)` は完全相対論NC-UPFとスピノル状態数の検査を行い、SOC用DFTモデルを作ります。`prepare_dft` はこの入力を配置するだけで計算を始めません。`prepare_wse2_soc` のサンプルは明示した構造ファイルを使います。これはCalculatorのSOC力、バンド、`tddft_from_atoms` のSOC対応を完了したという意味ではありません。現在の小規模検査は[Ni](NICKEL_SOC.md)を参照してください。


SOC入力には全元素の完全相対論NC-UPFを用意します。状態数は省略時に電子数の2倍で、電子温度は`electronic_temperature_k=300`のように指定できます。ASEの初期磁気モーメント・初期電荷の指定は引き継ぎません。SOC入力を作れることと、ASE CalculatorのSOC力や`tddft_from_atoms`のSOC対応は別です。
擬ポテンシャルは構造変換時・計算時には自動取得しません。必要なら計算前に[専用取得手順](PSEUDOPOTENTIALS.md)で用意します。
