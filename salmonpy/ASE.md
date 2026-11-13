# ASEとの独立したバンド連携

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
  エネルギー基準・換算係数の記録です。`native/` と `result/` も従来どおり保存します。
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
- PZ/KYのASE Calculator・力・固定セル最適化は独立した [ASE_CALCULATOR.md](ASE_CALCULATOR.md) に対応。応力・MPI起動は未実装。
  構造からのGS計算と全エネルギー取得は別途 [DFT.md](DFT.md) に対応。
- ASE 3.29.0とmacOS arm64上で実GSの結合を検証。斜交セルは変換と本体経路処理の
  単体検証であり、斜交セルの実GS全体の検証ではありません。

仕様の参照先：
https://docs.ase-lib.org/ase/dft/kpoints.html
https://docs.ase-lib.org/_modules/ase/spectrum/band_structure.html

## 構造からGSを用意する

`dft_from_atoms` でASE非依存のDFTモデルを作れます。詳細は [DFT.md](DFT.md)。
