# GS配列の可搬アーカイブ（GS/2）

`salmonpy.gs_archive` は、呼び出し側が用意したGS配列をZarr v3互換形式へ保存する部品です。`save_ground_state` と `GroundStateArchive.open` を使います。既存の `GroundStateResult`（`salmonpy.gs/1`）とは別の入口です。SALMONのネイティブバイナリの自動変換、SCFの収束認定、GS/2からのTDDFT初期化はまだありません。

## 合成配列での例

以下は配列の受け渡し例で、実材料の収束した基底状態ではありません。

```python
import numpy as np
from salmonpy.gs_archive import save_ground_state, GroundStateArchive

archive = save_ground_state(
    'new-gs-array-archive',
    wavefunctions=np.ones((1, 2, 2, 2, 2), complex) / 8,
    density=np.full((2, 2, 2), 2 / 64),
    eigenvalues=[[-0.5, 0.2]], occupations=[[2., 0.]],
    kpoints=[[0., 0., 0.]], weights=[1.],
    cell_bohr=np.eye(3) * 4, fractional_positions=[[0., 0., 0.]],
    atomic_numbers=[2], electrons=2., energy_hartree=-1., xc='pz',
    pseudopotentials=[dict(atomic_number=2, filename='He.example',
                          sha256='0' * 64, source='synthetic example; no potential supplied')],
    orbital_block=1,
)
a = GroundStateArchive.open(archive.directory)
u = a.wavefunctions(0, 0, 1)  # 第0 k点の軌道[0,1)、shape=(1,z,y,x)
occupations = a.array('occupations')
```

配列とメタデータの読み出し結果は不変です。書き込み中に入力配列を書き換えたり、読み込み中にアーカイブを編集したりしないでください。保存先は新規ディレクトリのみで、失敗時は作成中のディレクトリを除去します。`manifest.json` を最後に書くため、中断した未完成品は読み込めません。

## 数値と座標の規約

初版は周期系・非スピン分極のみです。波動関数はセル周期部分 `u[k,band,z,y,x]`（complex128）で、Bloch関数は `psi=exp(+i*k_cartesian・r)*u` です。

`cell_bohr` の各行が格子ベクトル a1,a2,a3、逆格子の各行は `B=2*pi*inverse(A).T` です。`kpoints[k]` はBの係数（reduced reciprocal coordinates）であり、Cartesian座標ではありません。実空間は原点0で、重複する周期端点を含まない格子です。

`r(ix,iy,iz)=(ix/nx)*a1+(iy/ny)*a2+(iz/nz)*a3`、`ix=0,...,nx-1`（y,zも同様）。`grid_zyx` と配列の末尾軸は `(nz,ny,nx)` です。原子座標は格子ベクトルの係数、波動関数は bohr^(-3/2)、密度は bohr^(-3)、固有値・全エネルギーは hartree です。

k点重みは非負で和1、占有数は各軌道0〜2、電子数は `sum_k,b weight[k]*occupation[k,b]` と一致させます。密度は非負、セルの体積は正です。これらと全配列の有限値は検査しますが、波動関数の規格化・直交性、密度との整合、SCF収束は認定しません。

## 保存形式

ルートの `manifest.json` はschema `salmonpy.gs/2`、物理メタデータ、配列のレイアウト、全データ・ZarrメタデータのSHA256を持ちます。`zarr.json` はZarr v3のgroupです。

配列は `wavefunction_real`、`wavefunction_imag`、`density`、`eigenvalues`、`occupations`、`kpoints`、`weights`、`cell_bohr`、`fractional_positions`、`atomic_numbers` です。原子番号のみint64、他はfloat64。無圧縮のbytes codec・little endian・C順序です。複素数を実部と虚部に分け、複素型の拡張に依存しません。

波動関数のチャンクは `(1, orbital_block, nz, ny, nx)` で、k点1個・軌道ブロック単位、実空間は分割しません。固有値・占有数は `(1, orbital_block)`。他の配列は全体で1チャンクです。最後の軌道ブロックも完全なチャンク長を持ち、配列外は0で埋めます。`orbital_block` の既定16はI/O粒度で、メモリ上限ではありません。軌道数以下に丸めます。

擬ポテンシャルは元素ごとのファイル名・SHA256・入手元の参照のみで、本体を保存・自動取得しません。同じ原子番号で異なる擬ポテンシャルを混在させる形式は未対応です。

## 読み込みと予算

`open` は全ファイルの寸法・ハッシュ・有限値を検査しますが、全波動関数を一括確保しません。チャンク単位のバッファと物理メタデータ用配列を使います。`wavefunctions(k,start,stop)` は必要な軌道チャンクだけを読み、その都度ハッシュを再検査します。余分なファイル・ディレクトリ、シンボリックリンク、未対応のレイアウトは拒否します。

固定のマシン依存上限は設けません。`open(...,max_bytes=...)` はアーカイブ全ファイルの合計、`array(...,max_bytes=...)` と `wavefunctions(...,max_bytes=...)` は返す配列のバイト数を制限します。これらはプロセスRSS上限ではなく、チャンク・検査・不変コピー用の追加バッファが必要です。SHA256は破損検出であり、配布元の認証ではありません。

## Fortranとの接続範囲

同梱 `include/salmonpy_gs_stream.f90` はlittle-endianのIEEE real64環境で、1次元の実数バッファをチャンクとして読み書きできます。`gs_read_chunk(path,values,status)` と `gs_write_chunk(path,values,status)` を使い、status=0が成功、1が形式/I/O、3が確保失敗、4が非有限値です。読み込み失敗時は呼び出し側の配列を保持し、書き込みは既存ファイルを上書きしません。big-endian環境は拒否します。Fujitsuでのコンパイルは未検証です。

Pythonチャンク `(1,block,nz,ny,nx)` のC順序は、Fortranの `(nx,ny,nz,block)` と同じx最速のバイト順です。呼び出し側で平坦バッファの長さを `nx*ny*nz*block` にして渡します。JSON生成・ハッシュ照合・配列の再構成はこの小モジュールには含めません。

既存SALMONの `wfn.bin` はstream形式の分岐がありますが、実数/複素数・MPI分割・スピン等の配置確認が必要です。`info.bin` と `occupation.bin` は逐次非書式レコードです。実際のGS/1→GS/2変換とGS/2→SALMON初期化には、次段階で監査済みランタイム側の明示的アダプターが必要です。GS/2をそのまま既存のrestartディレクトリへ渡さないでください。

形式の根拠: [Zarr core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html)。NumPy以外の実行時依存は増やしていません。独立したZarr実装との互換性は合成配列で検証します。

保存入力がcomplex128でなければ、NumPy変換で波動関数全体のcomplex128コピーを作る場合があります。1軌道だけの読み出しでも軌道チャンク全体を読むため、max_bytesはその一時バッファを制限しません。Fortran helperはint8の生バイトI/Oを使い、実数のバイト順を変えるランタイム設定の影響を避けます。IEEE real64と8bitのファイル記憶単位を検査します。
