# GS配列の可搬アーカイブ（GS/2）

`salmonpy.gs_archive.save_ground_state`はGS配列をZarr v3互換形式で保存し、`GroundStateArchive.open`で読みます。`GroundStateResult`の`salmonpy.gs/1`とは別の保存形式です。SCFの再実行や収束の再認定は行いません。

| 入力・用途 | 使用する入口と条件 |
|---|---|
| GS配列を保存・読む | `save_ground_state`／`GroundStateArchive.open`。メタデータと配列の整合性を検査 |
| Γ点のGS/1から書き出す・TDDFTへ渡す | `GammaRestartConverter`。下記の直列・単一軌道ファイルなどの条件を確認 |
| 複数k点のGS/1から書き出す・TDDFTへ渡す | `SerialRestartConverter`。対応する変換器・計算プログラム・擬ポテンシャルを明示 |
| 保持中のGSから書き出す・GS/2を読み込む | `NativeGSState.export_archive`／`NativeGSSession.load_archive`。保持GSガイドの通常k点メッシュの条件を確認 |

各入口のスピン・MPI・ファイル形式の制約は以下の節と[NATIVE_GS.md](NATIVE_GS.md)で確認してください。複数k点対応はΓ点専用の変換器にも適用される、という意味ではありません。

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

既存SALMONの `wfn.bin` はstream形式の分岐がありますが、実数/複素数・MPI分割・スピン等の配置確認が必要です。`info.bin` と `occupation.bin` は逐次非書式レコードです。GS/1→GS/2書き出しは下記の厳密なΓ点条件で対応し、複数k点・MPI分割の変換は未対応です。GS/2をそのまま既存のrestartディレクトリへ渡さず、次の専用変換器を使います。

形式の根拠: [Zarr core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html)。NumPy以外の実行時依存は増やしていません。独立したZarr実装との互換性は合成配列で検証します。

保存入力がcomplex128でなければ、NumPy変換で波動関数全体のcomplex128コピーを作る場合があります。1軌道だけの読み出しでも軌道チャンク全体を読むため、max_bytesはその一時バッファを制限しません。Fortran helperはint8の生バイトI/Oを使い、実数のバイト順を変える計算プログラム設定の影響を避けます。IEEE real64と8bitのファイル記憶単位を検査します。

## Γ点GS/2からTDDFTへ

初版は周期系・Γ点1個（k=0、重み1）・非スピン分極・固定イオン・CPU直列計算プログラムです。格子は各軸8点以上、電子数は偶数、XCは既存DFT入口のpz／libxc_pz／libxc_pbeです。Libxcには対応する計算プログラムが必要です。アーカイブの密度・固有値は計算処理再開レコードにはせず、波動関数と占有数からSALMONが密度・ポテンシャルを再構成します。収束状態かどうかは呼び出し側が保証してください。

専用変換器は、使用する監査済み直列計算プログラムのbuild.jsonと残っているソース・CMakeビルド情報を指定して明示的に構築します。初版はGNUコンパイラ・Unix Makefilesのビルド情報が対象です。計算プログラムの実際のFortran_FLAGSとコンパイラを引き継ぎ、対象実行ファイルのSHA256へ固定します。元の波動関数read/write手続きを外部ソースから変更せずに切り出した合成試験で、作成ファイルが一致することを確認します。ソース・全オブジェクトの禁止関数監査を省略しません。

```text
python tools/build_gs_restart_converter.py RUNTIME_BUILD_JSON NEW_CONVERTER_DIRECTORY
```

wheel内のスクリプトは `salmonpy/_runtime_source/tools/build_gs_restart_converter.py` にあります。importや変換時に自動ビルドしません。記録と同じGNUコンパイラ／実行時I/O設定が必要です。MPI版計算プログラム、別のSHA256の実行ファイル、Fujitsuのコンパイラ設定は現版では拒否します。

```python
from salmonpy.gs_native import GammaRestartConverter
from salmonpy.gs_archive import GroundStateArchive
from salmonpy.dft import Pseudopotential
from salmonpy.fdtd import SalmonExecutable

converter = GammaRestartConverter.from_build('NEW_CONVERTER_DIRECTORY/build.json')
engine = SalmonExecutable('AUDITED_SALMON_EXECUTABLE', converter.target_sha256)
archive = GroundStateArchive.open('existing-gs-array-archive')
# species/local_channelと実ファイルは利用者が指定。SHA256はアーカイブと照合されます。
potential = Pseudopotential.from_file('Si_rps.dat', atomic_number=14, local_channel=2)
with archive.tddft_initial_state(converter=converter, executable=engine,
        pseudopotentials=[potential], max_steps=1000) as initial:
    with initial.new_state(dt=0.02) as state:
        current = state.step([0., 0., 0.001])
```

`write_bundle(archive, NEW_DIRECTORY, executable=..., pseudopotentials=...)` は独立したrt.inp・擬ポテンシャル・計算処理再開レコードを作るだけです。TDDFTプロセスはnew_stateで初めて起動します。変換中だけ実部・虚部の一時ストリームを保持し、完了後は削除します。元のアーカイブは変更せず、上書きしません。

GNUのレコードI/O設定GFORTRAN_CONVERT_UNITも記録し、後の起動・復元時に変更されていれば拒否します。unitごとの変換指定は未対応です。変換器と計算プログラムのビルド記録は利用者が信頼できるものを指定してください。ハッシュは配布者認証ではありません。

合成の直交規格化軌道で、元のSALMON読取処理と変換器の一致、および実計算プログラムのTDDFT初期化・独立状態・復元後電流の一致を検証します。実材料のSCF収束や物理量の精度を認定したものではありません。

擬ポテンシャルのhashはファイルと元素を照合するもので、元GSのHamiltonianの一致まで認定しません。local_channel、擬ポテンシャルとXCの整合、元GSの物理的な正当性は利用者が確認してください。変換中は実部・虚部の一時ストリーム合計16×格子点数×軌道数バイトと、同量の計算処理複素出力が必要です。一時ストリームは完了後に削除し、波動関数のhash照合は1MiBずつ読みます。

## 既存のΓ点GS/1を書き出す

```python
from salmonpy.dft import GroundStateResult
from salmonpy.gs_native import GammaRestartConverter

converter = GammaRestartConverter.from_build('NEW_CONVERTER_DIRECTORY/build.json')
gs = GroundStateResult.open('existing-native-gs-result')
archive = gs.export_archive('new-portable-gs', converter=converter, orbital_block=16)
```

書き出し側も同じGNUコンパイラ・対象計算プログラムに固定した専用プログラムです。古いimport専用変換器では利用できず、新しい変換器を明示的にビルドしてください。以前の変換器・固定バイナリを上書きしません。

対象はSALMONpyのDFTが生成したGS/1、Γ点1個、非スピン分極、直列計算プログラム、実数または複素数の単一wfn.binです。元の全ファイルhash、実行ファイルSHA256、info.binの軌道数・k点数・反復数・ランク数、占有数と電子数、格子・構造・XC・擬ポテンシャルとgs.inpの物理設定を照合します。手編集した入力や別の入力書式は推測せず拒否します。密度混合履歴を使わない初期推定波動関数による再開入力も、既知の形式のみ対応します。

固有値は `native/gs_eigen.data` の単一k点・単一スピン・原子単位の表から取得し、表の占有数とバイナリの占有数も一致させます。ファイル欠損・eV表・非有限値・不一致は拒否し、0や推定値では埋めません。

密度は`rho(r)=sum_b occ[b]*abs(u_b(r))**2`を軌道順に加算して再構成します。元GSの混合密度をそのまま保存した値ではありません。全エネルギーと固有値は元の報告値で、再構成密度から再計算しません。波動関数と占有数の値は保持しますが、密度混合・反復履歴を含む完全なSCFチェックポイントではありません。

`metadata['provenance']`に再構成密度の由来、元GSマニフェスト・計算プログラム・書き出しプログラムのhash、実数/複素数の別、局所チャネル、レコードI/O設定を記録します。この任意フィールドがない既存GS/2も読み込めます。書き出したアーカイブを再入力するときは局所チャネルも自動照合します。hashは元GSの真正性や物理的妥当性を認証するものではありません。

元GSは変更せず、検証した3つの計算処理レコードを私有一時領域へコピーします。書き出しプログラムは1軌道ずつ処理し、Python側は分割実数ストリームと複素数の作業ファイルをメモリマップして、全軌道を一括でヒープ確保しません。メモリマップのページキャッシュはOS管理でRSS上限ではなく、一時ディスク領域も必要です。完了・失敗の両方で一時領域とマップを解放し、失敗時は作成中の出力のみを削除します。

実数・複素数の合成計算処理レコードで、元SALMONの読み書きとの値の一致と、GS/2経由で再入力したTDDFTの電流・保存復元の一致を検証します。合成試験用のGS/1完了メタデータはAPI検査用であり、SCFを収束させた実材料の結果ではありません。

## 複数k点：明示的な直列変換器

`SerialRestartConverter.from_build('NEW_BUILD/build.json')`で、複数k点まで検証した新しいGNU変換器を明示指定します。既存の`GammaRestartConverter`はΓ点限定のままです。古い変換器の動作も変更せず、新しい変換器のビルドは別ディレクトリで行います。同じ`GroundStateResult.export_archive`と`GroundStateArchive.tddft_initial_state`に渡せます。

対象は座標軸に沿った正の対角セル、周期系・非スピン分極・固定イオン・CPU直列1ランクです。書き出しはSALMONpy DFT生成の通常メッシュに限定し、native/gs_k.dataの直交・逆格子係数形式を必須にします。逆格子、k点の順序と座標、等重みを入力kmeshと照合します（座標atol1e-14、重みrtol1e-14）。非直交セルのCartesian形式や、手指定file_kwを使った既存GSの書き出しは拒否します。読み込み側は、明示的な逆格子係数と非負・和1の重みを受け取り、不均一重みも扱います。

計算処理ではk点を外側・軌道を内側に保持し、占有数はFortranの(軌道,k点)順、GS/2は(k,軌道,z,y,x)順です。密度は`sum_k,b weights[k]*occupations[k,b]*abs(u[k,b])**2`をk点・軌道順に加算して再構成します。重み付き電子数と16×格子点数×軌道数×k点数のサイズを検査し、添字・容量計算は64ビットを使います。計算処理の波動関数作業領域は1軌道分、占有数と重みの領域はk点数に応じて増えます。

再入力はkpoints.datをfile_kwで明示し、そのhashもreceiptに記録・照合します。レコードI/O環境の固定、元データの保護、独立ワーカーの条件は従来どおりです。provenanceには任意のkpoint_source='native/gs_k.data'が追加されます。

波動関数と占有数は計算処理のビットを保持します。一方、k点テキストの桁数制限とfile_kwでの重み再規格化があるため、一般のk点・重み・電流のビット一致は保証しません。元GSと再入力の電流比較はrtol1e-12、atol1e-14で判定します。実計算プログラムの合成試験では2点・8点と不均一重み2点がビット一致し、3点・6点の差は約1e-18でした。同一の再入力状態の復元後はビット一致しました。これはSCF収束や実材料の物理精度の認定ではありません。

```python
from salmonpy.gs_native import SerialRestartConverter
converter = SerialRestartConverter.from_build('NEW_MULTIK_CONVERTER/build.json')
archive = gs.export_archive('new-portable-multik-gs', converter=converter)
initial = archive.tddft_initial_state(converter=converter, executable=engine,
                                    pseudopotentials=potentials, max_steps=1000)
```

複数k点でも、`state.save_checkpoint('new-history')`と、新しく作った同一モデル・固定dtの状態での`state.restore_checkpoint('new-history')`を使えます。GS/2アーカイブと対応する変換器・計算プログラム・擬ポテンシャルは別途必要です。チェックポイントは履歴だけを持ち、GS配列を複製しません。合成試験では2・3・6・8点と不均一重み2点について、6ステップ保存後に別Pythonプロセスで復元し、その後20ステップがビット一致しました。異なる入力履歴の2状態を交互に進めた結果、片方を閉じた後の継続も、独立実行とビット一致しました。これは各状態が専用ワーカープロセスを持つ経路の検証であり、計算処理TDDFTの同一プロセス内での複数コンテキスト対応を示すものではありません。

## 保持計算処理 GSからの保存

`NativeGSState.export_archive` で通常k点メッシュ・DFT設定を保持するGS（全ランク呼び出し、rootが保存）を保存できます。
対象範囲・来歴・検証の限界は [NATIVE_GS.md](NATIVE_GS.md) を参照してください。
完全なSCFチェックポイントや計算処理 GSハンドルの復元ではありません。

## 計算処理保持GSへの初期化

`NativeGSSession.load_archive` はrootでGS/2を一度変換し、新しい計算処理所有者へ
SCFなしで読み込みます。対応する変換器・計算プログラム・擬ポテンシャルを明示します。
[NATIVE_GS.md](NATIVE_GS.md) と examples/portable_gs.py を参照してください。
配列の初期化であり、保存したエネルギー・固有値・混合履歴・途中のSCF状態の復元ではありません。

GS/2経由の保持GSは `salmonpy.retained-gs-reexport/1` の来歴で再保存できます。
計算処理の軌道・占有数と再構成密度を、元データから引き継ぐエネルギー・固有値・
物質・k点情報と区別します。今回の収束は不明のままです。保持GSガイドを参照してください。
従来の来歴は引き続き読め、新しい来歴には対応版のリーダーが必要です。

TB-mBJはルートの `xc_parameters: {"cval": null}`（自動）または固定係数を保存します。計算処理変換はこの情報の欠落を拒否します。他の自作XCでも有限JSON値のパラメータ辞書を保持できますが、内蔵計算処理変換は未対応パラメータを黙って捨てません。
