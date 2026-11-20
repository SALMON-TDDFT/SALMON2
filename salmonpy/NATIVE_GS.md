# メモリ上に保持するGS

> 現行ビルド配置（2026-11-19）：通常XCと高速PR版HSEは同じ
> runtime/sources/solverからビルドします。以下のstandard／通常版はAPIの
> 能力プロファイルを指し、別の保守ソースではありません。既定はHSE無効です。

[English](NATIVE_GS.en.md) · Topic: `retained_native_gs`

通常版ランタイムでは、`salmonpy.native_gs` はDFTのGSを同じプロセス内に複数保持します。
各状態は入力・波動関数・密度・占有数・固有値・力・XC・MPI資源を個別に所有します。
異なるセル・格子・k点・電子数・原子配置・XCの計算を作成した後も、先の状態を参照できます。
Pythonの通常の `DFT.run` は従来の別プロセス実行を維持します。

## 使い方

```python
from types import SimpleNamespace
from mpi4py import MPI
from salmonpy.dft import SCFSettings
from salmonpy.native_gs import NativeGSBackend

# dft_a / dft_b は既存の不変な DFT 設定。仕事用フォルダは全ランクから見える場所に作る。
context = SimpleNamespace(solver_comm=MPI.COMM_WORLD, closed=False, _dependents=0)
backend = NativeGSBackend(build_record, work_root,
                          decomposition=(1, MPI.COMM_WORLD.size, 1, 1, 1))
with backend.open(context) as session:
    a = session.run(dft_a, scf=SCFSettings(), material_id="film-a")
    b = session.run(dft_b, scf=SCFSettings(), material_id="film-b")
    energy_a = a.energy_hartree
    density_a = a.array("density")
    eigen_b = b.array("eigenvalues")
    a.close()
    energy_b = b.energy_hartree  # a の終了後も有効。
# 全状態を終了し、native接続だけを解除。MPI_Finalizeはホストの責任。
```

`decomposition` は `(k, orbital, x, y, z)` で、積を群のランク数に合わせます。
既存の分割済みコミュニケータを `context.solver_comm` に渡しても使えます。
同じ群の全ランクが、同じ順序・入力・物質ID・ハンドルで呼び出してください。
メインスレッドで順次実行します。LD/SBEなどの軽い部品は引き続き個別に利用できます。

`session.run` はSCFを実行します。`state.info` はエネルギー、残差、許容値、反復数、
SALMON本体が記録した収束フラグ、原子数、バンド数を返します。SCF未収束は警告のみで、
最後の反復状態を保持します。`material_id` と不変の `calculation` も状態に残ります。

## 配列とMPI

`state.array(quantity)` は `NativeGSArray(values, lower_indices, quantity, unit)` を返します。

| quantity | 内容 | 単位 |
|---|---|---|
| density | 価電子の密度（スピン合計） | bohr^-3 |
| occupations | 占有数 | 無次元 |
| eigenvalues | 一電子固有値 | Hartree |
| forces | 原子の力 | Hartree/Bohr |
| orbitals | 実数または複素軌道を複素数で返す | native軌道振幅の規約 |

`values` は読み取り専用のコピーです。native状態を変更せず、ハンドル終了後も使えます。
軸と添字はSALMONのネイティブ配列順で、`lower_indices` はFortranの下限です。
密度・軌道は各ランクの局所配列で、軌道にはハローを含む場合があります。
軌道の軸は `(x,y,z,spin,orbital,k,image)`。占有数・固有値は `(orbital,k,spin)`、
力は `(component,atom)` で、後者3種類はnativeで複製される配列です。
全ランクの配列を自動連結せず、複製を二重に足さないでください。

保持状態と、読み出したコピーの分だけメモリを使います。大きな軌道は必要時だけ読み出します。
`session.resources()` は状態数、生存中の子通信数、通信管理スロットの高水位を返します。
一方をcloseしても他方は保持され、全状態を閉じれば生存中の子通信数はゼロへ戻ります。
スロットの高水位が残ることは、生存資源のリークとは区別します。

## 対応範囲と導入

Python入口は周期・全k点・非スピン偏極・固定イオンのPZ、Libxc-PZ/PBEと、
機能宣言のある高速PR版ランタイムでのHSE06です。HSEはk点方向のMPI分割を使います。
SOC／+U／DC、GSの完全なSCF保存／再開、
スレッド並行は対応外です。逐次のGS／TDDFT切替は後述の入口を使います。1プロセスでnative電子計算の接続は1生存期間です。
ハンドルの参照・配列取得・closeも群内の集合呼出しです。全状態を閉じてからホストMPIを終了します。

`sp_gs_protocol()==1` を持つ、明示的にビルド・監査した標準ランタイムを指定します。
古いビルド済みランタイムにはこの入口がない場合があります。importでの取得・ビルドは行いません。
`tools/build_embedded_tddft.py` で構築し、準備済みソースと全アプリケーションオブジェクトを監査します。
MacではNumPyのBLASと衝突しない名前解決と、MPI COMMONシンボルの共有を維持してリンクします。
Fortranホストの開発プローブも `mpi_link_command` を使い、この規約に揃えます。
Linuxの既存リンク方式は維持します。今回の実測はApple Silicon／GNU／OpenMPIで、他環境は別検証です。

仕事用フォルダ・入力・擬ポテンシャルが全ランクで同じ内容に見えることを、native呼出し前に検査します。
一般のFortran STOP、メモリ不足、MPI／プロセス障害の復帰は保証しません。
正常計算の更新順を維持し、未収束を理由に停止する制限は追加していません。

同梱例：`mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT`。
これは粗い格子・4反復の所有確認用です。物性の収束確認を済ませた入力ではありません。
GS再開バンドルは自動出力しません。通常のSALMON診断テキストと進捗出力は生成します。
公開native契約は [standard_gs.h](native_adapters/standard_gs.h) を参照してください。

## 既存GSの読み込みとメモリ上のTDDFT接続

`session.load(directory, input_file="rt.inp", material_id="silicon")` は、対応する
入力と既存GSを読み込みます。SCFを回さず、新たなGS再開データも書きません。
全ランクから変更されていない元ファイルを参照できる必要があります。整合性検査は
旧形式バイナリの物理的な由来を証明しないため、入力との対応は利用者が確認します。
読み込み時の `info.converged` は `None` です。residual/tolerance のゼロは未測定を示し、
収束したという意味ではありません。読み込みでは力を計算しないため、forcesの取得は拒否します。

```python
with backend.open(context) as session:
    gs = session.load(existing_bundle, material_id="silicon")
    model = gs.tddft_model(dt=0.003, max_steps=100)
    gs.close()
    first = model.new_state("first")
    second = model.new_state("second")
    j = first.step([0.0001, 0, 0])  # transverse、A/c → 電荷電流j
    first.checkpoint(checkpoint_path)  # ランクごとに異なる保存先を指定
    # 新しい状態で checkpoint(..., restore=True) とすれば時間発展の再計算なしで復元
    first.close()
    second.close()
    model.close()
```

モデルが初期軌道・占有数・入力値のコピーを所有するため、GSを閉じた後も状態を作れます。
複数状態は同じ初期GSから独立に進みます。`session.run` が返すGSも同じ
`tddft_model` を使えます。GS読み込み・計算とRTステップは同一セッション内で順次
切り替えられます。セッションを閉じると、RT状態・モデルを先に閉じて接続を解放します。
現在はGSのMPI群・分割を引き継ぐ標準周期系・全k点・非スピン分極・固定イオンの経路です。
占有RTバンド数を超える軌道並列数は拒否します。自動再分割やスレッド同時実行には未対応です。
`sp_tddft_memory_gs_protocol()==1` と `sp_gs_import` を持つランタイムを明示的に使います。

```sh
mpiexec -n 2 python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_WORK
```

## 保持GSからの明示的なポータブル保存

`session.run` で作った通常k点メッシュの状態は、サブグループ全ランクから同じ呼び出しでGS/2形式へ保存できます。

```python
from salmonpy.gs_archive import GroundStateArchive

archive = gs.export_archive("new-portable-gs", orbital_block=16)
gs.close()
if archive is not None:  # サブグループのrootのみ。他ランクはNone
    archive = GroundStateArchive.open("new-portable-gs")
    u = archive.wavefunctions(0)  # k=0、(band,z,y,x)、読み取り専用
```

新しいnative ABIやコンパイルは不要です。波動関数・占有数・固有値を保持し、
密度の実領域の添字を使って波動関数のハローを除外します。保存する密度は
占有数と波動関数からバンド順に再構成し、元の混合密度は保存しません。
物質ID・GS識別子・ランタイムのハッシュ・記録された収束フラグ・反復数・
擬ポテンシャルのハッシュと局所チャネルを来歴として残します。
擬ポテンシャル本体は同梱しません。SCF未収束は警告のみで最後の状態を保存します。
既存の保存先は上書きせず、書き込み失敗時は新規作成中の不完全な保存先を除去します。

現段階の対象は、DFT設定を保持する周期・非スピン分極・固定イオン・通常k点メッシュ・
MPIサブグループの状態です。Γ点と複数k点に対応します。既存ファイルから読み込んだGSは不変DFT設定を
持たないため未対応です。GS/2経由の読み込み状態は後述の情報を保持して再保存できます。完全なnative SCF状態の復元は残作業です。
これは保存機能の範囲で、通常のGS/TDDFT計算の制限は増やしません。
非立方格子・非対称ハロー・実数／複素数相当の合成native配列で確認しました。
実GSでの書き出しと、その後のTDDFT結果の一致は未測定です。
今回は新たなSCF計算も実GSの収集も行っていません。

これは配列の交換形式で、完全なSCF再開ではありません。混合履歴・力・ソルバー資源・
元の混合密度は保存しません。`GroundStateArchive.open` は配列を読み込みますが、
native GSハンドルを復元しません。nativeへの取り込みは
[GS_ARCHIVE.md](GS_ARCHIVE.md) の明示的な対応コンバーターとランタイム・
擬ポテンシャルの条件に従います。
サンプルの `python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT --export-archives`
は、明示的に2つの粗い2k点GSを計算してrootで保存します（このデモは1または2ランク）。
4反復は所有権の確認用で物理的な収束を保証しません。今回の検証ではこのSCFサンプルは
実行していません。保存時はハローを含む全波動関数のコピー、再構成密度、チャンク用の
メモリが必要です。`orbital_block` は保存単位で、ピークメモリの上限ではありません。

k点はSALMONの `lattice.f90:init_kvector` に従い、0始まりのiについて
`(i+0.5)/n-0.5`、x最速・次にy・zの順、重み1/nkとします。
これは不変の `DFT.input_text` が生成するシフトなし・非Γ中心の全k点メッシュです。
任意のk点表や対称性による縮約点を推測する処理ではありません。
密度はnativeのk・バンド順で `重み[k]*占有数[k,b]*abs(u[k,b])**2` を足し合わせます。
合成の2・3・6・8点で各k・バンド・座標・重み・電子数を確認し、Γ点の全保存ファイルは
b66b97fdとバイト単位で一致しました。任意／縮約k点の保持状態・native GS
ハンドルの復元は残作業です。nativeへ取り込むコンバーターには別のランタイム・セル等の
条件があり、保存できることが全てのnative読込との互換性を保証するわけではありません。

### MPI保存の所有権と失敗時の動作

サブグループ全ランクが同じ引数で `export_archive` を呼びます。
k点・軌道・空間と混合分割は、ランク番号を推測せずnativeの添字と実領域から組み立てます。
複製された占有数・固有値の一致、全k点・バンド・格子の重複なしの被覆を検査します。
rootだけがGS/2を書き込み、`GroundStateArchive` を返し、他ランクは `None` を返します。
保存先の書き込み権限はrootにだけ必要です。他ランクで後から読む場合は、そのファイルが
見える場所で明示的に開いてください。`gs.close()` は引き続き全ランクで呼びます。

rootはランクごとに1タイルずつ受信し、全ランクの送信配列を同時には保持しません。
ただし全体の複素波動関数、ローカル配列のコピー、再構成密度、保存用バッファは必要です。
任意のメモリ上限は設けません。実際の確保失敗・書き込み失敗は全ランクに通知し、
元のGS状態は変更しません。新規の不完全な保存先は除去し、既存出力は保持します。
この保存だけでは、MPI分割を変更したnative状態の復元はできません。
MPI通信自体の故障やプロセス消失からの復帰は対象外です。

実MPI4ランク・合成配列で、k点／軌道／空間／混合分割（不均等バンド分割を含む）、
直列版との全保存ファイル一致、範囲・複製・被覆・メタデータ・設定・確保・書込の失敗、
元配列の保持と再試行を確認しました。SALMONの実GS保存やSCF/RT物理の検証ではありません。
`mpiexec -n4 python tools/test_retained_gs_archive_mpi.py NEW_OUTPUT` で実行できます。
インストール後も `salmonpy/_runtime_source/tools` に同梱し、native GSやランタイムを読みません。

### 既存GS読み込み時の参照入力ファイル

`session.load` は、擬ポテンシャル・原子配置ファイル・`file_kw` のk点表を、
相対サブディレクトリごとコピーして保持します。各MPIランクで内容を検証し、
メモリ上のGS→TDDFTモデルにも引き継ぎます。引用符内のファイル名の `!` は保持し、
実際のコメントだけを除去します。必要な参照ファイルが欠けた場合はnative初期化前に
拒否し、既存GS・TDDFT状態は変更しません。必須ファイルを読めない場合は初期化できず、
数値の未収束を警告だけにする方針とは区別します。

既存GSの共通比較ランナーに `--external-k` を追加しました。

```sh
mpiexec -n 4 python tools/test_memory_gs_mpi.py BUILD EXISTING_GS NEW_WORK --size 2 --split orbital --external-k
```

既存GSと入力の作業コピーを作り、通常メッシュを同じ外部k点表に置き換え、元の
ファイル経路と比較します。SCFは実行せず、参照GSも新規作成しません。元の入力は
変更しません。作業出力から収集するのは診断の要約だけとしてください。
この任意ランナーは明示シフト・Γ中心指定のない通常メッシュを対象としますが、
`session.load` の通常機能に新しい制限を加えるものではありません。
PZの2点、Libxc-PBE／TB-mBJの8点について、内部MPI2、独立状態と兄弟状態、
再計算なしの直接再開、RT進行中の参照ファイル欠落の拒否を検証しました。
インストール後も `salmonpy/_runtime_source/tools` に同梱します。
読み込んだGSの `calculation` は引き続きNoneです。正規化した物質モデルの紐付け・
旧形式読み込みGSの物質情報の完全な紐付け・完全なGS/SCF状態の復元・元GSの物理的な由来は別の残作業です。
ファイルの整合性や時系列の一致は、物理的なSCF収束を証明しません。

## 保持した宣言入力のメタデータ

任意の `metadata` extra（f90nml）を管理ノード等で事前に導入してください。
ネット接続できない計算ノードで暗黙に取得しません。必須依存はNumPyだけのままで、
通常のimport・run・loadではf90nmlを読み込みません。

```python
# 既存セッション内。input_metadata自体はMPI集団操作ではありません。
metadata = gs.input_metadata()
# closeは従来どおり全ランクで呼びます。
gs.close()
record = metadata.as_dict()  # nativeハンドルを持たないJSON互換のコピー
cell = metadata.settings["cell_bohr"]
```

`GSInputMetadata` は不変で、`as_dict` は編集可能な独立コピーを返します。
利用者の物質ID、保持した入力・参照ファイルのSHA256、source_idを残します。
source_idは宣言された入力バイトの識別子で、物理的な同一物質・バイナリの由来・
チェックポイント互換性を証明しません。書式やファイル名でも変化し得ます。
native側のidentityや利用者の物質IDは変更しません。
`gs.input_metadata()` はGSを閉じた後も、保持したテキスト・参照データだけから
呼べます。元ファイルの読み直し、native呼び出し、MPI通信、再計算は行いません。

明示されたau/A_eV_fsのセルをbohrへ、対応するインライン／外部ファイルの
Cartesian／相対原子座標を分率座標へ整理します。原子種の順序、指定された格子数・幅、
通常kメッシュ条件または外部k点表の座標・入力重み、XC条件、種の添字ごとの
擬ポテンシャルのハッシュ・局所チャネルを返します。nativeで規格化された重みや
省略時の既定値は推測しません。単位不明・曖昧・未対応の形状はNoneとnotesで表します。
これは宣言入力の記述で、実効nativeモデルの問い合わせではありません。

読み込みGSの収束・物理的由来は不明のままです。runで作ったGSでは、既に取得した
収束フラグを記録します。説明の読み出しでSCFを回したり、GSを再認証したり、
未計算の力を公開したり、読み込みGSのcalculationを設定したりはしません。
native復元や読み込みGSのポータブル保存を許可するものでもありません。
通常の計算・読み込みに新たな制限は加えません。任意の説明処理で依存不足・
入力解釈エラーが起きてもGS状態は変えず、呼び出し側で計算とは別に扱えます。

自作ホストでも `salmonpy.gs_input.describe_gs_input(text,
assets=相対名からbytesへの対応, material_id="film")` を使えます。
ランタイムやファイルアクセスは不要です。
`python -m salmonpy.examples.memory_gs BUILD EXISTING_GS NEW_WORK --describe-input`
は、GSを閉じた後にrootの記述を表示します。rootの説明処理が失敗した場合は、次の
集団操作へ進む前に全ランクで合意します。共通比較ランナーも `--describe-input` に
対応します。実際に説明を呼ぶランクにはf90nmlが必要で、importやヘルプではMPIを
初期化しません。

## ポータブル配列から保持GSを初期化する

`session.load_archive` でGS/2を既存nativeローダーへ接続できます。
事前に検証済みのGammaRestartConverterまたはSerialRestartConverterと、その対応する
直列SalmonExecutableを用意します。暗黙のコンパイル・取得・SCFは行いません。
全MPIランクが同じ条件で呼び、rootだけがアーカイブを読み変換します。その後の新規
変換先ディレクトリは、サブグループ全ランクから参照できる場所に置きます。

```python
with backend.open(context) as session:
    gs = session.load_archive(
        "portable-gs", "new-converted-bundle",
        converter=converter, executable=paired_serial_engine,
        pseudopotentials=potentials, material_id="silicon")
    origin = gs.archive_origin  # 不変のハッシュ・由来情報。close後も保持
    model = gs.tddft_model(dt=0.003, max_steps=100)
    gs.close()
    state = model.new_state("response")
    j = state.step([0.0001, 0, 0])  # transverse A/c → 電荷電流
    state.close()
    model.close()
```

material_idを省略すると、保持GSからの保存来歴に物質IDがあれば継承し、なければ
"material"とします。明示したIDを優先し、元のIDや過去の宣言収束状態はarchive_originで
今回の読み込み収束状態(None)と区別して保持します。アーカイブのmanifest、変換器、
対応直列実行ファイル、共有ランタイムのハッシュを残します。input_metadataも利用でき、
calculationはNone、力は未計算のままです。元の保存データ・以前の状態は上書きしません。

新しいnative所有者へ軌道・占有数を初期化する機能です。保存済みの密度・固有値・
エネルギーをnative再開レコードとして復元するものではなく、密度・ポテンシャルを
初期化時に再構成します。混合履歴・ソルバーの内部資源・途中のGS反復は復元しません。
変換エラーはnative呼び出し前に全ランクで合意し、新規の不完全な変換先を除去します。
変換済みの完全なディレクトリは、その後native受付が失敗した場合も利用者用に保持します。
既存保存先は上書きしません。native STOP・MPI通信故障・プロセス消失の扱いは従来どおりです。

変換器の対応条件は引き継ぎます。周期・非スピン分極・固定イオン・偶数電子の閉殻、
各方向8格子以上、対応する共有ランタイムでのPZ／Libxc-PZ／PBEが対象です。
Gamma版は厳密Γ点、SerialRestartConverterの複数k点は正の軸に沿った直交セルが対象です。
変換はrootで直列に行い、その後のnative初期化には対応するk点／軌道／空間MPI分割を使えます。
任意のRT/SCFチェックポイントの再分割ではありません。対応する直列実行ファイルは、
ここでは形式を結び付ける識別子で、実際のソルバーとして起動しません。

記録されたFortranレコードマーカー・整数幅の不一致は変換前に拒否します。
コンパイラー・ソースの照合情報が不明、またはその他のフラグが異なる場合は警告だけで、
明示的に選んだnativeローダーへ進みます。古い記録にないコンパイラーのハッシュを
推測で補いません。未検証のコンパイラー・環境での一般的なnative I/O安全性は保証しません。
レコード変換の環境条件はランク間で一致させ、変換済みの記録・内容を読み込み前に検査します。

`salmonpy.examples.portable_gs.current_from_archive(session, ...,
A_steps=..., dt=...)` はdt,2dt,...でのA/cからjを返し、この計算の所有者を解放します。
繰り返す場合は同じセッションを保持し、毎回新しい変換先を指定します。変換済みデータは
利用者が管理します。ネット接続できない計算ノードでは、管理ノードで変換器を準備し、
検証済みの配置記録と実行ファイルを事前に移してください。

Γ点の配列読み込みも等価な明示Γ点表を使い、native DFTで複素軌道を確保します。
既定の実数軌道分岐へ入らず、ゼロでない虚部を保持します。生成k点表のハッシュを
変換記録へ含め、この経路をarchive_originにも記録します。複素Γ点の合成データで、
軌道・空間MPI分割とも軌道の完全一致、独立／ファイル経路の電流、直接再開が通りました。

## GS/2から読み込んだGSを再保存する

`load_archive` が返す生存中のGS所有者は、同じ集団API `export_archive` で再保存できます。
旧形式の `load` だけで読み込んだ状態は完全な物質モデルが不明なため、部分的な入力説明を
使って推測で保存可能にすることはありません。

```python
with backend.open(context) as session:
    gs = session.load_archive("source-gs", "new-native-bundle",
        converter=converter, executable=paired_serial_engine,
        pseudopotentials=potentials, material_id="silicon")
    saved = gs.export_archive("new-gs-array-archive")  # rootのみarchive、他はNone
    # GSと既存のRTモデルは変化しません。
    gs.close()
```

読み込み時に保持するのは小さなメタデータ・固有値・k点・構造の配列と不変の擬ポテンシャルです。
元の波動関数全体や生きたアーカイブオブジェクトはキャッシュしません。GS所有者が生存していれば、
元アーカイブを削除した後も再保存できます。軌道・占有数はnativeから読み、ハローを除外します。
密度は元のk点重みを使って再構成します。構造・格子・XC・擬ポテンシャル参照、k点座標・入力重み、
参照固有値・エネルギーは元アーカイブから引き継ぎ、新たに計算した値とは扱いません。

来歴の新しい種類 `salmonpy.retained-gs-reexport/1` に、energy_source／eigenvalue_source／
sampling_sourceがsource_archiveであること、元manifestのSHA、引き継いだ収束の宣言を記録します。
今回のconvergedはNoneのままです。参照値からSCF収束やnative内部状態の復元を主張しません。
この来歴形式を読むには対応版が必要で、既存GS/2やDFT計算由来の来歴形式は変更しません。
calculationはNone、力は未計算のままです。複素配列は実部・虚部を演算せずコピーし、ゼロの符号も
保持します。NaN／Infinityは全ランクで拒否します。

非一様な重み・不規則k点も、通常メッシュと推測せずそのまま保存します。既存DFT由来の保存の
演算順は維持します。サンプルの `current_from_archive(..., reexport_to="new-array-archive")` で、
RT状態を作る前の保存も選択できます。再読込は新しいモデルの初期化で、既存RTチェックポイントを
別モデル・別MPI分割へ移す機能ではありません。完全なSCF混合履歴・内部資源・反復途中の再開は別の残作業です。

### MPI全ランクでの配列取得

`gs.array(quantity)` はGSのMPI群の全ランクから、同じ対応済みの量を表す
文字列で呼びます。不正な型・未対応の量・インポート時に計算していない力・
ランク間で異なる量の指定は、native配列取得前に全ランクで `InputError` に
なります。指定を直して再試行でき、保持中のGSと動作中のTDDFT状態は変わりません。
ランクの呼び忘れやAPI呼出順の不一致からの復旧を保証するものではありません。
正常入力の数値と取得可能な量は変更していません。

### SCF続行のための状態分離

新しくビルドする標準ランタイムでは、収束判定用の直前密度・局所ポテンシャル、
電子数規格化、混合なし反復のカウンター、終了時の反復数・残差・収束フラグ・
イオン更新フラグを、各GSが所有する混合状態に保持します。既存の入口は
新しいSCF呼出しごとに従来どおり初期化し、演算順と収束判定方針を維持します。
GSファイル・GS/2読込み・RTチェックポイントの形式は変更していません。
これは内部状態分離の段階で、SCF続行の公開APIと完全保存はまだありません。
直前場の配列2本は、GSをリセット／解放するまで保持されます。今回の変更は
標準版のみで、HSEのソースプロファイルは変更していません。

## 高速版で既存GSを読み込む

上のSCF例は通常版ランタイムを対象とします。gs_abiプロトコル1を明示した高速版
ビルドでも、SCFなしの既存GS読み込みを使えます。

```python
context = SimpleNamespace(solver_comm=MPI.COMM_WORLD, closed=False, _dependents=0)
backend = NativeGSBackend(fast_build_record, shared_work_root,
                         decomposition=(MPI.COMM_WORLD.size, 1, 1, 1, 1))
with backend.open(context) as session:
    a = session.load(bundle_a, material_id="film-a")
    b = session.load(bundle_b, material_id="film-b")
    density = a.array("density")  # 元の状態から独立した読み取り専用コピー
    occupations = b.array("occupations")
    a.close()                     # bは引き続き参照できる
```

機能宣言のない古い高速ライブラリは自動で採用・再構築しません。入力は既存の
波動関数と互換である必要があります。全k点・縮約k点の初期化でモデル別に
対称性作業領域を切り替え、必要なsym.datも引き継ぎます。読み込みでは不要な
RT因数学習と予備SCFだけを無効にし、XC・HSE遮蔽係数は保持します。
呼び出しは主スレッドでMPI群内の全ランクが同じ順に行います。

高速版では密度・占有数・軌道を取得できます。infoのenergy_hartree・residual・
tolerance・convergedは未計算を示すNoneで、物理的なゼロではありません。
固有値・力は計算していません。import専用記録のsession.runと、読み込み状態の計算済みアーカイブ出力は利用できません。tddft_modelは
下記の追加機能を持つランタイムで利用できます。通常版の既存機能は維持します。読み込みだけで収束を保証しません。

既存サンプルに読み込みモードを追加しました。
`mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD NEW_OUTPUT --import-bundle GS_A --second-bundle GS_B`
各状態のmaterial_idを保持し、出力JSONの未計算値はnullとします。合成入力での
独立性検証であり、実物質の物理精度の検証ではありません。


### 保持GSから高速PR版HSEの応答へ

新しくビルドした高速版で `gs_abi.memory_rt == 1` と
`sp_fast_hse_memory_gs_protocol() == 1` を確認できれば、
`gs.tddft_model(dt=..., max_steps=...)` を使えます。HSE06・原子単位・
周期スカラー系の互換GS、k点のみの内部MPI 1／2ランクが対象です。
通常版の既存GS→TDDFTは維持し、PZ／TB-mBJ／Libxcの通常RTは、下記のHSE無効プロファイルで利用できます。

```python
with backend.open(context) as session:
    gs = session.load(bundle_a, material_id="film-a")
    model = gs.tddft_model(dt=0.001, max_steps=32)
    gs.close()  # RTモデルが初期軌道・占有数を独立に所有する
    a = model.new_state("point-a")
    b = model.new_state("point-b")
    try:
        j_a = a.step([0.0, 0.0, 1.e-5])
        j_b = b.step([0.0, 0.0, -1.e-5])
    finally:
        a.close()
        b.close()
        model.close()
```

接続時にSCF、GS再開ファイルの複製・出力、時間発展の再実行はしません。
静的ポテンシャル・対称性・RT資源は従来どおり初期化します。初期軌道と占有数は
モデルが値として所有するため、その分のメモリは必要です。ゼロコピーではありません。
元のGSを閉じて元データを削除した後も新しい状態を作れます。別GSの逐次読み込みでも
既存RTの設定を維持します。全RT状態を閉じるまでsessionとホストMPI群を保持します。
モデルを閉じると新規作成は停止し、既存の子状態は利用できます。主スレッド・全ランク
一致・逐次呼び出しの契約は変わりません。

AはCartesian A/c、戻り値は電荷電流密度jです。全k点ではx／y／z・回転する外場を
使えます。縮約時は実際の外場が保持対称性に適合する必要があり、例のz方向は
検証したC2z入力に合わせたものです。直接保存／復元は既存の同一ランタイム・MPI分割の
契約に従います。任意のRT／SCF途中状態を初期化する入口ではありません。

既存サンプルの `--import-bundle` に `--current-steps 16` を追加できます。
`tools/test_memory_gs_mpi.py BUILD FIXTURE NEW_OUTPUT --size 2 --fast-hse` は
4ランクで従来ファイル経路と比較します。合成入力のTaylor／学習付きPT-CN・全k／縮約kで
電流・エネルギー、状態独立性、GS終了・元データ削除後の新規作成、GS逐次読み込み、
再実行なしの直接再開を検査します。実装の検証であり、物理的収束の保証ではありません。


### 高速版ソースの通常XC接続（HSE無効）

同じ高速PR版の準備ソースから、通常XC用のCPUホストMPIライブラリを構築できます。

```sh
python tools/prepare_hse_rt_state.py NEW_SOURCE --explicit-routing
python tools/build_embedded_tddft.py NEW_SOURCE NEW_BUILD --jobs 4 \
  --cmake-arg=-DUSE_HSE=OFF --cmake-arg=-DUSE_LIBXC=ON
```

MPI・コンパイラ・BLAS・Libxcの場所は環境に合わせて明示します。HSE用のFFTW依存は不要です。
ビルド記録はsalmonpy.embedded-tddft/1を維持し、source_profile="fast-pr"、
ordinary_rt_abi=1、gs_abi/1を宣言します。新しいHSE無効ビルドではmodesにrunも含みます。既存のNativeGSBackendの
session.load→gs.tddft_model、またはNativeTDDFTBackendのファイル／レジストリ経路を使います。
PZ・TB-mBJ・Libxc-PZ／PBE、周期全k点・スカラー・固定イオン・middlepoint・A/c駆動が対象です。
物質IDを継承し、電荷電流密度jを返します。

HSE無効構成ではsession.runでSCFを実行し、エネルギー・固有値・力を取得できます。
SCF未収束は警告とし、最終反復の状態を保持します。session.loadしたGSの未計算メタデータはNoneです。
計算済みGSは明示的なexport_archiveに接続します。今回の検証は出力配列の組立てまでを
メモリ上で確認し、GSアーカイブ自体は保存していません。GS終了後の独立状態作成と
直接RT再開にGSファイルの複製・時間発展の再実行は不要です。初期配列のコピー分は
メモリを使います。全RT状態を閉じるまでsessionとMPI群を保持します。別GSの逐次読込みでは
RTの入力・タイマーを退避／再開します。主スレッド・全ランク一致・一接続寿命の条件は維持します。

通常TDDFTの接続コードは標準版と共有し、コンパイル時に高速版GSの値所有シードへ接続します。
高速版では差分格子とXCの作業領域が状態所有なので、旧来の共有領域切替は行いません。
ファイル入力ではAPIが指定したdt・最大ステップ・MPI分割を反映してから入力スナップショットを
確定し、標準版と同じ順序にします。高速HSEと標準版の既存経路は維持します。
通常XCとHSEを同じGS session内のレジストリに混在させる入口と、旧ソースツリーの削除は
未実装です。HSE有効構成の新しいGS計算入口は後述します。最終的なソース統合は残っています。

retained_gsサンプルの--import-bundle／--current-steps 16も利用できます。
合成入力でk点・軌道・空間の内部MPI2、直接再開、異なるXCのGS逐次読込みを検証します。
従来バイナリとの比較ではPZ・TB-mBJの完全一致と、Libxcの微小差（今回最大約5.5e-20au）を
区別します。再開・独立状態の比較は完全一致です。実物質の収束を保証する検査ではありません。

通常XCのGS接続検証では、PZ／Libxc-PZ／PBE、軌道・空間MPI2で独立実行と配列が完全一致しました。
変更前の高速版とのPZ／PBEエネルギー・テキスト出力も一致しました。4反復の未収束PBEで
旧標準版との差は約3.4e-7 Hartree残り、標準版との同一性や物理的収束を保証しません。

計算済みGS→TDDFTも、同じtddft_model入口を使います。例えばHSE無効・Libxc有効の
ビルド記録で、次のサンプルは4反復のGSを作り、最初のGSを閉じて16ステップ進めます。

```sh
mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD_JSON NEW_OUTPUT --current-steps 16
```

この小さい例は物理的収束の例ではありません。SCF未収束は警告だけです。波動関数や
GS再開ファイルは保存せず、結果JSONにGSの数値要約と電荷電流を記録します。
物質IDを引き継ぎ、x/y/zを含むA/cで独立した2状態を駆動します。
既存GSを読む場合の--import-bundleも引き続き使用できます。

統合検証では、PZ／Libxc-PZとLibxc-PBE、軌道・空間MPI2で、GSから作ったモデルを
保持したまま別のGSを計算します。GS終了後の電流は独立実行と完全一致し、
16ステップの途中で直接保存・復元した状態も、時間発展の再実行なしで一致しました。

### 高速PR版HSEの計算済みGS

新しい高速HSEビルドはgs_abi.modesにrunを宣言します。以前のimport専用ビルドに
計算機能を推定して追加しません。DFT(..., xc="hse06")を作り、
NativeGSBackendのsession.runに渡します。MPI分割は(size,1,1,1,1)、
周期・全k点・非スピン偏極・固定イオンです。既存のGS→TDDFTと同じく、
元のGSをcloseした後もモデルは独立した初期配列を保持します。

```sh
mpiexec -n 2 python -m salmonpy.examples.retained_gs BUILD_JSON NEW_OUTPUT --xc hse06 --current-steps 16
```

新しいHSE入力はexx_pre_scf_threshold=0を明示し、指定したSCF反復を最初からHSEで行います。
PBE予備SCFの進捗をHSE収束と誤認しません。未収束は警告とし最終反復の状態を返します。
HSEのバッチDFT.runは出力レポートが未接続のため、この入口へ案内します。
通常HSEへのフォールバックはありません。力の物理的収束・精度の検証は別途必要です。

macOSのHSEライブラリは通常TDDFTと同じリンク処理を使い、NumPyのBLASと
指定したBLASの複素数ABIを混同しません。機能追加には新しい明示ビルドを使い、
既存の固定ランタイムは書き換えません。

## 再開用配列を保存せず、エネルギー・力を得る

バッチDFT.runは既定で従来どおり再開データを保存します。エネルギーと力だけ必要なら
実行時にsave_restart=Falseを指定します（モデル設定は変わりません）。

```python
gs = calculation.run(engine, "new-output", scf=SCFSettings(max_iterations=4),
                     save_restart=False)
print(gs.energy_hartree, gs.forces_hartree_per_bohr, gs.restart_saved)
```

対象はバッチ入口のPZ／Libxc-PZ／PBEです。4反復は接続確認用で、物理的な収束を保証しません。
未収束でも警告とともに結果を返します。restart_savedは実際に再開用3ファイルが揃ったかを
示し、省略した場合はFalseです。エネルギー・力の読み込みは可能ですが、バンド初期化、
TDDFT初期化、warm start、GSアーカイブ出力には再開用配列を保存した結果が必要です。
save_restartを記録していない従来の結果は保存ありとして扱います。
ランタイムが保存省略を無視しても計算失敗にはせず、警告して実際の保存状態を記録します。

高速版ソースもgs_status.dat（SALMONPY_GS_V1、Hartree）と
gs_forces.dat（SALMONPY_FORCES_V1、Hartree/Bohr）を出します。SCFや力の計算順は変えません。
HSEのnative入口とCLIもこのレポートを出しますが、DFT.runのHSE選択は
高速PR版の来歴確認が未接続のため、引き続きNativeGSBackendを使います。

Siサンプルは--no-restartと--max-iterationsに対応します。
`python -m salmonpy.examples.dft_si EXECUTABLE SHA256 PSEUDO NEW_OUTPUT --no-restart --max-iterations 4`

### 一本化に向けた高速版NLCCの移植

高速ソースに、通常版で検証していたNLCC付きLibxcのエネルギー重み付けと内殻の力を移しました。
Libxcが使う価電子＋内殻密度と、エネルギー密度の重みを一致させます。力は各物質が持つ
原子種・境界条件・XCポテンシャルと、実際にサンプルした内殻表の微分から計算します。
HSE用のプロジェクター微分、状態所有、既存のスピン処理は保持します。

V1の打切りとV2の全サポートを分け、V2は斜交セル・セルをまたぐ長い内殻・周期像を
同じ補間規約で扱います。広いV2内殻の力積分には補償和を使います。IEEE検査はスカラーで、
使用する表の範囲を検査し、V2の未使用パディングは参照しません。

検証は与えた密度のエネルギー微分と力、OMPおよびMPIの分割比較です。
実材料の収束を保証せず、ASE CalculatorのNLCC力に必要なランタイム検証記録を
自動的に付与するものでもありません。従来ランタイムの保護は維持します。
NLCCなしの既存計算は数値回帰で確認します。旧ソースの削除・既定ビルドの切替は
残る経路の確認後に行い、この移植だけで一本化完了とはしません。

## SCF再開に向けた進捗の保持（2026-11-19）

共通の高速版ソースでも、前回密度・前回局所ポテンシャル・電子数による規格化・
収束フラグ・反復数・残差・混合前の収束回数・原子更新フラグを、混合状態が所有します。
状態をコピーすると作業配列も独立します。計算式と更新順は変更しません。
現在のSCF入口は新しい呼び出しごとにこの進捗を初期化します。これは保存に必要な
所有権の整備であり、途中のSCFを再初期化なしに続ける公開APIや、完全な保存形式の
完成ではありません。混合履歴・CG・各ポテンシャル・SCFキャッシュを含む継続は次の作業です。
開発順はSCF保存・再開→共有状態・異常復帰・native ABI→NLCC／単位整合性です。
高速HSEのノルム診断と物理的な格子収束キャンペーンは、利用者の指示で保留しています。

## SCFの途中停止・メモリ内保存／再開

新しく明示ビルドしたCPUランタイムでは、周期・非スピン分極・
固定イオンのPZ／Libxc-PZ/PBE（明示的な対称性あり／なし）と高速PR版HSE06（全k点、または後述の縮約k点）でSCF protocol1を利用できます。従来のrun()は
変更せず、stepsを指定した場合だけ反復境界で止めます。

```python
from salmonpy.dft import SCFSettings
with backend.open(context) as session:
    with session.run(calculation, scf=SCFSettings(max_iterations=100),
                     material_id="film-A", steps=5) as state:
        with state.snapshot() as saved:
            state.continue_scf(steps=10)
            state.restore(saved)  # 初期化・過去反復の再実行なし
            state.continue_scf(steps=20)
        print(state.scf_status)
```

総反復上限・許容値・物質入力は固定し、stepsは呼ぶたびに指定する反復予算です。
途中で収束すれば終了します。上限で未収束なら警告して結果を保持します。
意図的な一時停止はSCF失敗ではありません。途中停止中の力は未評価であり、
SCF区間の完了後に取得できます。

NativeSCFSnapshotはFortranメモリに独立した数値コピーを保持します。対象は
軌道・CG、占有、混合履歴と率・収束進捗、密度・スピン・各ポテンシャル、
Poisson／逆格子／XC作業領域、初期化後の実効入力です。GS→TDDFT用の元の
宣言入力とは分けて保持します。継続時には初期化・初期混合・収束作業領域の
リセットを飛ばします。復元にHamiltonian計算やSCFの再実行はありません。
scf_statusのinitializationsとexecuted_iterationsで回数を確認できます。
復元後に意図的に分岐を再計算した分は、実行回数に加算されます。

全操作は内部MPI群の全ランクが同じ順序でメインスレッドから呼びます。
スナップショットは同じ生存GS状態・同じビルド・MPI分割に属します。
使用後は明示的にcloseしてください。GS状態のcloseでも全コピーを解放します。
コピー分のメモリを使います。今回はディスク保存・pickle・プロセス終了後の再開・
MPI分割変更には対応していません。テキスト出力位置は巻き戻しません。
不正／ランク不一致の選択は数値状態の変更前に拒否しますが、一般のOOM・
初期化・プロセス故障まで原子的に回復する保証はありません。
SOC・+U・DC・イオン運動のSCF継続はまだ対象外です。高速PR版HSE06には後述のk点MPIを使います。
旧ランタイムで従来のrun/loadは使えます。新APIは対応する明示ビルドが必要で、
不足時に自動コンパイルはしません。

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT
```

この小規模例はGS再開ファイルの出力を無効にし、数値配列をRAM内で比較して、
スカラーの検証記録だけ保存します。GS波動関数・密度データセットは収集しません。

### 高速PR版HSE06のSCF途中保存

通常XCと同じ `run(steps=...)`、`snapshot()`／`restore()`、
`continue_scf(steps=...)` を、CPU・全k点の高速PR版HSE06でも使えます。
Gamma点と内部k点MPIに対応し、`decomposition=(comm.size, 1, 1, 1, 1)` と
HSE計算対応の明示ビルドを選びます。HSE無効の通常ランタイムでは使えず、
低速HSEへのフォールバックはありません。

混合履歴に加えて、状態所有のACE／減衰演算子、source・action・占有数の
キャッシュ、factor／PT-CN履歴の配列とフラグ、交換エネルギー、支持領域の
判定、実際の混合反復齢、PBE事前SCFの段階・連続判定回数・開始反復を保持します。
FFT／backendカーネルは生存資源として残し、数値部分のコピー中だけ一時的に
移動します。復元はカーネルの複製・再初期化・交換演算を行いません。
物質状態のポインタは保存データではなく元の生存コンテキストへ戻します。
復元後にPBE→HSEの本来の切替点へ達した際の演算子構築は、連続実行でも
必要な将来の計算であり、復元時の再計算ではありません。

HSEへ移行済みなら、反復上限で未収束でも警告して最終反復を保持します。
未収束を収束済みとは表示しません。PBE事前段階が終わらない結果はHSE結果とは
扱えません。小さな試験の緩い閾値は段階切替の検査用で、物理精度の保証ではありません。
保存は同じ生存状態・ビルド・MPI分割でのメモリ内保存です。縮約k点は後述し、SOC/+U/DC、
イオン運動、ディスク保存／プロセス再開は今回の拡張に含みません。

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume HSE_BUILD_JSON NEW_OUTPUT --xc hse06
```

密度・占有数・軌道・固有値・力・メタデータが、新旧ランタイムの連続実行と
一致することを検証しています。PBE→HSE切替の前後での復元、異なる物質を
挟んだ独立性、未収束時の警告と保持、GSから続けた16ステップのTDDFT電流も
確認しました。実装の回帰検証であり物理的収束の検証ではありません。
GSの波動関数・密度データは収集していません。

### 共通の保存／復元APIへの接続

`state.checkpoint_capability` は生存状態用のnative-scf-memory/1を宣言します。
`state.checkpoint_target()` は保存、`state.checkpoint_target(snapshot)` は直接復元も
接続します。CHECKPOINTING.mdの `execute_checkpoint` を使えます。
scf_resumeサンプルもこの共通入口を使用します。保存データは各ランクのFortranメモリに
残り、接続の作成だけでは保存やSCF反復は走りません。新しい対応ランタイムは独立した
数値候補の準備にも対応します。restore_checkpointで全体の準備成功を確認した後、
各候補をpublishしてから継続します（CHECKPOINTING.md参照）。旧ランタイムは
保護付き準備を拒否し、暗黙の直接復元へのフォールバックはありません。
複数部品の直接復元は引き続き原子的ではありません。

### 明示的な対称操作と縮約k点HSEのSCF

`DFT(..., symmetry=SymmetryOperations(operations, directions="yyn"))` で
分率座標の回転と並進を渡せます。operationsの形は(n,3,4)で、最初の3列が回転、
4列目が分率並進です。入力はコピーして不変に保持します。対称性の自動検出、セルの
標準化、ネットワーク参照は行いません。高速PR版は既存の分率符号付き置換行列の範囲で、
セル・格子・原子・k点の適合性も検査します。利用者が結晶と、後のRTでは実際の外場を
保存する部分群を選びます。全体の偏光をz方向に限定する機能ではありません。

```python
from dataclasses import replace
import numpy as np
from salmonpy.dft import SymmetryOperations
ops = np.zeros((2, 3, 4))
ops[:, :, :3] = np.eye(3)
ops[1, :2, :2] = -np.eye(2)  # z軸まわりのC2。恒等操作も含む。
calculation = replace(calculation, symmetry=SymmetryOperations(ops, directions="yyn"))
```

`input_assets` は生成したsym.datと擬ポテンシャルを含みます。native／バッチとも同じ
保持済みのバイト列を書き、MPIでは実行前に資産の内容を照合します。対称性のハッシュを
物質・warm startの識別に含め、GS→HSE TDDFTにも資産を引き継ぎます。
既定のsymmetry=Noneでは従来の入力本文と計算経路を維持します。
縮約k点のSCF継続はCPU高速PR版HSE・全q点・内部k点MPI1／2に対応し、全k点と
同じ停止／継続／snapshot／直接復元／保護付き数値候補のAPIを使います。
PZ／Libxc-PZ／Libxc-PBEの明示的な対称性付きSCFにも同じAPIを使えます。
内部k点MPI2と空間MPI2で、物質を交互に進めた場合の独立性・保護付き復元・連続実行との一致を確認しました。
対称性・密度対称化のテーブルは各GS
コンテキストが所有し、継続時に切り替えます。固定の初期化テーブルは元のモデルに残り、
別のソルバーを初期化する方式ではありません。

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume HSE_BUILD_JSON NEW_OUTPUT --xc hse06 --symmetry c2z --protected
```

C2で縮約したSCFの復元、異なるセル・格子の物質を挟む実行で、軌道・占有数・エネルギー・
密度と続く16ステップのz駆動RTが旧実装の連続実行と完全一致しています。
対称力の意図的な単位修正だけは、後述の旧基底を復元して比較します。これは実装の検証で、物理精度の認定ではありません。
別の固定格子Si試験では全8点／縮約6点をSCF収束させ、エネルギーは1e-6 Hartree以内、
密度は相対1e-4・絶対1e-7 bohr^-3以内でした。
比較から既存の単位不整合も見つかりました。力をCartesian座標から直接格子基底へ変換して
対称化した後、平均をCartesianへ戻していませんでした。恒等操作・C2・一般セルの
実コード試験で旧処理の失敗を確認し、Fortranで `Amat * average(Ainv * F)` を返すよう
修正しました。**対称化した力は、正しいCartesianのHartree/bohrへ意図的に変わります。**
対称性なしの力、エネルギー、密度、固定イオンでの電流は変更しません。
旧対称力との回帰比較だけは旧基底からの復元を行いますが、アプリ側の自動変換ではありません。

修正後の同じ固定格子・SCF収束済みSiでは、全8点／縮約6点の力も1e-6 Hartree/bohr以内で
一致し、最大差は約1.15e-7でした。これは対称性等価性と単位の検証で、格子・k点収束、
エネルギー微分、構造最適化までの物理精度を認定するものではありません。
NLCC／ASEの検証ゲートは解除せず、過去の出力も上書きしません。

新しい `sp_gs_cartesian_force_protocol()==1` が座標基底の修正を宣言します。
宣言のない旧ランタイムの対称力は、警告して値を保持し、Hartree/bohrと誤表示せず
基底未確認と表示します。Cartesianの対称力には対応する修正版を明示的に選んでください。
暗黙のビルドやPythonでの数値変換は行いません。縮約GSアーカイブ、ディスク再開、
MPI分割変更は引き続き別の残作業です。

開発用の小さな検査は `python tools/test_symmetry_force_units.py NEW_OUTPUT` です。
Fortranコンパイラを明示的に使い、実際の対称化モジュールを監査・ビルドして調べます。
新しいGS計算や、ビルド済みランタイムを通常利用するための前提ではありません。

### SCF途中状態のディスク保存・別プロセスでの再開

対応する新規ランタイムでは `state.save_scf(new_directory)` で反復境界を保存し、
別プロセスの `session.run(calculation, scf=scf, material_id=..., restart=directory)`
で読み込めます。これは収束済みGSの `export_archive()` による再利用とは別の機能です。
同じランタイムのSHA256、物質ID、入力・擬ポテンシャル・対称性資産、MPIランク数と
分割が必要です。入力条件を変えて続けるwarm startや、MPI再分割の形式ではありません。

```python
# 保存するプロセス
with session.run(calculation, scf=scf, material_id="Si", steps=3) as state:
    state.save_scf(checkpoint_directory)  # 新しい保存先。既存の保存は上書きしない。
# 別プロセス。同じcalculation/scfとMPI分割を指定する。
with session.run(calculation, scf=scf, material_id="Si", restart=checkpoint_directory) as state:
    assert state.scf_status["executed_iterations"] == 0
    state.continue_scf(steps=5)
```

Fortranが各rankの数値配列、密度混合履歴・混合率・反復位置、初期無混合の進捗、
XC段階、高速HSEのACE・履歴・占有・キャッシュをstream形式で保存します。
Pythonは版付きmanifestとファイルのSHA256を扱い、全rankの保存成功後にmanifestを
書きます。不完全な保存は再開用として受け付けません。保存先は計算結果領域を指定します。

読み込みはまず全rankの入力・ファイルを検査します。通信、格子、擬ポテンシャル、
Ewald、FFT／Libxcなど実行資源を準備し、保存配列を独立候補へ読み、全rankの成功を
確認して公開します。初期波動関数・密度の生成、混合の初期化、初回エネルギー評価、
SCF反復は通りません。`executed_iterations` は今回のプロセスで実際に実行した回数です。
保存前の実行回数はmanifestに残します。FFTハンドル、ポインタ、MPIコミュニケータ番号は
保存せず、実行側の資源を使用します。HSEのFFTカーネルは必要な次反復で準備されます。

検証ではPZ／Libxc-PZ/PBEと高速PR版HSE、k点MPI2／通常XCの空間MPI2、
C2縮約k点、PBE事前段階からHSEへの切替で、3反復保存＋別プロセス5反復が連続8反復と
完全一致しました。rankファイル破損は状態作成前に拒否し、native読込失敗でも既存状態を
保持します。これは実装検証で、物理収束の認定ではありません。未収束は警告して保持します。

既存のサンプルで `--save-scf NEW_CHECKPOINT` を指定すると停止位置を保存できます。
別の実行で `--restart CHECKPOINT` を指定してください。物質設定、`--xc`、対称性、
`--iterations`、MPI分割は揃えます。保存済みのSCF位置はファイルから復元します。
形式は `salmonpy.native-scf-disk/1`、ABIは `sp_gs_scf_disk_protocol()==1` です。
別ランタイム・MPI分割変更・可搬GSアーカイブへの変換は、この形式の対象外です。

### 残差の反復履歴と未収束時の図

`save_scf()` は数値チェックポイントに加え、`diagnostics/residuals.csv`、
`diagnostics/residuals.png`、`diagnostics/report.json` を出します。PNGは残差の変化を拡大した図と、残差・
収束閾値を対数軸で重ねた図を並べ、未収束ならその旨を表示します。CSVの列は反復番号、残差、
閾値、PBE事前段階フラグです。初期の無混合反復は残差を評価しないため、この表には含めません。
残差の定義は指定したSCF収束指標です。小さな残差だけで物理的収束を保証するものではありません。

`state.array("scf_history")` は同じ4行×反復数の読み取り専用配列です。
履歴はFortranが保持し、メモリsnapshot・ディスク保存・再開で引き継ぐため、再開前後を
1本の図にできます。保存を伴わず図だけ必要なら `state.write_scf_report(new_directory)` を
呼びます。再開後にこの関数を呼ぶと、それまでの全履歴を出力できます。

描画には任意依存のMatplotlibを使います。導入されていない場合や描画が失敗した場合も、
警告だけにしてチェックポイントとCSVを保持します。描画失敗をSCF計算失敗と表示しません。

### 共通APIからのディスク保存と保護付き復元

`state.scf_disk_capability` は保存形式・実際のMPI分割・物質／入力／ランタイム識別を
宣言します。ファイルは読みません。`state.checkpoint_target(directory=path)` で
`execute_checkpoint()` と `restore_checkpoint()` に接続できます。
メモリ用の引数snapshotとdirectoryは同時には指定しません。

```python
from salmonpy.checkpointing import execute_checkpoint, restore_checkpoint, CheckpointState
# a,bは同じsolver_comm上にある独立したSCF状態。各rankが同じ順で呼ぶ。
groups = [{"name": "scf", "ranks": list(range(comm.size)), "point_ids": ["A", "B"]}]
targets = {"A": a.checkpoint_target(directory=checkpoint_a),
           "B": b.checkpoint_target(directory=checkpoint_b)}
saved = execute_checkpoint(comm, groups, states=targets)
# 必要な計算を進めた後、保存した位置へ戻す。
live = CheckpointState(states={"A": a, "B": b})
prepared = restore_checkpoint(comm, groups, states=targets,
                              saved=saved["compatibility"], live=live)
for candidate in prepared["states"].values():
    candidate.publish()  # 全候補を公開してから計算を再開する。
```

ファイルのSHA256・入力を照合し、各状態の独立した数値候補へ読み込んでから、共通APIが
全体の成功を確認します。途中失敗では候補を破棄し、元の状態はそのまま継続できます。
失敗後に再試行するときはcheckpoint_target()を作り直します。保存形式と読み書きは
バックエンド側に任せる、既存の共通契約を使います。

1状態だけなら `state.restore_scf(directory)` も使えます。複数状態への直接呼出しは
全体として原子的ではないため、まとめて保護する場合は上のrestore_checkpoint()を使います。
保存・復元はSCF反復を実行せず、通信・FFT・Libxc資源を共有コピーしません。
公開時の致命的なメモリ不足・プロセス消失まで巻き戻せるという保証ではありません。

サンプルは `--disk --save-scf NEW_CHECKPOINT --protected` で共通APIのディスク経路を
実行します。--diskを省いた既存のメモリ経路は維持します。

### 復元候補の後片付け

native SCFの準備済み候補は `candidate.close()`（discardと同じ）とwith文に対応します。
publish前にwith文を抜けると候補だけを破棄し、元の数値状態を保持します。
publish済み／破棄済み候補のcloseは何もしません。close後のpublishはできません。

Pythonが保持するGS状態をcloseすると、その状態が所有する未公開候補も先に破棄します。
メモリsnapshotをcloseすると、そのsnapshotから作った未公開候補を先に破棄します。
これにより、候補準備後に利用者コードの例外でwith文を抜けても、状態とMPI資源を閉じられます。
別のsnapshotや外部Cホストに属する候補は、このsnapshotの終了では破棄しません。

各MPIランクは同じ順でclose/publish/discardを呼びます。ガベージコレクションから
MPI操作を実行する機能ではありません。C ABIで直接準備した候補はそのホストの責任で
破棄してください。Pythonは外部ホストの予約を勝手に解放せず、native側の拒否を保持します。

### k点MPI数を変えて再開する

新しい対応ランタイムでは、PZ／Libxc-PZ／Libxc-PBEと高速PR版HSE06の周期・非スピン分極・固定イオンSCFを、
k点だけで並列化した配置間で再開できます。確認済みは1→2ランクと2→1ランクです。
通常のrestartは従来どおり同じMPI分割を要求し、変更時だけ明示的に指定します。

```python
# 保存時と同じcalculation/scf/material_id、同じランタイムを使う。
# 新しいMPI群はdecomposition=(comm.size,1,1,1,1)で開く。
with session.run(calculation, scf=scf, material_id="Si",
                 restart=checkpoint_directory, repartition=True) as state:
    assert state.scf_status["executed_iterations"] == 0
    state.continue_scf(steps=5)
```

新規保存のmanifestにはrankごとのglobal k範囲を記録します。入力のうちk並列数だけを
変更し、セル・格子・k点集合・電子数・XC・混合条件・閾値・反復上限・資産は同じにします。
保存形式、ランタイムSHA256、物質ID、旧k範囲の重複・範囲間の欠落、全rankファイルのSHA256を
作成前に検査します。Fortranが必要な旧rankファイルを読み、波動関数とCG作業配列の
k部分を再配置します。密度・混合履歴・Poisson／XC数値状態・残差履歴は保持し、
継続用の並列設定を新しい実行資源に結び直します。新しい担当k範囲がすべて埋まることは
公開前にnative側でも検査します。保存元のファイルは変更しません。

読込直後の波動関数・密度・占有・固有値・残差履歴は保存値と完全一致しました。
続く5反復は、保存元と同じ初期状態からの連続8反復とrtol1e-9／atol1e-10で比較しています。
既定のGaussian初期化では乱数seedに担当開始k点が含まれ、MPI分割によって初期軌道が
変わる場合があります。初期化自体は変更せず、別の初期状態を比較基準にしません。これは実装の検証であり、格子・k点などの物理収束を認定するものではありません。

高速PR版HSEでは、密なACE因子・source/actionキャッシュ・k点別factor/PT-CN履歴も
Fortranで再配置します。全k点／縮約k点とPBE準備段階の保存で1↔2を検証しています。
HSE空間分割は後述の専用入口を使います。packed/distributed ACEと
異なるランタイムへの移行は対象外です。
共通checkpoint_target／保護付き復元は、引き続き同一MPI分割で使います。
必要なABIはsp_gs_scf_k_repartition_protocol()の通常XC版1、HSE版2です。
暗黙にはビルドしません。

既存サンプルでは、保存側に `--k-parallel --save-scf NEW_CHECKPOINT`、
異なるランク数の再開側に `--restart CHECKPOINT --repartition` を指定できます。
通常XCで両側とも同じ2×2×2 k格子を使い、その他の入力も揃えます。

### 通常XCの軌道並列数を変えるSCF再開

同じ `run(restart=..., repartition=True)` で、純粋な軌道並列配置にも対応します。
両側の `decomposition=(1, comm.size, 1, 1, 1)` と、追加ABI
`sp_gs_scf_o_repartition_protocol()` の版1または2を使います。新規保存はglobal軌道範囲を
記録します。これがない旧保存も同一配置では読めますが、軌道配置変更には使えません。
物理入力・資産・ランタイムは同一に保ちます。PZ／Libxc-PZ／PBEの1↔2ランクで、
復元直後の軌道部分の完全一致、SCF再計算ゼロ、継続結果の既定rtol1e-9／atol1e-10
での一致を確認しました。Fortranが軌道とCG作業配列を再配置し、密度・混合・XC・
Poisson・残差履歴を保持します。実軌道・複素軌道ともnative側で扱い、
Pythonのローダーでは数値配列を並べ替えません。

サンプルの保存は `--orbital-parallel --save-scf NEW_CHECKPOINT`、再開は
新しいランク数で `--orbital-parallel --restart CHECKPOINT --repartition` です。
軌道配置は明示的な対称性指定がなければ固定Gamma格子を使います。高速HSEの軌道並列は後述の専用入口を使います。
通常XCの空間分割は後述の入口を使います。k並列と軌道並列の相互変更、
通常XCの混合配置は後述の専用入口を使い、別ランタイムへの移行は未対応です。
通常の復元と共通の保護付き復元は、従来どおり同一MPI配置で使います。

### 高速PR版HSEのGamma軌道並列SCF：同じ配置への直接再開

`sp_gs_hse_orbital_protocol()` の版1または2 を宣言する新しい明示ビルド版では、
`decomposition=(1, comm.size, 1, 1, 1)` とGamma点のみの高速PR版HSE06 GSを
保持できます。既存Gamma交換経路に従い、CPU・周期・非スピン分極・固定イオン・
直交セル・整数の全占有スピン対を対象とし、対称性・局在化・ペア省略は使いません。
既存のk点並列HSE経路はそのまま使えます。

`save_scf()` と `run(restart=...)` は同じランタイム・入力・MPI分割で使います。
軌道／CG、密度／混合履歴、ACEに加え、Gamma交換の進行・gauge・source作業状態も
各状態が所有し、そのまま保存します。バックエンド関数ポインタと通信IDは保存せず、
通信資源は復元先に結び直します。保存状態の再構築にSCF反復やACE再計算は使いません。

既存サンプルは `--xc hse06 --orbital-parallel --save-scf NEW_CHECKPOINT` で保存し、
同じランク数で `--xc hse06 --orbital-parallel --restart CHECKPOINT` として再開します。
ランク数を変える場合は後述のprotocol2入口を使います。複数k点の軌道HSE、
HSE空間GSは後述の専用入口を使います。混合配置、軌道並列GSからRTへの受け渡しは未対応です。2ランクの別プロセス再開と
PBE準備段階の保存は実行の検証であり、物理的な格子・k点収束を認定するものではありません。

### 高速PR版HSEの軌道並列数を変える

保存側と再開側の両方に、固定設定 `hse_orbital=True` を指定します。1ランクでも
同じGamma軌道交換の計算経路を選び、既定のk点並列経路と区別するための指定です。
HSE-orbital protocol2と、配置変更用orbital-repartition protocol2が必要です。
既存の既定経路・同一配置への再開は維持し、物理入力とランタイムは揃えます。

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, comm.size, 1, 1, 1), hse_orbital=True)
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortranが保存したglobal軌道番号に沿って、密なACE因子とsource/action/work、
交換用source・previousの列を再配置します。共有保持するgauge・占有・フェーズ・
密度混合履歴は維持し、通信資源を復元先へ結び直します。保存状態を作り直すための
ACE評価やSCF反復は行いません。保存元のファイルは変更しません。
別プロセスの1↔2、PBE準備段階の1↔2、2原子・4バンドの不均等な2↔3を
検証しました。軌道の復元直後は完全一致、初期SCF実行数は0、その後の5反復は
保存元と同じ初期状態からの連続8反復と既定rtol1e-9／atol1e-10で一致しました。

サンプルは両側に `--xc hse06 --orbital-parallel` を指定し、保存側に
`--save-scf`、異なるランク数の再開側に `--restart CHECKPOINT --repartition` を
付けます。1ランクの保存元にも明示的な軌道経路の指定が必要です。packed/distributed
ACE、タイル形式gauge、局在疎表現、複数k点の軌道HSE、空間／混合配置、別ランタイム
への移行と軌道並列GS→RTは、この密なGamma経路の対象外です。
共通の保護付き復元は従来どおり同一MPI配置で使います。

### 通常XCの空間並列数を変えるSCF再開

`sp_gs_scf_r_repartition_protocol()` の版1または2 の対応ランタイムでは、周期・非スピン分極・
固定イオンのPZ／Libxc-PZ／PBEを、純粋な空間分割間で直接再開できます。各ホストを
`decomposition=(1, 1, nx, ny, nz)`、`nx*ny*nz == comm.size` で開きます。
並列設定だけを変え、ランタイム・物質ID・物理入力・資産・SCF上限・混合条件は揃えます。

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, 1, comm.size, 1, 1))
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

新しいmanifestは各rankの格子内部範囲 `grid_bounds` を3組の始点・終点として記録し、
範囲・重複・欠落・入力／ランタイムの一致・全ファイルのハッシュを検査します。
Fortranがglobal格子座標に従って密度・ポテンシャル・軌道／CGの内部値、Fourier密度、
全ての密度／ポテンシャル混合履歴を再配置します。周期境界の通信領域は保存した
内部値をコピーして埋め、軌道を時間発展させません。混合の世代・残差・占有・
スカラー状態も保持します。比較するのは格子内部の値で、古い作業用haloの一致は要求しません。

既存のSCF実行ゼロのbootstrapで、新分割の静的格子・変換・通信資源を準備します。
FFTの部分変換バッファとXCの再利用作業領域は履歴ではないため、新しい配置のものを
使います。保存値の再構築に初期軌道生成・密度評価・固有状態計算・SCF反復は行いません。
格子内部の完全な被覆と全rankの成功を確認してから公開し、途中読込失敗では元の状態を保持します。

実軌道／複素軌道の1↔2、x/y/z方向の変更、16³格子の2↔4、
2×2×1→1×1×2分割を検証しています。内部値の復元直後は完全一致し、
続く5反復は保存元と同じ初期状態からの連続8反復とrtol1e-9／atol1e-10で一致します。
これは実装の検証であり、物理収束の認定ではありません。サンプルは保存側に
`--spatial-parallel --save-scf NEW_CHECKPOINT`、ランク数を変えた再開側に
`--spatial-parallel --restart CHECKPOINT --repartition` を指定します。
grid_boundsがない旧保存も、同一配置では復元できます。HSEの空間分割は後述の専用入口を使います。
通常XCの混合配置変更は後述の専用入口を使い、別ランタイムへの可搬化は未対応です。
共通の保護付き復元は引き続き同じMPI分割で使います。

### 高速PR版HSEの空間並列数を変えるSCF再開

保存側・再開側の両方に、固定設定 `hse_spatial=True` を指定します。1ランクでも
同じCPU・全支持Gamma交換経路を選ぶ指定で、hse_orbitalとは同時に使いません。
HSE-spatial compute protocol1とspatial-repartition protocol2が必要です。
物理入力・ランタイム・資産・物質ID・SCF設定は揃えます。

```python
backend = NativeGSBackend(build_record, work_root,
    decomposition=(1, 1, comm.size, 1, 1), hse_spatial=True)
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortranが保存したglobal格子座標とx方向が最速の配列順序に従って、密なACE因子、
交換用source/action/work、spatial source/previousの格子行を再配置します。
共有保持するgauge・占有・フェーズ・混合履歴は引き継ぎ、通信資源を結び直します。
保存状態の再構築にACE評価・固有状態／密度計算・SCF反復は行いません。
新分割の静的資源と履歴を持たないFFT／XC作業領域は通常XCの空間再開と同じ扱いです。
保存元データは変更しません。

HSEの別プロセス1↔2、PBE準備段階の1↔2、16³格子の2↔4、
2×2×1→1×1×2分割で成功しました。格子内部の復元直後は完全一致し、
SCF実行数0、続く5反復は保存元と同じ初期状態からの連続8反復と
既定rtol1e-9／atol1e-10で一致します。途中読込失敗でも元の状態を保持し、
復旧後の継続が可能です。これは実装の検証であり、物理収束の認定ではありません。

サンプルは両側に `--xc hse06 --spatial-parallel` を指定し、保存側に
`--save-scf NEW_CHECKPOINT`、ランク数を変えた再開側に
`--restart CHECKPOINT --repartition` を付けます。既定のk／軌道経路と
従来CLIの挙動は維持します。対象はCPU・周期・非スピン分極・直交セル・固定イオン・
全占有の整数スピン対・Gamma点・全支持・密なACE・タイルなしgaugeです。
packed/distributed ACE、タイルgauge、局在疎表現、複数k点／混合配置、
別ランタイム移行と空間HSE GS→RTは対象外です。

### 通常XCのk点・軌道・空間を組み合わせた混合並列での再開

`sp_gs_scf_m_repartition_protocol()==1` の対応ランタイムでは、PZ／Libxc-PZ／PBEの
純粋分割から混合分割、混合分割同士を、同じ明示的な
`run(restart=..., repartition=True)` で再開できます。decompositionは
`(k, orbital, x, y, z)` の順で、積を新しいホストMPI群のランク数と一致させます。
物理的な格子・k点集合・電子数／バンド数・SCF／混合条件・物質ID・資産・
ランタイムは維持し、k点数をランク数から決め直さないでください。

```python
# 保存側の例：(2, 2, 2, 1, 1)、8ランク。
# 再開側の例：(1, 2, 1, 2, 1)、4ランク。
backend = NativeGSBackend(build_record, work_root,
                         decomposition=(1, 2, 1, 2, 1))
with backend.open(context) as session:
    with session.run(calculation, scf=scf, material_id="Si",
                     restart=checkpoint, repartition=True) as state:
        assert state.scf_status["executed_iterations"] == 0
        state.continue_scf(steps=5)
```

Fortranが各保存rankをglobal x/y/z・軌道・k点の5次元領域として扱い、
grid_bounds／orbital_bounds／k_boundsから軌道とCG配列を再配置します。
周期haloも内部値からコピーします。密度・ポテンシャル・Fourier密度・混合履歴は
元からk点／軌道群の間で複製されているため、各空間領域から一つだけ取り込み、
重複して加算しません。軌道とスカラー場を別々に被覆検査し、全rank成功後に公開します。
静的資源と履歴を持たないFFT／XC作業領域は空間再開と同じ扱いで、
保存状態を作り直すための密度評価・固有状態計算・SCF反復は行いません。
通常の同一配置への復元は維持します。

実軌道／複素軌道の2→8、8→4、2原子系の4→8、混合4→2、
同じ4ランクでの分割変更を検証しています。内部の軌道・密度の復元値は完全一致し、
SCF再計算0回、続く5反復は保存元と同じ初期状態からの連続8反復と
既定rtol1e-9／atol1e-10で一致します。途中ファイル読込失敗でも元の配列と
保存元を維持します。独立した混合状態と既存ランタイムとの通常計算の一致も確認します。
物理収束の認定とは区別してください。

既存サンプルは両側に `--decomposition K O X Y Z` を指定し、
物理的なk格子は固定2×2×2を使います。保存側は `--save-scf NEW_CHECKPOINT`、
新しいランク数・分割の再開側は `--restart CHECKPOINT --repartition` です。
HSEの混合分割移行・異なるランタイムへの可搬化・共通保護付き復元でのMPI分割変更は
この入口の対象外です。
