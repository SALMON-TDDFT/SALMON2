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

新しく明示ビルドしたCPUランタイムでは、周期・非スピン分極・対称性なし・
固定イオンのPZ／Libxc-PZ/PBEと高速PR版HSE06でSCF protocol1を利用できます。従来のrun()は
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
保存は同じ生存状態・ビルド・MPI分割でのメモリ内保存です。対称性、SOC/+U/DC、
イオン運動、ディスク保存／プロセス再開は今回の拡張に含みません。

```sh
mpiexec -n 2 python -m salmonpy.examples.scf_resume HSE_BUILD_JSON NEW_OUTPUT --xc hse06
```

密度・占有数・軌道・固有値・力・メタデータが、新旧ランタイムの連続実行と
一致することを検証しています。PBE→HSE切替の前後での復元、異なる物質を
挟んだ独立性、未収束時の警告と保持、GSから続けた16ステップのTDDFT電流も
確認しました。実装の回帰検証であり物理的収束の検証ではありません。
GSの波動関数・密度データは収集していません。
