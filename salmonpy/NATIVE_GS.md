# メモリ上に保持する標準GS

[English](NATIVE_GS.en.md) · Topic: `retained_native_gs`

`salmonpy.native_gs` は、標準DFTのGSを同じプロセス内に複数保持します。
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

Python入口は既存DFT設定と同じ周期・全k点・非スピン偏極・固定イオンのPZ、Libxc-PZ/PBEです。
高速PR-HSE、SOC／+U／DC、保存／再開、
スレッド並行や同時GS／TDDFTは今回の対応外です。1プロセスでnative電子計算の接続は1生存期間です。
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
