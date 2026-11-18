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
持たないため未対応で、native GSハンドルの復元も残作業です。
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
