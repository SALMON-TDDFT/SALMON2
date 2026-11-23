# オフライン計算ノード向け入力の準備

[English](PREPARATION.en.md)

施設が許可する管理・ログインノード、または手元PCで準備します。
`salmonpy.preparation` は明示した**ローカルファイル**を集め、登録済み検証関数で
物質設定を確認し、照合用manifestを作ります。取得・導入・ビルド・ジョブ投入は行いません。

## Pythonからの利用

```python
from salmonpy.preparation import prepare_inputs, verify_inputs

# recipeとregistryはCONFIG.mdに従って用意する。
# 検証関数は元のcwdでなく、渡された作業用base_directoryを基準に参照を検証する。
bundle = prepare_inputs(
    recipe, "new-job-inputs", registry=registry, base_directory="downloaded",
    files={"structures/cell.cif": "my-supercell.cif"},
    provenance={"structures/cell.cif": {
        "source": "user-provided supercell", "license": "user-managed"}}
)
verify_inputs(bundle)
```

`files` は「転送先の相対パス→手元のファイルパス」です。手元の相対パスは
`base_directory`が基準で、絶対パスも指定できますが自動記録はしません。
物質パラメータには転送先の相対パスを指定します。構造・擬ポテンシャル・初期状態と
その参照先を含め、必要なファイルをすべて列挙してください。ディレクトリ、symlink、
上位パスへの脱出、予約メタデータ名、大文字小文字の衝突、既存保存先は拒否します。
欠損・非通常ファイルは検証関数や出力作成に入る前にまとめて報告します。依存の再帰探索は行いません。

`provenance`はファイルごとの任意JSONメタデータで、省略時は出典不明と記録します。
分かる範囲で提供元ID・版・ライセンスを記録します。記載内容の真偽を認証する機能ではありません。
認証情報を入れないでください。今回のヘルパー自身はCIF/ASE構造の読込や格子幅の解決を行わず、
専用アダプターは末尾のリンクを参照してください。自作検証関数で明示的に行うことはできますが、その仕様の検証が必要です。
任意パラメータに隠れたファイル依存を収集側が判定することはできません。

出力は`materials.json`（元レシピ）、`resolved.json`（採用設定カタログ）、指定位置のコピー、
および他の全ファイルのSHA256・バイト数・出典を記録した`manifest.json`です。
収集側がsymlinkや元機種のパスを挿入することはありません。コピー対象は内容で、実行権限や
仮想環境の保存ではありません。wheelや計算プログラムZIPも明示的に同梱できますが、依存解決や
OS/Python/MPIの適合確認は別です。[RUNTIME.md](RUNTIME.md)に従ってください。

信頼する検証関数は一時作業ディレクトリで動きます。一時的な絶対パスを戻り値に保持せず、
元ファイルのパスに依存せず、コピーを書き換えず、計算起動・データ取得をしないでください。
戻り値は移設可能な相対参照と採用値にします。関数の副作用の隔離・巻き戻しは行いません。
検証後にだけ保存先を作り、失敗時は一時領域を除去し、既存出力を上書きしません。

## 転送後の確認

名前を保ってディレクトリ全体を転送します。準備時に`manifest.json`のSHA256も別途記録してください。
計算ノードでは次を実行します。

```sh
python -m salmonpy.preparation transferred-inputs --manifest-sha256 RECORDED_SHA256
```

Pythonでは`verify_inputs(directory, manifest_sha256=...)`です。検証はローカルファイルだけを読み、
登録計算処理をimport・呼出ししません。ファイル・ディレクトリの欠損/余剰、symlink、
不正manifest、ハッシュ・サイズ不一致を拒否します。戻り値は独立したmanifest辞書で、
CLIは整合性確認だけの短い結果を表示します。内部整合性だけを調べる場合はhash指定を省略できます。
manifest自体も改変された場合、その内部hashだけでは認証できません。別途記録したhashは信頼する経路で保持します。
OSエラーはそのまま、不正な内容は`InputError`として報告します。

成功は科学入力・依存の完全性、物理的妥当性、ソルバー実行成功を保証しません。
選択計算処理で読み込み・検証し、計算先で小規模確認してください。
出力先は入力一式とは別の新規ディレクトリにします。入力内へ結果を追加すると全件照合で拒否します。
任意のユーザーコードの通信を遮断する機能ではありませんが、検証器自身は通信しません。

## 実行可能なサンプル

ローカルの`coefficient.json`へ`{"coefficient": 2.0}`を書いて実行します。

```sh
python -m salmonpy.examples.prepare_inputs coefficient.json new-inputs
python -m salmonpy.preparation new-inputs
```

自作検証関数を明示登録して、合成のスカラー応答ファイルをまとめる例です。
入力準備と整合性の例であり、Si/DFT計算ではありません。

構造・格子を事前検証する明示登録例は[STRUCTURE.md](STRUCTURE.md)と`examples.prepare_structure`です。汎用収集器が構造を自動解釈するわけではありません。

## 計算と分離したコード・入力検査（任意）

ユーザーが施設の規則に従い、管理・ログインノードなどで独自に実行します。

```sh
python -m salmonpy.preparation input.json --code calculation.py --report check.json
```

準備済み入力ディレクトリなら`input.json`の代わりにそのパスを渡し、必要に応じて
`--manifest-sha256 RECORDED_SHA256`を指定します。Pythonでは`check_job(code, inputs, manifest_sha256=None)`です。

このコマンドは**成功しても検査して終了するだけ**です。渡した計算コードのimport・実行、import先の追跡、
材料プラグインの呼び出し、MPI起動、ビルド、ソルバー実行、ジョブ投入は行いません。
レポートは新規ファイルに保存し、既存のものを上書きしません。計算は別途、スケジューラを通じて割り当てられた計算ノードで実行します。
検査済みレポートは計算の必須条件ではありません。混在計算サンプルの`--input-check`と検査パイプラインは取り除きました。
検査と計算の混同を避けるため、従来のMPI事前検査サンプルも除去しました。各計算処理の通常の入力検証と異常処理は維持します。

### 対象を限定する

- Pythonコード：構文・コンパイル可否とファイルhashだけ。コンパイルしたコードを実行したり、import先を確認したりしません。
- JSONファイル：構文・重複キー・有限値だけ。科学パラメータのスキーマを推測しません。
- 準備済み入力一式：`verify_inputs`によるファイル一覧・サイズ・hashの照合。

レポートは`checks`、`errors`、`not_checked`を分けます。終了コード0／passedはこれらの検査だけが通った意味で、
管理ノードでの計算許可、実行成功、物理的妥当性を保証しません。コードがその入力を使うか、動的なパス、
任意の自作モデルの組合せ、計算処理 namelistの意味、計算プログラム・MPI・GSの互換性、収束・メモリ容量は保証しません。
成功後の自動修正・計算実行はありません。コードと入力の問題は独立に記録し、検出した問題があれば終了コード2を返します。

あらゆる組合せをこの検査器で覆うことは目指しません。今後の科学的検査は対応するスキーマ・計算処理と範囲を明示し、
未知の組合せは「禁止」ではなく「未確認」とします。既存の`preflight_inputs(..., registry=...)`は信頼する検証関数を明示実行する別の入口です。
`--code`と`--resolve-builtins`は併用せず、この検査経路から任意の検証コードが動くことを避けます。


## FDTDを準備だけで終了する（管理ノード向け）

`prepare_fdtd(calculation, new_directory, mpi=None)`は計算処理形式のinputfileと
材料ファイルを新しいディレクトリへ保存するだけです。実行ファイルは不要です。
サブプロセス、MPI初期化、ビルド、計算、ジョブ投入は行いません。戻り値はPathで、
実行メソッドはありません。検査に通った後も自動で計算へ進みません。

```python
from salmonpy.preparation import prepare_fdtd
from salmonpy.fdtd import FDTDMPI
from salmonpy.examples.prepare_fdtd import calculation  # 小さな接続確認用LD例

folder = prepare_fdtd(calculation(), "fdtd-input", mpi=FDTDMPI((2, 1, 1)))
# ここで終了。転送・投入・実行は別操作です。
```

```sh
python -m salmonpy.examples.prepare_fdtd NEW_DIRECTORY --mpi-grid 2 1 1
```

`salmonpy-prepared-fdtd.json`はstatus=prepared-not-run、MPI分割、ランク数、各入力の
SHA256とサイズを記録します。ランチャーや実行ファイルの環境依存パスは保存しません。
既存ディレクトリを上書きしません。これはSALMON用の入力ファイル群で、
prepare_inputsの物質アーカイブとは別形式です。verify_inputsをこの形式へ適用する
ものではありません。全入力の物理収束や、選んだ計算プログラムとの互換性の保証でもありません。

準備先を必要に応じて計算ノードへ転送し、施設のジョブスクリプトで、同じ機能を
持つ監査済みSALMONを明示的に起動してください。実行時は記録済みの実行ファイル
SHA256を照合し、ソース／オブジェクト監査を繰り返しません。計算先のディレクトリを
作業場所にし、stdinへinputfileを渡します。例えばSlurmの割当内での起動部分は
` srun --ntasks=2 "$SALMON_EXECUTABLE" < inputfile > stdout.log 2>&1 `です。
これはサイト依存の例で、管理ノードで直接実行する指示や自動投入ではありません。
キュー・課金先・時間・スレッド数などは施設の規則に従ってジョブスクリプトへ指定します。

## サイト固有のランチャー引数

FDTDMPIの既定rank_argumentsは('-n','{ranks}')で、従来の起動を維持します。
command(executable)は引数タプルを返すだけで、実行・ファイル確認をしません。

```python
from salmonpy.fdtd import FDTDMPI
mpi = FDTDMPI((2, 1, 1), launcher=("site-launch",),
              rank_arguments=("--tasks={ranks}",))
argv = mpi.command("/site/salmon")  # プレビューだけ
# ('site-launch', '--tasks=2', '/site/salmon')
```

rank_arguments=()ならランク数の引数を追加しません。割当から数を決めるサイト向けで、
実際の割当数は入力のMPI分割と一致させる必要があります。任意の引数を自動的に
推測・検査する機能ではありません。{ranks}以外は文字列のままで、シェル解釈や
環境変数展開を行いません。jsrunなどは、資源集合数・集合内タスク数とMPIランク数の
対応をサイトに合わせて指定してください。単にsrunへ名前を変えれば済むとは限りません。

SalmonExecutable.runは、直列Pythonが明示的に子プロセスを起動する従来入口です。
MPIランク内の二重起動拒否は残しています。Python自体をMPIアプリにする場合は、
COUPLED.mdのParallelContextと計算処理同じプロセス内での計算部品を集合的に使ってください。
DFT／TDDFTのバッチ投入経路をこの追加だけで共通化したとは扱わないでください。


## DFT／TDDFTも準備だけで終了する

`prepare_dft(model, new_directory, scf=None, save_restart=True, initial_guess=None)`は、DFT.input_textと
同じ入力・擬ポテンシャル・対称性ファイルを保存し、Pathを返します。SCFもプロセスも
起動しません。既定のsave_restart=Trueは将来実行するDFTでGSを保存・再利用する設定で、
準備時にGSを作る意味ではありません。Falseなら将来の再開データ出力を無効にします。
既存DFT.runとwarm startの数値処理は変えません。未収束時は結果を残し、残差確認のため
後続処理を一旦止めます。既定は新規入力で、initial_guessを明示した場合だけ初期推定を準備します。

```python
from salmonpy.preparation import prepare_dft
from salmonpy.dft import SCFSettings

folder = prepare_dft(model, "dft-input", scf=SCFSettings(max_iterations=300))
# modelはDFTまたはASE連携で作成済みのモデル。ここで終了します。
```

`prepare_tddft(rt_input_file, new_directory, files=None)`は、指定した計算処理 RT入力を
バイト単位で保持し、filesで明示したローカルファイルだけをコピーします。dt、波形、
ゲージ、並列分割、directory_read_dataなどを自動変更しません。GSや実行ファイルの
存在を管理ノードへ要求せず、GSを自動探索・読込・複製しません。

```python
from salmonpy.preparation import prepare_tddft

folder = prepare_tddft("rt.inp", "tddft-input",
                       files={"pseudo1_rps.dat": "Si_rps.dat"})
# GS参照はrt.inpに記述したまま。投入・計算は別操作です。
```

filesの名前は入力内の参照と合わせてください。依存ファイルを自動推定しません。
明示した大きなファイルもストリームコピー・ハッシュ計算し、GSを転送する場合は
利用者が対象ファイルを明示して管理します。入力内の既存GS参照は計算ノード側で
解決してください。設定、GS形式／XC／格子・擬ポテンシャルの整合性、実行計算プログラム、
サイトの並列条件や収束を認証する入口ではありません。

各ディレクトリにはsalmonpy-prepared-dft.json／salmonpy-prepared-tddft.jsonを保存し、
status=prepared-not-runとファイルのSHA256・サイズを記録します。これは
prepare_inputsの物質アーカイブとは別形式です。既存ディレクトリは上書きしません。

```sh
python -m salmonpy.examples.prepare_electrons dft Si_rps.dat NEW_DIRECTORY
python -m salmonpy.examples.prepare_electrons tddft rt.inp NEW_DIRECTORY --asset pseudo1_rps.dat Si_rps.dat
```

このサンプルも準備で終了します。runやsubmitの自動継続はありません。DFTの例は
小さなSi入力の作り方で、物理収束を保証しません。計算は割当済みジョブ内で、
同じ機能を持つ監査済みSALMONへgs.inpまたはrt.inpをstdinから渡して明示的に起動します。
実行ファイルの記録済みSHA256を照合し、ソース／オブジェクト監査を毎回繰り返しません。
HSEを使う場合はHSEのみで、prepare_dftはそのprovenance／capability検査を
置き換えません。従来のDFT.runからHSE batch実行を解禁する変更ではありません。

PythonのA→Jワーカーを使う場合も、準備とは別です。TDDFTInitialState.from_salmonは
既存GS束のsnapshot loaderで、new_stateはワーカーを起動します。管理ノードの
入力準備だけを意図した時に、これらを自動で続けて呼ばないでください。
caller-owned MPI内では既存の計算処理同じプロセス内での計算経路を使います。

### 初期推定を明示して準備だけ行う

`prepare_dft(model, "new-prepared", initial_guess=previous_result)` は、
DFT.runと同じ整合性検査と入力生成を使います。info.bin・occupation.bin・wfn.bin
の3ファイルだけをinitial_guess/へ複製し、公開前にコピーのハッシュを照合します。
元GSの記録とrequired_executable_sha256も残します。元GSは変更せず、実行ファイルを
読み込まず、SCFやジョブ投入にも進みません。計算ノードで記録された実行ファイルの
SHA256を選択・確認してから、別のジョブとして起動してください。

セル・格子・k点・元素・擬ポテンシャル・電子数・状態数・XCは一致が必要です。
原子位置は既存warm start同様に変更可能です。波動関数の初期推定であり、混合履歴を
含む厳密なSCF再開ではありません。他形式や自作バックエンド、HSEの生存状態間の
初期推定移行は別作業です。prepare_electronsの例には明示的な
`--initial-guess GS_DIRECTORY`を追加しました。例と同じ物質モデルのGSを指定します。
