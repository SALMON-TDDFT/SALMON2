# 自作物質モデルのC ABI（プレビュー1）

[English](C_ABI.en.md)

同じ計算処理物質モデルを、Pythonから差し込む場合とC/C++/Fortranの呼び出し元プログラムから
呼ぶ場合で同じ規約を使います。公開ヘッダは[include/salmonpy_material.h](include/salmonpy_material.h)、
Fortranの宣言は[include/salmonpy_material_iso.f90](include/salmonpy_material_iso.f90)です。
Python入口は`salmonpy.c_material.CMaterial`です。

対象は固定dt・CPU・直列の端点A駆動物質モデルです。ABI1はプレビューで、
SALMON本体やLD/SBEの既存の内部Cシンボルをすべて公開する宣言ではありません。
計算処理TDDFTの共有変数、一回限りの初期化、複数コンテキスト対応は今回変更しません。

## 共通のデータと所有権

- 入力はAc=A/c、出力は電荷電流密度j。実験室座標・原子単位・3成分のIEEE754 binary64。
- モデルと時間発展状態は別々の正のint64ハンドル。0は無効。ライブラリをまたいで渡しません。
- モデルは設定をコピーして所有し、各状態は全可変履歴を所有します。入力や出力のポインタを保持しません。
- 呼出しは直列化します。ABI1には内部MPI集団処理やスレッド安全性の保証を含めません。
- 初期時刻はtime0、nステップ後はtime0+n·dtです。初期Aを明示し、初期電流は物理値を返します。
- 入力portのlagと初期Aを使うかはテーブルに宣言し、PythonのTimeSpecにも保持します。

状態バッファは呼出し側が確保します。3成分は連続配置で、入力と出力は重ねません。
モデル設定は`configuration`のバイト列で、解釈はプラグインが定義・文書化します。
共通ABIがモデル固有のパラメータ型やファイル形式を推測することはありません。

## 問合せとコールバック

`salmonpy_material_query(version, table_capacity, &table)`でABI番号とテーブルサイズを
照合します。呼び出し元プログラムは返却status、abi_version、struct_size、必須コールバック、
backend名と時刻宣言を確認してから呼びます。これは同じプラットフォームのC ABIであり、
ファイルや通信の形式ではありません。表とbackend文字列はライブラリの寿命中有効です。

| コールバック | 役割 |
|---|---|
| `model_create` | 設定をコピーしてモデルを作る。失敗時の出力ハンドルは0 |
| `model_destroy` | 状態が残っていればBUSY、なければモデルを解放 |
| `state_create` | モデル・dt・time0・初期Aから独立状態を作る。失敗時は0 |
| `state_observe` | 時間発展させず、最後に成功した計算時刻と電流を返す |
| `state_step` | 次端点Aを受け、同じ端点の時刻と電流を返す |
| `state_destroy` | 状態の全資源を解放。破棄済みハンドルは拒否 |

すべてint32のstatusを返します。0は成功、1は引数、2はハンドル、3はメモリ、
4は数値、5は使用中、6はABI、7は内部エラーです。失敗時には最後に成功した計算状態を変更せず、
出力の時刻・電流は呼び出し元プログラムが無視します。構築途中の資源は失敗した作成コールバック自身が解放します。
診断文字列は任意の呼出し側バッファへ書き、短ければNUL終端して切り詰めます。
C++例外を境界の外へ出したり、Fortranのstop/error stopで呼び出し元プログラムを終了したりしません。

## Pythonから使う

```python
from salmonpy.c_material import CMaterial
from salmonpy.materials import check_material, bind_material
from salmonpy.crystal import respond

# configurationの形式、ライブラリとSHA256はプラグインの提供者が指定します。
with CMaterial(library_path, library_sha256, configuration_bytes) as model:
    report = check_material(model, drive)
    response = respond(model, drive)
    binding = bind_material(model, key='sample', identity=model.identity)
    # bindingを使うResponseSessionは、このwithの内側で作成・解放します。
```

CMaterialはライブラリのSHA256を照合してから明示的にロードし、設定をコピーします。
identityはライブラリSHA・設定SHA・backend名から作ります。通常のimportでは
ライブラリの読込み、コンパイル、MPI初期化を行いません。
Pythonでもモデル・状態の同時呼出しは避けてください。計算処理呼出しは別名の
ライブラリを含めてプロセス全体で直列化しますが、Python側の寿命管理をスレッド共有可能とする保証ではありません。

`model.close()`は新しい状態作成を止め、既存状態が閉じるまで計算処理モデルの解放を延期します。
既存状態は引き続き使えます。必ずwithかcloseで解放し、GCによる自動解放に頼りません。
状態の`close()`は繰返し可能です。C呼び出し元プログラムは状態を先に破棄してからモデルを破棄します。

非有限なPython入力は呼出し前に拒否します。計算処理の非0 statusは、コードを持つ
`CABIError`（RuntimeFailure派生）として返し、Pythonの最後に成功した計算時刻・電流も保持します。
計算処理側にも失敗時の原子性が必要です。成功を返しながら不正な時刻や電流を返した場合は、
状態を無効にして再利用を拒否します。壊れた計算処理状態をPythonで勝手に巻き戻すことはしません。

ライブラリは信頼できる実行コードとして扱います。SHAは発行元の署名や障害の隔離ではありません。
CMaterial外で先にロードしたり、ロード後のファイルを交換したりせず、版ごとに別の不変ファイルを使います。
ABI検査が失敗した場合も、ロード済みパス・ファイル・イメージとSHAの対応を固定します。
別名やOSローダーの再利用が別ビルドのSHAと混ざれば拒否します。

## C++・Fortran呼び出し元プログラムの実行例

明示的な開発用build入口です。既存保存先や既存バイナリは変更しません。

```sh
python tools/build_c_material_example.py NEW_BUILD --cc cc --cxx c++ --fc gfortran
```

インストール済みwheelでは、`salmonpy/_runtime_source/tools/build_c_material_example.py`にあります。
ヘッダ・ISO宣言・例のソースはwheel内にも含まれます。buildにはC/C++/Fortranコンパイラが必要です。
普通のwheel利用やimportでコンパイラを要求・自動起動しません。
この例のbuild手順はGNU互換Fortran・Unix系を対象とし、今回の実測はApple Silicon Macです。

例のモデル設定は、計算処理C double二つの`{gain, decay}`です。Pythonでは
`struct.pack('=dd', gain, decay)`、Fortranでは同じ連続double二つを渡します。
例はj_new=decay·j_old+gain·(A_new−A_old)という製造解であり、実材料の構成則ではありません。
Cプラグイン、C++呼び出し元プログラム、Fortran呼び出し元プログラムを別のソース・オブジェクトから構築し、
ソースと全アプリケーションオブジェクトを禁止ルーチン検査してからリンクします。
ソースのスナップショット・ログ・ハッシュを新規build先へ残します。

実測では両呼び出し元プログラムの独立状態、設定コピー、無効・破棄済みハンドル、BUSY、
NaN入力後の状態保持が通りました。Python参照との差は0、同じCプラグインを
MPI点分散へ載せた2/3ランクでも参照との差は0でした。物理検証や計算処理TDDFTの
再入可能性を検証したという意味ではありません。
