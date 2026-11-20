# オフライン計算ノード向け入力の準備

[English](PREPARATION.en.md) · Topic: `offline_inputs` · Paired revision: `2026-11-15.1`

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
専用アダプターは末尾のリンクを参照してください。自作検証関数で明示的に行うことはできますが、その契約の検証が必要です。
任意パラメータに隠れたファイル依存を収集側が判定することはできません。

出力は`materials.json`（元レシピ）、`resolved.json`（採用設定カタログ）、指定位置のコピー、
および他の全ファイルのSHA256・バイト数・出典を記録した`manifest.json`です。
収集側がsymlinkや元機種のパスを挿入することはありません。コピー対象は内容で、実行権限や
仮想環境の保存ではありません。wheelやランタイムZIPも明示的に同梱できますが、依存解決や
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
登録バックエンドをimport・呼出ししません。ファイル・ディレクトリの欠損/余剰、symlink、
不正manifest、ハッシュ・サイズ不一致を拒否します。戻り値は独立したmanifest辞書で、
CLIは整合性確認だけの短い結果を表示します。内部整合性だけを調べる場合はhash指定を省略できます。
manifest自体も改変された場合、その内部hashだけでは認証できません。別途記録したhashは信頼する経路で保持します。
OSエラーはそのまま、不正な内容は`InputError`として報告します。

成功は科学入力・依存の完全性、物理的妥当性、ソルバー実行成功を保証しません。
選択バックエンドで読み込み・検証し、計算先で小規模確認してください。
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

## ユーザーが任意で行う投入前検査と開始時のエラー記録

投入前検査は、ユーザーが施設の規則に従って管理・ログインノードなどで独自に実行する任意の手順です。
SALMONpyがそのノードへ接続して自動実行したり、検査レポートをジョブ投入の必須条件にしたり、ジョブを投入したりはしません。
この手順を省いても通常の計算は実行できます。ユーザーがファイルと材料設定を確認したい場合のコマンドです。

```sh
python -m salmonpy.preparation new-inputs --resolve-builtins --report preflight.json
```

`--resolve-builtins`は組み込み材料のレシピ用です。自作バックエンドではPythonから
`preflight_inputs(bundle, registry=your_explicit_registry)`を呼びます。
結果は`status`、エラーの段階・種類・メッセージ、manifestのhash、検査範囲を持つJSON化可能な辞書です。
検証関数は前述の「計算起動・ダウンロードをしない」契約を守ります。registryを省くとファイル整合性だけを確認します。
ランタイム・MPIの適合、native namelistの意味、GS互換性、収束、メモリ容量の保証ではありません。

CLIはJSONを表示し、成功は終了コード0、不正入力は2を返します。`--report`は入力一式の外側に新規ファイルを作り、
既存レポートは上書きしません。記録失敗もコマンドの失敗になります。投入したジョブと一緒にレポートを保管してください。
CLIの表示形式は変わりますが、`verify_inputs`自体の例外とmanifest辞書を返す契約は維持します。

転送後、計算開始前にも確認するジョブスクリプトの例です。

```sh
set -eu
python -m salmonpy.preparation transferred-inputs --manifest-sha256 RECORDED_SHA256 --report job-input-check.json
# 上のコマンドが成功した場合だけ、ここで計算を起動します。
```

各ランクから見えるファイルを確認する実行可能サンプルも用意しました。

```sh
mpiexec -n 2 python -m salmonpy.examples.mpi_input_preflight transferred-inputs job-ranks-check.json --manifest-sha256 RECORDED_SHA256
```

このサンプルは入力検査だけで、ソルバーを起動しません。全ランクの結果を集約してから成否を決め、同じ終了コードで終了します。
一部ランクのファイル欠損でも、ジョブ全体でソルバーの作成へ進みません。本番スクリプトでも状態作成の前に同じ流れを置きます。
不正入力は記録を見てジョブ外で修正し、再投入します。自動入力修正や対話的な再試行は必要ありません。
SCF未収束や助言用メモリ予算超過は従来どおり警告とし、新しい数値停止しきい値は設けません。

### 通常の混在計算への接続

`examples.mpi_mixed_native`に任意の`--input-check STAGED_INPUT_BUNDLE`と
`--input-manifest-sha256 RECORDED_SHA256`を追加しました。worldの全ランクがソルバー状態を作る前に
入力一式を照合し、新規出力先へ`input-preflight.json`を保存します。不正入力なら全ランクが終了コード2で終了します。
出力先の作成失敗もbarrierより前に共有し、他ランクを待たせません。既存の出力先は上書きしません。

開始時の確認はファイル整合性だけです。入力一式から計算引数を自動設定したり、別途指定する`--layers`、
ランタイム、GSがレシピと一致することを保証したりはしません。準備した参照と実行引数を揃え、材料設定の検査は投入前に行ってください。
それより後のランタイム・数値異常の扱いは、サンプルの既存の終了／トランザクション方針を維持します。

通常の混在サンプルは、`--input-check`を指定しない場合も、明示指定した層配置・モデル群・TDDFT材料カタログの
読み込み結果を各ランクで確認します。`OUTPUT/startup-inputs.json`へランク、入力の段階、パス、例外の種類と原因を残します。
読込・JSONエラーは全ランクで集約し、ソルバー状態を作る前に終了コード2で終了します。成功の記録はこの読み込みの成功だけを表し、
native検査や計算全体の成功を保証しません。管理ノードでの操作や事前レポートは不要です。その他のnative・ランタイム異常には別の方針が残ります。
