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
検査と計算の混同を避けるため、従来のMPI事前検査サンプルも除去しました。各バックエンドの通常の入力検証と異常処理は維持します。

### 対象を限定する

- Pythonコード：構文・コンパイル可否とファイルhashだけ。コンパイルしたコードを実行したり、import先を確認したりしません。
- JSONファイル：構文・重複キー・有限値だけ。科学パラメータのスキーマを推測しません。
- 準備済み入力一式：`verify_inputs`によるファイル一覧・サイズ・hashの照合。

レポートは`checks`、`errors`、`not_checked`を分けます。終了コード0／passedはこれらの検査だけが通った意味で、
管理ノードでの計算許可、実行成功、物理的妥当性を保証しません。コードがその入力を使うか、動的なパス、
任意の自作モデルの組合せ、native namelistの意味、ランタイム・MPI・GSの互換性、収束・メモリ容量は保証しません。
成功後の自動修正・計算実行はありません。コードと入力の問題は独立に記録し、検出した問題があれば終了コード2を返します。

あらゆる組合せをこの検査器で覆うことは目指しません。今後の科学的検査は対応するスキーマ・バックエンドと範囲を明示し、
未知の組合せは「禁止」ではなく「未確認」とします。既存の`preflight_inputs(..., registry=...)`は信頼する検証関数を明示実行する別の入口です。
`--code`と`--resolve-builtins`は併用せず、この検査経路から任意の検証コードが動くことを避けます。
