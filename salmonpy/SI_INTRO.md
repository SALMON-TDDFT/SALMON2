# Si入門：wheel導入から電荷電流まで

[English](SI_INTRO.en.md) · Topic: `si_intro` · Paired revision: `2026-11-14.2`

実際のSi基底状態をSCF計算し、指定したA/cパルスで時間発展させ、電荷電流を固定の
期待値と比較する入門例です。入力・Si擬ポテンシャル・期待値・出典はwheelに含めます。
計算本体は別のビルド済みランタイムを使います。この手順にコンパイルやソース取得は不要です。

## 1. 配布されたwheelを導入する

OS・CPU・Python版に対応するwheelを用意してください。この新しい例を含む
`salmonpy.examples.si_intro` が必要です。コミット `6f8b36c6` の既存wheelには
まだこの例が入っていません。ソースの更新だけではwheelは公開されないため、
更新済みwheelを配布者から受け取ってください。

macOS／Linuxでは次を実行します。

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --only-binary=:all: ./SALMONPY_WHEEL.whl
python -m salmonpy.doctor --json
python -m salmonpy.examples.si_intro --help
```

`SALMONPY_WHEEL.whl` は仮名です。実際のPython版・プラットフォーム名付きファイル名に
置き換えます。対応wheelがない場合にソースビルドへ進まないよう、`--only-binary=:all:`
を保持してください。この例のPython依存はNumPyのみです。Windows用計算本体はこの
入門例では提供しません。Windows版LD/SBE wheelの対応とTDDFTの対応は別です。

## 2. 対応するビルド済みランタイムを取得する

[SALMONpyランタイム試用版Release](https://github.com/otobetomohito/SALMONpy/releases/tag/runtime-2026.11.14)から取得します。[CI 43554502676](https://github.com/otobetomohito/SALMONpy/actions/runs/43554502676)で3環境が成功したコミット `8006261b9e48541b591e496a813b3f06ca5fc1b1` の配布物です。

| Machine / 環境 | ZIP |
|---|---|
| Apple Silicon Mac | [salmonpy-runtime-macos-15.zip](https://github.com/otobetomohito/SALMONpy/releases/download/runtime-2026.11.14/salmonpy-runtime-macos-15.zip) |
| Intel Mac | [salmonpy-runtime-macos-15-intel.zip](https://github.com/otobetomohito/SALMONpy/releases/download/runtime-2026.11.14/salmonpy-runtime-macos-15-intel.zip) |
| Linux x86_64 / glibc | [salmonpy-runtime-ubuntu-22.04.zip](https://github.com/otobetomohito/SALMONpy/releases/download/runtime-2026.11.14/salmonpy-runtime-ubuntu-22.04.zip) |

非公開リポジトリのため、アクセス権のあるアカウントでログインします。Release資産にはActions成果物の有効期限はありません。**外側ZIPをそのまま保持**し、自動展開されたフォルダは渡さないでください。`RUNTIME_ARTIFACT.zip` は実際のファイル名に置き換えます。単一ランク・PZ／TB-mBJ用でLibxcは含みません。準備時にOS互換性を検査します（arm64 Mac >=15.8.2、Intel Mac >=15.7.9、Linux glibc >=2.35）。ハッシュ検査を無効化したり未検証のZIPで代用したりしません。

## 3. 準備して計算する

ソースツリー外の新規保存先を選びます。準備では外側ZIPの固定ハッシュを照合し、通常の
ランタイム導入処理で内側ZIP・実行ファイル・依存ライブラリを検証します。固定入力と
擬ポテンシャルも手元にコピーします。暗黙のダウンロードはありません。

```sh
python -m salmonpy.examples.si_intro prepare new-si --runtime-artifact RUNTIME_ARTIFACT.zip
python -m salmonpy.examples.si_intro run new-si --threads 1 --memory-budget-gib 4
```

この例では**ジョブ全体で4 GiB**の許容量を記録し、ソルバー1スレッドで実行します。
4 GiBはプロセスごとの許容量ではありません。利用許可のある予算を指定してください。
`--threads` は1〜4で、他の数値ライブラリのスレッド数は例の中で制限します。
DFTとTDDFTは順番に実行します。メモリ予算は記録値であり、OSによるハード上限ではありません。
計測不能は計算失敗と分けて報告されますが、予算内の証明にはなりません。未測定の見積もりを
根拠に規模を拡大しないでください。

既存の保存先を上書きしません。再試行には別の新規保存先を使い、失敗ログを残してください。

## 4. 結果を確認する

`new-si/result/report.json` に結果を保存します。計算・期待値検証・任意の資源計測の
成否は別々に確認します。本体が正常終了しても期待値検証の失敗を成功扱いしません。
`result/gs/native/stdout.log` がGSログ、`result/rt` が応答アーカイブです。
`result/current.csv` と `result/current.npz` に時刻・A/c・jを保存します。CSV列は
`time_au, Ac_x_au, Ac_y_au, Ac_z_au, j_x_au, j_y_au, j_z_au` です。
出力は原子単位の**電荷電流密度 j**です。SALMONの電子物質電流はJ = −jであり、返された
電流を再度符号反転してはいけません。

`examples/si_intro_data/expected.json` に全21時刻の期待値を保存しています。
基準のエネルギーは約 **−31.7465390575 Hartree**、電流成分の絶対最大値は約
**1.47275×10⁻⁶ a.u.** です。自動照合の許容差はエネルギー絶対誤差1e−7 Hartree、
電流の各成分が `atol=1e-9`、`rtol=1e-5` です。SCF反復数は診断値で、完全一致の条件ではありません。
期待値は固定したmacOS arm64ランタイムで作成しています。他OSの成果物があることだけでは、
この入門例のOS間の数値一致を検証したことにはなりません。

設定は立方Si 8原子、セル長10.26 bohr、32電子・32状態、実空間格子16³、k点2³、PZです。
GS残差の基準は1e−10です。時間発展はdt=0.02 a.u.で20ステップのみで、短いz偏光パルスは
接続を学ぶための設定です。**格子・k点・時間刻みについて収束した光学物性データではありません。**
研究用途では格子・k点・状態数・時間刻み・伝播時間・パルス条件を別途検証してください。

## 同梱データの出典

`examples/si_intro_data/input.json` に固定条件、`provenance.json` に入手元と
ランタイム識別子・各ファイルのハッシュを保存します。`Si_rps.dat` はSALMONの固定コミット
`3322d0561c0cad5180679136b3eddbe8b909487c` の `samples/exercise_04_bulkSi_gs` から取得しています。
SHA256は `8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c` です。
元のApache-2.0ライセンスを `LICENSE.txt` として保持します。旧KY書式には生成時のXCメタデータが
ないため、固定されたPZ Si例という出典を根拠とし、生成条件を自動認定したとは扱いません。
GS→RTは同じ物理条件・擬ポテンシャルを引き継ぎます。この例は任意の外部GSの整合性判定器ではありません。
