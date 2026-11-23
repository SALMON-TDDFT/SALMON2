# DFT・TDDFT・FDTDの計算プログラムを用意する

[English](RUNTIME.en.md)

LD・SBEの最初のテストには追加の計算プログラムは不要です。対応wheelを導入し、
`python -m salmonpy.quickstart` を実行してください。
DFT・TDDFT・FDTD・バンド計算には、SALMONpyの使用する入口に対応した実行ファイルまたは共有ライブラリを追加します。

## 使い方に合うものを選ぶ

| 使い方 | 必要なもの | 案内 |
|---|---|---|
| LD・SBE・波形・構造・解析 | SALMONpy wheelと必要な任意依存 | [最初の手順](README.md) |
| 別プロセスでDFT・TDDFT・FDTD・バンド計算 | 対応する実行ファイル | このページ |
| 同じMPIジョブ内で電子計算とFDTDを結合 | 呼び出し元とMPIが一致する共有ライブラリ | [結合手順](COUPLED.md) |
| HSEの基底状態・時間発展 | HSEに対応する計算ライブラリ | [保持GS](NATIVE_GS.md)・[HSE](SHARED_HSE.md) |

対応するビルド済み配布物を優先し、ない場合だけ明示的に構築します。
import・入力検査だけで、取得・構築・MPI初期化・計算実行は始まりません。
wheelが使えることと、追加の計算プログラムが同じ環境で使えることは別に確認します。

## 既存の配布物を使う

OS・CPU・必要な機能・Pythonパッケージとの対応が確認された配布物を選びます。
一般のSALMON実行ファイルに、SALMONpy用の逐次呼出しや状態管理の入口があるとは限りません。

導入済みの配布物を検証する例です。

```python
from salmonpy.runtime import load_runtime
engine = load_runtime('my-runtime', required=['tddft'])
```

`engine` をDFT・TDDFT・バンドの `executable=` に渡します。
Libxc-PZ/PBEが必要な場合は `required=['tddft', 'libxc']` とし、Libxcのない配布物へ置き換えません。
ハッシュ一致・機能宣言だけで、すべての入力や実行時プロトコルの互換性を保証するものではありません。

信頼できる配布元から取得したローカルカタログとZIPを使う例です。

```sh
python -m salmonpy.runtime status --catalog downloaded/catalog.json --require tddft
python -m salmonpy.runtime install my-runtime --catalog downloaded/catalog.json --archive downloaded/RUNTIME_ZIP --require tddft
python -m salmonpy.runtime verify my-runtime --require tddft
```

`RUNTIME_ZIP` はカタログに載る実際の計算プログラムZIPです。保存先は新しいフォルダにします。
カタログと全ファイルのハッシュ、ZIPの展開先を検査します。ハッシュは配布元の署名ではないため、カタログの入手元を確認してください。
公開されたHTTPSカタログも使えますが、privateリポジトリの資産はアクセス権のある方法で事前に取得します。

Actions成果物は外側のZIPを展開し、中の `catalog.json` と計算プログラムZIPを上記へ渡します。
[旧Si入門](SI_INTRO.md)の専用 `prepare` は外側ZIPを受け取る別の手順で、混同しません。
Si例は2026-11-23に新しいApple Silicon用ZIPでローカル再検証済みですが、公開Releaseは未作成です。旧2026-11-14配布物は現行の電子数付き応答と合わず、現在の入口では受け付けません。

## 配布物を作る開発者向け

手動ワークフロー **Runtime bundles and TDDFT integration** の `scope=bundle` は、
Linux x86_64・Apple Silicon Mac・Intel Macで構築・移設後の検証を行い、成功した配布物を保存します。
`scope=worker` はLinux一環境の短いワーカー検査で、配布物は作りません。
設定があることと、特定の版の実行成功を区別してください。通常のwheel CIが検査するのはLD/SBEの配布物です。
Windows用の本体配布物は現行のワークフローでは提供しません。

## 対応品がない場合に構築する

ソース一式・ファイル目録・監査ツールはwheelにも含まれます。別のSALMONリポジトリは不要です。
コンパイラ、CMake、BLAS/LAPACK、必要な場合のMPI・Libxcは構築環境に準備します。
以下はGNU Fortran/C・LP64 OpenBLASを使う例で、コンパイラ名とライブラリの場所は環境に合わせます。

```sh
python -m salmonpy.runtime build NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

`--jobs` は構築時の並列数です。計算時の並列数・メモリ上限ではありません。
macOSの `/usr/bin/gcc` は通常Apple Clangなので、OpenMPを使う場合は実際のGNU Cコンパイラを指定します。
Windows用の本体ビルドは、この手順で検証済みとは説明しません。

- Libxc：対応するFortranインターフェースを用意し、`--cmake-arg=-DUSE_LIBXC=ON` と `--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX` を追加します。
- MPI：施設のMPIと一致する `mpifort`・`mpicc` を指定し、`--cmake-arg=-DUSE_MPI=ON` を追加します。MPI版の構築だけで、単一ランク用のPython入口が並列化されるわけではありません。
- HSE：`--cmake-arg=-DUSE_HSE=ON` を明示します。既定はOFFです。FFTWが必要で、用意されていなければ構築中に固定版を取得する経路があります。オフラインでは事前に用意し、`--cmake-arg=-DFFTW_INSTALLDIR=FFTW_PREFIX` を指定します。[詳細](../runtime/README.md)。

成功時は `build.json` に `built-and-audited`、実行ファイルの場所とSHA256を記録します。
失敗時も `build.log` と準備したソースを残します。既存の検証済みバイナリを上書きしません。
構築と監査の成功は、物性の収束や全機能の実行検証を意味しません。

```python
import json
from pathlib import Path
from salmonpy.fdtd import SalmonExecutable
record = json.loads(Path('NEW_BUILD/build.json').read_text())
if record['status'] != 'built-and-audited':
    raise RuntimeError('runtime build did not complete')
engine = SalmonExecutable(record['executable'], record['sha256'])
```

維持するソースは `runtime/sources/solver` の一系統です。変更はそこへ統合済みで、
過去の状態分離パッチを追加適用する手順は不要です。ソース・全アプリケーションオブジェクトの
禁止固有値ソルバー監査は構築時に維持します。変更していない検証済みプログラムの再利用時は、記録SHA256を照合します。

## インターネットに接続できない計算ノード

1. 計算先のOS・CPU・Python・MPI・施設モジュールを確認します。
2. 接続可能で許可された環境でwheelと依存、計算プログラム、入力・構造・擬ポテンシャル・初期状態を準備します。出典・ライセンス・ハッシュを記録し、認証情報を計算入力に含めません。
3. 相対参照を保って転送し、ローカルwheelは `pip --no-index --find-links`、計算プログラムはローカルのカタログとZIPで導入します。
4. 施設の計算資源へ投入して小規模例を確認します。管理ノードでの準備や検査が成功しても、そのまま計算へ移行しません。

不足するファイルは計算前に用意します。暗黙の取得・ソースビルドは行いません。
擬ポテンシャルはローカルに準備します。許可付きの取得は[専用の入口](PSEUDOPOTENTIALS.md)を使います。GS/2は参照・SHA256・出典を保持し、計算中の自動取得はしません。
[準備だけを行う入口](PREPARATION.md)・[構造ファイル](STRUCTURE.md)も参照してください。

## NLCCを含む力を使う

NLCC付き擬ポテンシャルの力には、バイナリと対象条件に結び付いた検証記録が必要です。
`load_runtime(..., required=['nlcc-forces'])` または検証したローカルビルドの
`load_validated_build(build_json, validation_json)` を使います。
検査・配布への収録方法は[ASE Calculator](ASE_CALCULATOR.md)を参照してください。
力の実装検証と、格子・k点・SCFの物理的収束は別に確認します。

## 同じジョブ内で検証結果を再利用する

多数の物質点を作る場合、最初の完全検証を明示的に再利用できます。

```python
import secrets
from salmonpy.runtime import prepare_runtime_job, load_runtime
validation = prepare_runtime_job(
    'my-runtime', 'new-job-validation.json',
    job_id='CURRENT_JOB_ID', job_token=secrets.token_urlsafe(32),
    required=['tddft', 'single-rank'],
)
engine = load_runtime('my-runtime', required=['tddft', 'single-rank'], validation=validation)
print(engine.runtime_validation_report)
```

計算・MPI初期化・取得・構築は始まりません。ジョブの代表が一度だけ準備し、必要なワーカーへ
そのジョブ内で配ります。秘密はログ・入力・チェックポイントへ書かず、次のジョブでは新しく作ります。
省略時は毎回完全検証です。対応する更新情報の高速経路はPOSIX用で、使えなければ警告して完全検証へ戻ります。
再利用時も実行ファイルのSHA256は照合します。付随ファイルの再ハッシュと木の走査を省けますが、
ファイル数に比例する更新情報の検査は残ります。信頼できない更新情報を返す共有ファイルシステムでの安全性を保証する仕組みではありません。

別プロセス型TDDFTではproviderとbackendの両方へ `runtime_validation=validation` を渡します。
対応していない自作コードへ、この機能を強制しません。共有ライブラリの読込みとは別の機能です。
`mpi_runtime_receipt --help` の例は、明示的な `--run` まで計算もMPI初期化も行いません。
[メモリと計測](MEMORY.md)・[MPI配置](DISTRIBUTED.md)を参照してください。

## 状態の管理・保存について

| 目的 | 現在のガイド |
|---|---|
| 同じ物質の独立したTDDFT状態 | [SHARED_TDDFT.md](SHARED_TDDFT.md) |
| 異なる物質を同じMPI群で順に使う | [TDDFT_REGISTRY.md](TDDFT_REGISTRY.md) |
| HSE・対称性 | [SHARED_HSE.md](SHARED_HSE.md)・[FAST_HSE_SYMMETRY.md](FAST_HSE_SYMMETRY.md) |
| 基底状態・SCF途中状態 | [NATIVE_GS.md](NATIVE_GS.md) |
| 保存・再開と復元失敗時の扱い | [CHECKPOINTING.md](CHECKPOINTING.md) |
| 自作プログラムからの呼出し | [API.md](API.md)・[C_ABI.md](C_ABI.md) |

同一プロセスでの状態切替は、任意の並行スレッド実行とは別です。
対応するXC、MPI分割、保存形式、復帰できる異常の範囲を各ガイドで確認してください。
