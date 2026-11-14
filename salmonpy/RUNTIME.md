# 計算ランタイムの準備

LD／SBEのwheelには共有ライブラリを同梱しています。DFT／TDDFT／FDTD／バンド計算には、リンク済み実行ファイルと依存ライブラリを使います。オブジェクトファイルだけでは利用者側でリンクが必要になります。import時のダウンロード・コンパイルはありません。

## まず既存品・配布品を使う

利用者のメモリ上限・CPU数・目的を確認し、OS・CPU種別と必要機能を調べます。検証済み実行ファイルがある場合は、記録済みSHA256で `SalmonExecutable` を作れます。管理された配布物は次で再検証します。

```python
from salmonpy.runtime import load_runtime
engine = load_runtime('my-runtime', required=['tddft', 'libxc'])
# DFT／TDDFT／バンドの入口へ executable=engine として渡す。
# FDTDは engine.run(calculation, 'new-output', threads=2)。
```

公開カタログを使う場合の手順です。`CATALOG_URL` は公開元から得たURLに置き換えます。

```sh
python -m salmonpy.runtime status --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime install my-runtime --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime verify my-runtime --require tddft libxc
```

**[独立SALMONpyの試用版Release](https://github.com/otobetomohito/SALMONpy/releases/tag/runtime-2026.11.14)を配布しています。** 同梱カタログには3環境を登録しました。非公開リポジトリなので、認証済みブラウザで `catalog.json` と対応するOS・CPU名の内側ZIPをダウンロードし、下記の `--archive` で導入します。URLによる匿名ダウンロードは利用できません。Libxc版はこのReleaseに含みません。候補が複数なら `--artifact ID` で指定します。対応品がなければエラーとbuild入口を返し、勝手にコンパイルしたりLibxcを外したりしません。保存先は新規ディレクトリです。HTTPS、ZIPと全ファイルのSHA256、展開パスを検査します。ハッシュは署名ではないので、カタログは信頼できる配布元から取得します。

## 現在のブランチから取得する

Actionsの **Full runtime bundles** はLinux x86_64・macOS arm64・Intel Mac用のビルドを設定しています。別の場所へ移した実行ファイルでTDDFT接続テストが成功した場合のみ配布物を保存します。設定済みと実行成功は区別してください。これらは単一MPIランク・PZ／TB-mBJ用で、Libxcを含みません。手元で検証したmacOS arm64のLibxc版は別プロファイルです。Windows本体ランタイムは未提供です（LD／SBEのWindows wheelとは別）。

現在のブランチの成功したActions実行から `salmonpy-runtime-...` を取得します。**Actionsの外側ZIPを展開**し、中のランタイムZIPと `catalog.json` を使います。

```sh
python -m salmonpy.runtime status --catalog downloaded/catalog.json --require tddft
python -m salmonpy.runtime install my-runtime --catalog downloaded/catalog.json --archive downloaded/RUNTIME_ZIP --require tddft
```

`RUNTIME_ZIP` は実際のファイル名です。ブランチの配布物には公開URLがないため `--archive` が必要です。利用者側のコンパイラやリンカーは不要です。ビルド環境より古いOSは保守的に対象外とします。計算入力・擬ポテンシャルは別途必要です。

## 対応品がない場合だけ明示的にビルドする

必要な開発環境と利用者の計算資源を確認してから実行します。

```sh
git clone https://github.com/SALMON-TDDFT/SALMON2.git UPSTREAM
python -m salmonpy.runtime build UPSTREAM NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

これはGNUコンパイラ用の例です。macOSではGNU Cコンパイラの実際のパスを指定します。通常 `/usr/bin/gcc` はApple Clangです。Git・CMake 3.22以上・GNU Fortran/C・LP64 OpenBLASは事前に必要です。システムへの自動インストールはしません。Windows本体の直接ビルドは未検証で、Linux／WSLや別途検証した環境を使います。`--jobs` はビルド並列数で、計算スレッド数やメモリ制限ではありません。

Libxcが必要なら互換性のあるFortran対応Libxcを準備し、`--cmake-arg=-DUSE_LIBXC=ON` と `--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX` を追加します。

固定版パッチ・共有カーネル参照・監査ツールはwheelにも含めるため、SALMONpyのGitチェックアウトは不要です。`build.log` と `build.json` にコマンド・成否・実行ファイルのハッシュを保存します。失敗してもログとソースを保持し、新しい場所で再試行します。成功はビルドとオブジェクト監査の完了を表し、科学的精度の検証ではありません。対象ワークフローの回帰検証後、記録したハッシュを使います。

```python
import json
from pathlib import Path
from salmonpy.fdtd import SalmonExecutable
record = json.loads(Path('NEW_BUILD/build.json').read_text())
if record['status'] != 'built-and-audited':
    raise RuntimeError('runtime build did not complete')
engine = SalmonExecutable(record['executable'], record['sha256'])
```

MPI版は計算機のMPI環境に合わせて明示的に構築します。DFT/バンドとプロセス分離TDDFTは単一ランクのままです。FDTDはFDTDMPIから既存の空間MPI分割を利用できます。FDTD_MATERIALS.mdを参照してください。今回の準備機能は、計算順序・独立状態・A→Jの境界を変更しません。

Experimental host-MPI TDDFT and native E/H coupling now have a separate entry: [COUPLED.md](COUPLED.md). One point per solver group and one native lifetime per process; do not infer unrestricted reentrancy or restart support.

## NLCC力の検証記録

`load_runtime(directory, required=['nlcc-forces'])`は、ハッシュで確認したバンドル内の
検証記録が実行ファイル・対応パッチ・CPU非スピン分極NCの範囲に一致した場合だけ、
CalculatorのNLCC力を有効にします。通常の実行ファイル指定や機能名だけでは有効になりません。
自分で検証したビルドには`load_validated_build(build_json, validation_json)`を使います。
詳しくは[Calculator](ASE_CALCULATOR.md)を参照してください。

開発者向けの明示的な手順は、`tools/run_nlcc_native_tests.py`でXC・コア力検査を実行し、
`tools/validate_nlcc_runtime.py BUILD_JSON XC_TEST_JSON CORE_TEST_JSON PBE_UPF NEW_CHECKS_DIRECTORY`
で固定Si8の実行検査と記録作成を行います。同梱時は`tools/package_salmon_runtime.py`へ
`--nlcc-validation CHECKS/validation.json --nlcc-potential PBE_UPF`を追加します。
同梱後のバイナリでも実行検査を行い、新しいSHA256へ結び付けます。import時にビルドや計算はしません。
これは実装回帰の記録であり、実材料の格子・k点・SCF・力の収束確認は別途必要です。

## 実験用：RT状態が時間設定を所有する版

既定ランタイムを変更せず、外部に準備した対応版ソースへ`_runtime_source/tools/prepare_rt_step_config.py BASE NEW_DIRECTORY`を明示適用できます。開発チェックアウトでは`python tools/prepare_rt_step_config.py BASE NEW_DIRECTORY`です。新しい外部ディレクトリだけを作り、元ソース・既存実行ファイルは変更しません。失敗した準備は状態記録とともに残り、実行承認されたソースではありません。

この追加パッチはRT状態にdt・Taylor次数・伝播方式・場更新設定を持たせる最初の分離です。PZの合成入力で比較する実験版であり、同一プロセス内の複数材料対応を完成させたものではありません。孤立系のcalcVbox、イオン運動、単一スケールMaxwellなどには共有設定参照が残ります。FDTDとの結合は引き続き共通dtが必須です。

HSEの基準には利用者のPR予定高速版を用います。今回のPZ検証はHSEの速度・精度評価ではなく、このパッチをHSEのPRブランチへ自動適用しません。HSEのACE因子・予測履歴・作業領域の状態分離と検証が別途必要です。ビルドは明示的な依存ライブラリと対応コンパイラを使い、禁止固有値ドライバーのソース・全アプリケーションobject検査とCMake検査を維持します。
