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
python -m salmonpy.runtime build NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

これはGNUコンパイラ用の例です。macOSではGNU Cコンパイラの実際のパスを指定します。通常 `/usr/bin/gcc` はApple Clangです。CMake 3.22以上・GNU Fortran/C・LP64 OpenBLASは事前に必要です。システムへの自動インストールはしません。Windows本体の直接ビルドは未検証で、Linux／WSLや別途検証した環境を使います。`--jobs` はビルド並列数で、計算スレッド数やメモリ制限ではありません。

Libxcが必要なら互換性のあるFortran対応Libxcを準備し、`--cmake-arg=-DUSE_LIBXC=ON` と `--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX` を追加します。

固定した本体ソース一式・ファイルのハッシュ一覧・監査ツールはwheelにも含めるため、SALMONpyのGitチェックアウトは不要です。`build.log` と `build.json` にコマンド・成否・実行ファイルのハッシュを保存します。失敗してもログとソースを保持し、新しい場所で再試行します。成功はビルドとオブジェクト監査の完了を表し、科学的精度の検証ではありません。対象ワークフローの回帰検証後、記録したハッシュを使います。

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

## 実験用：高速HSE版のRT状態分離

HSEは利用者のPR予定版`c9e8ce431f43eefd4a88998134adb4796f8ec837`を別プロファイルとして固定します。`python tools/prepare_hse_rt_state.py NEW_DIRECTORY`で、同梱ソースを新しい外部ディレクトリへコピーし、検証した完成形のファイルを配置します。`--reference-only`は比較用の元コードを作ります。両方とも必須のソース・object検査をCMakeの構成・増分ビルドへ保持し、PRチェックアウトは変更しません。既定のSALMONpyランタイムは置き換えません。

この段階では、ACE因子、受理済み予測履歴、PT-CN履歴とdt、Taylor試行、交換作用キャッシュを`s_exx_rt_state`へまとめています。Taylor/PT-CNの寿命管理は状態を引数に取り、異常なdtや段階の指定を変更前に拒否します。所有する確保済み配列のバイト数も診断できます（履歴配列のdescriptor等を含み、RSSではありません）。

既存のHSEアダプターはまだ標準状態1つを使います。幾何・kernel・Wannier、公開legacyフラグ、入力設定などの共有状態が残るため、同一プロセス内の複数材料・同時実行に対応したという意味ではありません。交換演算の高速経路は維持します。合成試験の数値比較は、実材料の精度や性能改善の認定とは区別します。FDTD結合は共通dtのままです。

## 実験用：HSE状態を交換作用へ明示的に渡す経路

前節の準備に`--explicit-routing`を加えると、変更前後のファイルのhashを照合したうえで経路分離済みの完成形を配置します。`--reference-only`とは併用できません。既定の準備・既存の呼び出し方は維持します。

交換・履歴処理の24手続きは、指定されたRT状態を内部呼び出しにも渡します。各呼び出しのローカル参照を使い、プロセス全体のactive状態を切り替えたり全配列をコピーしたりしません。公開legacyフラグと交換エネルギーはRT状態が所有し、旧名は標準状態だけを参照する互換用aliasです。aliasの関連先を利用者が変更することは対応範囲外です。

systemへ結合しない従来の状態単独入口では、対象は共通の交換設定を使う直列CPU・scalar dense ACE作用と寿命管理です。非標準状態によるrefresh、直接/full交換、packed・分散・spin経路は未対応で、共有資源へ到達する前に拒否します。`exx_refresh`と`exx_add_action`で`rt_status`を指定した場合、この能力不一致は非ゼロを返して状態/出力を保ちます。省略時のnativeエラーはプロセス停止になり得ます。一般的な失敗原子性や公開C ABIへの適合ではありません。

本体へ接続する開発用入口は `exx_bind_system(system, owner, status)` です。
Hamiltonian、XC/全エネルギー、Taylorの試行状態、PT-CNのfreeze/履歴、再開情報が
systemに結び付いた所有者を使います。高速PR版のkernel、ACE、試行配列、履歴を
所有者ごとに保持し、プロセス全体のactiveポインタ切替や汎用kernelへの置換はしません。
`SALMONPY_HSE_CONTEXT=1` を指定した明示的な検証実行では、本体のmain_tddftも
この入口を使います。既存バイナリやPython TDDFT APIは変更しません。

対象は直列CPU・周期・非スピン分極HSE・全実空間格子・全支持交換で、プロセスの
入力設定は共通かつ不変とします。結合済み所有者のrefresh/actionは空間・軌道・k点MPI、
Wannier、対称性、spin経路を拒否します。同じ格子点数/k点数で格子間隔の違う状態は
所有できますが、呼び出しの途中でglobalなXC・mesh・伝播法・factor設定を切り替える
使い方は未対応です。設定不変は呼出側の前提であり、設定全体を独立化したわけではありません。
FDTDの共通dtは維持します。

ownerは全呼び出しが終わるまで生存させ、ownerや結合済みsystemのコピー・代入を
しないでください。同時呼び出しは不可です。終了時は `exx_unbind_system(system)`、
続いて `rt_state_clear(owner)` を呼び、所有FFT計画を解放します。使用中のownerへの
二重bindは拒否し、再利用前にclearを要求します。内部Fortranの開発用入口であり、
安定したC/Pythonの複数TDDFT APIではありません。native内部の数値・確保エラーは
プロセス停止になり得ます。FFT計画の確保量が不明なため、資源付きownerに対する
`rt_state_bytes` はstatus=2・bytes=-1を返します。使用量0と解釈しません。

```sh
python tools/test_hse_explicit_state.py HSE_BUILD_RECORD NEW_OUTPUT --native-resources --threads 2
python tools/test_hse_rt_state_runtime.py HSE_BUILD_RECORD REFERENCE_RECORD FIXTURES NEW_OUTPUT --native-context --threads 2
```

前者はセルの違う状態のACE/直接交換を交互実行し、別プロセス参照と比較して、破棄・
再構築による干渉も検査します。後者は本体TDDFTを32ケースの合成kmesh/dt/伝播法/factor設定で
動かし、電流・エネルギー・予測/solver診断を比較します。HSE資源の独立性と本体の数値回帰を
検証するもので、TDDFT全体の複数状態の交互実行や実材料の収束精度を保証しません。
単独コンパイル用の `test_hse_rt_state.py` はRT状態のみの親プロファイル用であり、
本体接続済みプロファイルには上の全objectを使う検査を用います。

## 実験用：独立したHSE資源の部品

`salmonpy/native_adapters/hse_resource_context.f90` はSALMONpy固有のソースとして、
このリポジトリとwheelに同梱します。前節の固定した高速PR版のkernelモジュールへ
リンクします。本体ソース一式も同梱し、コンパイルは明示的に行います。
段階別の変更ファイルの完成形は runtime/overrides に保持します。

`hse_context` はFFT kernelと受理済みRT/ACEデータを所有します。`hse_open` で
直交した直方体格子、均一な全k点メッシュ、格子間隔、遮蔽・混合率を固定し、
`hse_refresh` に規格直交した複素占有軌道を渡し、`hse_apply` で混合率を掛けた
ACE作用を得ます。配列順序は(grid, orbital, k)です。占有はスピンごとに1で、
スピンの倍率は含めません。電荷電流を返す入口ではなく交換計算の部品です。
格子体積を含む軌道の重なりが単位行列に対して1e-8以内であることを検査します。
CPUのstrided FFTに限定し、汎用kernel作用を使います。最適化済み本体の全経路、
性能、本体TDDFTとの等価性を保証するものではありません。

失敗した構築も含めて `hse_close` を明示的に呼びます。繰り返しのcloseは可能です。
コンテキストのコピー・代入、メモリの直列化、同時呼び出しは不可です。
コピー可能なRT履歴と異なり、FFT計画のポインタを所有します。
statusは0=成功、1=入力・寿命の誤り、2=アダプターで捕捉した確保失敗、
3=kernel/ACE/数値エラーです。戻り値で報告するrefresh失敗は受理済みACEと軌道を、
apply失敗は出力を保持します。FFT作業領域は変わり得ます。既存kernel内部の
確保・FFT失敗によるプロセス停止までは回復しないため、一般的な失敗原子性を
備えた公開ABIではありません。

```sh
python tools/test_hse_explicit_state.py HSE_BUILD_RECORD NEW_OUTPUT --resource-context --threads 2
```

ソースと全リンク対象のアプリケーションobjectを検査してから、異なる格子・k点数・
遮蔽・混合率の状態で交互refresh、独立したkernel/ACE参照との一致、不正入力時の
状態保持、両方の破棄順序を試験します。実材料の収束検証ではなく合成回帰試験です。
この部品はまだ `exx_refresh` やPython TDDFT状態へ接続していません。
Wannier、空間分割・MPI、対称性、spin、予測、全エネルギー・電流、最適化経路の
接続は次の段階です。FDTDの共通dtは維持します。

## ソースの所有と更新

標準ソルバーと高速PR版HSEの本体ソース一式を、元のライセンス・由来情報とともに
リポジトリとwheelに保持します。準備時は全ファイルのハッシュを検証してコピーし、
段階別の完成形について変更前後のハッシュを確認して配置します。Gitやパッチ適用は不要です。
旧形式の UPSTREAM OUTPUT 引数も受け付けますが、使用するのは同梱した固定ソースであり、
指定された外部チェックアウトは変更しません。HSEの保守対象は高速版のみです。標準版の
USE_HSE=ON は明示的に拒否し、高速版のreference-onlyは比較専用とします。

source_tree_sha256 は同梱ツリーの一覧、replacements_sha256 は段階別の完成形一覧を
識別します。旧patch_sha256は既存バイナリの検証記録との対応・由来情報としてのみ残し、
新しいソースの検証済み証明には使いません。NLCCの読込はsource-tree識別子が存在する場合、
未知の値を追加で拒否します。今後の数値変更には改めて検証が必要です。
ソース・一覧・ハッシュは同じコミットで更新し、共有カーネルの一致と影響するプロファイルの
ビルド・回帰を確認します。高速HSEの本体接続は上記の所有者経路を使い、独立した汎用資源部品は診断用として別に保持します。
