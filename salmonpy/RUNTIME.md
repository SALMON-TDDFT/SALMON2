# 計算ランタイムの準備

LD／SBEのwheelには共有ライブラリを同梱しています。DFT／TDDFT／FDTD／バンド計算には、リンク済み実行ファイルと依存ライブラリを使います。オブジェクトファイルだけでは利用者側でリンクが必要になります。import時のダウンロード・コンパイルはありません。

## インターネットに接続できない計算ノード

計算ノードから外部へ接続できることを前提にしません。取得と計算を分け、
接続可能な手元端末や、施設が許可する準備環境で必要なものを揃えてから転送します。
ログインノードも通信可能とは限らないため、施設の運用に従ってください。

1. **計算先の環境を先に確認**：OS・CPU・Python版、必要機能、MPI実装と施設の
   モジュール環境を記録します。手元Mac用のwheelやランタイムをLinuxノードへ転送しても使えません。
2. **必要な依存をまとめて準備**：対応するSALMONpyとNumPy、利用する場合のASE等の
   依存wheel、必要なランタイムZIPと信頼するカタログ、構造・スーパーセルファイル、
   擬ポテンシャル、入力、初期状態を取得します。参照先の外部ファイルも含め、版・出典・
   ライセンス・SHA256を記録します。必要なMPI/BLAS/Libxc等が配布物に含まれるか、
   施設環境で提供されるかも確認します。認証情報は計算入力に埋め込みません。
3. **転送後にローカルで検証**：ファイルのハッシュを照合し、相対参照を保って配置します。
   Pythonの導入はローカルwheel一式を使い、パッケージ索引へのアクセスを無効にします。
   ランタイムはローカルの `--catalog` と `--archive` を両方指定して導入し、
   `verify` と小規模例を実際の計算ノード上で実行します。ログインノードで通るだけでは十分ではありません。
4. **不足はジョブ投入前に解消**：データベースURLを実行時に参照せず、取得済みの構造ファイルを
   指定します。ファイルや依存が不足した場合は具体的に報告し、暗黙のダウンロード・
   ソースビルドへ切り替えません。対応ランタイムがない場合のビルドは、依存を揃えた
   施設対応環境で明示的に行います。

この案内は転送用一式の自動収集・依存解決機能ではありません。自作バックエンドも、
実行中の通信や取得を必要としないよう準備してください。SALMONpyのimportは通信しませんが、
利用者のコールバックによる通信を禁止する仕組みではありません。

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
変更済みの完成形は runtime/sources/solver に統合しています。

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


## 実験用：状態ごとの固定dtを持つTDDFT全体の独立状態

`salmonpy/native_adapters/tddft_context.f90` に内部Fortran開発用の `rt_open`、
`rt_step`、`rt_close` を用意します。本体の初期化・時間発展ルーチンの引数である
波動関数、密度、Hartree/XC/局所ポテンシャル、伝播履歴、擬ポテンシャル作業配列、
XCハンドル、HSE所有者、出力ハンドルを各状態が保持します。同じrestartから各状態を
別々に初期化し、初期化済み状態をコピーしません。stepは端点A/cを受け取り、電荷電流密度j
（native物質電流の逆符号）と全エネルギーを原子単位で返します。エネルギーは毎ステップ計算します。

事前にプロセス共通のnative入力と直列ランタイムを初期化します。具体例は
`tests/native/test_tddft_contexts.f90`です。安定したPython/C ABIではありません。
対象は周期・非スピン分極の高速HSE、直列CPU/FFTE、固定イオン、transverse、初期A=0、
dt以外は共通かつ不変のglobal設定です。通常/learned Taylor、直接Taylor、PT-CNを検証対象とします。
MPI、OpenACC、FFTW Poissonビルド、第2波形、restart、空間分布・射影・分解電流などの
追加出力は拒否します。状態名は英数字・underscore・hyphenに限定し、既存出力を上書きしません。

格子/k点数・バンド/電子数・XC・遮蔽・温度・伝播/予測法・factor設定を保存し、
stepや別状態の初期化前に変更を拒否します。これは保護であり、salmon_global全体の分離ではありません。
入力・初期化設定には共有が残り、global変数の任意の変更は未対応です。
異なるXC・格子・k点を同一プロセス内で併存させる作業は残ります。
内部入口 `rt_open(..., delta_t=...)` で状態ごとの有限・正のdtを固定します。省略時は入力dtを取得します。
初期化、Taylor係数、PT-CN、電場評価、出力時刻は状態のdtを参照します。
途中のdt変更、適応刻み、FDTDサブサイクルは未対応です。
状態をコピー・同時実行せず、rt_closeで明示的に解放します。不正Aや検査対象設定の変更は
進める前に拒否しますが、native内部の確保・数値エラーはプロセス停止になり得ます。RSS上限はありません。

```sh
python tools/test_tddft_contexts.py HSE_BUILD_RECORD SYNTHETIC_FIXTURES NEW_OUTPUT --reference FAST_REFERENCE_RECORD --threads 2
```

本体のTDDFTステップで異なるA波形の2状態を交互実行し、別プロセス参照、逆順実行、
片方を破棄・再生成しながら他方を継続する場合を比較します。電流・エネルギーだけでなく、
軌道・密度・ポテンシャル・HSE履歴も比較します。返す電流の符号とエネルギーは、
元の本体が対応するAcos2パルスでも直接照合します。元の本体は任意波形ファイルを拒否するため、
その制限を外す変更はしていません。バイナリtraceは合成試験の診断用で、restartではありません。
ローカルでの実装回帰であり、物理収束、異なる格子/XC/k点、MPI、新しいPython APIの保証ではありません。

検査ツールに `--different-dt` を付けるとdt=0.001/0.0015を、そのdtを入力した別プロセス参照と比較します。
交互実行時には旧global dtを意図的に変更し、逆順・破棄/再生成も検査します。
ローカル検証の記録は `docs/tddft-state-settings-validation-2026-11-15.md` に残します。


### 現在の高速PR版での共有カーネルと場のランタイム

新しいexplicit-routing高速PR版は、通常版とSBE・FDTD・LDのソースを共有し、
準備時にnative/wheel側のカーネルとのバイト一致を検査します。同じソースをHSE有効／
無効でビルドでき、無効側はHSE用FFTWをリンクしません。

NativeEHFieldとHSE物質群にはNEW_BUILD/build.jsonを使います。高速版で場を使う
記録にはfield_abiのprotocol1／salmonpy.eh-centered-current/1宣言が必要です。
MPIバッチFDTDのSalmonExecutable.from_mpi_buildには別の
NEW_BUILD/executable-build.jsonを渡します。それぞれ実体のSHA256を記録します。
旧embeddedの場記録も引き続き使え、import時のビルド・通信は増やしません。

場・物質はホスト管理の別MPI群で実行します。標準GS／TDDFTのソース統合と、
同一プロセスで場・電子状態を同時所有する対応は未完です。既定の準備入口と
配布済みバイナリは切り替えていません。

## 開発用：XC作業配列とステンシル表の所有

explicit高速HSEプロファイルでは、XC格子配列とPZ作業配列を各`s_dft_system%xc_work`が所有します。
`calc_xc(..., workspace=...)`は指定したPZ作業配列を再利用し、省略時は呼び出し内だけで確保します。
形状・スピン数の変更では対象の配列だけを再確保します。XC設定・Libxcハンドルはコピーしません。
周期境界のmodulo表と最適化ステンシルの可否フラグは各`s_stencil`が所有し、
Hamiltonianとステンシルの処理へ明示的に渡します。状態の破棄で配列も解放されます。

`tools/test_owned_workspaces.py HSE_BUILD_RECORD NEW_OUTPUT --threads 2`で、
元のPZルーチンとの比較、異なる2形状、形状・スピン数の再確保、表の内容、他方の解放後の継続を検査します。
これは作業配列の検査であり、異なる格子のTDDFT全体を検証したものではありません。
入力・初期化設定などに共有が残るため、TDDFT全体の共通格子/k点/XC制約は維持します。
Apple Siliconのローカル検査ではAVX/SVE最適化カーネル、MPI、GPU、同時呼び出しは検証しません。
記録は`docs/tddft-workspaces-validation-2026-11-15.md`です。


## 開発用：状態ごとのOpenMP分割方法

ステンシル・電流・力のOpenMP分割フラグを各`s_stencil`へ保存し、計算中は共有フラグを参照しません。
内部入口`rt_open(..., openmp_mode=...)`は`auto`・`orbital`・`rgrid`を受け付け、
3カーネルに共通の指定を状態作成時に固定します。省略時は従来の個別入力と自動選択の閾値を保ちます。
不正値は初期化前に拒否します。これはループの分割方法であり、スレッド数の変更ではありません。
同時呼び出しは引き続き未対応です。

全体検査に`--different-openmp`を付けると、Aを軌道方向、Bを実空間方向に分割し、
各々の別プロセス参照・逆順・破棄/再生成と比較します。交互実行時は旧入力文字列も変更します。
`--different-dt`との併用が可能です。共通格子/k点/XC制約は維持します。
力の分割方法も接続しましたが、固定イオンの電流回帰を力の数値検証とは扱いません。
記録は`docs/tddft-openmp-validation-2026-11-15.md`です。


## 開発用：明示的な作成時設定

`structures:s_rt_run_config`は`delta_t`・`max_steps`・`stencil_mode`・`current_mode`・
`force_mode`を持つ値オブジェクトです。native入力を初期化した後、`rt_input_config(options)`で
従来入力を取得し、必要な項目を設定して`rt_open(state, name, status, config=options)`へ渡します。
値をコピーするため、作成後に元のoptionsを変更しても状態は変わりません。

dtは有限・正、容量は1以上かつ既定整数の最大値未満、分割指定はauto/orbital/rgridです。
従来のdelta_t/openmp_modeキーワードとの同時指定は拒否し、不正値は初期化前に返します。
容量上限の検査は端点整数のオーバーフロー対策であり、メモリ見積もり・確保成功の保証ではありません。

RT配列と最終ステップ判定は状態の容量を使い、上限超過は更新前に拒否します。
`--different-capacity`ではA=12/B=6ステップを、異なるdt・分割方法、別プロセス参照、
逆順・破棄/再生成、元設定の容量変更、旧global nt変更と組み合わせて検査します。
一方が上限に達しても他方は継続できます。格子/k点/XCなど残る入力は共通かつ不変のままです。
汎用入力オブジェクトや公開Pythonバックエンドの完成とは扱いません。
記録は`docs/tddft-config-validation-2026-11-15.md`です。


## 開発用：明示的な格子境界の初期化

`structures:s_rgrid_config`に周期性（0または3）、指定点数、格子間隔、整合確認用の長さを保持します。
`grid_config_from_input(lengths, config)`で従来入力を取得し、Maxwellでは既存のal_em/dl_em規則を保ちます。
`initialization_sub:init_grid_whole_config(lengths, spacing, config, grid, status)`は
明示的に渡した値だけから境界・halo範囲・添字配列を生成します。旧init_grid_wholeは互換ラッパーとして残します。

座標配列とMPI分割表は従来どおり変更せず、後続の初期化処理が担当します。
完全な格子・ソルバーのコンストラクターではありません。不正な設定はstatus=1、確保失敗はstatus=2で返し、
成功するまで渡された格子を変更しません。整数範囲の検査はメモリ予算の保証ではありません。
旧ラッパーの失敗は従来同様にプロセス停止となります。

`tools/test_grid_config.py HSE_BUILD_RECORD NEW_OUTPUT`で、偶数/奇数の周期・孤立境界、
点数/間隔指定、異なる2形状、旧globalの変更、Maxwell整合性、不正呼び出し後の状態保持を検査します。
下位の孤立格子検査を孤立系TDDFTの追加対応とは扱いません。
セル・k点・擬ポテンシャル・XC/HSE初期化には共有入力が残り、TDDFT全体の共通格子/k点/XC制約は維持します。
記録は`docs/grid-config-validation-2026-11-15.md`です。


## 開発用：明示的なセル形状の決定

`initialization_sub:resolve_cell(config, geometry, status)`は値オブジェクト
`s_cell_config`（長さ・ベクトル・間隔・点数）だけを読み、`s_cell_geometry`
（基本ベクトル・長さ・間隔・直交性）を返します。直交セルでは明示的な長さ、対角ベクトル、
点数×間隔の優先順位を保ち、斜交セルでは元のベクトルとノルム計算を使います。
明示的な格子間隔は長さ÷点数より優先します。不正・非有限・非正・体積ゼロのセルは
結果を変更せず非ゼロstatusで返します。規格化した行列の行列式で特異セルを検査しますが、
ほぼ特異なセルの数値条件までは保証しません。

DFT/RT初期化もこのルーチンを使います。互換層には後続処理向けのglobal alの正規化が残ります。
セル/k点/擬ポテンシャル/XCの共有入力の完全分離や、異なるセルのTDDFT全体への対応ではありません。
`tools/test_cell_config.py HSE_BUILD_RECORD NEW_OUTPUT`で優先順位・斜交セル・独立性・
失敗時の保持を検査します。記録は`docs/cell-config-validation-2026-11-15.md`です。


### k点・擬ポテンシャル・XC入力の所有

高速HSE開発プロファイルでは、内部Fortranの入力を次の値型へ分離しています。

- `s_kpoint_config` → `lattice:init_kvector_config(system, unfold, config)`。メッシュ・シフト・k点ファイル・unfold・対称操作をコピーし、実効メッシュを`system%kmesh`に保持します。HSEもこのメッシュと対称操作を参照します。
- `s_pp_input_config` → `salmon_pp:read_pslfile(system, pp, config)`。元素・ファイル・局所チャネル・マスクなどを`pp%input_config`へコピーします。検出形式とパーサーの作業フラグも`pp`が所有します。
- `s_xc_config` → `salmon_xc:init_xc_config(functional, config, status)`。汎関数名・スピン・係数・HSE遮蔽・rVV10設定を保持します。初期化済みLibxcハンドルはコピーせず、再初期化前に`finalize_xc`を呼びます。`calc_xc(..., pre_scf_active=...)`は評価段階を明示し、省略時は本計算です。

従来の入口は入力を取得する互換層として残します。明示経路は共有入力を一時的に切り替えません。
これは内部の初期化・評価部品の分離です。SOC/+U専用の擬ポテンシャル処理、SCF段階・アルゴリズム設定、MPIとファイル単位には共有が残ります。異なるk点・XC・物質を持つTDDFT全体の同居はまだ許可しません。内部入力は検証済み設定を前提とし、I/O失敗などの全失敗原子性は保証しません。

検査ツールは`tools/test_kpoint_config.py`、`tools/test_pp_input_config.py`、`tools/test_xc_input_config.py`です（引数は`HSE_BUILD_RECORD NEW_OUTPUT`）。検証範囲と結果は`docs/shared-input-validation-2026-11-15.md`を参照してください。


### スカラー擬ポテンシャルの格子配置

`prep_pp_sub:init_ps(..., config=...)`は内部値型`s_pp_grid_config`をコピーし、`ppg`に保持します。
周期性・Poisson/FFT選択・診断設定は作成時に固定し、省略時の従来入口も最初の取得後は保持します。
原子数・原子種対応は`system`、元素数は読み込んだ擬ポテンシャル配列から取得します。
原子数・種対応・格子形状・格子配置設定・radial配列構成が変わる更新は拒否します。
読み込んだradial数値表は変更せず同じものを再利用してください。検査対象は設定と配列構成で、
毎ステップの表内容のハッシュ計算・全コピーは行いません。原子位置の変更には従来の
`dealloc_init_ps` → `init_ps`の順序を使い、同じ設定で新規作成した結果との一致を検査します。
キャッシュ・補間・射影子の通信表も、共有の原子数ではなく所有配列のサイズを使います。
NLCCの格子配置は、`system%kion`と読み込み時に保持した周期性を参照します。

明示的に独立な入口はスカラー・+Uなしです。従来のSOC/+U入口は残しますが、専用モジュールの
共有状態を分離したとは説明しません。MPI通信資源・グローバルなタイマー・異物質TDDFT全体の
制約も残ります。検査は`tools/test_pp_grid_config.py HSE_BUILD_RECORD NEW_OUTPUT`です。


### XC評価フェーズの独立性

高速HSEの内部入口`exx_set_phase(system, pre_scf)`は、そのsystemの評価フェーズを固定します。
`exx_pre_scf(system)`で取得し、HSEの有効判定・熱的SCF補助処理・半局所XC・rVV10は同じフェーズを使います。
RT初期化は、共有フラグにかかわらず本計算フェーズを所有します。別の計算がpre-SCFに入っても、
作成済みRT状態や新しく作るRT状態の評価を切り替えません。
フェーズ変更時はbind済みの独立HSE ownerだけの交換資源・予測履歴を破棄します。未bindの状態は
作成前のフェーズ指定に使えますが、既存の共有交換キャッシュを転用してはいけません。
pre-SCFから本計算へ戻す明示的な利用者は、現在の軌道から`exx_refresh`を実行してから交換作用を使います。

GS/DCのpre-SCF制御は、閾値・必要連続回数・現在の連続回数・開始反復を各systemに保持します。
GS、軌道計算、DC断片／全体系エネルギー、hybridメタデータは同じフェーズを参照します。
局所判定→集団通信→回数更新→診断→切替→ポテンシャル再構築の順序は従来どおりです。
SCF呼び出しごとに回数と開始反復をリセットしますが、完了したpre-SCFを再開させません。

フェーズ未所有のsystemと引数なし`exx_enabled()`には互換の共有フラグ参照が残ります。
DCの通信・物質情報交換、適応交換制御、従来の交換キャッシュは共有です。
GS/DC全体の独立化や任意の並行実行対応とは説明しません。

`tools/test_pre_scf_state.py HSE_BUILD_RECORD NEW_OUTPUT`で制御状態の独立性を、
`tools/test_pre_scf_runtime.py HSE_BUILD_RECORD ORIGINAL_FAST_HSE_RECORD NEW_OUTPUT`で小さいGS計算4条件を検査します。
後者の閾値は意図的に緩く、切替とスカラー出力の回帰検査です。物理収束の検証ではありません。
MPI有効ビルドの記録を渡し、`--dc-ranks 2`または`--dc-ranks 4`を加えると、2断片を各1または2空間ランクで比較します。
`--mpiexec PATH`でランチャーを指定できます。空間分割HSEの要件に従いrestart出力を無効にし、各ランクのOpenMP／BLASを1スレッドに制限します。
ローカルでDC/MPI計8条件が元の高速HSEと一致しました。複数ソルバーの同時独立実行を保証する検査ではありません。
記録は`docs/xc-phase-validation-2026-11-15.md`と`docs/pre-scf-validation-2026-11-15.md`です。

### SCF交換制御の所有先

内部の`exx_set_scf_flags(system, ...)`、`exx_adaptive_is_ready(system)`、`exx_is_frozen(system)`、
`exx_support_has_changed(system)`は、systemに結び付いたHSE資源を参照します。
GSの適応交換準備、一時的な交換固定、支持領域の変更通知もこの所有先へ振り分けます。
`exx_check_localization(system=system)`は同じ所有先を検査します。集団通信や固定／復帰の順序は変えません。
未結合のsystemは従来の所有先を使います。完全に独立したGSソルバーの構築、残る設定変数の分離、同時呼び出し対応とは区別してください。

### GSが所有するHSE資源

GS入口はsystem初期化後に局所的なHSE資源を結び付けます。対象は既存の独立資源判定で認める範囲に限定します。
通常の`theory='dft'`、構造最適化なし、周期・非スピン偏極・非スピノルのHSE06、CPU上の全支持領域交換、
Wannier・対称性縮約なし、負のtemperature設定（整数占有経路）、電子計算の各分割が1ランクの場合です。
DC、分散／有限温度／局在化交換などは従来の所有先を保ちます。物理モデルの対応範囲を広げる変更ではありません。

選択時は`EXX_GS_OWNER explicit`を出力し、結果・checkpointの出力完了後に結合解除と資源解放を行います。
`tools/test_pre_scf_runtime.py ... --require-gs-owner`は変更版でこの経路が選ばれたことも検査します。
GS全状態・設定の独立化は未完で、同時実行用GS APIは提供していません。

GSの同一入力による連続実行検査は`tools/test_gs_lifecycle.py HSE_BUILD_RECORD ORIGINAL_FAST_HSE_RECORD NEW_OUTPUT`です。
直列プローブを明示的にコンパイルし、1プロセスで2回GSを実行して元コードの別プロセス実行と比較します。
pre-SCFの4条件でスカラー結果と、波動関数・密度混合履歴を含む全生成restartファイルがバイト一致しました。
比較後はrestart実体を削除し、ハッシュとログを残します。小さい同一入力の連続利用の検査であり、
異なる入力の連続利用、並行GSコンテキスト、TDDFTへのcheckpoint読み込みは保証しません。
記録は`docs/gs-lifecycle-validation-2026-11-15.md`です。

### GS状態の保持（内部Fortran）

`salmonpy_gs_context`に`s_gs_context`、`gs_run(ctx,status)`、`gs_close(ctx)`を用意しています。
従来のGS作業変数（軌道、密度・ポテンシャル、占有・力、混合履歴、擬ポテンシャル、XCハンドル、対象構成のHSE資源）を保持します。
**実引数のコンテキストは`TARGET`付きで宣言し、アドレスを変えず、実行中の状態をコピーしないでください。**
破棄前に明示的にcloseします。従来の`main_dft`は実行・解放を行うラッパーとして残します。

生存中のコンテキストを再実行すると変更前にstatus=1で拒否し、正常実行は0を返します。closeは複数回呼べます。
nativeエラーによるプロセス停止は残ります。入力・並列・タイマーは共有のため、安定公開ABI、明示入力GS API、並行実行対応ではありません。
複数状態の保持は同一入力のスカラーHSEで検証します。Aを保持してBを実行・解放し、電子配列、占有・力、混合履歴、
一部のHSEキャッシュが変わらないことを検査します。全メモリの網羅比較ではありません。

### GSの明示SCF設定

`gs_input_config(value)`は共有入力のSCF反復上限・収束指標・収束閾値・pre-SCF閾値・必要連続回数を`s_gs_run_config`へ取り込みます。
`gs_run(ctx,status,config=value)`は初期化前に検証・コピーします。省略時も同じ入力を一度取り込む互換経路です。
`max_iterations`は正かつ整数最大値未満、両閾値は有限・非負、`pre_scf_steps`は正が必要です。
不正設定はコンテキストを有効化する前にstatus=1で拒否します。

SCF反復上限と収束判定（適応交換の準備判定閾値を含む）は保持した設定を使います。
渡した値や該当する共有入力を後から変えても、選択済み設定は変わりません。
`convergence_metric`はrho_dne/norm_rho/norm_rho_dng/norm_pot/norm_pot_dngの5方式です。
指標値の演算・規格化と表示も保持した設定に従います。閾値は選んだ指標に対応する内部原子単位で渡します。
軌道の占有数を決める入力、構造・物質などの設定は未分離です。異なる入力全般の独立化や並行実行を保証しません。

明示設定でpre-SCFを有効にする場合も、従来の条件（新規・静的hybrid、restart/MD/構造最適化なし、
DCでは正の温度、checkpoint間隔・終了時刻指定なし、密度混合と密度収束指標）を入口で再検査します。
従来実行入口でも不正な実行設定は非ゼロ終了として報告し、計算を飛ばして正常終了させません。

ローカルファイルをまとめるPython入口は[PREPARATION.md](PREPARATION.md)です。対応環境・依存の自動解決やインストールは行いません。

### GSの状態別混合設定

内部の`s_gs_run_config%mixing`（`s_mixing_config`）には`method`、`mixrate`、`alpha_mb`、
`beta_p`、`nmemory_p`、`auto_mixing`、`update_ratio`、`initial_no_mix_steps`を保持します。
`gs_input_config`で取り込み、必要なら値を変更して`gs_run`へ渡します。初期化前にコピーし、
実行中の係数・前回残差・密度／ポテンシャル履歴は各状態の`s_mixing`が別々に保持します。
自動半減とpre-SCF等での係数リセットは、この状態の初期係数を使います。

`method`はsimple/broyden/pulay/simple_potential/simple_dm、`auto_mixing`はy/n。
係数とupdate_ratioは有限値、Pulay履歴長は1〜21、無混合反復数は0以上・整数最大値未満です。
simple_dmは従来どおりspin-orbit有効時のみ。履歴の格納容量21と演算順は変更していません。
これは内部Fortran設定であり、新しいPython APIや一般の並行GS実行を提供するものではありません。
legacy初期化は入力をその時点で取り込みます。別入口のDFT-MDの自動混合判定は共有入力のままです。

検証: `tools/test_gs_lifecycle.py ... --mixing pulay --auto-mixing --distinct-mixing`。
同じ物質で混合条件だけ異なる二状態を、元の実行ファイルの各条件と比較します。
詳細と未検証範囲は`docs/gs-mixing-config-validation-2026-11-16.md`を参照してください。

### GS収束指標の電子数規格化

`s_gs_run_config%normalization`は`nelec`と`nelec_spin(2)`を値として保持します。
`gs_input_config`で共有入力から取り込み、`gs_run`が状態へコピーします。
非スピン分極ではnelec、スピン分極ではnelec_spinの和が正ならその和、そうでなければnelecを使う従来規則を維持します。
DCでは従来どおり状態内の`dc%elec_num_tot`が優先されます。
非DCのrho_dneで分母が0以下になる設定は初期化前にstatus=1で拒否します。
その他の収束指標の演算・格子規格化は変更しません。

これは収束判定の規格化入力です。**物質の電子数・電荷・占有数を設定する入口ではありません。**
通常は取り込んだ値をそのまま使います。物質入力の全面独立化やスピン分極/DCの実測検証済みとは案内しません。
`tools/test_gs_lifecycle.py ... --different-normalization`は、診断用に第二状態だけ分母を2倍にして、
その残差だけ半分になり、電子状態出力は変わらないことを検査します。
検証範囲は`docs/gs-normalization-config-validation-2026-11-16.md`に記録します。

### 状態が保持する占有数入力

内部`s_dft_system%occupation_input`はnelec、nelec_spin、nstate、spin、spinorbit、temperature、redistribution_delayを初期化時にコピーします。初期充填・通常ne2mu・通常GS再分配条件・ギャップの占有バンド番号はその値を使います。未初期化systemには従来入力の互換経路があります。温度は内部原子単位です。

生存中の値は変更しません。gs_runの新しい占有数上書きAPIではありません。DC専用ソルバー・並列配置やHSEの他の設定等には共有入力が残り、任意の電子数・温度の同時GS実行を保証しません。

検証: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT`、`tools/test_gs_occupation_runtime.py BUILD_RECORD REFERENCE_RECORD NEW_OUTPUT`。前回の未コミットファイルは再起動後に見つからず、会話から再構成しました。復旧後の証拠は`docs/gs-occupation-recovery-validation-2026-11-16.md`です。

### HSEの状態別占有数入力

HSE交換の電子数・温度・状態数（nelec/temperature/nstate）は、同じ`occupation_input`を参照します。初期化前の従来経路は維持し、演算式・更新順は変更しません。

独立性の検証範囲は、直列CPU・周期境界・非スピン分極・整数占有・全支持の高速HSEです。電子数4/2の二状態について、共有入力を書き換えた単独・交互実行、解放・再作成、作用と交換エネルギーを検査します。有限温度GSは従来版との一致を検査しますが、有限温度/Wannierの複数状態独立性は未検証です。nstateのWannier分岐も新たな独立性検証の対象ではありません。

外側のTDDFT設定一致ガードは維持します。HSE以外の共有入力が残るため、電子数や温度が異なる任意のTDDFT状態を同時に作れる新APIとは案内しません。生存中の占有数設定は変更しません。
検査は`tools/test_hse_explicit_state.py BUILD_RECORD NEW_OUTPUT --native-resources`（必要なら`--threads 2`）。検証記録は`docs/hse-occupation-inputs-validation-2026-11-16.md`です。

### TDDFTの電子数診断

`subdip`の電子数チェックも、状態の`occupation_input.nelec/nelec_spin`を使います。2スピンで正のスピン電子数合計があればそれを優先し、その他はnelecを使う従来規則を維持します。上下限の厳密な不等号、出力、停止処理は変更しません。未初期化状態は従来入力を参照します。
これは広い範囲の異常検知（指定電子数の0.1倍から10倍の間）であり、精密な電子数保存や物理的収束を保証する検査ではありません。共有I/O・イオン情報や外側TDDFTガードは残ります。
検査: `tools/test_occupation_inputs.py BUILD_RECORD NEW_OUTPUT --dip`。検証範囲は`docs/dip-occupation-inputs-validation-2026-11-16.md`を参照してください。

### 高速PR版HSEの直接保存／復元（内部Fortran入口）

`native_adapters/tddft_context.f90` の `rt_checkpoint(ctx,path,restore,status)` は高速PR版HSEの数値状態を直接保存します。保存は `restore=.false.`、復元は `.true.`。新しいプロセスで同じ入力・GSから `rt_open` した直後に復元し、`rt_accepted_step(ctx)+1` から残りを進めます。保存時刻までの時間発展は繰り返しません。静的な交換演算子・FFT資源は新しいプロセスで作り直し、生のFFTハンドルは保存しません。

保存対象は波動関数の両バッファ、密度・局所／Hartreeポテンシャルと履歴、RT履歴、場・電流・エネルギー基準・受理済みステップ、ACE、学習因子／PT-CN履歴、交換計算のキャッシュです。復元前に入力・GS・実行ファイル・アダプター・保存ファイルのハッシュを呼び出し側で確認してください。同じバイナリと設定専用のnative形式で、汎用アーカイブではありません。設定不一致は拒否し、途中まで復元して失敗した状態は以後の `rt_step` を拒否します。

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python tools/test_hse_direct_checkpoint.py FAST_HSE_BUILD_RECORD SYNTHETIC_FIXTURES NEW_OUTPUT
```

ツールは高速PR版の監査済みオブジェクトからアダプターを明示コンパイルし、別プロセスで連続・保存・直接再開を比較します。2／3 k点、Taylor（ACE／学習ACE／高速PR版内full-action）と学習PT-CNを6／7ステップで分割し、残り6／5ステップの電流・全エネルギー・内部トレースが完全一致しました。1／2 OpenMPスレッドで検証しています。通常版HSEや比較用ベースラインは実行エンジンとして使いません。

現段階は周期・スカラー・固定占有・固定イオン・全支持領域のserial CPU高速HSE内部Fortran入口です。安定したnative ABI、GPU、対称性・軌道／実空間分散ACE・縮約q点などへの対応を示すものではありません。内部k点MPIとFDTD結合は下記の専用入口を使います。合成軌道の数値回帰であり、実材料の物理収束は別途確認します。既存nativeテキスト出力の過去のファイル位置は再開しません。

保存は交換資源が初期化された受理済みステップで行います。初期化直後のステップ0保存は検証対象外です。

高速PR版HSEの外側MPI／FDTD接続には、明示ビルドした専用ワーカーと `FastHSEResponseBackend` を使います。このワーカーモードでは、外側MPIの各物質点から別のserialプロセスを呼びます。構築・結合・直接再開はCOUPLED.mdを参照してください。HSE内部MPIへの対応ではありません。

内部k点MPIには `tools/build_fast_hse_mpi.py` と `salmonpy.fast_hse_mpi.NativeFastHSEBackend` を使います。物質群は1／2ランク、軌道／実空間分割なし、メインスレッドで1プロセス1回のnative寿命に限定します。MPIの管理はホストが行い、呼出中はcwdを一時変更します。同じ分割のランク別直接保存でHSE数値履歴を復元し、静的FFT資源は再構築します。4ランクのFDTD/SBE/HSE例と失敗時の制約はCOUPLED.mdを参照してください。

高速版では通常の初期化入力の保存領域とRT出力管理も共通ソースから構築します。
アダプター側の同一出力ソースは削除しました。MPI・直列ビルドはソルバーの
構築済み出力管理を使い、ABIが適合する旧準備済みビルドでは同じ共通ソースを明示的に構築します。
HSE固有設定の完全な保存は従来どおり維持します。通常のGS保持・埋め込みTDDFTの
高速版への移行は、まだ完了していません。

高速HSEのモデル登録でも、通常版と同じモデル別パーサー作業領域を使います。
一時ファイル接続・namelist診断・座標入力の識別・共通入力と既定値の保存領域を
読み込みごとに分離しました。高速HSEの既定値・単位変換・検査は維持し、
HSE固有設定の完全な保存には従来のhse_model_inputsを引き続き使います。
パーサーの共通valuesだけを完全なHSEモデルとして扱うことはできません。
従来のバッチ読み込みではMDとの互換性のため名前付き一時ファイルを維持します。
今回の変更はパーサー移行であり、GS保持・通常の埋め込みTDDFTの統合は残ります。

高速版の内部GSコンテキストに、既存の互換なGS波動関数を読み込む専用経路
`gs_run(..., import_only=.true.)`を追加しました。従来の初期化を通し、SCF・
PBE予備計算・GSチェックポイント出力を行いません。各状態がXC・HSE・MPI資源を
保持し、一方の破棄は他方に影響しません。共通gs_seedモジュールを使う
`gs_export_seed`で軌道・占有数・共通入力・格子と分割の情報を独立した値として
取り出せます。この値は完全なHSEモデルや可搬チェックポイントではありません。
PythonのNativeGSBackend.loadも、gs_abi/1を明示した高速版ビルド記録に接続しました。
密度・占有数・軌道と物質IDを保持し、未計算のエネルギー・残差・許容誤差・収束は
Noneとして返します。通常版のSCF経路は維持します。高速版のsession.run、固有値・力、
アーカイブ出力は未接続です。GSから高速PR版HSEのTDDFTモデルを作る入口は
下記の追加機能を宣言したランタイムで利用できます。
読み込みは収束・エネルギー・力を保証しません。MD・DC・構造最適化の読み込みと
実物質での収束検証は、今回の検証範囲に含めていません。

高速版の追加機能 `gs_abi.memory_rt=1` は、ビルド時にコンパイル済みプロトコルを照合して宣言します。保持GSから高速PR版HSEへ、GS再開ファイル・SCF・時間発展の再実行なしで接続できます。NATIVE_GS.md参照。高速版session.run・固有値／力・アーカイブ書き出しは未接続です。通常XCのRTは下記のHSE無効プロファイルで利用できます。タイマー初期化は通常版と共有し、GS／RT逐次接続で二重確保しません。

高速版ソースもbuild_embedded_tddft.py・USE_HSE=OFFで通常XC用ホストライブラリを構築できます。source_profile=fast-pr、ordinary_rt_abi=1、読み込み専用gs_abi/1を宣言します。接続コード・worker・MPI bridgeは標準版と共有し、差分格子・XC作業領域は高速版の状態所有を使います。PZ／TB-mBJ／Libxc-PZ／PBEの保持GS接続・独立状態・直接RT再開はNATIVE_GS.md参照。GS計算と出力の移植は残っています。

## ソース配置と既定ビルドの統合（2026-11-19）

通常XCと高速PR版HSEは同じ `runtime/sources/solver` を使用します。
既定ビルドはHSE無効、高速HSEは同じソースから明示的に有効化します。
旧ソースは開発用履歴アーカイブに保存し、wheel/sdistには含めません。
配置変更で数値kernelは変更せず、高速版固有の初期化設定はモデルごとに保持します。

高速PR版HSEには`fftw3.f03`を備えたFFTW3も必要です。軽量LD/SBE wheel自体の依存追加ではありません。
明示的なソースビルドは導入済みFFTWを優先します（`--cmake-arg=-DFFTW_INSTALLDIR=YOUR_FFTW_PREFIX`）。
利用可能なものがない場合は現行HSEビルダーが固定版FFTW3.3.10をダウンロード・構築します。
オフライン環境では事前に用意し、設定ログの`HSE: using installed FFTW`を確認してください。
import時の自動構築・ダウンロードは行わず、通常のUSE_HSE=OFFにはこのHSE用FFTWは不要です。


## ジョブ単位の検証受領書を明示的に使う

同じジョブ内では、最初の完全検証を明示的に再利用できます。

```python
import secrets
from salmonpy.runtime import prepare_runtime_job, load_runtime
validation = prepare_runtime_job(
    "managed-runtime", "new-job-validation.json",
    job_id="CURRENT_JOB_ID", job_token=secrets.token_urlsafe(32),
    required=["tddft", "single-rank"],
)
engine = load_runtime("managed-runtime", required=["tddft", "single-rank"], validation=validation)
print(engine.runtime_validation_report)  # 秘密トークンは含みません。
```

どちらの呼び出しもSALMON・MPI・取得・ビルドを起動しません。ジョブ内で最初に一度
準備し、物質点ごとには作り直しません。受領書はランタイムの外の新しいファイルに
保存し、既存ファイルを上書きしません。ランタイム内へ出力せず、ジョブ中は内容を
変更しないでください。validation=を省略すれば、従来どおり毎回完全検証します。
環境変数や大域キャッシュによる暗黙の切り替えはありません。

代表ランクが割り当てごとに新しい秘密を生成し、RuntimeValidationを呼び出し側の
通信（例：comm.bcast）で同じジョブのワーカーへ配ります。実パス、目録、ファイルと
ディレクトリの更新情報、ジョブIDをHMAC-SHA256で結び付けます。
秘密そのものは受領書に保存せず、repr(validation)でも隠します。
コンテキストの全フィールドをログに出したり、物質入力・チェックポイントのJSONや
コマンドラインに秘密を書かないでください。秘密の有効範囲は一つのジョブです。
次の投入では新しく生成します。job_idには呼び出し側が現在の割り当てIDを指定します。
全スケジューラを自動推定する機能ではありません。長時間ジョブを止める一律の期限は
設けません。配布元の署名や科学的正しさの証明ではなく、信頼するランタイムを
ジョブ内で再利用する仕組みです。配布元・ビルドの信頼確認は従来どおり必要です。

再利用時も、実行ファイル本体は必ず読み、記録したSHA256と照合します。
小さなruntime.jsonのハッシュと受領書の署名を確認し、ファイル・ディレクトリの
更新情報を前後に照合します。付随ライブラリ・データのハッシュ用読み込みと、
ディレクトリ木の再走査を省きます。同じユーザーの悪意ある変更や、古い・粗い更新情報を
返すファイルシステムに対する防御保証ではありません。ファイル数に比例するstatは
残り、実行ファイルが容量の大半なら削減幅も小さくなります。実行ファイルのI/Oは
無くしません。

受領書が古い・読めない、ジョブや秘密が異なる、配置や更新情報が変わった場合は、
警告して通常の完全検証へ戻ります。正しいランタイムなら計算を続けられます。
実データの破損は通常の整合性検査で拒否します。自動で受領書を書き換えません。
必要機能と、実行ファイルに結び付いたNLCC力の検証条件は維持します。
作成時は全ファイルを検証します。更新情報を使う高速経路は現状POSIXのみで、
それ以外では警告して完全検証を維持します。共有FSの更新情報の整合性やノード間の
識別情報は利用施設で確認してください。ローカルMacでの成功は多ノード保証ではありません。

プロセス分離型TDDFTでは、TDDFTMaterialProvider(...)とTDDFTResponseBackend(...)の
両方へruntime_validation=validationを明示します。物質の物理的な識別情報や
電子状態へ秘密を入れません。この指定を持たない設定済みカタログの経路は完全検証を
維持します。埋め込みbuild.json・nativeライブラリの読み込みは今回の対象外です。
自作バックエンドは受領書に対応する必要がなく、非対応を理由に拒否しません。

```sh
# 表示のみ。MPI初期化・ファイル作成・計算を行いません。
python -m salmonpy.examples.mpi_runtime_receipt RUNTIME NEW_RECEIPT --job-id CURRENT_JOB_ID
# 割り当て済み資源内で、施設のランチャーに置き換えて実行。SALMONは起動しません。
mpiexec -n 4 python -m salmonpy.examples.mpi_runtime_receipt RUNTIME NEW_RECEIPT --job-id CURRENT_JOB_ID --loads 3 --report NEW_REPORT.json --run
```

サンプル自身はmpiexecを起動せず、ランチャーを入れ子にしません。
rank0で完全検証後、同じジョブ内でだけコンテキストを配ります。各ランクが明示的に
再利用し、ハッシュ対象バイト数・更新情報の照合数・木の走査数・時間・呼び出し元MPIの
生存を記録します。秘密は記録しません。検証成功後に物理計算を自動開始することも
ありません。投入前の独立した検査や利用者独自モジュールに受領書を強制しません。

ローカルの合成試験（18ファイル、約58MiB、うち実行しないデータ用の代用品2MiB）では、
再利用時のハッシュ対象は実行ファイル扱いの2MiBと小さな目録になり、木の走査はゼロです。
4 MPIランクで各3回の読み込みを確認し、呼び出し元MPIも維持しました。
この合成ファイルは計算用ランタイムとして実行するものではありません。
FEFS・多ノードの性能、GPU、物理量の再現を検証したという意味ではありません。


## 全体通信を最小限にする

ジョブ全体の通信は、必要な初期化・異常時の合意・明示的に要求された最終集約に
絞ります。受領書によるランタイム再利用は各ランクで独立に行います。
検証専用mpi_runtime_receiptは、初期bcast一回とrank0への最後のgather一回だけです。
allgather、追加の生存確認allreduce、Barrierは使いません。最後のgatherでMPIの
生存も確認できるため、全ランクへの診断結果の複製を省きます。
サンプルは検証・報告後に終了し、ソルバーを起動しません。

NativeEHField.advanceのPythonラッパーでは、担当場グループ内の検査を六回から
二回へ統合しました。最初に入力エラーと電流、次にnativeの状態コード・Eの有限性・
候補Aの有限性・積分エラーを合わせて交換し、全メンバーが同じ判断をします。
native内部の物理的な通信をこの回数に含めてはいません。符号付きゼロの一致、
nativeに入る前の拒否、復帰可能な異常の状態保持、合意後の確定は維持します。
数値更新の式や順序は変更しません。

物質内部の総和や近傍交換、必要な異常合意は残します。
ResponseSessionのAは所有点分だけを担当代表へScattervで配り、その物質群の中だけに
配布します。場グループには応答サービスの入力を送り返さず、単一ランク群への追加の
配布も省きます。場パケットも担当点の辞書だけをscatterし、値・単位・位置・時刻を
保持します。Jは従来から所有点分を代表へ回収しており、場グループへの返却と
コールバック・native更新の順序は維持します。

場の代表がplanの先頭でなくても、点IDが交差していて担当点数が違っていても対応します。
初期化時の対応表でJを全点順へ戻し、更新ごとの点検索を省きます。
communication_reportのreceived_A_rows_by_groupに各群の受信行数を記録します。
logical_A_bytes_per_stepは論理的なデータ量で、実ネットワーク転送量ではありません。
全体通信を追加せず、異常・有限性の合意は残しています。

保存互換性のためprivateなresponse._Aの形は全点のままですが、物質側では担当行だけを
現在値として更新し、場の代表だけが応答サービスの全点A履歴を持ちます。
非担当行は複製された観測データではありません。全体の観測には場の代表の出力を
使ってください。全行を保存した旧チェックポイントも読めます。
nativeの場A/Eと物理状態をゼロで代用する変更ではありません。
native場の担当J配信はCOUPLED.mdの専用節を参照してください。外側の連続する入力合意もCOUPLED.mdの専用節に従い統合しており、
混在計算全体の通信が最適化済みという意味ではありません。
