# 高速PR版HSEと対称性

[English](FAST_HSE_SYMMETRY.en.md) · Topic: `fast_hse_symmetry` · Revision: `2026-11-18.8`

HSEの開発・検証は高速PR版を本筋とし、代表k点・star再構成・orbit再利用・ACEと対称性を
中心に進めます。通常版HSEを代替エンジンにはしません。標準SALMONの対称性縮約とは演算が
異なるため、標準側に合わせて高速PR版のアルゴリズムを変更しません。

| 層 | 現在の範囲 |
|---|---|
| 高速PR版の交換コア | 対称性マップ・位相キャッシュ・FFT資源を独立所有。周期・非スピン分極・CPU・full-q・内部kランク1／2を検証 |
| 対称性／orbit／ACEの比較 | 縮約4代表k点と独立に列挙した全8k点。C2z／mirror-x、異なるセル・遮蔽係数・軌道 |
| 内部Fortran TDDFT | 明示opt-in、1プロセス1初期化／内部k1／2、factor-off Taylor。直接／候補再開を検証 |
| Python／native-C／FDTD | 新プロトコルの明示opt-in。group1／2・factor-off Taylorの保護付き混在を検証。既定は拒否 |
| 保存復元 | 内部交換コアV2と明示的内部TDDFTの直接／候補再開を検証。限定した公開の保護付き混在再開を検証 |
| 対称性付き内部k MPI2 | 専用capabilityを持つ新ライブラリで明示対応。全／縮約比較と直接再開を検証 |

## 独立した交換コアの所有

内部の `exx_bind_system` 経路で、限定した内部kランク1／2の対称性付きコアを扱えるように
しました。各systemがk点設定を持ち、各HSE ownerが対称性マップ・位相・FFT資源を所有します。
全格子検査にはそのsystemのngrid・primitiveセル・格子幅を使い、共有num_rgridには依存しません。
対応するrt_status入口では、格子次元／点数の不一致や未対応owner分割をrefresh／actionの変更前に
拒否します。それ以外の不正な対称操作・原子配置・重みはnative停止の対象であり、一般的な
数値異常復帰を保証するものではありません。新しいコア範囲はexx_factor_mode='off'とTaylor／ACE
または比較用full-actionに限定し、学習／PT-CNの対称性は更新前に拒否します。

セル・物質・対称群を変更する際は新しいownerを作り、生存中の所有設定を書き換えません。
軌道の更新は既存のrefresh契約で行います。FFTプラン・コミュニケータを所有する状態を
単純コピーしないでください。GSの自動binding方針は変更していません。

ローカル試験では、C2zとmirror-xの2状態に異なるセル・遮蔽係数・軌道を持たせ、交互の更新・
作用、格子／分割の拒否、解放／再作成を検証しました。full-actionとACE、通常のstar展開と
高速orbit再利用が、OMP1／2・MPI1／2プロセスで通っています。MPIプロセスごとに設定と呼び出し回数も変え、出力の有限性と一致を検査します。小規模な交換コア回帰であり、
収束した実材料や縮約k点TDDFTの検証ではありません。native試験は
`_runtime_source/tests/native/test_hse_symmetry_resources.f90` に同梱します。
対応する監査済み高速PR版モジュール／ライブラリから明示的にコンパイルします。

## 状態所有の場の検査

内部Fortran入口 `symmetry_check_field(map, A, status)` と
`exx_validate_vector_potential(system, A, status)` で、一様なA/cが保持したCartesianの
対称操作と両立するか検査します。共有の対称性設定は参照しません。statusは0が適合、
1が不正／未準備、2が場による対称群の破れです。非ゼロAは最大成分で規格化してから
方向の許容誤差1e-10を適用するため、振幅が微小でも禁止方向は通りません。
ゼロAにも有効なマップが必要です。この入口は状態変更・通信を行いません。

TDDFTアダプターでは数値状態のバックアップ／時間発展より前に検査し、既存の集団通信で
失敗を共有します。独立したC2z／mirror-x、傾いたセル、微小／巨大な場、不正値、設定不整合、
終了済みownerを検証しました。場とマップの境界の検査であり、対称性情報全体やnative初期化の
保証ではありません。通常worker／native-Cのopenと結合保存復元の禁止は維持し、無制限な対称性付き時間発展・
結合再開への対応とは扱いません。明示opt-inの範囲は後節を参照してください。明示ビルド用プローブは
`_runtime_source/tests/native/test_hse_symmetry_field.f90` です。

## 内部交換コアの対称性checkpoint

既存の内部 `hse_checkpoint_write/read/capture/install` で、対応済みの周期・非スピン分極・
CPU・full-q・内部kランク1／2の交換コアについて、値だけの対称性マップ／キャッシュを扱えます。
非縮約ownerは従来のV1をバイト一致で維持します。対称性ownerは同一バイナリ用V2として、
保持したマップ・kernel識別情報・source／ACE／履歴値・任意の全位相キャッシュを保存します。
このraw streamは呼び出し側でアーカイブ・モデルを認証する前提で、独立した可搬形式ではありません。

復元先には一致する独立ownerと演算子を用意します。演算子構築は行いますが、受理済み時間ステップの
再生は行いません。読込時には、静的なマップ値と信頼済み配列寸法を照合してから、範囲を制限した
キャッシュを候補へ読みます。検査だけ・候補作成だけでは現状態は変わりません。途中で切れた保存、
異なるマップ／位相分割、不正キャッシュは拒否します。反映は値配列の移動のみで、FFT／C／MPI資源は
複製・置換しません。

nativeプローブ `_runtime_source/tests/native/test_hse_symmetry_checkpoint.f90` で、full／ACE、
star／orbits、キャッシュ有無、OMP1／2、MPI1／2プロセスの独立COMM_SELFコアを検証します。
ownerの作り直しを含め、復元直後の作用・エネルギーと、その後のrefresh／actionのcheckpoint
バイト列が連続実行と一致します。これは交換コアの保存復元です。別枠の内部TDDFT検証は次節を参照してください。
既定のPython／native-C TDDFTの拒否は維持します。限定した公開opt-inは後節を参照してください。学習／PT-CNの
対称性、縮約内部k MPI、対称性付き混在計算の失敗・再開は残作業です。

## 明示的な内部TDDFTの検証

内部Fortranの `rt_open(..., allow_symmetry=.true.)` で、交換コアとTDDFT全体のcheckpointを
接続します。周期・非スピン分極・CPU・full-q・内部kランク1／2・factor-offのTaylor／ACEと、
zを保つ対称群に限定します。full-actionは比較検証用のみです。初期化は1プロセス1回で、
独立したTDDFTの実行には新しいプロセスを用意します。

密度／電流の対称化情報はまだプロセス共通・不変です。フラグ・操作行列・格子・逆行列を所有入力と
照合し、step／current／checkpointと、復元候補の反映直前に検査します。不一致では候補を反映せず、
設定を戻せば反映を再試行できます。異なる物質の同一プロセス内共存や、不正な初期化資材からの
制御された復帰を保証するものではありません。

実Aの所有検査は、場が対称性を破る場合に内部status34を返します。バックアップ／時間発展前の
拒否なので、内部呼び出し側はAを直して再試行できます。このstatusは後述の新しい公開プロトコルで宣言します。旧バイナリから推定しません。
通常のworker／旧C入口は対称性を拒否し、新しい明示opt-inだけが限定範囲を有効にします。

nativeの `test_hse_symmetry_tddft.f90` と、補助ツール
`prepare_fast_hse_symmetry_fixture.py`／`test_fast_hse_symmetry_tddft.py` を
`_runtime_source/tests/native`／`_runtime_source/tools` に同梱します。対応するMPI高速PR版の
モジュール／ライブラリへ明示的にコンパイルし、ハーネスへ
`PROBE PSEUDO NEW_OUTPUT --ranks 2 --threads 2` を渡します。各プロセスは独立したCOMM_SELFへ
接続し、縮約内部k MPI2ではありません。8^3格子・全8k点と縮約4k点のC2z／mirror-x合成状態で、
電荷電流／エネルギーを比較し、奇数／偶数の直接・候補再開で内部トレースと最終checkpoint全体の
バイト列が一致することを検証します。受理済みステップは再生しません。GS／光学収束の実計算例ではなく、
初期状態作成でダウンロード・SCF・ビルドは行いません。

対称性付き保存では、保存されない重なりキャッシュの無効フラグを保存用コピー内でtrueにそろえます。計算中のフラグと非縮約V1は変更しません。

## 公開native／Pythonと混在opt-in

`NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)` を明示します。
既定Falseは従来の `sp_fast_hse_init` を選びます。Trueには `sp_fast_hse_init_symmetry` と
`sp_fast_hse_symmetry_protocol()==1` が必要です。新initの引数は従来と同じ
`(const char *path, double dt, int capacity, double J[3])` です。opt-inをbackendの識別情報に
含め、内部kランク1／2を許可します。入力操作行列は変更しません。未対応・未知のプロトコルは
native寿命を取得する前に拒否します。学習／PT-CNの対称性は未対応、full-actionは
比較検証用のみです。

有限かつ全ランク一致のstatus34と専用プロトコルが揃った場合だけ、Pythonは受理状態を保つ
`InputError` を返します。Aを修正してから再試行します。群は入射偏光だけでなく、実際のFDTD局所場と
適合しなければなりません。Aの投影・微小成分の切捨て・別エンジンへの自動切替は行いません。
全体の巻戻しには `CoupledStepper.step_protected`／`--protect-steps` を使い、
`StepRecoveryError.original_preserved` がTrueかつ解放成功の場合だけ再試行します。

通常の `mpi_mixed_native` 例は `--hse-symmetry` と、等価な全k点初期状態を明示した
`--hse-reference-bundle FULL_K_BUNDLE` を受け付けます。従来のシリアル高速PR版と比較するための
参照で、任意のGSを自動変換しません。`--hse-ranks 1` または `--hse-ranks 2` を指定します。配置・LD・SBE・半局所TDDFTと
直接保存再開の入口は従来どおりです。5／6ランクで、field1／2・場内2種類LD・SBE1・PBE-TDDFT2・
対称性HSE1を検証しました。故障・保護付き復元と新規直接再開で、受理済みステップの再生なしに
軌跡を維持します。合成初期状態の接続試験であり、物質・場の適合性と物理収束は別途確認します。

## 極小勾配XCの安定評価

Libxc5.2.3のHSE06は、WPBEh(0) − 0.25 WPBEh(ω) ＋ PBE相関を組み合わせます。
WPBEh(0)のE1項とlog項を別々に微分すると、ほぼゼロの勾配で大きな項の相殺が起きます。
元のphase0合成状態では全8／縮約4 k点の電流差約5.7e-10になり、刻みの縮小でも消えませんでした。

現在の高速PR用完全ソースでは、縮約勾配s≤1e-4で同じ式をu=s²で評価し、相殺を解析的に
まとめてから微分します。捨てるE1級数の項は倍精度以下です。遮蔽交換のuが丸め誤差以下の
場合は解析的二次極限を使い、相関は元のLibxcで評価します。範囲外の点のLibxc出力は変更しません。
密度・勾配・Aのfloor／平滑化、許容誤差の緩和、通常HSEへの変更は行いません。
スピンのスケーリングと固定coreの密度・勾配も検証します。

元のphase0を変えず、32ステップで全／縮約k点の電流・エネルギー一致、FDTD結合、
故障復帰・保護付き復元・新規プロセスでの直接再開を検証しました。
通常のPBEは配列のバイト一致、高速HSE Taylor／PT-CNと混在は丸め誤差内の差を確認しました。
合成状態の数値検証であり、実材料の物理収束の保証ではありません。

明示的な検証入口（ビルド・ダウンロードは行いません）:

```sh
python tools/test_hse_low_gradient.py BUILD_RECORD COMPILED_NATIVE_PROBE NEW_OUTPUT
mpiexec -n 2 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD_RECORD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT
# 元のFDTD局所場を使う場合: --drive ORIGINAL_ARRAYS_NPZ --point-index 1
```

native probeは `tests/native/test_hse_low_gradient.f90` を、対応する新しい
`hybrid_semilocal.mod` とライブラリに対して明示的にコンパイルします。
80桁エネルギー差分、密度・遮蔽係数・切替境界・スピン・coreと通常勾配のバイト一致を調べます。
ライセンスと式の由来は `runtime/NOTICE.md`。元ソース・既存バイナリは変更しません。
全k参照ワーカーとMPIライブラリの**両方**を新ソースから明示的に構築してください。
対称性protocol1だけでは、この数値修正の有無は識別できません。

`probe_hse_low_gradient.py BUILD_RECORD NEW_JSON` は今後も生のLibxcを診断します。
wrapperを直しても `low_gradient_sensitivity_detected=True` になり得ます。
この診断をwrapperの失敗・合格と取り違えず、新しい上記検証とソース／バイナリのハッシュを残します。
過去のphase.23は接続用fixtureとして残しますが、今回の修正確認は元のphase0です。

次はより広い初期化・並列分割・実材料での検証です。高速PR版のみの方針を維持します。
メモリ予算は警告のみ、分子の作業は保留です。


## 縮約k点の内部MPI（2ランク）

内部2ランクでは `sp_fast_hse_symmetry_kmpi_protocol()==1` が追加で必要です。
旧ライブラリはnative寿命の取得前に拒否します。CPU・非スピン分極・周期・full-q・
factor-off Taylorの範囲を維持し、各ランクが連続した代表k点とそのstar位相を所有します。
空の担当範囲、実空間／軌道分割、学習／PT-CN対称性は今回の対象外です。

保存時は全k点を各ランクへ複製せず、担当範囲の位相キャッシュを保存します。
復元先で独立に初期化した範囲・MPIランク数・担当と照合し、不一致を拒否します。
**保存時と同じMPI分割で再開してください。分割変更を伴う再開には対応しません。**
状態を保つstatus34はMPI内で最大値を共有し、ランク数倍しません。

```sh
mpiexec -n 4 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD_RECORD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT --group-size 2
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD_RECORD REDUCED_K_BUNDLE NEW_SAVE --mode save --group-size 2 --checkpoint CHECKPOINT
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD_RECORD REDUCED_K_BUNDLE NEW_RESUME --mode resume --group-size 2 --checkpoint CHECKPOINT
```

二つの独立群（各2ランク）で異なる遮蔽係数を使い、7→32ステップの電流・エネルギーと
最終checkpointバイトを連続実行と照合します。受理済み時間ステップの再計算はありません。
初期化に必要な演算子の構築は行います。元のphase0・弱い局所場で全8／縮約4 k点を
OMP1／2で比較しました。通常の混在例も `--hse-ranks 2` を受け付けます。
旧ライブラリの拒否検査は `--mode deny-old-k2`（旧ビルド専用）です。
従来のCOMM_SELF検証の記述は履歴として残します。物理収束の保証ではありません。


## 初期化前の対称性ファイル検査

Python側でnativeの読込み・寿命取得・作業コピー作成より前に `sym.dat` を読み取り検査します。
標準形式は1操作につき3行、各行4個の有限実数です。D指数、コンマ、`2*0` の反復指定と
行末の `!` コメントを扱います。入力を書き換えたり、不足値を補ったりしません。
`yn_symmetry` の無効軸による選択をFortranと同じ順で適用し、保持されたfast-PR用の
符号付き置換行列、恒等操作、格子並進を法とする重複、群の閉性を検査します。

一部のランクだけに不正入力があっても、既存の群内合意で全員が `InputError` を返します。
native寿命は未取得なので、元入力を修正し、物質providerを作り直して再試行できます。
保持済みproviderの入力を後から変更すると既存のハッシュ検査に拒否されます。
続く幾何事前検査は下節を参照してください。k-star・擬ポテンシャル内容・再開バイナリの
互換性と、その他のnative初期化・停止経路は残作業です。

検証用には、前節のsaveコマンドへ `--admission-fault truncated`、`nonfinite` または
`nonclosed` を追加します。各群の最後のランクだけを壊し、全員の拒否と修正後の開始を検証します。
通常入力の計算式と結果は維持します。初期化異常全体が復帰可能になったとは解釈しません。

軸指定は小文字 `y`／`n` の1〜3文字です。不正な指定を黙って補正せず、事前に拒否します。


## セル・原子配置・格子の事前検査

`NativeFastHSEBackend.create_state` は、対称性あり／なしの周期高速PR入力について、
native読込み・寿命取得・作業コピーの前に有効セル・格子・原子配置も読み取り検査します。
ASE/f90nmlの追加依存はありません。各ランクが同じcollective入口へ参加する必要があります。
検査で拒否されたら元入力を直し、新しいproviderを用意して再試行します。

- セル：`al`、対角`al_vec1..3`、`num_rgrid*dl`、ベクトル指定の優先順位をSALMONと合わせ、有限・正の長さ／体積、逆行列を検査。現native接続は軸整列直交セルに限定し、非直交・軸入替は事前拒否します。
- 格子：正の整数点数／幅とnative整数範囲を検査。実際の点数はSALMONと同じ `int((length+1e-10)/spacing)`、幅指定のセル整合許容差も同じ。ceilや入力の補正はしません。
- 原子：inlineまたはbundle内の明示ファイルのCartesian／縮約座標を1つだけ指定。原子数、有限座標、1..nelemの種別IDを検査。既存の高速PR接続と同じく入力単位はauのみです。
- 対称性あり：セルのCartesian等長性、全格子の整数写像（角で誤差と整数範囲を評価し、格子全体の配列は作らない）、同一種別の原子の一対一対応を検査します。

数値設定は実数／整数literal、D指数、反復、正の添字・範囲指定を扱います。
null項やその他の未対応書式は推測せず拒否します。座標は標準の「label x y z species」の5項形式。
ファイル名の大文字小文字を保持し、外部へ抜ける参照は拒否します。入力は書き換えません。

検証ツールの `--admission-fault cell|cell-metric|cell-orientation|atoms|grid|grid-translation` で、各MPI群の1ランクを
壊し、全員の拒否・寿命未取得・修正後開始を確認できます。正常な電流・エネルギーと直接再開も維持します。
これは入力の構造検査です。GSとの同一性、k-star、擬ポテンシャル内容、再開バイナリ、未対応の
native初期化を検証したことにはなりません。物理収束も別途確認してください。

軸を入れ替えた右手系／左手系の合成入力で現ランタイムのACE初期化停止を確認したため、
これらを黙って変換せず対応外として拒否します。一般セルのnative初期化は残作業です。
