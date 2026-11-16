# SALMONpy: AI向け利用・実装ガイド

LD単独利用のPython試作APIを提供しています。実際の入口は `LorentzDrude` の `respond`、`new_state`、`permittivity` です。構築方法・単位・時間位置は README.md と examples/ld_response.py を参照してください。FDTDには外部実行ファイルを使うバッチ入口があります。WORKFLOWS.mdを参照してください。結晶は共通A駆動インターフェースを持ち、内蔵SBE速度ゲージと初期状態ローダーを提供します。TDDFTはPZ／TB-mBJ／Libxc-PZ／Libxc-PBE・中点伝播のプロセス分離型接続を提供します。存在しないAPIを例示しないでください。

## 物質応答の境界

結晶の公開インターフェースはA駆動のみです。TDDFTとSBEは、内部状態と新旧Aを受けて時間発展しJを返す交換可能な物質応答として設計します。結晶にE駆動の選択肢を追加しないでください。SBEは既存の演算順を保った共有Fortranカーネルに接続しています。

E,Hを扱うFDTDとA駆動の結晶の間には結合処理を置きます。単位、ゲージ、時間位置、符号、座標回転、Jの体積規格化は現行マルチスケール実装に準拠して確定し、アダプターで明示します。Aの積分履歴は結合側、波動関数・密度行列は物質側が保持します。

ドルーデ＝ローレンツはE駆動の構成則として独立させます。結晶用A駆動APIに合わせるための数式変更はしません。

## 今回利用できるFortran境界

`native/lorentz_drude.f90` の `s_ld_state` がP、J、旧J、旧E、係数、次の電場更新用の補助量Qを所有します。公開関数は `ld_allocate`、`ld_coefficients`、`ld_store_old`、`ld_update`、`ld_finalize` です。具体的な配列形状・呼び方は `tests/native/test_ld.f90` と、外部に準備したSALMONソースの `src/maxwell/fdtd_eh.f90` を参照してください。

`ld_coefficients` の累積先は呼び出し前にゼロ初期化します。`ld_store_old` と `ld_update` は既存演算順を保つため分割E成分を受け取ります。`rj*_fdtd_ld` は物理電流の単純和ではなく補助量Qです。

`fdtd_eh_kernel.f90` は差分演算 `eh_fd` を提供します。`fdtd_eh` は引き続き時間ループ、境界、通信、観測、再開ファイルを調整します。実験的な単一ステップ・結合入口は `salmonpy.coupled.NativeEHField`。COUPLED.mdの範囲に限定します。

## 変更してはいけない更新順

旧J/E保存 → E更新（LD係数を含む）→ 光源 → 旧Q注入 → E通信 → LDのJ/P/Q更新 → 線形応答 → 旧H保存 → H更新 → H通信 → TTM → 観測・保存。

LDは電流を返すだけでは分離できません。E更新係数への寄与も必要です。係数を固定したまま時間刻みを変更しないでください。独立計算には別々の状態インスタンスを使います。

## 検証と禁止事項

ビルド前に `python tools/check_forbidden_eigensolvers.py --source-root native` を実行します。構成・増分ビルド時の検査を外しません。新規実行ファイルは全アプリケーションオブジェクトを検査し、再利用時は記録されたSHA256と照合します。禁止された固有値ドライバーへのフォールバックはありません。PZHEEVDの検証には固有値、残差、直交性、ワークスペース容量を含めます。

IEEE判定は各要素へのスカラー呼び出しに限定します。GSデータを収集しません。計算結果は利用者が指定した書き込み可能な出力先へ保存し、コードツリーに格納しません。

変更後は独立LDの旧演算との比較とFDTD全体の比較を行います。単一スレッドのビット一致と並列の数値許容差検証を区別します。再開互換性や未検証の計算機での動作を推定だけで保証しないでください。

## 入力の寿命

INPUT_CONTRACT.md に従います。材料係数はモデル作成時、dtと初期条件は状態作成時、Eは各stepで指定します。`material.config` と `state.config` は確定値の読み取り専用レコードです。設定変更は新しいモデル・状態の作成で行い、private属性を書き換えません。`respond` は毎回新規状態です。設定レコードはチェックポイントではありません。

## 配布

本体ランタイムの準備は RUNTIME.md / RUNTIME.en.md と manual.json の `runtime` タスクに従います。既存の監査済みランタイムを再利用し、次にOS・CPU・必要機能に合う配布物を探します。対応品がない場合だけ明示的な `python -m salmonpy.runtime build` を案内します。importやinstall失敗を理由に自動ビルドしません。独立SALMONpyの非公開Releaseと期限付きActions配布物を区別し、公開URLを作り話で補わないでください。Libxcを必要とする計算には `--require libxc` を付け、PZ版に勝手に代替しません。取得・検証した `load_runtime` の戻り値を既存の `executable=` に渡します。

wheelの構築と修復・OS別検証は PACKAGING.md に従います。検証済み環境とCI設定のみの環境を混同しません。wheel利用者にはFortranコンパイラを要求せず、ソース構築時のみ要求します。

構造・波形・解析の単独入口は WORKFLOWS.md に記載しています。FDTDの一ステップAPIはCOUPLED.mdの条件のみ例示できます。周期構造、斜め入射、未対応のTDDFT計算方式を対応済みとして例示しないでください。共通Acos2パルスとFDTD側のCEP規約の差はアダプターで0.25周期補正します。


## SBEの接続

`salmonpy.sbe.SBEInitialState.from_salmon` は原子単位・単一スピンの
`_k.data`, `_eigen.data`, `_tm.data` を読みます。`bands` はファイルの全バンド数、
`active_bands` は任意の伝播バンド数です。占有数は現行SBEと同じ閉殻規約です。
`new_state(dt=..., initial_A=...)` ごとに独立した密度行列を生成します。
結晶はA/c駆動のみ。中点Aで四次Taylor伝播、新しい端点AでJ評価の順を保持します。
補正次数と非局所補正は初期モデルで固定し、dtは状態作成時に固定します。
SBEの長さゲージ・バイナリ再開は未対応です。四次Taylor法なので厳密な
ユニタリー性や自動誤差制御を保証しません。APIは刻み幅を自動変更しません。


## TDDFTの接続と寿命

`salmonpy.tddft.TDDFTInitialState.from_salmon` は自己完結した入力ディレクトリを
複製し、ハッシュ付きの初期スナップショットを所有します。初期GSのバイナリ解釈は
SALMONの既存ローダーが行います。保存用バンドルは1組を状態間で共有し、各 `new_state` は専用プロセスと独立した作業コピーを所有し、
波動関数、Vh履歴、予測子修正子、XCや通信の共有変数を他の計算と共有しません。
`with` で解放してください。共通 `respond` は成功・失敗とも自動的に解放します。
PZ／TB-mBJ／Libxc-PZ／Libxc-PBE、周期系、固定イオン、非スピン分極、中点伝播、CPU一MPIランクが対象です。
実行ファイルは検証済み配布物を取得するか、明示的に構築・監査し、SHA256を指定します。wheelだけでTDDFTが動くとは
説明しないでください。dt、最大ステップ数、初期Aは状態作成時に固定します。
共通の出力は電荷電流密度 j です（2026-11-14 から）。SALMONの電子物質電流 J = −j は状態の内部に保持し、`state.matter_current` で参照できます。符号を二重に反転しないでください。
異常終了時は受理済みA履歴を再実行して復元し、復元失敗時は状態を閉じます。
実GSによる検証と合成波動関数による接続回帰を区別してください。

Libxc-PZ/PBEには `USE_LIBXC=ON` の外部SALMONが必要です。Libxc全汎関数の対応ではありません。GSとRTで汎関数・擬ポテンシャルを一致させてください。現ローダーはGS生成条件との整合を自動判定しません。PBEの検証はヘッダーでPBE・NC・スカラー相対論・NLCCありを確認したSi UPFを使用しました。詳細は docs/libxc-validation-2026-11-13.md を参照してください。

## 後処理用アーカイブ

OUTPUTS.md に従います。`crystal.respond(..., output_directory=...)`、
`TDDFTState.finish(directory)`、`outputs.OutputArchive` と `outputs.save_arrays` が
実在する入口です。TDDFTのfinishは状態を閉じ、初期ファイルハッシュと選択した
本体出力を残します。大きい分布は初期入力のanalysisで選択・間引きます。
波動関数チェックポイントとprojection/unfoldingは未対応です。NPZアーカイブを
restartファイルとして渡さないでください。別ノード通信やコア固定は未実装・未検証です。

## バンド図

BANDS.mdを参照してください。`bands.BandStructure.from_salmon` と `.plot` が入口です。
固有値の単位、実際のCartesian k、セルの選択、余剰点、収束対象バンド、エネルギー基準を
明示してください。経路上の価電子帯上端を全域のVBMと呼ばないでください。

バンドの計算入口は `BandInitialState.from_salmon` → `run(BandPath, ...)`。
初期束の外の新規保存先、対応した監査済み実行ファイルを使い、未収束を成功扱いしません。

ASE連携は任意の独立モジュール `salmonpy.ase`。ASE.mdを参照。
`run_bands` は構造を照合し、ASEのk点列をそのまま計算します。Calculator・力は独立したASE_CALCULATOR.mdの対応範囲を参照してください。

GS作成は `salmonpy.dft.DFT`、ASE変換は `salmonpy.ase.dft_from_atoms`。DFT.mdに従います。ASE入口の電子数はKY/NC-UPFの価電子数から中性系として決定・照合します。直接DFTでは電子数を明示します。軌道数・格子条件・k点は必須で、収束保証のある自動既定値はありません。
`salmonpy.ase.tddft_from_atoms` はGSを保存して既存TDDFT初期状態を返します。ASE.en.mdとexamples/ase_tddft.pyも参照してください。結晶はA/c駆動、戻り値は電荷電流j。状態は独立で、ランタイムの自動取得・コンパイルはしません。
収束済み `GroundStateResult` からバンド／TDDFT初期状態を作れます。PZ／Libxc-PZ／Libxc-PBE・KY／NC-UPFのCalculatorと本体の力はASE_CALCULATOR.mdを参照してください。

`salmonpy.ase_calculator.SalmonCalculator` は任意依存の独立モジュール。energy/forcesのみ、PZ／Libxc-PZ／Libxc-PBEとKY／NC-UPF・単一ランクが対象です。Libxc対応ランタイムと擬ポテンシャルの整合を確認します。
NLCC付きUPFの力は検証記録なしではCalculatorが計算前に拒否します。信頼できるバンドルをload_runtime(..., required=['nlcc-forces'])、または監査済みビルドと検証記録をload_validated_buildで読み込んだ場合のみCPU・非スピン分極NC・PZ/Libxc-PZ/PBEの力を有効にします。SalmonExecutable(path,sha256)だけでは無効です。記録はバイナリ・既知パッチ・微分検査・固定Si8実行検査に結び付け、移動後にも実行検査します。互換パッチのXC重み・全密度勾配・明示的コア力の修正は旧バイナリに自動適用されません。実装回帰と物理的な収束を区別し、メタGGA・OpenACC・スピン分極Libxcの力を検証済みと説明しないでください。
`warm_start=True` は直前の互換なGSの波動関数だけを再利用し、既定Falseです。セル・計算条件が変われば新規SCF、元データが破損すればエラーです。旧原子位置・混合履歴は再利用しません。
結晶TDDFTの基本はtransverse（trans_longi=tr）。Longitudinalは明示的な場の選択です。孤立系・non-collinear・Ehrenfestは初版の目標として残しますが、今回のASE APIの対応済み機能には含めません。
粗い格子のエネルギー有限差分はegg-box効果を含むため正解と仮定せず、剛体並進と格子収束を調べます。

## 初心者への最初の案内

manual.jsonの `beginner_intake` を必ず適用します。新しい計算を選ぶ前に、未回答の
「使ってよい総メモリ」「使ってよいCPU並列数」「やりたいこと」を短くまとめて聞きます。
既に回答された条件は再質問せず、「不明」も受け付けます。装置の全メモリ・全コア数を
自動的に利用許可と解釈しません。目的に応じて構造規模・精度・所要時間の希望を補います。

最初に、今できること／追加準備が必要なこと／未対応と、予算内で始められる小規模例を示します。
計算規模の上限は実測・暫定見積もり・不明を区別します。未測定の原子数上限などを断言しません。
必要なら予算に収まる小規模試行でピークメモリ・時間・保存量を測定し、条件を記録して更新します。
見積もりはメモリ強制制限ではありません。自動容量推定器・ハードメモリ制限は現状未実装です。

プロセスとスレッドを合計したCPU使用量、全子プロセスを含むメモリ、保存容量を扱います。
メモリの余裕は許可された上限の内側に確保します。現在のPython起動が単一MPIランクである
機能を、任意のMPI並列数で実行可能と案内しません。収まる計算規模と物理的な収束精度は別です。

実行可能な用途別Pythonサンプルと必要な入力は [examples/README.md](examples/README.md) を参照してください。AI検索の各taskの `sample_files` からも辿れます。

## 初見の環境診断と検索

`python -m salmonpy.doctor --json` でライブラリの読み込み可否、追加依存の所在、
任意指定した外部実行ファイルのハッシュ照合を確認できます。計算・コンパイル・導入は行いません。
外部SALMONの `hash_matches` を機能対応や実行成功と読み替えないでください。
`manual.search` は常に `getting_started` と `task_catalog` を返します。一致しない検索語では
目的別一覧を示して目的を確認し、対応不能とは即断しません。

本リポジトリはSALMONpy専用です。同梱本体ソースはruntime/sources、段階別完成形はruntime/overridesです。ビルド用srcは `runtime/README.md` の手順でリポジトリ外に準備します。

## 並列経路を選ぶ

LD中心の光伝播には `fdtd_mpi` タスクを選び、既存MPI FDTDの各ランク内で多数のLD状態を管理します。
起動元は直列Pythonです。MPI版の完成ビルド記録が必要で、配布用直列ランタイムとは区別します。
所定のAに対する独自応答は `distributed_response`、二点TDDFTは `distributed_tddft` を選びます。
これらはPythonを一度MPI起動し、標準例では一物質群に一物質・一マクロ点を割り当てます。
内部MPIのTDDFTとE/Hフィードバックには `native_coupled` タスクとCOUPLED.mdを使います。

## 場の部品（j→A）

`salmonpy.maxwell` の場の部品は `step(j_n) → A_{n+1}` です。`Transverse`／`Longitudinal`／`ThinFilm`／`Multiscale1D` は SALMON の tr／lo／2d／1次元多スケールと同じ式・順序・初期条件です。部品間のインターフェースでは電荷電流 j を受け渡し、物質の内部は SALMON の物質電流 J を保持します（j = −J）。crystal プロトコルは j を返します。J_matter を宣言した部品だけ `maxwell.couple` が変換します。場の部品に J を直接渡さないでください。t=0 の電流は場に戻しません。1次元では x が伝播方向です。座標回転、背景誘電率のある1次元媒質、2次元・3次元の Weyl 伝播、再開は未対応として扱い、対応済みと例示しないでください。時刻配置は `salmonpy.timing`（`salmonpy.time/1`）で宣言し、`couple()` が結合前に検査します。新しい部品には `timing` を実装し、検査を回避する変換を書かないでください。詳細は FIELDS.md を参照してください。
これは一群一点・一プロセス一回のネイティブ寿命に限定した実験的APIです。
既存FDTDバッチと応答サービスを並べるだけで結合したと説明しません。
ネイティブ結合のE時刻n·dtと電子時刻(n+1/2)·dtを混同しません。再開、回転、サブサイクルは未対応です。
実行前に不足する資源予算を確認し、例のランク・スレッド・状態数と、未測定のメモリ使用量を明示します。
機能ごとの実行例・出力形式は examples/README.md / README.en.md にあります。

## 次元選択の推奨

電子応答を結合する計算では、TDDFTは負荷が大きいため、まず1次元の伝播・層構造を提案します。
3次元を推奨するのはSBEまたはユーザー作成のカレント計算コードを使う場合です。
3次元TDDFTが不可能とは説明せず、必要性と状態数・メモリ・時間の予算を確認します。
ここでの次元は場・マクロ配置で、微視的な電子状態の次元ではありません。
現行40³・二点サンプルは接続回帰用です。自動1次元化、未実装SBE結合アダプター、
周期境界を対応済みとして案内せず、実在する入口の制約を確認してください。LD単独FDTDは別の経路です。

## 資源計測と計算の成否を分ける

導入・計算・数値検証・資源計測の結果を別々に記録します。任意の計測器や
`time -l`、sysctl、監視APIの失敗だけで、成功した計算を失敗と判定しません。
計算プロセス自身の終了コード・ログ・出力検証を根拠に判定し、非ゼロ終了や
数値検証の失敗を計測問題として除外しないでください。ラッパーの終了コードしか
分からない場合は計算成否を不明とし、既存出力を保存して、新しい出力先で再試行します。

取得できないピークメモリ等は unavailable / not_measured と null で記録し、0や推定値で
代用しません。計測対象（親プロセスのみ／子を含むジョブ全体）・単位・取得方法を残します。
計測に失敗しても予算内という保証はできません。ハード上限が必要な実行では、有効な
制限手段を確認するまで実行・拡大しません。初回と環境修正後の結果は分けて残します。

wheelからの利用テストはソース側の `docs/wheel-user-test.md` / `.en.md` と
`tools/test_wheel_user_workflow.py` を使います。指定wheelを新規venvに導入し、
LD・合成SBE・波形・形状・保存の公式例を実行します。ソースビルドへ自動切替せず、
wheelハッシュと実際の読込元を記録します。これはスモークテストであり、先行Notebookの
追加の解析解・保存則検証や、外部SALMON・MPIの検証を含みません。

## 公開契約とメモリ上の状態保存（0.1.0a5）

`public_api` と `snapshots` タスクを参照。公開名は `salmonpy.api.PUBLIC_API`、
版は `salmonpy.__version__`。例外は `salmonpy.errors` の親と分類を使い、
既存の組込み例外による捕捉も維持します。OS・依存・コールバックの例外がすべて
SalmonpyErrorになるとは説明しません。詳細はAPI.md/API.en.md。

プロセス分離TDDFT・LD・SBE・Pythonの4種類の場が `snapshot()`／`restore(token)` に対応します。
同じモデルオブジェクト、同じdt・開始時刻と場の容量・駆動条件を要求します。
保存はメモリ上の不変コピーで、失敗したrestoreは状態を変えません。保存数と状態サイズの
分だけ追加メモリが必要です。トークンは編集・シリアライズせず、Zarr GSや出力NPZと混同しません。
手動結合では同じ受理済み時刻で場と物質を両方保存します。MPI全体の原子的復元、可変dt、
埋め込みTDDFTやネイティブE/Hの復元は未対応です。プロセス分離TDDFTのディスク再開は下記の専用APIで対応し、
出力アーカイブやネイティブ波動関数checkpointとは区別してください。
examples/snapshot_restore.pyは実行可能な例です。結晶の基本はtransverse・A/c駆動です。

プロセス分離TDDFTのrestoreは別ワーカーでA履歴を再生し、各端点電流をrtol=1e-12・atol=1e-14で照合後に切替えます。失敗時は元ワーカーを保持します。同一初期モデル・dt・初期A・容量・実行ファイル・スレッド数を要求し、復元コストは履歴長に比例、検証中は2ワーカー分の資源が必要です。SNAPSHOTS.mdを参照。

## TDDFTのディスク再開

TDDFTState.save_checkpoint(directory)とrestore_checkpoint(directory, max_bytes=None)を使います。
新規の保存先のみ、Zarr v3の無圧縮little-endian履歴とSHA256 manifestです。GS・擬ポテンシャル本体は含みません。
利用者が同じ初期バンドル・実行ファイル・dt・初期A・capacity・threadsを用意し、別Pythonプロセスから再開できます。
全履歴再生の費用は残り、物質電流を内部保存し公開値は電荷電流です。既定上限は受信モデルの最大ステップ数×48バイトと形式の管理量から算出し、max_bytesは利用者の追加予算です。全RSS上限ではありません。
破損・不一致・再生失敗では現状態を保持します。場やMPIを含む一括再開は未対応です。SNAPSHOTS.md/enを参照してください。

## 自作物質モデル

materials.check_material(initial, drive)で利用者指定A波形について独立性・再現性・有限値・端点時刻・入力所有権を検査します。
materials.bind_material(initial, key=..., identity=...)は同じnew_state/step(A)モデルをMaterialBindingへ変換します。明示TimeSpecが必要です。
MPIは1ランク/物質群の計算点分散であり、単一の直列モデルの内部並列化ではありません。identityは作者が版・全設定・初期データを識別して指定します。
宣言されたJ_matterのみjへ変換し、物理精度やfailure atomicityを証明しません。本番でチェッカーを暗黙実行しません。
native TDDFTは既存のmanaged provider、LDは既存E駆動adapterを使います。MATERIALS.md/enとexamples/custom_material.pyを参照してください。

## 自作ネイティブ物質のC ABIプレビュー1

c_material.CMaterial(library, sha256, configuration_bytes)は公開ヘッダinclude/salmonpy_material.hのABI1を読みます。
モデル/状態は独立int64 handles、固定dtの端点Ac→電荷電流j、明示TimeSpec、statusコード、失敗時原子性が契約です。
CMaterialは共通respond/check_material/bind_materialで使え、モデルclose後も既存状態が閉じるまでnativeモデル解放を延期します。
ライブラリは信頼できる実行コードで、SHA照合は障害隔離ではありません。ロード済みファイルは置換せず版ごとに別ファイル、全native呼出しはprocess-wide直列化、Pythonモデル/状態の同時呼出しは未対応です。
C++/Fortranホスト例は明示build、ソース/全object監査を保持。Cプラグインの製造解検証であり、LD/SBE永続handle化や複数native TDDFT contextを実装済みと説明しません。
C_ABI.md/enを参照。MPI collective C ABI、波動関数checkpoint、GPUは別作業です。

## SBE永続ハンドル

SBEInitialStateのモデル配列を一度ネイティブへ準備し、状態は独立したaccepted/trial密度を持ちます。
通常stepはA/J3成分だけ渡し、density/snapshot時だけ全密度をexportします。close/context managerと補助GC解放に対応。
初期モデルはdata-onlyのasdict/pickle/deepcopyを維持し、nativecacheは外部weak-key登録で移送しません。
shared velocity kernels/演算順序/既存OMP外側k分割は不変。常駐メモリはモデルnativecopyと各状態のrho2枚を必要とします。
明示サイズoverflow/allocateエラーは拒否しますが、未変更kernelの自動work領域OOMを回復可能とは説明しません。
既存固定nativebinaryを上書きせず、新規build/installで検証できます。SBE.md/enを参照。LD/TDDFTのhandle移行を実装済みとしません。

## LDの状態保持

LD.md／LD.en.mdに従います。係数と極ごとの履歴は状態ごとに保持し、stepではEと総P/Jだけを受け渡します。close／withで解放し、respondとMPI応答サービスも所有状態を解放します。固定のメモリ上限はありません。snapshot／restoreは同一モデル・固定dtでのみ利用し、状態そのもののpickleは案内しません。共有LDカーネルと物理的P/Jの規約は不変です。

TDDFTの保存用GSは状態の参照で保持し、初期モデルclose後も既存状態の再生に使えます。初期モデルclose後のnew_stateは拒否し、最後の状態解放後に保存用も削除します。N状態でN+1組、restore検証中は追加の作業用コピーが必要です。論理コピー量と実ファイルシステム上の物理占有量は区別します。

## GS/2配列アーカイブ

`salmonpy.gs_archive.save_ground_state`／`GroundStateArchive.open` は明示的に渡した周期・非スピン分極のGS配列のZarr v3入出力です。GS_ARCHIVE.md／GS_ARCHIVE.en.mdに従い、k,band,z,y,xのセル周期u、reduced逆格子k、原点0・端点なしの格子を使います。SALMONpy DFT生成のΓ点GS/1はGroundStateResult.export_archiveでGS/2へ書き出せます。SCF収束の再認定は行いません。GS/2からのTDDFTはΓ点・直列GNUランタイムに固定した監査済みGammaRestartConverterを明示指定した場合だけ利用できます。擬ポテンシャルは参照だけで自動取得しません。Fortran stream helperはチャンクの実数I/Oだけで、完全なSALMONアダプターではありません。

Γ点GS変換ではgs_native.GammaRestartConverterを使い、対象実行ファイルSHA256・実コンパイラ/Fortran_FLAGS・擬ポテンシャルhash・GFORTRAN_CONVERT_UNITを照合します。MPI／多k／別ABIの汎用変換とは説明しません。合成軌道での実ランタイム接続試験と物理的なGS検証を区別します。

GS書き出しは新しいexporter付きGammaRestartConverterを明示指定します。実数/複素数の直列Γ点のみ、原子単位のgs_eigen.dataとネイティブ占有数の一致を必須にします。密度は軌道と占有数から再構成し、混合履歴の保存と呼びません。provenanceに由来・元GS等のhash・局所チャネルを記録し、再入力時の局所チャネルも検査します。実GSを収集せず合成データで接続試験することと物理検証を区別します。

複数k点は新しいSerialRestartConverterを明示指定し、座標軸に沿う直交セル・直列・非スピン分極のみです。GammaRestartConverterはΓ限定を維持します。元GS書き出しは検証済み通常kmeshとgs_k.dataを照合し、再入力はhash付きkpoints.datをfile_kwへ渡します。密度・電子数はk点重みを含めます。波動関数/占有数と、k点テキスト・重み再規格化による電流丸め差を区別し、一般メッシュの電流比較はrtol1e-12/atol1e-14を使います。

## 本体ソースの管理（2026-11-15）

本体ソースはruntime/sourcesに同梱し、段階別完成形はruntime/overridesで管理します。外部SALMONのcloneや互換パッチ適用を案内せず、`python -m salmonpy.runtime build NEW_BUILD --jobs N`を明示的ビルド入口にします。HSEは高速PR版のみ。標準版でUSE_HSEを有効にして代替しません。元コミット・ライセンス・全ファイルのハッシュを保持し、import時はビルドしません。旧patch_sha256は由来と過去記録の互換用識別子で、新しいソースの認証ではありません。新規数値変更はsource-tree/replacementsのハッシュと対応する検証を更新します。

## 高速HSEの本体接続（開発用）

RUNTIME.md/enのnative owner経路に従います。`--explicit-routing`で準備した高速版は、
明示的な検証時に`SALMONPY_HSE_CONTEXT=1`でTDDFT本体へ独立HSE ownerを結合します。
FFT・ACE・履歴の交互実行検査と、本体電流/エネルギーの32ケース比較を区別します。
共通かつ不変のglobal設定、直列・非スピン分極・全支持の範囲です。TDDFT全体の複数状態・
異なるglobal XC/dt・MPI/Wannierの完全な独立性や、Pythonからの新規HSE APIを案内しません。
owner/結合済みsystemはコピーせず、unbind後clearで解放します。資源付き状態のbytes=-1・
status=2は計測不可です。GitHub Actionsは利用枠節約のため、新規起動せずローカル検証を優先します。

## TDDFT全体の独立状態の開発用検査

RUNTIME.md/enの「状態ごとの固定dtを持つTDDFT全体の独立状態」に従います。tddft_context.f90は
nativeの全step引数を所有し、dt以外は同じ設定・同じGSから独立初期化した2状態の交互実行を検査します。
異なるA波形、逆順、片方の破棄/再生成、別プロセス参照との配列比較、本体Acos2参照を使います。
rt_stepの返すjは電荷電流、エネルギーは毎ステップ、単位は原子単位です。
内部Fortranの開発用入口であり、既存Python TDDFTの制約を解除したとは説明しません。
XC作業配列・modulo表は各状態へ移しましたが、global入力・初期化設定の共有は残ります。異なるXC/grid/k、MPI、同時呼出しは
未対応です。保存した基本設定の変更や非有限Aは進める前に拒否しますが、全native失敗の原子性や
RSS上限は保証しません。traceは合成検査専用で、checkpointではありません。

内部Fortranのrt_open(..., delta_t=...)では状態ごとの固定dtを指定できます。有限・正の値のみ、
省略時は入力dtを取得します。異なるdtの別プロセス参照と旧global dtを変更した交互実行を比較します。
適応刻み・途中のdt変更・FDTDサブサイクル・公開Python APIの拡張とは説明しません。


XC作業配列はsystem%xc_work、周期境界表・最適化可否はstencilが所有します。
calc_xcのworkspace省略時は呼び出し内で確保します。下位の異形状PZ/配列検査を、異なる格子/XCの
TDDFT全体の対応と読み替えないでください。共通格子/k点/XCの保護は維持しています。
AVX/SVE最適化経路・MPI・GPUはこのローカル検証の対象外です。


OpenMPのステンシル/電流/力の分割フラグは状態ごとに固定します。内部rt_openのopenmp_modeは
省略時に従来入力、指定時はauto/orbital/rgridです。プロセスのスレッド数設定・MPIランク分割・
同時実行対応とは別です。異なる分割方法の交互実行試験を、異なる格子/XCや力の検証と説明しません。


内部Fortranの作成時設定はs_rt_run_configです。rt_input_configで従来入力から値を取得し、
rt_open(config=...)へ渡します。dt・最大ステップ数・3種類のOpenMP分割指定をコピーして固定します。
configと旧delta_t/openmp_modeの混在は拒否します。容量の整数検査をメモリ保証と扱わず、
格子/XCなど入力全体が独立化されたとは説明しません。上限到達後のstepは更新前に拒否します。


s_rgrid_configとinit_grid_whole_configは、実空間格子の境界・halo・添字配列のみを初期化する内部入口です。
座標配列やMPI分割表は後続処理で扱います。明示的な設定での異形状・孤立境界の検査を、
異なる格子/XCや孤立系TDDFT全体の対応と読み替えないでください。不正値は既存格子を変更せず返します。


resolve_cellはs_cell_configからセル形状と格子間隔を決める内部ルーチンです。結果は値として独立し、
不正値・体積ゼロでは変更しません。native互換層のglobal al正規化とk/pseudo/XCの依存は残るため、
セル設定全体の独立化や異なるセルのTDDFT対応済みとは案内しません。


内部Fortranのk点・radial擬ポテンシャル・XC入力はs_kpoint_config / s_pp_input_config / s_xc_configへ分離しました。
RUNTIME.md/enの「入力の所有」に従い、値設定はコピー、初期化済みLibxc/HSEハンドルはコピー禁止です。
calc_xcのpre_scf_active省略は本計算を意味します。互換ラッパーは既存SCF段階を明示的に渡します。
下位部品の異設定・共有入力変更テストを、異物質TDDFT全体・MPI・同時呼出し対応と説明しないでください。
SOC/+U専用処理・SCF/アルゴリズム・MPI/ファイル単位の共有は残り、共通grid/k/XCガードを維持します。
検査範囲はdocs/shared-input-validation-2026-11-15.mdで確認し、未検証PP形式や物理収束を保証しません。


スカラー擬ポテンシャルの格子配置はs_pp_grid_configをppgが保持し、原子数・原子種はsystemから取得します。
init_ps(config=...)の独立経路はscalar/+Uなし。作成後の設定・原子数・種対応変更は拒否し、位置変更だけ従来順序で更新します。
SOC/+U専用状態、MPI資源、タイマー、異物質TDDFT全体の共有解消済みとは案内しません。
検証はdocs/pp-grid-validation-2026-11-15.mdとtools/test_pp_grid_config.pyを参照してください。

PP格子を更新する間、読み込んだradial数値表と格子形状は不変です。表内容は毎stepで再hashせず、設定・構成のみ検査します。


内部HSEのexx_set_phase(system, pre_scf)とexx_pre_scf(system)で評価フェーズを明示所有できます。
RT初期化は本計算フェーズを所有し、共有pre-SCFフラグ変更の影響を受けません。
明示的にpre-SCFから戻す際は現軌道でexx_refresh後に交換作用を使います。
GS/DCのpre-SCF閾値・連続回数・必要回数・開始反復も各systemが所有し、集団通信と切替・再構築順は維持します。
DC通信・物質情報交換、適応交換制御、従来交換キャッシュは共有のため、全SCF制御の独立化・並行実行対応とは案内しません。
小さいGSの4条件比較とDC/MPIの8条件比較は切替・数値の回帰検査で、物理収束や複数ソルバーの同時独立実行の保証ではありません。
DC比較はMPI有効ビルド記録を渡し、tools/test_pre_scf_runtime.pyに--dc-ranks 2または4を指定します。
空間分割HSEではwrite_gs_restart_data='no'が必要です。MacではMPI CラッパーがApple Clangを選ぶ場合があり、
GNUの-fopenmpを使う構成ではCコンパイラもGNUに合わせます。MPIソケット起動の制限は計算失敗と区別します。
検査はtools/test_pre_scf_state.pyとtools/test_pre_scf_runtime.py、記録はdocs/pre-scf-validation-2026-11-15.mdを参照します。

評価フェーズの変更はbind済み独立HSE ownerの資源・履歴を無効化します。共有交換キャッシュの転用はしません。

GSの適応交換準備・交換固定・支持領域変更フラグと局在化検査はsystemに結び付いたHSE ownerへ振り分けます。
未結合なら従来の共有ownerを使います。完全なGS独立化・残る設定の分離・同時呼び出し対応とは案内しません。
内部APIはexx_set_scf_flagsと3つの問い合わせ関数、exx_check_localization(system=system)です。

通常GSの限定構成ではHSE資源をGS内部で生成・結合し、結果/checkpoint出力後に解除・解放します。
判定はexx_bind_gs_systemで、status=0は結合成功、2は対象外で従来owner維持、1は結合衝突です。
対象は構造最適化なしの通常DFT、周期・スカラーHSE06、全支持領域CPU交換、Wannier/対称性縮約なし、
整数占有（temperature<0）、各電子分割1ランクです。DC/有限温度/空間分割は対象外です。
ログEXX_GS_OWNER explicitで選択を確認できます。GS全体の独立化や並行実行対応とは説明しません。

GSの同一入力・同一プロセス2回連続実行と元コードの別プロセス実行は、
tools/test_gs_lifecycle.pyでスカラー結果と全生成restartファイルを比較します。
4条件で一致していますが、異なる入力の併存、同時実行、TDDFTへのrestart読込の検証とは区別してください。
比較後にrestart実体は削除し、ハッシュを保持します。

内部GSコンテキストsalmonpy_gs_contextは、実引数をTARGET付きで宣言して安定したアドレスを保ちます。
生存中のs_gs_contextをコピーせず、破棄前にgs_closeを呼びます。gs_runの再実行はstatus=1で拒否、正常終了は0です。
入力・並列・タイマーなどの共有は残るため、明示入力のGS API完成や並行実行対応とは案内しません。
複数状態保持の検証は同一入力のスカラーHSEに限定します。

内部GSの収束指標を含む5つのSCF設定はgs_input_configでs_gs_run_configに取り込み、gs_run(config=value)へ渡します。
max_iterationsは正かつ整数最大値未満、threshold/pre_scf_thresholdは有限・非負、pre_scf_stepsは正です。
不正値は初期化前にstatus=1、正常設定はコピーして使います。省略時は共有入力を一度取り込みます。
convergence_metricはrho_dne/norm_rho/norm_rho_dng/norm_pot/norm_pot_dng。指標と閾値の内部単位を一致させます。
GS混合はs_gs_run_config%mixingにmethod/mixrate/alpha_mb/beta_p/nmemory_p/auto_mixing/update_ratio/initial_no_mix_stepsを保持します。
初期化前にコピーし、適応係数と履歴は各状態で独立です。係数は有限、nmemory_pは1〜21、無混合回数は0以上・整数最大値未満。
methodはsimple/broyden/pulay/simple_potential/simple_dm、auto_mixingはy/n。simple_dmにはspin-orbitが必要です。
内部Fortranのみ。占有数用電子数・構造・物質入力の分離完了、DFT-MD設定の独立化や並行GS対応とは説明しません。
同じ物質で異なる混合条件の二状態を元SALMONと比較します。物理収束の保証ではありません。

明示pre-SCF設定も入力読込と同じ適用条件を再検査します。restart/MD/最適化、checkpoint指定や密度以外の収束などを
設定構造体で迂回しないでください。不正な実行設定は内部status=1、従来実行入口では非ゼロ終了になります。

## 物質IDと設定カタログ

`material_config`タスクとCONFIG.md/enに従います。実装済みは`salmonpy.materials/1`の
設定検証・不変化だけ。明示登録したPython resolverを使い、JSONからimportや計算起動をしません。
初期状態は参照IDとsource/prepareの形だけを検証し、読み込み・物理整合の検証済みとは扱いません。
構造と幅の独立アダプターは後述のstructure_gridを参照。全計算の`salmonpy.simulation/1`、MPI配置、
出力/restartへのID継承は後続です。スーパーセルは別の構造ファイルを参照する設計です。
MateriAppsを構造取得APIと呼ばず、実際の提供元と出典を記録します。数値例を物性データと扱いません。

## オフライン計算ノードへの案内

計算ノードのインターネット接続を前提にしません。未確認なら接続可否と準備・転送先を確認し、
既知の施設制約は再質問しません。RUNTIME.md/enの手順で取得と実行を分け、計算先に合う
wheelと依存、ローカルカタログとランタイムZIP、構造/スーパーセル・擬ポテンシャル・入力・
初期状態と全参照データを事前準備します。版・出典・hashを記録し、転送後の照合と実計算ノードでの
小規模確認を案内します。実行時のDB検索・URL取得、暗黙ダウンロードや自動ビルドを提案しません。
不足物を列挙して準備段階へ戻します。ログインノードの通信可否も推測しません。
オフライン一式の自動収集・依存解決、ユーザーcallbackの通信遮断は未実装です。

明示的なローカルファイル収集は`offline_inputs`タスクの`salmonpy.preparation.prepare_inputs`、
転送後の照合は`verify_inputs`です。PREPARATION.md/enに従い、登録済みresolverと
コピー先相対名→元ファイルの一覧を指定します。汎用収集器が自動で依存探索・構造解釈・格子幅解決・科学的妥当性検証を行うとは説明しません。検証器はbackendを呼ばず通信しません。
manifestのhashを別経路で保存し、出力は入力一式の外へ置きます。

## 構造ファイルと格子幅

`structure_grid`タスクとSTRUCTURE.md/enを参照。StructureInputはASE Atomsまたは
明示形式cif/extxyz・非負frame番号でローカル構造を読み、独立コピーと出典hashを保持します。
CIF対称性展開はASE読込に含まれますが、追加スーパーセル化・標準化・折り返しはしません。
grid_from_spacingは正の軸方向直交セルに限定し、幅の単位はangstrom/bohrを明示。
ceilと最低点数/明示倍数で切り上げ、実際の幅と点数を保存します。MPI数から変更しません。
ASEのdft_from_atoms/tddft_from_atomsはgridまたはgrid_spacing＋grid_spacing_unitを選択し、
両指定を拒否。直接DFT/Calculatorは従来どおり点数を渡します。格子生成を収束・資源保証と扱いません。
examples.prepare_structureが既存オフライン収集APIへの明示登録例です。自動backend登録、
擬ポテンシャル選択、実DFTの準備完了とは説明しません。混合・部分占有は代表元素で計算せず拒否します。直接ASE幅指定は採用点数のみ保持し、希望値等を保存する場合はヘルパーのレポートを記録します。

## 組み込み物質の設定接続

`builtin_materials`タスクとBUILTIN_MATERIALS.md/enに従います。builtin_material_registryを
明示作成し、resolve/load_materials→material_from_catalog→new_stateへ接続します。
LD/SBEに加えDFT設定接続とGS→TDDFTの橋渡しは下記を参照。LDはE、SBEはA/c駆動で、符号や更新順は既存と同じです。
SBEは既存の3種類のSALMONテキストをhash付きで読み、変更時は再解決を要求します。
モデル/状態にmaterial_idとdefinition_digestを保持しますが、point/run ID、MPI割当、
出力/restartの自動ID継承やinitial_statesレシピの選択は未対応。record()を結果と明示保存します。
状態作成で導入済みnativeが必要です。登録/解決/モデル作成にnative起動・通信・コンパイルはありません。

## DFT物質設定からTDDFTへ

`dft_material_config`タスクとDFT_CONFIG.md/enに従います。salmonpy.dft/dftを明示登録し、
ローカル構造・PP・XC・states・kmesh・点数/幅から既存DFTを作成します。UPF元素は原子番号と照合し、
KY元素やPPのXC由来まで自動認定したとは説明しません。固定入力と実行時SCF/threads/runtimeは分離します。
ConfiguredMaterial.run→ConfiguredGroundState.tddft_initial_state→new_stateでIDと定義hashを継承。
既存GS検証とtransverse・A/c→j、独立workerを維持します。warm startは結果ラッパーを展開して既存検査を通します。
任意GSのJSON選択・point/run ID・ディスクID自動保存・MPI接続は後続です。configured_siは短い接続検査です。

## 保存GSをIDで選ぶ

`saved_gs_selection`とINITIAL_STATES.md/enに従い、tddft_from_catalogを明示呼出します。
GS/1のsource相対パスとground-state.jsonのSHA256、DFT material_idを指定。全nativeファイルと
構造/PP/ライセンスを管理ノードで収集・転送し、verify_inputsで照合します。GS生成時と同じ
実行ファイルhashを要求し、機種間restart互換性やSCF再実行は暗黙に提供しません。
material/initial_state IDと各digestはモデル・状態に保持し、record()を結果に保存します。
既存の「initial_states未対応」は汎用factoryと任意レシピに関する制約です。この入口は
完了したSALMONpy GS/1のみ対応し、prepare、GS/2、point/run IDやMPI自動割当は対象外です。
サンプルのAはA/c、返す電流は電荷電流密度j = -J_SALMONです。収束・メモリ保証とは区別します。

## 直列の点・run対応表

`point_runs`とPOINT_RUNS.md/enを参照。PointRunに順序付きpoint_id→準備済みモデルを明示し、
一つのrun/pointに一状態だけ作成します。物質と初期状態のID/hashを保持し、save_arraysで
配列のpoint_axesを明示して保存、open_point_arraysで対応表を検査します。LDはE、結晶はA/c。
自作モデルにもrecord/new_stateで接続できますが、隠れた状態共有・実行能力は認定しません。
既存ガイドのpoint/run未対応は従来API・MPI・restart自動継承の制約です。この入口は直列出力
に限定し、時刻同期、点の自動並べ替え、checkpoint再開やMPI割当を提供しません。

## MPIの利用者ID

`mpi_point_identity`とMPI_IDENTITY.md/enに従い、ResponseSessionまたは
CoupledStepper.from_multiscaleへ同一point_record/point_mapを全ランクから明示指定します。
整数通信番号は既存のまま、利用者の文字列point_idへ明示対応します。run_idをランクごとに生成しません。
prepare/create前に全ランクの一致を検査し、自群のSolverContext.point_recordsでIDを参照できます。
session.save_arraysは集団呼出しで、引数は場leaderだけです。open_point_arraysで対応表も検証します。
metadataは宣言値の整合性検査でありproviderの実体・物理・restart互換性の認定ではありません。
従来ガイドのMPI ID未対応は省略時の自動継承やrestartについての制約で、この入口は明示的な継承です。

## 記録付きモデルのMPI binding

`recorded_material_binding`とRECORDED_MATERIALS.md/enを参照。設定済みSBE/自作端点Aモデルを
bind_recorded_materialで接続し、response_identityでIDを手で写さず点記録を生成します。
物質群は1ランク。全ランクで同じrun_id/point_mapを指定し、モデルをセッションより長く保持します。
recordは準備と状態生成で再照合しますが、虚偽のhash宣言は検出しません。LDのE変換は明示の専用経路、
TDDFTは管理ランタイム専用providerを使います。設定済みTDDFTの専用factoryは後続です。

## 設定済みTDDFT/LD専用binding

`configured_native_binding`とCONFIGURED_NATIVE.md/enを参照。bind_configured_tddftは
管理runtimeを明示しsnapshotを借用します。binding/prepareではコピーもworker起動もしません。
bind_configured_ld_endpointsは既存のA端点差分→E変換を明示的に選ぶ入口です。
response_identityでIDを継承し、物質群は各1ランク、TDDFTは独立子worker。準備モデルはセッション解放後にclose。
群内native MPIや半陰的E/H結合・資源上限・物理収束の保証へ拡張解釈しません。

## 点ID付きTDDFT checkpoint

`point_checkpoint`とPOINT_CHECKPOINTS.md/enを参照。PointStateのsave_checkpoint/restore_checkpointは
既存TDDFT履歴に厳密なrun/point/material/initial ID照合を重ねます。再開時は同じrun_idを明示。
別点/runへの分岐やMPI全体の再開は未対応。初期GS/runtimeは別途検証して用意し、archiveから自動取得しません。
復元は履歴再計算で最大2worker、max_bytesは内部archive予算でRSS保証ではありません。
LD/SBEや自作backendのディスクserializerとして案内しません。生TDDFT archiveはPointStateへ直接渡せません。

## 共通A/E/Hローカル入力

`field_input`とFIELD_INPUT.md/enを参照。FieldPacketは任意A/E/HをFieldValueで持ち、
Noneとゼロを区別。TimeSpec・単位・ゲージ・位置をFieldContractで宣言し、FieldAdapter.step_fieldsへ渡します。
AはA/c、HはBではありません。結晶はAのまま、LDはE。長さゲージはE必須ですが、native孤立系の接続完了ではありません。
既存stepへの橋はlegacy_quantityを明示し、A/c:auまたはE:auの端点値だけを渡します。
自作step_fieldsには宣言成分だけ渡し、半ステップや単位・位置を検査します。補間・微分・単位変換はしません。
MPIのA/E/H転送は後述の明示的なpacket経路を使います。磁気物理の対応はbackendごとに確認します。状態は借用、backend例外でadapterは無効です。

## 孤立系E-workerの現在の範囲

FIELD_INPUT.mdの検証用E1プロトコルは、ソースから明示的にビルドしたstandard runtimeで
PZ・固定イオン・中点伝播の双極子を比較する開発用接続です。結晶のTDDFTInitialStateにEを
渡したり、検証用driverを公開初期状態ローダーとして案内したりしないでください。
公開入口はISOLATED_TDDFT.mdのload_isolated_runtimeとIsolatedTDDFTInitialStateです。
ローカル検証記録に結び付いたバイナリに限定し、メモリ上のsnapshot/restoreを提供します。
既存配布runtimeの自動対応、電流密度への暗黙変換、孤立系ローダーのMPI応答サービスへの直結は未対応です。孤立系のsave_checkpoint/restore_checkpointは専用E/dipoleスキーマで履歴再実行を行います。GSとruntimeは別途指定し、同一容量・dt・位置・版などを照合します。周期系形式と混同しないでください。

分子双極子からFDTDへの一方向接続はDIPOLE_RADIATION.mdを参照してください。電荷双極子、
FDTDセル体積、配向・個数を明示し、最初の区間電流をゼロにする保持時間を隠さないでください。
双方向の自己場帰還や静的双極子の静電場を計算済みとは案内しないでください。

## MPI場packet
ResponseSessionのfield_contractsとstep_fieldsで明示的なA/E/Hを転送できます。
FIELD_INPUT.md/enとmpi_field_inputs.pyを参照。戻り値は電荷電流密度jです。
孤立系の双極子を自動で規格化したり、場の時空間補間を推測しないでください。
既存step(A)との混在・MPI全体のrollback・native結合driverへの自動接続は未対応です。

場の補間にはfield_sampling.SampledField/sample_packetを使えます。
利用者が選ぶ空間の凸重みとexact/linear時間方針を明示し、範囲外へ外挿しません。
成分の時計・単位を保ち、Hを推定で補いません。native結合は後述の明示された時刻配置に限定します。

NativeEHField/MultiscaleFieldのsample_fields()で1ステップ後からA/E標本を取得できます。
nativeではA=(n+1/2)dt、E=n*dt、1DではA=n*dt、E=(n-1)dtです。
Hをゼロで補わないこと。sample_fields自体は履歴を蓄積せず、結合更新もしません。

H対応native ABIではsample_magnetic_field()を全fieldランクから明示的に呼べます。
raw Hは(n+1/2)dt、保存用の時間平均Hはn*dtです。両者を混同しないでください。
旧ABIのH欠落を0で補わず、対応runtimeを明示ビルドするかHを必要としない経路を選びます。

PacketCoupledStepperは宣言済みの直接標本に限りpacket往復結合を提供します。
NativeEHではA/Hとjがoffset1/2、Eが0。1DではA/jが0、Eがlag1です。
t0/clock_offsetは0、位置・単位は厳密一致。一般補間や既存電子providerの自動変換はしません。
FIELD_INPUTの表とmpi_packet_coupling例を参照。双極子をjとみなさないでください。

既存A駆動bindingのpacket化はpacketize_A_bindingを明示的に呼びます。
field_phaseはendpoint（1D）またはyee_half_step（native E/H）を必須指定。
元の端点A/j宣言、ゼロ初期Aを要求し、電荷電流の符号を再反転しません。
絶対時刻依存の自作モデルへ機械的に適用せず、時刻対応を確認します。
response.timing_recordsに元/変換後の宣言が残ります。FIELD_INPUT.md/enと
mpi_packet_couplingの--adapt-a、--runtime/--bundleを参照してください。

packetize_A_bindingはleader_only・一物質群1ランクに限定し、collective物質backendを拒否します。

分子への電場選択はdipole_source.select_molecular_fieldで方針を明示します。
provided_driveまたはenvironment_scattered_onlyを使い、対象自身だけの真空参照と
配置・観測・源履歴の一致を宣言します。宣言は計算条件の独立検証ではありません。
DIPOLE_RADIATION.md/enとmolecular_field_selection例を参照してください。
純粋な標本選択であり、native/孤立系TDDFTの双方向結合・放射反作用の完成とは説明しません。

分子のnative/reference帰還は開発用振動子で順序検証済みですが、保持したEによる
結合仕事の差が1次で減るだけで、結合全体の保存則は未達です。振動子単体の約2次精度と
混同せず、孤立系TDDFTの双方向対応を宣言しません。新しい検証スキルの一般判定も、
SALMON固有の更新順・固定dt・資源予算やユーザー所有の物理収束判断を上書きしません。

内部GSの収束規格化はs_gs_run_config%normalization（nelec/nelec_spin）に取り込み、状態へコピーします。
非分極はnelec、分極は正のnelec_spin和を優先し、DCは状態内elec_num_totを優先します。
非DCのrho_dneの分母0以下は初期化前にstatus=1。通常は取り込んだ値を変更しません。
これは電子数・電荷・占有数を変更するPython APIではありません。物質入力の全面独立化とは案内しません。
二倍の分母で残差だけ半分にする診断と、通常条件で元SALMONとの一致確認を区別します。

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

## LD・SBE・TDDFTの混在に向けた入口

`mpi_mixed_multiscale`は、同じ1次元場にLDと合成SBEを配置する接続サンプルです。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
```
場1ランク、物質点2ランクで動かし、非ゼロの応答と独立した直列更新との一致を検査します。点ID10/20と物質キー、A/c・電荷電流j、通信時間を保存します。LDは端点A差分からEを得ます。native E/Hの半陰的LD更新と同じ方式ではありません。

4ランクにして`--runtime RUNTIME --bundle GS_BUNDLE --input rt.inp`を指定すると、管理済み直列TDDFTワーカーを第三点に接続する入口を使えます。Si入門系の実GS（PZ、16³格子、2×2×2 k点）を用いた4ランク・100ステップの三種類混在実行は、直列計算と完全一致しました。実GSと対応ランタイムを明示的に用意してください。これは管理済み直列ワーカーの検証で、TDDFT群内部のMPI並列とは別の経路です。光学物性の収束は保証しません。記録: `docs/mixed-three-material-validation-2026-11-16.md`。

最終構成は、FDTDランク内の多数のLD状態と、別MPI群のSBE/TDDFT状態の混在です。この小規模サンプルの「LD一セルに一ランク」を本番の配置方針にしないでください。実材料の収束や三次元混在対応を保証するサンプルではありません。

## native E/H内のLDと、SBE・内部MPI TDDFT

`mpi_mixed_native`は、場1ランクが二種類のLD領域を持ち、別の1ランクがSBE、別の2ランクが1点のTDDFTを担当する接続試験です。外側4ランクでTDDFTの軌道分割を内部2ランクにします。
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT
```
`--tddft-ranks 1`なら外側3ランクです。`--without-ld`はLD領域を真空にする対照条件です。BUILD_JSONは同じMPI環境で明示ビルドした`tools/build_embedded_tddft.py`の記録、SERIAL_RUNTIMEは応答再計算の比較用です。自動ビルド・ダウンロードは行いません。

LDはnative E/Hの既存の係数・補助量・更新順を使います。電子状態はA/c駆動を維持し、受け渡しは電荷電流jです。点ID10はSBE、20はTDDFT、LD領域ID1/2は場の材料格子上で管理します。電子点の周辺はLDと重複させません。半ステップずれの場時計と物質の局所時計を区別します。

この40³格子・電子2点の小さな三次元配置は接続検査用です。実用TDDFTは引き続き一次元伝播を推奨します。巨大な多数LD、分散した複数FDTDランク、全計算再開、物理収束の検証とは区別してください。ネイティブ異常は回復可能ではなく、サンプルはMPIジョブを停止します。

検証済み: 128ステップ、TDDFT内部1/2ランクの場・電流比較、直列応答再計算、LDなし対照、LD単独の従来一致、wheel導入後実行。範囲と誤差は`docs/native-mixed-material-validation-2026-11-16.md`を参照してください。

### FDTD群も並列化する

場をy方向に2分割する例（場2＋SBE1＋TDDFT1＝外側4ランク）:
```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD_JSON SERIAL_RUNTIME GS_BUNDLE NEW_OUTPUT --field-grid 1 2 1 --tddft-ranks 1 --compare REFERENCE_OUTPUT
```
REFERENCE_OUTPUTには、同じサンプルを場1ランク・TDDFT1ランク（外側3ランク）で実行した出力を指定します。`--compare`を省略すれば単独実行です。比較は同じステップ数・dt・材料設定・点順を要求し、A/c・E・jの履歴を確認します。必要な外側ランク数は`product(field-grid)+1+tddft-ranks`です。

x/y/z各方向の2分割は、128ステップで場1ランクの結果と完全一致しました。y/z分割では両LD材料が両方の場区画に分かれます。場2ランクとTDDFT内部2ランクを同時に使う5ランク構成は、今回の4コア枠では未検証です。大規模性能・メモリ上限の保証ではありません。記録: `docs/distributed-field-mixed-validation-2026-11-16.md`。

## 混在層の構成入口

配置には `salmonpy.layers.assemble_native_layers` と `examples/mixed_layers.json` を使用する。物理設定を配置JSONへ重複させず、準備済みの物質IDとモデル／MaterialBindingを渡す。Bohr単位、厚み切り上げ、全電子セルの独立点化、必要ランク数と実際の厚みを利用者へ示す。ランク予算はメモリ制限ではない。LDと電子点の禁止された近接を真空で黙って変更しない。構成関数はMPI・GS計算・ダウンロードを開始しない。自作コードも明示的なMaterialBindingで接続する。詳細と実行例はCOUPLED.md／COUPLED.en.md。

## native混在サンプルの直接保存／復元

`examples/mpi_mixed_native.py` の `--direct-save-at N` と `--direct-resume PATH` を使う。入力・nativeライブラリ・GS・分割・全ステップ数は同じに保ち、新しい出力先へ復元する。場/PML/LD/SBE/周期PZ／TB-mBJ／Libxc-PZ・PBE中点TDDFTの内部状態と積算A・直前電流・時計を直接復元する。`restored_step` と `native_steps_this_job` で過去の時間発展を再計算していないことを確認する。参照用の独立シリアル計算は検証のため別途実行される。保存は同じ環境・MPI構成の同梱サンプル用であり、汎用バックエンド保存やHSE・未検証XCの直接復元を推定しない。詳細はCOUPLED日英ガイド。

TB-mBJ（自動／固定cval）とLibxc-PBEは、場1＋SBE1＋TDDFT内部軌道MPI2の4ランクでも、直接保存63→残り65ステップの再開が連続128ステップとビット単位で一致した。分割を同じに保つこと。場2＋SBE1＋TDDFT2は5ランク必要で、4コア枠では未検証。

HSEの実行エンジンは高速PR版のみ。通常版HSEや比較用ベースラインを利用者向けエンジンに追加しない。高速HSEの直接保存／復元はACE・因子履歴・専用伝播器・native資源を含む対応が必要で、中点伝播のXC許可リストへ名前だけを追加しない。

高速PR版HSEは `native_adapters/tddft_context.f90` の内部serial Fortran入口で直接保存／復元に対応した。`rt_checkpoint(ctx,path,restore,status)` と `rt_accepted_step(ctx)` を使い、新しいプロセスで同じ元GS・入力から初期化した後に内部数値状態を読み、残りだけ進める。ACE・因子／PT-CN履歴とFFT資源の再構築を含む。`tools/test_hse_direct_checkpoint.py` は1／2スレッド、2／3 k点、4伝播設定、6／7分割の32ケースで電流・エネルギー・内部トレースが完全一致した。これはhost-MPI/FDTD結合の高速HSE対応や安定したPython/C APIを意味しない。通常HSEは追加しない。詳細はRUNTIME日英ガイド。

高速PR版HSEの外側MPI／FDTD接続は `fast_hse.FastHSEMaterialProvider` と `FastHSEResponseBackend` を使う。専用serialワーカーを明示ビルドし、各HSE点は1ランク群にする。HSE内部MPIではない。場1／2ランク＋SBE1＋HSE1と場内LDを結合し、共有保存先へ直接保存／復元する。`mpi_mixed_native --fast-hse-build ... --layers mixed_hse_layers.json` の手順はCOUPLED日英ガイド。初期電流は実測値を返し、差し引かない。ゼロ場サンプルでは1e-12 au以下を要求。入力の高速伝播器と学習条件を守り、通常HSE・Gammaの未対応学習モードへフォールバックしない。全ランクの保存検査後にHSE LOADを行い、counter/current/Aを照合する。

内部k点MPIの高速PR版HSEは `fast_hse_mpi.NativeFastHSEBackend` を明示ビルドしたホストMPIライブラリへ接続する。物質群1／2ランク、軌道／実空間分割なし、1群1点、メインスレッド・1プロセス1寿命。場1＋SBE1＋HSE2と場内LDの4ランク例は `--fast-hse-mpi-build ... --hse-ranks 2`。同じ分割・バイナリ・共有保存先でランク別直接保存／復元し、過去ステップの再計算をしない。静的初期化と独立serial参照比較は別途行う。全shardの検査後に反映するが、後段I/O失敗や群間での原子的巻戻しは保証しない。失敗後はジョブを中止する。入力、制約、計測範囲はCOUPLED日英ガイドに従う。

層構成では任意の `backend.supported_group_sizes` を検査する。高速HSEは(1,2)。対応ランク数の宣言は空でない正の整数の重複なしtuple。未宣言の自作バックエンドへ未確認の能力を推定しない。数値・再開の能力検査は別途必要。

対応ランク数の任意宣言はResponseSessionの直接接続でも物質準備前に集合検査する。宣言を全群で照合し、tupleの順序差は同一視する。Aパケットラッパーも宣言を保持する。
