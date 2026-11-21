# FDTDの独自材料と界面

[English](FDTD_MATERIALS.en.md) · Topic: `fdtd_materials` · Paired revision: `2026-11-21.1`

検証済みv2ランタイムで、`python -m salmonpy.examples.fdtd_layers MANAGED_RUNTIME NEW_DIRECTORY`
を実行します。1プロセス・1スレッドの人工的な有限サイズ二層の例であり、収束した無限多層膜の光学計算ではありません。

等方的な分散材料は `LorentzDrude(plasma_frequency=..., strength=[...], resonance=[...], damping=[...])`
で定義します。すべて原子単位で、周波数・減衰率は角周波数です。材料ごとに極数を変えられ、共鳴周波数0がDrude極です。
定数の等方媒質には `Dielectric(epsilon=..., mu=..., conductivity=...)` を使います。
別コードの係数を移す前に[LDの規約](README.md)を確認してください。係数は利用者が検証したものを指定します。
実験値への明示的なフィッティングはLD/OPTICAL_DATAを参照してください。収束保証はありません。サンプル係数は人工的なものです。

材料をタプルとして `FDTD` に渡し、形状の材料ID 1は最初の材料、2は次の材料に対応させます。0は真空です。
直方体と楕円体を組み合わせられます。IDの変更時は材料タプルと形状のIDを同時に変更してください。
形状が重なる格子では後の形状が優先されますが、材料IDの大小に物理的な優先順位はありません。
生成される材料IDのcubeを確認し、真空PMLから材料を離して、格子幅と時間刻みの収束を調べてください。
現在の上限は200材料、各LD材料100極、1000形状です。

既定の界面には、異種分散モデルを平均化しない、座標で決まる階段状の割り当てを使います。

- 電場の辺の両端が非真空なら、座標が小さい側のセルの材料モデル全体を使用します。どちらかが真空なら、従来の真空界面の係数を維持します。
- 磁場の位置では、固定した空間順序で最初の非真空セルを選びます。Hxでは `(i,j,k)`、`(i,j+1,k)`、`(i,j,k+1)`、`(i,j+1,k+1)` の順、HyとHzも対応する横断面で同様です。すべて真空なら真空係数です。
- ネイティブMPIの初期化時に材料格子の面・辺・角を通信します。E/HとLDの時間更新順は変更しません。

検証対象の単一材料・真空の結果は従来と一致します。異種材料界面では旧来のID優先処理から結果が変わります。
サブピクセル精度、回転に対する不変性、収束次数の改善を保証する方法ではありません。
観測点の従来形式のログに出る材料パラメーターは隣接する**セル**の値で、補間される電場の辺の有効係数ではありません。

複数材料の実行では旧v1ランタイムを拒否します。v2を明示的に取得・ビルドしてください。importでの自動ビルドはありません。
`salmonpy-run.json` の `interface_assignment` に方式を記録します。各 `run` は独立した計算状態を作ります。
材料係数や時間刻みを変えるときは新しい計算を作成します。

公開バッチAPIの対象はLDと定数誘電体です。任意の電流コールバック、異方性テンソル、表形式の分散、周期境界、斜め入射、
E/H FDTD内でのTDDFT/SBE材料は未対応です。Pythonから既存MPI FDTDを起動する入口は下記を参照してください。
TDDFT/SBEのA→J入口および既存A形式のマルチスケール計算は別の機能です。


## MPI FDTDの各ランク内で多数のLDを管理

FDTDは元々SALMONでMPI並列化されています。今回のAPIは既存の空間分割を指定するもので、新しいMaxwell並列化ではありません。
各FDTDランクが担当する電磁場格子とLDのP/J/Q履歴を保持し、材料IDごとに係数を選びます。
一つのランクに複数種類の材料と多数のLDセルを置け、LDセルごとの専用MPI群は不要です。
既存の半陰的なE/LD/H更新順を維持し、別のA端点駆動LDサービスはこの計算には使いません。

```python
from salmonpy.fdtd import FDTDMPI, SalmonExecutable
solver = SalmonExecutable.from_mpi_build("MPI_BUILD_DIRECTORY/build.json")
result = solver.run(calculation, "NEW_OUTPUT", mpi=FDTDMPI((2, 2, 1)), threads=1)
```

mpiexecや初期化済みmpi4pyの内部ではなく、新しい直列Pythonプロセスから実行します。
FDTDMPIの3整数が空間分割、その積がMPIランク数です。既定の `launcher=("mpiexec",)` はホスト用の引数タプルへ変更でき、シェル展開は行いません。
各方向に最低2セル/ランクが必要です。MPI起動の対象は現在macOS/Linuxです。

`python -m salmonpy.runtime build` にUSE_MPI=ONとホストのMPI C/Fortranコンパイラを指定して明示的に構築します。
RUNTIME.mdに従い、ビルド並列数はメモリ・CPU予算で選んでください。完成時のビルド記録とバイナリのハッシュでMPI対応を確認します。
対応情報のない古い記録は再ビルドが必要で、CMakeCacheの書き換えだけではMPI対応になりません。配布済み直列ランタイムはMPI版ではありません。
別の手段で監査済みの実行ファイルを用意した場合、低水準の `SalmonExecutable(...,mpi_enabled=True)` は呼び出し側による明示的な対応表明です。
バイナリの自動判定ではありません。

`salmonpy-run.json` に起動引数、ランク分割・数、スレッド数、完了・失敗・中断を記録します。
時間上限・中断時はローカルの起動プロセス群を停止しますが、遠隔ノードの後始末はMPI/スケジューラの起動機能にも依存します。
出力先はジョブごとに分けます。threadsはOpenMP設定で、別設定のBLASスレッドがあればCPU予算へ含めてください。

```sh
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_DIRECTORY --grid 2 2 1
```

例は1・2・3極の人工的なLD三材料を使い、材料セル数を表示します。メモリ使用量は未測定です。本計算の規模決定には資源見積もりと収束確認が必要です。

TDDFTなどの電子応答をE/Hへ返すには、別の[ネイティブ結合入口](COUPLED.md)を使います。上記のバッチAPIの材料リストへTDDFTを渡す方式ではありません。

## 横方向周期のスラブ

`Boundaries(('pml','pml','periodic','periodic','periodic','periodic'))`でx方向をPML、y・z方向を周期境界にできます。周期指定は同じ軸の両面へ設定し、周期位相0で接続します。座標と材料IDプレビューは従来の中心原点を維持します。Fortranの場用通信が各軸の境界を使い、電子状態計算の周期設定を変更しません。PMLに接する材料への既存検査も維持します。旧ランタイムではなく、この軸別境界に対応したランタイムを用意してください。

`python -m salmonpy.examples.fdtd_periodic_slab RUNTIME NEW_OUTPUT`は誘電体、`--ld`付きはLDスラブです。Fresnelの複素透過率との比較は明示的な`tools/test_fdtd_periodic.py AUDITED_MPI_BUILD_JSON NEW_RESULTS`で行います。計算はFortran、フーリエ和と解析式は検証用Cで評価します。256×8×8、1／2ランクの検証で、3周波数の相対複素誤差は誘電体で最大約0.0554%、LDで約0.0375%でした。横方向の観測点と1／2ランクの結果も一致しています。これは指定した合成スラブの検証であり、任意材料や全条件の精度保証ではありません。現在の共通格子初期化には細いMPI分割での最小領域制約が残るため、この例は横8格子を使います。

背景誘電率は `LorentzDrude(..., epsilon_infinity=2.25)` で指定できます。既定値1の結果は維持します。非単位背景でもimpulseの時間発展は継続し、旧規格化スペクトルだけを警告して省略します。PlaneWaveの背景を含む伝播、CFL、P/JとDの区別は [LD.md](LD.md) を参照してください。

フィットしたLDのvalid_frequency_rangeも保持します。帯域外のキャリアや広帯域impulseは警告だけで継続します。パルス全体の帯域確認はLD.check_excitationへ利用者が用意したEスペクトルを渡してください。

## 全モード共通の等方テンソル設定

`model.tensor_config` / `config.tensor_config` は独立1成分のLDTensorModelConfigを返します。
epsilon_ij=epsilon*delta_ijです。LorentzDrude.from_config(tensor)、FDTDのmaterials、
LDMaterialProviderで共通に受け取ります。三次元E/Hに限定する方針は撤回しています。
permittivity_tensor(omega)は読み取り専用の複素(...,3,3)配列です。旧scalar入口も保持します。
FDTD.tensor_materialsと実行記録で表現を確認できます。native E/H FDTDは等方・独立対角・相反非対角テンソルに対応しています。単独／MPI物質応答の独立対角・非対角更新は下記の実装を参照してください。酸化物と仮定・帯域はOPTICAL_DATA.mdを参照してください。

## 対角テンソルFDTDと波長板

格子xyzに主軸が沿う `TensorLorentzDrude` をFDTDのmaterialsへ渡せます。
各軸の背景・極重みから、EとLD電流の半陰的係数をFortranが作ります。
P/J/Qの更新順とYee格子の物質界面割当は従来どおりです。
相反非対角の更新も下記のとおり接続済みです。非対角成分を黙って0に置き換えません。

対角入力は `media_type='ld-diagonal'` と `salmonpy-ld-diagonal.dat` の組です。
`input_text()` だけを保存して終えず、`input_assets()` も保存してください。
通常の `SalmonExecutable.run` と `NativeEHField` は資産を自動配置し、来歴にはSHA256を記録します。
資産がない旧ランタイムへ新入力を渡さず、明示ビルドまたは対応配布を利用します。
scalarのepsilon_em/f_ldは先頭軸の互換情報で、完全な対角値はこの資産にあります。
元の等方入力は資産を持たず、従来と同じ入力・演算経路です。

```sh
# 既定はプレビューのみ。管理ノードでこの確認後に計算が自動開始されることはない。
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --plate quarter
# 計算資源を確保した環境で明示的に開始
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --plate quarter --run
# 半波長板、または空間・時間の細分化
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON OTHER_DIRECTORY --plate half --run
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON REFINED_DIRECTORY --spacing 4 --run
```

サンプルはn_x=1.5、n_y=2の人工的な無損失板です。実物質の光学精度検証ではありません。
`analysis.waveplate_diagnostic` は、全パルスと反射の減衰が入ったE時系列から
指定角周波数の複素振幅と正規化StokesをFortranで計算します。
比較先は真空中の平行平板・垂直入射・等振幅同位相の入射で、界面Fresnel反射と
多重反射を含みます。単なる幾何学的位相差との比較ではありません。
S3=2 Im(Ex conj(Ey))/S0、フーリエ規約exp(+i omega t)を明示します。
一つの周波数での検査であり、パルス全帯域や実物質の収束保証ではありません。
非単位背景／対角LDのimpulseでは、伝播とraw出力を続け、従来のscalar光学規格化は省略します。

phasorは台形則によるフーリエ積分で、単位はE(au)×時間(au)です。パルスのピーク振幅ではなく、偏光の評価には規格化した比を使います。

## 相反非対角native E/H FDTD

実験室系のSPD背景行列epsilon_infinityと、極ごとのPSD強度行列を持つ
TensorLorentzDrudeを渡します。共鳴・減衰は極の行列成分間で共通のスカラーです。
反対称な磁気光学応答は対象外です。native E/Hバックエンドの一次元相当・二次元相当・
三次元の構造で同じ材料表現を使います。別のmaxwell.pyのAのみの場へ異方性背景を
追加する機能ではありません。非ゼロBloch位相は新しい演算子では未対応です
（Python FDTDにはその指定口はありません）。

完全テンソルはmedia_type='ld-tensor'とsalmonpy-ld-tensor.datを使い、バッチと
NativeEHFieldの両方で自動配置します。対応ランタイムと資産SHA256を確認します。
独立な6成分を保持し、位置の異なるYee電場を、対になる辺→節点の平均とその転置平均で
結びます。物質が所有する辺へ射影し、別物質の極履歴は混ぜません。
対角前処理付き共役勾配法でEを解き、同じ極演算子をJに加えてから従来のP・Qを更新します。
数値処理はFortranです。等方・対角入力は従来の演算経路を保ちます。

FDTD(..., tensor_rtol=1e-10, tensor_maxiter=200)で連立更新の許容誤差と回数を
指定します。有限値のまま回数に達した場合は警告して継続します。残差を確認して調整してください。
非有限演算や非正の曲率の場合は試行を拒否し、nativeホストでは元の状態へ戻して再試行できます。
PEC固定電場は作用範囲から除きます。完全テンソルの直接チェックポイントは材料行列と
更新条件を照合し、不一致は反映前に拒否します。従来の等方・対角の保存形式は維持します。
場チェックポイントのMPI再分割を追加したものではありません。

```sh
# プレビューのみ。ソルバーも出力先も開始しない
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --angle-degrees 30
# 資源を確保した環境で明示実行。入射の主軸2成分は等振幅・同位相
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --angle-degrees 30 --run
```

waveplate_diagnostic(..., angle=pi/6)は実験室系Eを受け、phasorとStokesを実験室xyで
返します。位相差の報告は主軸成分間です。angleはラジアン、サンプルのCLIは度です。
パルスと反射の減衰を含む時間長が必要です。人工的なFresnel平行平板の検証であり、
実物質の精度保証ではありません。

ローカル回帰検査（CIでの必須検査にはしません）:
```sh
mpiexec -n 2 python tools/test_eh_deep_recovery_mpi.py BUILD_JSON NEW_OUTPUT --baseline --resume-midpoint --tensor
```
既存のnative復元検査を共用します。ソースの作業コピーと、確保した計算資源が必要です。


## LD専用：中心周波数に合わせた界面平滑化

界面係数の準備を任意で有効にします。

```python
from salmonpy.subpixel import SubpixelSmoothing
from dataclasses import replace
smoothed = replace(calculation, smoothing=SubpixelSmoothing())
# PlaneWaveではsource.pulse.omegaを参照角周波数として使います。
# Impulseには中心周波数が無いためSubpixelSmoothing(reference_omega=...)を指定します。
```

対象は明示的なLD材料を含むnative E/H FDTDです。一次元的・二次元的・三次元の
配置で利用できます。SBE・TDDFTや利用者のA→j状態には適用しません。
少なくとも片側がLDである界面だけを平滑化し、相手は通常の誘電体でも構いません。
LDが無い場合や磁性・導電率付き背景が未対応の場合は警告し、元の係数を保持します。
native応答点の近傍は界面平滑化から除外します。電子密度や電流の規格化、
Maxwellとの結合契約を変更する処理ではありません。

専用Fortranモジュールが直方体・楕円体の局所平面の占有率と法線を求め、
入力の中心角周波数に対応する受動的テンソル係数を初期化します。
複素界面場の補正行列を固定し、各実背景・極強度を変換します。
元の極の共鳴周波数と減衰は保ち、参照周波数では複素実効誘電率に一致します。
法線方向の調和平均・接線方向の算術平均も含みます。
複素定数誘電率や負の極強度で代用しません。差分curl、ADE漸化式、E/J/P/Q/Hの
更新順は変更せず、既存のテンソル演算部が空間依存係数を扱います。

分散媒質に対する一つの周波数での近似です。Meepと同一の実装、全帯域での
二次収束を保証するものではありません。広帯域応答と格子・時間刻みを確認してください。
角・異なる複数界面・参照周波数で特異な混合は、該当箇所で警告し階段状の係数を
保持します。材料IDの大小とは無関係にLD同士の界面も扱います。
形状が重なる場合の後勝ち規則は同じです。真空PML・PECを維持し、材料はPMLから離します。

salmonpy-subpixel.dat（版2）はrun()/NativeEHFieldが配置し、モデル・保存データの
識別情報に含めます。バッチ入力を自分で書き出す場合はinput_text()だけでなく
input_assets()も出力してください。対応する監査済みランタイムが必要です。
不足・破損した係数ファイルは黙って無視しません。
salmonpy-subpixel-summary.jsonに参照周波数、全体の極スロット数、除外した応答点数、
平滑化・局所代替の標本数、整合誤差、背景からの保守的な時間刻み目安を記録します。
時間刻み目安の超過は警告だけで継続します。実際の安定性は別途確認してください。
この説明ファイルの書き込みが失敗しても、警告して計算を継続します。

メモリは各材料の有効極スロット総数×担当格子数に比例し、セルごとのテンソルを
保持します。同じ共鳴周波数・減衰を持つ極も現状は統合しません。
多数のLDモデルではnativeの配列容量の見積もりを確認してください。
直接保存は係数も含み、不一致モデルの復元を拒否します。
この機能によってFDTDのMPI再分割再開に対応するわけではありません。

```sh
# 表示のみ。計算も出力先の作成も行いません。
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --subpixel --angle-degrees 30
# 割り当て済みの計算資源内で明示的に実行
python -m salmonpy.examples.fdtd_waveplate BUILD_JSON NEW_DIRECTORY --subpixel --angle-degrees 30 --run
```

このサンプルは極強度ゼロのLD背景を用い、形状と偏光を物質フィットの精度から
分けて確認します。分散のある実験には実際の分散LDモデルを設定してください。
既定のsmoothing=Noneでは従来の階段状配置を使用します。
