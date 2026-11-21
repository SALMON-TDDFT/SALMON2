# Lorentz–Drude 応答

`salmonpy.LorentzDrude` は等方的で受動的な極の誘電率と、電場 E に対する物理的な分極 P・電流 J を返します。原子単位・Gaussian 単位系を使い、周波数は角周波数、時間は原子時間です。eV や fs は変換してから渡します。極の強度は無次元です。

```python
import numpy as np
from salmonpy import LorentzDrude

material = LorentzDrude(plasma_frequency=0.5, strength=[0.8, 0.2],
                       resonance=[0.0, 0.3], damping=[0.01, 0.02])
E = np.zeros((101, 3))
E[:, 0] = 1e-3
result = material.respond(E, dt=0.1)
with material.new_state(dt=0.1, initial_field=E[0]) as state:
    P, J = state.step(E[1])
    saved = state.snapshot()
    state.step(E[2])
    state.restore(saved)
    P, J = state.step(E[2])
```

`respond` の初期サンプルは P=J=0、時刻0です。`step` は次の時刻の E を受け取り、保存した直前の E と固定 dt で台形則を適用します。結果は極について総和した P/J であり、FDTD の補助量 Q ではありません。状態ごとの履歴は独立です。時間刻みを細かくして収束を確認してください。

連続モデルは `dP/dt=J`、`dJ/dt=f*wp**2*E/(4*pi)-gamma*J-w0**2*P` です。
`material.permittivity(omega)` は時間規約 exp(-i omega t) の
`epsilon_infinity+sum(f*wp**2/(w0**2-omega**2-i*gamma*omega))` を返します。
実軸上の極は拒否します。瞬間的な背景誘電率は `epsilon_infinity`（既定1）で指定でき、異方性・活性媒質・自動単位変換はありません。制約付きフィッティングは下記の独立した入口を使います。

## 状態の保持と解放

係数は状態の作成時に一度計算し、極ごとの P/J とともに計算処理側に保持します。通常の `step` は3成分の E を渡して3成分ずつの P/J を受け取り、全極の配列は書き出しません。`current` と `polarization` は受理済みの総和のコピーを返します。保存・復元時だけ全極の履歴を転送します。検証済みの入力でも計算が非有限値になれば、その呼び出しを受理せず、時刻・履歴を保持します。一括計算も呼び出し全体で受理・拒否します。

直接作った状態は `with` または `close()` で解放します。`respond` とMPI応答サービスは自動で解放し、閉じ忘れた状態もGCで解放します。閉じた状態のステップと保存・復元は拒否し、最後の時刻・総和は参照できます。単一の状態を同時に複数スレッドから操作しないでください。計算処理の管理呼び出しは直列化します。

固定のメモリ上限や状態数上限は設けません。配列寸法の整数範囲と確保失敗を検査します。保持する主要な計算処理配列は1状態につき `192*n_poles+48` バイト程度で、管理領域・Python側・入出力配列は別です。多数の待機状態では常駐メモリが増えるため、必要な期間だけ保持してください。極ごとの係数は dt が状態ごとに異なり得るため状態ごとに持ちます。

`material` はデータだけのオブジェクトで、使用後もpickleできます。進行中の状態そのもののpickleは使わず、同一モデル・同一設定の [snapshot／restore](SNAPSHOTS.md) を使います。公開のC物質ABIにLDのプラグインを追加する変更ではありません。ライブラリのビルドは明示的に行い、import時には行いません。

## 瞬時背景 epsilon_infinity

```python
material = LorentzDrude(plasma_frequency=0.5, strength=[1.0],
    resonance=[0.3], damping=[0.02], epsilon_infinity=2.25)
response = material.respond(E, dt=0.1)
D = material.displacement(E, response.polarization)
u_background = material.background_energy(E, H)
```

`epsilon_infinity` は有限の正数です。`material.config.epsilon_infinity` で確認でき、従来の `background_permittivity` 属性も同じ値を返します。物質JSONのparametersにも同名で指定できます。省略時は1で、既存の分散極の更新式・順序を変更しません。

P・Jは分散極のみです。D=epsilon_infinity*E+4*pi*Pを `displacement` で求められます。`background_energy` はmu=1で (epsilon_infinity*E^2+H^2)/(8*pi) を返し、分散極の蓄積エネルギーや散逸は含みません。E・P・Hは同じ場所・時刻へ揃えた(...,3)配列を渡します。これらと連続誘電率の計算はFortran側です。背景は履歴を持つ架空の極ではなく、状態設定と保存／復元の識別に含まれます。

FDTDには既存のepsilon_emを通して伝え、半陰的E更新の分母と分子に背景を含めます。LD極自体のP/J更新は共通です。真空/PMLを含む保守的CFLを使い、背景が1より小さければさらに厳しくします。遅い物質があるだけでは真空CFLを緩和しません。

背景が1以外のimpulseは、従来の分散極電流へのkickとして時間発展を続け、生の極応答と観測データを保存します。背景1を仮定した従来の規格化スペクトルは警告して省略します。光学応答にはPlaneWaveと明示した参照・規格化を使ってください。A→JのLDアダプターも返すのは極電流です。背景は結合する場の側へ設定する必要があり、非単位背景では警告します。`dielectric_function` は全電荷電流を前提とするため、分散極Jだけから背景を含む全誘電率を得たと解釈しないでください。LDモデル全体の連続誘電率には `material.permittivity` を使えます。

検証は合成LDです。MPI1/2で背景2.25の平板とFresnel参照、分散ゼロ極限と通常誘電体の一致、異なる背景の二物質を含むPythonから直接呼ぶFDTDの直接保存／復元を確認しています。実材料の光学精度の保証ではありません。

## 制約付きフィッティングと来歴

```python
from salmonpy.ld import fit_lorentz_drude, LorentzDrude
config, report = fit_lorentz_drude(omega, epsilon, poles=2, drude=True,
    epsilon_infinity=True, units="atomic_gaussian", optimizer="scipy", rtol=0.05,
    source_file="my_optical_data.txt", reference="利用者が指定する出典")
if config is not None:
    material = LorentzDrude.from_config(config)
print(report.as_dict())
# 最適化の収束と精度判定は別: report.optimizer_converged / report.accepted
```

入力は非負で昇順の角周波数（Drudeを含む場合は正の周波数のみ）と複素epsilon、または `nk=...` の(N,2)配列（n,k）です。非磁性mu=1を仮定し、exp(-i omega t)の吸収規約n+i*kを使います。n/kからepsilonへの変換もFortranで行います。磁性媒質では利用者がepsilonを求め、複素epsilonとして渡してください。単位の自動変換はしません。`dielectric_function` の結果を第1引数へ直接渡すと、有効な周波数だけを選び、除外数・窓・時間条件などを来歴へ保持します。テンソルを自動平均せず、偏光・成分ごとに指定します。実験データは利用者が用意し、自動取得や光学定数表の同梱はしません。

`poles` はDrudeを含む総数です。`drude=True` は最初の共鳴周波数を0に固定します。`epsilon_infinity=True` は正の背景を最適化し、Falseは1固定、正の実数はその値に固定します。`initial=LDModelConfig` またはLorentzDrudeで初期推定を渡せます。返すモデルはplasma_frequency=1に規格化し、極重みをstrengthへまとめます。非負のstrength・resonance・dampingと正の背景を保つ投影付きLM法をFortran側で実行します。最適化する背景の下限は1e-9です。

目的関数は `sum(abs(delta_epsilon/max(1,abs(epsilon)))**2)`、reportの相対RMSは非重み付き複素残差のノルムをデータのノルムで割った値です。`rtol`（既定5%）は後者の判定であり、普遍的な物理精度基準ではありません。局所解・非一意性・不足データの可能性があるため、極の数・初期値・帯域を変えて利用者が比較します。`starts`（既定4）は決定的な複数初期値、max_iterations（既定300）は各初期値の上限です。不十分でも処理全体を停止せず、警告とpassiveな候補を返します。有限な候補が作れない場合だけconfig=Noneです。候補を採用するかはreportを見て決めてください。既存の時間発展や追加ジョブは起動しません。

FitReportには適用帯域、相対/重み付きRMS、選択した極数・背景条件、初期値・反復設定・ソフトウェア版、計算ルーチン実測時間、出典、参照ファイルSHA256、入力配列SHA256を保持します。source_fileは参照ファイルの記録であり、自動解析したと偽らず、実際に使った配列は別のハッシュで識別します。as_dict()は独立したJSON保存可能な辞書です。recommended_dtは活性極・減衰・プラズマの最速スケールに対する0.1倍の精度目安で、安定性上限ではありません。場の格子/CFLを知らないため、stability_marginは未確定です。

適用帯域はconfig.valid_frequency_rangeへ引き継ぎ、物質JSON/provider/FDTDでも保持します。FDTDでは帯域外キャリアや広帯域impulseを警告します。キャリアが帯域内でも有限時間パルスの裾まで検証したわけではありません。`material.check_excitation(spectrum_of_E, rtol=0.01)` で利用者が与えたEスペクトルの帯域外成分を調べられます。この割合は等重みFFTビンの振幅二乗比で、吸収エネルギー比ではありません。警告だけで計算は継続します。

サンプル: `python -m salmonpy.examples.ld_fit --poles 1 2 3`。指定した極数だけを合成データで比較し、実材料の精度や大域最適解は保証しません。

CC0の実験データと保存済み候補の一覧は [OPTICAL_DATA.md](OPTICAL_DATA.md) を参照してください。

### 任意依存のSciPy最適化

`salmonpy[fit]` を明示的に導入し、`optimizer="scipy"` を指定します。
SciPyの境界制約付き `least_squares(method="trf", x_scale="jac")` が探索を担当し、
複素残差・初期値・最終誤差・時間刻みの目安はFortranが計算します。
最適化と有限差分ヤコビアンはSciPyが担当します。
互換性のため省略時は `optimizer="native"` です。SciPyがなければ導入方法を示す
エラーを返し、勝手にインストールしたり別の最適化へ切り替えたりしません。
保存済み係数を使う通常計算にSciPyは不要です。

SciPyではmax_iterationsは各初期値の関数評価上限max_nfevです。有限差分用の
追加呼出しは別です。版・全初期値の終了理由・選択初期値・許容値・optimalityを
報告します。終了条件の達成と、大域最適解・物理精度の保証は別です。

## 全モード共通の等方テンソル設定

`model.tensor_config` / `config.tensor_config` は独立1成分のLDTensorModelConfigを返します。
epsilon_ij=epsilon*delta_ijです。LorentzDrude.from_config(tensor)、FDTDのmaterials、
LDMaterialProviderで共通に受け取ります。三次元E/Hに限定する方針は撤回しています。
permittivity_tensor(omega)は読み取り専用の複素(...,3,3)配列です。旧scalar入口も保持します。
FDTD.tensor_materialsと実行記録で表現を確認できます。Pythonから直接呼ぶFDTDは等方・独立対角・相反非対角テンソルに対応しています。単独／MPI物質応答の独立対角・非対角更新は下記の実装を参照してください。酸化物と仮定・帯域はOPTICAL_DATA.mdを参照してください。

## 独立対角・非対角の相反テンソルLD

`TensorLorentzDrude` は極ごとに対称3×3の振動子強度を受け取ります。
背景epsilon_infinityは対称正定値行列、各strengthは正半定値行列です。
個々の非対角値は負でも構いません。共鳴・減衰は極ごとの非負スカラーです。
返す応答はepsilon_inf + sum(wp**2*F_q/(w_q**2-omega**2-i*gamma_q*omega))です。
検査・周波数応答・P/J時間発展・D/背景エネルギーはFortranが計算します。

```python
from salmonpy.ld import TensorLorentzDrude, LorentzDrude

model = TensorLorentzDrude(plasma_frequency=1.,
    strength=[[[.12,-.04,.01],[-.04,.20,.02],[.01,.02,.30]]],
    resonance=[.4], damping=[.03],
    epsilon_infinity=[[2.,0.,0.],[0.,3.,0.],[0.,0.,4.]])
epsilon = model.permittivity_tensor([.1,.2])  # (2,3,3), lab xyz
with model.new_state(dt=.01) as state:
    state.step([.001,0.,0.])  # 非対角強度によりJy/Jzも応答する
    saved = state.snapshot()
    state.step([0.,.001,0.])
    state.restore(saved)
copy = LorentzDrude.from_config(model.config)  # 独立したモデル再構成
# LDMaterialProvider(model.config)でMPI物質群へ渡せる
```

scalarのpermittivity()へ自動縮約しません。成分を利用者が明示的に選んでください。
極のP/Jと受理済みEは状態ごとに保持し、snapshot/restoreは再計算せず復元します。
数値異常を拒否したstepでは元の受理状態を維持します。同じ状態の同時スレッド実行は保証しません。
A端点アダプターのJは分散極のみです。瞬時背景の行列epsilon_infは場側が扱います。

対応範囲は相反・受動・実対称な極強度です。非相反／磁気光学的な反対称応答は未対応です。
一般の自動回転APIは未提供です。実験室系の行列を指定します。テンソルフィットは下記を参照してください。
単独／MPI物質応答とPythonから直接呼ぶFDTDの非対角電場更新へ接続済みです。
非対角成分を捨てず、対になるYee補間と連立電場更新で扱います。
材料表現は空間次元に依存しません。相反・受動な対称テンソルに限定し、反対称な磁気光学応答は別の実装が必要です。
サンプルは `python -m salmonpy.examples.ld_response --anisotropic` です。


3方向TDDFT応答からのテンソル抽出とフィッティング

`analysis.dielectric_tensor` に、独立した3回の計算の電荷電流と全電場を
`(3,N,3)`（計算、同時刻サンプル、実験室xyz成分）で渡します。各周波数で
`sigma=J E^-1`、`epsilon=I+4*pi*i*sigma/omega` を求めます。
SALMONの物質電流は呼び出し側で符号を反転してください。AからEへの変換や
単位換算は自動では行いません。DC・弱励起・独立性のない入力は未評価とし、
ゼロでは埋めません。線形性、時間長・窓、物理的収束は別途確認が必要です。

`fit_lorentz_drude_tensor` は任意依存のSciPyと、利用者が指定した初期テンソル
モデルを使います。背景SPDテンソルは初期値に固定し、共通の極周波数・減衰と
PSD強度 `F=L L^T` をフィットします。行列と残差はFortran、初期行列の分解は
コンパイル済みLAPACK、最適化はSciPyです。反対称成分も含む測定9成分を残差に
残すため、一般の磁気光学応答まで再現できるという意味ではありません。
精度合格と最適化の停止条件は別々に報告します。誤差超過は警告と候補の返却です。
初期極の強度が完全にゼロだと停留し得るので、非ゼロの初期強度を指定してください。

この抽出・フィット経路は合成応答で検証しました。実材料TDDFTの物理的収束を
新たに検証したものではありません。相反非対角Pythonから直接呼ぶFDTDへの接続はFDTD_MATERIALS.mdを参照してください。
フィットした候補でも、精度・帯域・時間刻み・格子の確認は必要です。

```python
from salmonpy.analysis import dielectric_tensor
from salmonpy.ld import TensorLorentzDrude, fit_lorentz_drude_tensor

# time: (N,), charge_current and total_E: (3,N,3).
# First axis selects three independent impulse runs; last axis is lab xyz.
epsilon = dielectric_tensor(time, charge_current, total_E,
    current_convention="charge", field_kind="total",
    units="atomic_gaussian", window="none",
    metadata={"runs": ["impulse-x", "impulse-y", "impulse-z"]})
initial = TensorLorentzDrude(plasma_frequency=1.,
    strength=[[[.3,.02,0.],[.02,.4,0.],[0.,0.,.5]]],
    resonance=[.7], damping=[.08],
    epsilon_infinity=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
model, report = fit_lorentz_drude_tensor(epsilon, initial=initial,
    units="atomic_gaussian", rtol=.05)
# Inspect report before use; this does not start an FDTD calculation.
```
