# CC0光学定数データとLD候補

RefractiveIndex.INFOがCC0として公開する数値データの小さな固定スナップショットです。
ネット接続なしで一覧・元のn/k・測定条件・出典・LD候補を参照できます。
読み込み時にダウンロード、フィット、コンパイル、MPI初期化、SALMON計算は行いません。
波長からの単位変換とn/kからepsilonへの変換は同梱Fortranカーネルで実行します。

```python
from salmonpy.optical_materials import list_materials, load_material

for row in list_materials():
    print(row['id'], row['recommended_fit'], row['fit_band_eV_actual'])

data = load_material('au-johnson-1972')
print(data.metadata['doi'], data.metadata['conditions'])
# 全元データ: data.omega (原子単位), data.energy_eV, data.nk, data.epsilon
model, report = data.model()  # 保存済み係数を読むだけ。再フィットしない
print(report['relative_rms_error'], report['optimizer_converged'])
state = model.new_state(dt=0.01)
state.close()
# 別候補を明示選択する場合: data.model(poles=8)
```

## 初期収録データ

| ID | 元の測定資料 | 条件の概要 | 候補極数 | epsilonの相対RMS |
|---|---|---|---:|---:|
| au-johnson-1972 | Johnson・Christy (1972) | 室温 | 6 | 4.11% |
| ag-johnson-1972 | Johnson・Christy (1972) | 室温 | 12（要再調整） | 2.74% |
| cu-johnson-1972 | Johnson・Christy (1972) | 室温 | 6 | 4.60% |
| al-mcpeak-2015 | McPeakほか (2015) | Siから剥離した蒸着Al膜、150 Å/s | 4 | 0.162% |
| si-aspnes-1983 | Aspnes・Studna (1983) | 室温、<111>、n型2.3×10^14 cm^-3 | 6 | 4.65% |

元データの波長範囲全体を保持しています。初期のフィット要求は1.5–6 eVで、
実際に含まれた測定点の端はcatalog/listのfit_band_eV_actualに記録します。
波長の丸めと離散的な測定点のため、要求端点と実際の適用帯域は一致しません。
表の誤差は指定帯域の複素epsilon全体に対する非重み付き相対RMSです。
各周波数の最大誤差や反射率・共鳴周波数の誤差を5%以内と保証するものではありません。

いずれも残差と内部時間スケールによるスクリーニング候補です。
現在の選択候補はSciPyの最適化終了条件に達しています。
データに対する残差が小さいことと、最適化の収束・物理精度・大域最適解は別です。
元の測定条件を別の薄膜、温度、結晶、表面粗さへ無条件に移用しないでください。

## 選択と警告

初期の金属・Siは全66試行（従来Fortran22件、SciPy標準初期値22件、既存係数も初期値に含めたSciPy22件）を保持します。
極数ごとに重み付き目的関数の小さいSciPy候補を選び、他の試行はalternative_trialsに保存します。
相対RMSが5%以下、recommended_dt*omega_maxが0.01以上、測定点の重み付き局所誤差が10%以下から最少極数を選びます。
安定性・物理精度の保証でも、他候補の利用禁止でもありません。
Agは全体RMS2.74%まで改善しましたが局所誤差26.7%で、needs_refinementとして既定推奨から外しています。
Al4は極端な内部時間尺度を解消してRMS0.162%です。旧候補もalternative_trialsに保持します。
Siは最少極数の6極（4.65%）を既定とし、8極（3.53%）も明示指定できます。
acceptedは全体RMSの判定だけなのでcatalog_qualityも確認してください。
明示的な極数指定、精度不足、未収束、帯域外励起は警告で知らせます。
計算可能な候補を一律に停止しません。場のCFLと時間刻みの収束は別途確認します。

n+i*kからepsilonへの変換は非磁性mu=1、時間規約exp(-i omega t)です。
P/Jは分散極のみ、瞬時背景はepsilon_infinityです。[LD.md](LD.md)を参照してください。
材料を一度読み込み、同じmodelから独立した状態を作って再利用してください。

## 元データと権利・再現

- データベース: https://github.com/polyanskiy/refractiveindex.info-database
- 固定コミット: e6d15e364f70de1cad8a96a20f75689f2fba8b21
- 公式配布条件: https://refractiveindex.info/about （CC0-1.0）
- データベース引用: M. N. Polyanskiy, Scientific Data 11, 94 (2024), https://doi.org/10.1038/s41597-023-02898-2
- 個別論文DOI: https://doi.org/10.1103/PhysRevB.6.4370 、https://doi.org/10.1021/ph5004237 、https://doi.org/10.1103/PhysRevB.27.985

optical_data/rawに元YAMLをバイト同一で保持し、CC0-1.0.txtに元のライセンス、
SOURCE_README.mdに上流README、catalog.jsonにパス・URL・コミット・SHA256・測定条件を保持します。
元論文PDF・図版・ウェブサイト本文は取り込みません。コードのライセンスとデータのCC0を区別します。
元データと派生係数も区別し、読み込み時に元ファイルと選択係数のチェックサムを確認します。

全試行の再計算は明示操作です。導入済みwheelではSDK内のtoolsから実行できます。

```bash
python tools/fit_optical_materials.py NEW_RESULTS --recorded-settings --optimizer scipy
```

--recorded-settings --optimizer scipyはAgの追加10/12極も含む各候補の保存済み設定を使います。
新しい極数・帯域を調べる場合は、この指定を外して--polesと--energy-rangeを明示します。新たな成果を同梱する場合だけ
--export-catalogを明示します。ネットからの再取得や既存の計算の実行はしません。

再フィットには任意依存 `salmonpy[fit]` が必要です。保存済みモデルの読み込みには不要です。

## 酸化物薄膜の追加

10データセットを保持します。酸化物は単結晶の一般値ではなく、次の試料条件のscalar n,k表です。

| ID | 試料 | 既定極数 | 全体相対RMS |
|---|---|---:|---:|
| sio2-rodriguez-2016 | 300 °C基板上への反応性電子ビーム蒸着膜 | 1 | 0.0430% |
| al2o3-zhukovsky-2015 | Si上20 nm ALD膜、成膜温度120 °C | 1 | 0.0332% |
| tio2-siefke-2016 | 350 nm ALD膜、広帯域測定と分散モデル | 8 | 4.47% |
| ta2o5-rodriguez-2016 | 300 °C基板上への反応性電子ビーム蒸着膜 | 2 | 2.39% |
| hfo2-franta-2015 | c-Si上112.7 nm蒸着膜、非晶質／微細多結晶、室温 | 1 | 1.29% |

要求帯域は1.5–6 eVです。Al2O3表は短波長端のため実際の上限が約5.88 eVであり、
各entryのfit_band_eV_actualが有効範囲です。表の波長全域をこのモデルで保証しません。
TiO2の8/10極は残差基準を満たしますが最適化上限に達したため警告を保持します。
多くの元表は測定を分散モデル／Kramers–Kronig解析で再構成した光学定数であり、
未加工の実験点とは区別します。吸収係数がないnのみの資料にk=0を補うことはしていません。
旧66試行と酸化物22試行を合わせて88試行を保持します。

## テンソル表現（全モード共通）

誘電率は相対誘電率epsilon[i,j]、座標は実験室xyz、時間規約exp(-i omega t)です。
`data.epsilon_tensor` は(N,3,3)、`model.permittivity_tensor(omega)` は
omega.shape+(3,3)を返します。複素吸収成分は保持され、戻り値は読み取り専用です。
scalar n,k資料で方向依存性が不明なら独立1成分のみを保存し、epsilon_ij=epsilon*delta_ij
として展開します。非対角は0で、対角3方向は同じ値です。異方性がないという実証ではなく、
方向情報のない資料に対する等方モデル仮定です。tensor_representation等に記録します。

```python
from salmonpy.optical_materials import load_material
from salmonpy.ld import LorentzDrude

data = load_material("sio2-rodriguez-2016")
measured_tensor = data.epsilon_tensor
model, report = data.model()
tensor = model.tensor_config  # LDTensorModelConfig: independent_components == 1
response = LorentzDrude.from_config(tensor)  # standalone / A-endpoint adapter
# FDTD(geometry, materials=(tensor,), dt=..., steps=...)
# LDMaterialProvider(tensor) is also accepted for MPI response groups.
```

FDTDはこの圧縮された等方テンソルを既存係数へ損失なく渡し、更新順・数式を保ちます。
A端点LDアダプターでも共通設定を使えますが、返すJは分散極のみで、瞬時背景は場側で別途扱います。
旧scalar APIは互換用に残します。全モード共通という方針は、全モードで一般異方性が完成した
という意味ではありません。native FDTDは格子主軸に沿った独立対角成分に対応しました。非対角電場更新、軸の自動回転、一般テンソルフィットは未対応です。単独／MPI物質応答のテンソル時間発展は対応しています。成分を黙って捨てる変換はしません。
