# 長波長の光学応答補正

`salmonpy.analysis.long_wavelength_correction` は [Otobeほか、PRB 85, 045134、付録A4](https://doi.org/10.1103/PhysRevB.85.045134) に基づく絶縁体の低周波処理です。フィット、伝導率変換、補正は観測用Fortranカーネル内で完結します。TDDFT時間発展・測定電流・Maxwellへのフィードバックは変わりません。

```python
from salmonpy.analysis import long_wavelength_correction
# observed は optical_conductivity / dielectric_tensor などの OpticalSpectrum。
# 周波数は原子単位の角周波数。以下の値は材料の推奨値ではありません。
result = long_wavelength_correction(
    observed, fit_range=(0.02, 0.035), cutoff=0.02,
    system="gapped_insulator", units="atomic_gaussian",
)
raw = result.raw_dielectric
corrected = result.dielectric
report = result.report
```

物理的な吸収端より低い信頼できるフィット帯域を指定してください。Fortranで Re ε=c0+c2*ω² を求め、cutoff未満の Im σ=-ω*(Re ε-I)/(4π) を構成します。電荷電流と exp(+iωt) の規約に従い、Iはテンソルの対角成分だけに入ります。SALMONpyの伝導率へ光速を再度掛けません。

既定の `absorption="preserve"` は論文の実誘電率の置換に合わせ、Re σとIm εを保持します。吸収側の低周波アーティファクトの修復は含みません。明示的な `absorption="analytic"` は Re σ=ω²*(d0+d2*ω²) もフィットする追加の絶縁体近似です。この場合だけモデルによるDC外挿を許します。既定ではDCの有効性を変えず、静的な実誘電率の推定値を別途報告します。

通常のDC伝導を持たないギャップ系に限定した任意適用です。金属、Hall応答、光励起自由キャリアには使わないでください。宣言だけでは適用条件を検証できません。ギャップや帯域の自動推定、実験値への固定、バンドギャップ補正、全周波数でのKramers–Kronig整合性の保証はありません。論文の1 eVや12.6を共通の既定値にはしません。

補正前後の独立した読み取り専用スペクトル、マスク、適用箇所、係数、残差、出典を返します。欠けたテンソル成分は補完しません。フィット不足は警告し元データを保持します。cutoff以上は変更しません。実装と符号は合成データで検証しており、新しい実材料TDDFTの検証を行ったという意味ではありません。
