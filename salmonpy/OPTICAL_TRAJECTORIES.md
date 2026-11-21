# 保存済み光学軌跡の再利用

追加のSALMON計算・GS生成・通信・MPI起動は行いません。新しい変換、重み済み電流の
集約、差の評価は観測用Fortranカーネル内で処理します。Pythonは配列の選択、
条件の組み合わせ、記録を担当し、元データを保持します。

## 時間長・窓関数の感度

```python
from salmonpy.analysis import optical_sensitivity, dielectric_function
result = optical_sensitivity(
    time, charge_current_x, total_E_x,
    frequencies=[0.0, 0.02, 0.04, 0.06],
    sample_counts=[2000, 4000], windows=["none", "cubic"],
    current_convention="charge", field_kind="total", units="atomic_gaussian",
)
reference = result.spectra[result.reference_index]
epsilon = dielectric_function(reference)
report = result.report
```

電流・全電場は同時刻の一成分、周波数は原子単位の角周波数です。sample_countsは
同じ軌跡の先頭からの標本数を指定し、開始時刻を変えません。
全条件を明示した共通周波数で評価し、補間・ゼロ詰め・隠れた再計算をしません。
比較の基準は最長区間と最初に指定した窓です。共通の有効周波数で相対RMS差を
報告し、基準がゼロ・比較点なしの場合は不明とします。収束合格とはしません。
読み取り専用の独立した配列と、標本数・時間長・dt・窓・名目周波数分解能
2π/(N*dt)・Nyquist・マスク・入力内容のSHA256を返します。
短い軌跡から周波数点だけを増やしても情報は増えません。Nyquist超過や弱い駆動を
マスクします。DC伝導率が有効でもDC誘電率は未定義の場合があります。

u=n/(N-1)として、none=1、hann=(1-cos(2πu))/2、cubic=1-3u²+2u³です。
片側のcubic窓は先頭のimpulseを保持します。Hannは先頭のimpulseを消すため、
駆動スペクトルがなくなる場合があります。変換はdt*sum(window*x*exp(+iωt))で、
利得補正や基線除去をせず、既存spectrum()の規格化に従います。開始絶対時刻の
位相も保持します。時間発展の変更ではありません。短縮するとパルスや減衰途中を
切ることがあり、小さな差だけでは物理精度・線形性・k点・格子・総和則を保証できません。
自動合否判定はしません。参照コードの標本数・周波数は材料の推奨値ではありません。

計算量は標本数×成分数×周波数数で、FFTとは異なります。光学帯域から調べたい
周波数を選んでください。大きな計算量には警告しますが停止しません。
物質フィットやOPTICAL_CORRECTIONの前に時間長・窓依存を確認してください。

## 重み済みk点寄与

```python
from salmonpy.outputs import read_decomposed_current
from salmonpy.analysis import k_current_contributions
raw = read_decomposed_current("Si_current_decomposed.data")
k = k_current_contributions(
    raw, dt=0.02, time0=0.0, current_convention="matter", units="atomic_gaussian",
    frequencies=[0.0, 0.02, 0.04], window="cubic",
    total_current=simultaneous_charge_current_xyz,
)
# k.currentは(sample,k,xyz)、k.totalは(sample,xyz)
# 周波数指定時のk.spectrum.valuesは(frequency,k,xyz)
```

SALMONの分解電流には占有数とk点重みがすでに掛かっています。再度掛けないでください。
Fortranは渡されたバンド・スピン成分をk点IDごとに足すだけです。
元ファイルの符号を明示し、matterならJ=-jを一度反転し、chargeなら保持します。
各k点が独立に持つ誘電率ではなく、電荷電流の寄与です。
打ち消しにより相対RMS寄与が1を超えても構いません。確率ではありません。
不足するバンド・スピン・k点は推定しません。IDだけでは完全性を証明できません。

時刻はnativeのステップ番号×runのdt+time0です。dtには観測間隔ではなく積分の
時間刻みを指定してください。任意のtotal_currentは元ファイルの符号とは無関係に、
同じ時刻の物理的な電荷電流です。利用者が明示的に時刻を揃え、補間・ずれの推定は
行いません。RMS差・相対RMS差・最大差を返し、基準ゼロなら相対差は不明とし絶対差を
残します。不一致を係数変更で修復しません。k点サンプリングの収束試験ではありません。
ファイルの単位が原子単位でない場合は先に変換してください。
ヘッダーの単位不明の場合も呼び出し側の明示的な単位宣言が必要です。

## 保存ファイルのサンプル

NPZのtime/current/Eに原子単位の時刻・電荷電流・全電場を入れます。
サンプルは一つの軸を選び、任意のk点解析を別に行います。

```sh
# 表示のみ。計算や出力先の作成は行いません。
python -m salmonpy.examples.optical_postprocess trajectory.npz NEW_OUTPUT --frequency-range 0 0.1 64
# 明示的な後処理のみ。ソルバーは起動しません。
python -m salmonpy.examples.optical_postprocess trajectory.npz NEW_OUTPUT --frequency-range 0 0.1 64 --windows none cubic --run
```

--decomposed-current FILEには--decomposition-dtと--decomposition-conventionも指定します。
NPZとの時刻合わせや異なる標本同士の全電流比較を自動で行いません。
重み済み電流の打ち消しなどの合成データで実装を検証しています。
実材料精度・GPU/Fujitsu・分散実行を新たに検証したという意味ではありません。
