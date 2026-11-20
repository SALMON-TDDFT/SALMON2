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
実軸上の極は拒否します。瞬間的な背景誘電率は `epsilon_infinity`（既定1）で指定でき、異方性・活性媒質・フィッティング・自動単位変換はありません。

## 状態の保持と解放

係数は状態の作成時に一度計算し、極ごとの P/J とともにネイティブ側に保持します。通常の `step` は3成分の E を渡して3成分ずつの P/J を受け取り、全極の配列は書き出しません。`current` と `polarization` は受理済みの総和のコピーを返します。保存・復元時だけ全極の履歴を転送します。検証済みの入力でも計算が非有限値になれば、その呼び出しを受理せず、時刻・履歴を保持します。一括計算も呼び出し全体で受理・拒否します。

直接作った状態は `with` または `close()` で解放します。`respond` とMPI応答サービスは自動で解放し、閉じ忘れた状態もGCで解放します。閉じた状態のステップと保存・復元は拒否し、最後の時刻・総和は参照できます。単一の状態を同時に複数スレッドから操作しないでください。ネイティブの管理呼び出しは直列化します。

固定のメモリ上限や状態数上限は設けません。配列寸法の整数範囲と確保失敗を検査します。保持する主要なネイティブ配列は1状態につき `192*n_poles+48` バイト程度で、管理領域・Python側・入出力配列は別です。多数の待機状態では常駐メモリが増えるため、必要な期間だけ保持してください。極ごとの係数は dt が状態ごとに異なり得るため状態ごとに持ちます。

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

検証は合成LDです。MPI1/2で背景2.25の平板とFresnel参照、分散ゼロ極限と通常誘電体の一致、異なる背景の二物質を含むnative FDTDの直接保存／復元を確認しています。実材料の光学精度の保証ではありません。
