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
`1+sum(f*wp**2/(w0**2-omega**2-i*gamma*omega))` を返します。
実軸上の極は拒否します。瞬間的な背景誘電率は1に固定され、異方性・活性媒質・フィッティング・自動単位変換はありません。

## 状態の保持と解放

係数は状態の作成時に一度計算し、極ごとの P/J とともにネイティブ側に保持します。通常の `step` は3成分の E を渡して3成分ずつの P/J を受け取り、全極の配列は書き出しません。`current` と `polarization` は受理済みの総和のコピーを返します。保存・復元時だけ全極の履歴を転送します。検証済みの入力でも計算が非有限値になれば、その呼び出しを受理せず、時刻・履歴を保持します。一括計算も呼び出し全体で受理・拒否します。

直接作った状態は `with` または `close()` で解放します。`respond` とMPI応答サービスは自動で解放し、閉じ忘れた状態もGCで解放します。閉じた状態のステップと保存・復元は拒否し、最後の時刻・総和は参照できます。単一の状態を同時に複数スレッドから操作しないでください。ネイティブの管理呼び出しは直列化します。

固定のメモリ上限や状態数上限は設けません。配列寸法の整数範囲と確保失敗を検査します。保持する主要なネイティブ配列は1状態につき `192*n_poles+48` バイト程度で、管理領域・Python側・入出力配列は別です。多数の待機状態では常駐メモリが増えるため、必要な期間だけ保持してください。極ごとの係数は dt が状態ごとに異なり得るため状態ごとに持ちます。

`material` はデータだけのオブジェクトで、使用後もpickleできます。進行中の状態そのもののpickleは使わず、同一モデル・同一設定の [snapshot／restore](SNAPSHOTS.md) を使います。公開のC物質ABIにLDのプラグインを追加する変更ではありません。ライブラリのビルドは明示的に行い、import時には行いません。
