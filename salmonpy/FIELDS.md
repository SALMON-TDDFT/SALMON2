# 場の部品（j→A）と物質との結合

[English](FIELDS.en.md) · Topic: `macro-field` · Paired revision: `2026-11-14.2`

`salmonpy.maxwell` は、光と物質の結合ループを場の側から閉じる部品です。
物質の部品（TDDFT、SBE、利用者のモデル）は A を受けて電流を返します（A→J）。
場の部品は電流を受けて次の時刻の A を返します（j→A）。`couple` は両者を交互に呼びます。

```
j_0 = 物質(A_0)  →  A_1 = 場(j_0)  →  j_1 = 物質(A_1)  →  A_2 = 場(j_1)  → …
```

## 規約

- 単位は原子単位、座標は実験室系です。A は SALMON の `Ac = A/c` で、E = −dAc/dt です。
- **インターフェースで受け渡す電流は電流（電荷電流密度 j）です。**
  物質の部品は内部で SALMON の電子物質電流 J を持ち、j = −J の関係にあります。
  `crystal` プロトコル（SBE と TDDFT の状態を含む）は j を返します。J_matter を出力すると宣言した部品だけ、`couple` が1か所で j = −J に変換します。
  場の部品を直接使う場合は、`step` に j を渡してください。
- 時刻：`state.step(j_n)` は時刻 n·dt の電流を受け取り、(n+1)·dt の A を返します（下の「時刻配置の宣言」で機械的に検査します）。
  物質側は、新しい端点の A で1ステップ進めてその時刻の電流を返します。
  これは SALMON の `time_evolution_step` と `calc_emfields` の順序と同じです。
- SALMON と同様に、**t=0 の電流は場に戻しません**。最初の `step` は外場だけで決まります。
- 設定は不変です。dt と最大ステップ数は `new_state(dt=..., steps=...)` で固定します。
  `state.config` で確定値を確認できます。
- 不正な入力や非有限の結果では状態を変えません（失敗は原子的）。

## 部品

| 部品 | SALMON での対応 | 内容 |
|---|---|---|
| `Transverse(external)` | `trans_longi='tr'` | A = A_ext。電流は場に影響しない |
| `Longitudinal(external)` | `trans_longi='lo'` | バルクの縦波配置。A_ind(n+1) = 2A_ind(n) − A_ind(n−1) + 4π j(n) dt² |
| `ThinFilm(thickness, external, epsilon=(ε1, ε2))` | `trans_longi='2d'` | 誘電体に挟まれた薄膜への垂直入射。返す A は膜内の透過場 |
| `Multiscale1D(cells, spacing, vacuum, incident, boundaries)` | 多スケール計算、`fdtddim='1d'` | Weyl ゲージの1次元伝播。セル 1..cells に物質点を置く |

`external` と `incident` には、`AcosPulse`、その列、`FieldSamples(kind='A')`（一様な場のみ。n·dt の格子上）、
または時刻配列を受けて (N,3) を返す関数を渡せます。
一様な場でNステップ進めるには、t=0からN·dtまでのN+1サンプルが必要です。

### Multiscale1D の詳細

- セル i の座標は x = i·spacing です。`vacuum=(左, 右)` は真空セルの数です。
- 入射波 A_inc(t − x/c) は左境界から入ります。初期の場は SALMON の `incident()` と同じ方法で作ります
  （左側を −3dt, −2dt, −dt で充填したあと、1回更新して t=0 の場にする）。
- 境界は `'abc'`（1次 Mur）、`'pec'`、`'periodic'` です。入射波を指定する場合、左境界は必ず `'abc'` にします。非ABC境界は `incident=None` のときだけ使えます。
- x は伝播方向なので、入射波に x 成分があるとエラーになります。
- 安定条件は c·dt ≤ spacing です。c·dt = spacing なら真空での数値分散がありません。
- `state.grid_A` は格子全体の A、`state.grid_E` は SALMON と同じ中心差分の E（時刻 (n−1)·dt）です。

## 時刻配置の宣言（`salmonpy.time/1`）

各部品は `timing` 属性で、自分の量をどの時刻に置くかを宣言します（`salmonpy.timing.TimeSpec`）。
`couple()` は結合の前に `salmonpy.timing.check_loop` で検査し、食い違いがあれば拒否します。
暗黙の半ステップずれや、黙った変換は行いません。

- 部品全体：`dt`、`t0`、`clock_offset`（部品の n 番目の時刻 = t0 + (n + clock_offset)·dt。0 か ½）
- 入力ポート：`quantity`、`offset`（0 か ½）、`sampling`（`point`／`interval_average`／`centered_difference`）、
  `lag`（出力 n+1 が使う最新入力が何ステップ前か）、`uses_t0`
- 出力ポート：`quantity`、`offset`、`sampling`、`at_t0`（`physical` か `zero`）
- 参考情報：`history`、`integrator`

検査規則：(1) dt と t0 が一致する。(2) 接続するポートの quantity・実効 offset・sampling が一致する
（変換は J_matter → j の符号反転だけ）。(3) 結合ループの lag の合計が 1 以上。
(4) 入力側が t=0 を使うなら、出力側の t=0 は `physical`。両側の宣言は `result.timing` と
出力アーカイブの `timing` に残ります。

| 部品 | 入力 | 出力 |
|---|---|---|
| 場の部品（4種） | j：offset 0、point、lag 1、t=0 は使わない | A：offset 0、point、physical |
| `SBEState`、`TDDFTState` | A：offset 0、point、lag 0 | j：offset 0、point、physical |
| E/H Yee の場（`coupled.NativeEHField`） | j：offset ½、point、lag 1、t=0 は使わない | A：offset ½、point、physical |
| E/H 結合での物質（応答サービスの状態） | A：offset 0、point、lag 0。電子系の時計（+½）上 | j：offset 0（+½） |
| 宣言のない利用者の部品 | `crystal` プロトコルの既定（上と同じ） | 同左 |

宣言のない物質の部品には `crystal` の既定（出力 j）を当て、警告を出します。E/H Yee 結合の `CoupledStepper` は、
全点の宣言を MPI で集め、電子系の時計（場の時計 + ½）に置いて `NativeEHField.timing` と照合します。結果は `stepper.timing` に残ります。

## 使い方

```python
import numpy as np
from salmonpy.fields import AcosPulse
from salmonpy.maxwell import Longitudinal, Multiscale1D, SPEED_OF_LIGHT, couple

pulse = AcosPulse(1e-3, .2, 60., polarization=(0., 0., 1.))
result = couple(material, Longitudinal(pulse), dt=.02, steps=4000)
# result.A: (steps+1, 1, 3)、result.current: 電荷電流密度 (steps+1, 1, 3)

layer = Multiscale1D(cells=3, spacing=SPEED_OF_LIGHT*.02, vacuum=(50, 50), incident=pulse)
result = couple(material, layer, dt=.02, steps=4000, output_directory='new-run')
```

`material` は `AInitialState`（`SBEInitialState` や `TDDFTInitialState` など）です。
1つ渡すと点ごとに独立な状態を作ります。点ごとに別の物質を使う場合は、点の数だけの列を渡します。
物質の状態は、成功・失敗にかかわらず `couple` が閉じます。

実行例：`python -m salmonpy.examples.field_coupling NEW_OUTPUT_DIRECTORY`

## 注意

- 陽的な結合です。速度ゲージの電流には A に比例する反磁性項があるため、
  おおよそ ω_p·dt < 2（ω_p² = 4πn_e）という安定性の上限があります。
- 座標回転（SALMON 多スケールの `rmat_ms`）、背景誘電率のある1次元媒質、
  2次元・3次元の Weyl 伝播、再開（チェックポイント）は未対応です。
- E/H の Yee 格子で結合する経路（`salmonpy.coupled`、[COUPLED.md](COUPLED.md)）は別の入口です。
  時刻配置が異なります（E は n·dt、電子系は (n+½)·dt）。

## 検証

`salmonpy/tests/test_maxwell.py` で次のことを確認しています。

- 縦波のプラズマ振動の周期が、leapfrog の離散分散式と一致する
- 薄膜の透過係数が、解析解 −iω/(−iω+γ) と一致する
- 1次元の真空伝播（c·dt = h）が解析解と相対 1×10⁻⁵ 以内で一致し、反射が残らない
- 1セル厚の1次元計算と薄膜モデルが一致する
- 同梱の SBE カーネルで `Transverse` と結合した結果が `crystal.respond` とビット一致する
- `Longitudinal` と結合した結果が SALMON の漸化式とビット一致する

SALMON 本体との直接比較は `tools/test_longitudinal_salmon.py` で行います（外部 SALMON が必要）。
同じ Si の合成初期状態を、SALMON の `trans_longi='lo'`／`'2d'`／1次元多スケール（1点）と、
SALMONpy の TDDFT ワーカー（`tr`）＋場の部品で計算します。400 ステップで A と電流の差は
相対 3×10⁻¹⁵ 以下で、SALMON の出力桁と一致しました。1次元では格子全体の A も一致します。
記録は `docs/maxwell-field-validation-2026-11-14.md` にあります。
