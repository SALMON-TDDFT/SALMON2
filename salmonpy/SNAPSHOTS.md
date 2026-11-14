# 軽い部品の状態保存と復元

[English](SNAPSHOTS.en.md) · Topic: `snapshots` · Paired revision: `2026-11-14.1`

`state.snapshot()` は不変の `StateSnapshot` トークンを返します。
`state.restore(token)` は受理済み状態を置き換え、`None`を返します。
同一プロセス・固定dt・メモリ上の操作です。時間発展、ファイル入出力、MPI通信、
ソルバーのコンパイル・読込みは行いません。

## 保存する履歴

| 状態 | 保存する時間発展の履歴 |
|---|---|
| `LDState` | 極ごとのP/J、直前のE、ステップ数 |
| `SBEState` | 密度行列、直前のA/c、内部の物質電流、ステップ数 |
| `Transverse`の状態 | A、ステップ数 |
| `Longitudinal`の状態 | A、2時刻分の誘起A、ステップ数 |
| `ThinFilm`の状態 | A、反射A、ステップ数 |
| `Multiscale1D`の状態 | マクロ点A、境界セルを含む全格子の旧・現・新A、格子電流、ステップ数 |

スキーマ・パッケージ版・時刻・互換条件・データ整合性も記録します。
元のモデルへの参照を保持するため、不要なトークンは破棄してください。
時間発展配列を複製するので、必要メモリは状態サイズと保存数に応じて増えます。
圧縮やメモリ上限の強制はありません。

## 互換条件と独立性

元の状態、または**同じモデルオブジェクト**から作った状態に復元できます。
dt・開始時刻が一致することが必要です。LDは初期Eの設定も一致させます。
場は最大ステップ数、形状・境界、サンプル済みの外場・入射波も一致させます。
同じ係数でも別々に構築したモデルは、この初版では拒否します。
関数で駆動波形を渡す場合、新規の場状態を作る際に評価します。同じ関数・モデルでも
サンプル値が異なれば復元を拒否します。

検証と配列確保をすべて終えてから状態を置き換えます。不正な復元は `SnapshotError` を
送出し、受理済み履歴を保持します。保存内容は不変のバイト列で、復元先には新しい配列を作ります。
一方の状態のstep／restoreは、別の状態や保存済みトークンを書き換えません。
同じ状態オブジェクトへの複数スレッドからの同時呼出しには対応しません。

```python
import numpy as np
from salmonpy.ld import LorentzDrude

model = LorentzDrude(plasma_frequency=1., strength=[1.],
                     resonance=[0.], damping=[0.1])
a = model.new_state(dt=0.02)
a.step([0.1, 0., 0.])
saved = a.snapshot()
expected_p, expected_j = a.step([0.2, 0., 0.])
a.restore(saved)
p, j = a.step([0.2, 0., 0.])
np.testing.assert_array_equal(j, expected_j)

b = model.new_state(dt=0.02)
b.restore(saved)             # saved.timeで独立した分岐を作る
b.step([-0.1, 0., 0.])      # aとsavedは変化しない
```

`python -m salmonpy.examples.snapshot_restore` でLD・SBE・4種類の場を実行できます。
導入済みLD/SBEカーネルが必要ですが、SALMON本体は不要です。
例は同じ環境での再実行の完全一致を検査します。異なる環境でのビット一致や
物理的収束を保証するものではなく、科学計算の回帰比較では従来の許容差を使います。

## 対象外

引き継ぎ指示書のZarr GSアーカイブ、出力NPZ、永続的な再開ファイル、異なる版・プロセス間の形式、
C ABI、可変dtの実装ではありません。トークンをシリアライズ・編集したり、出力アーカイブを
restoreへ渡したりしません。公開メタデータ `kind`、`schema`、`version`、`step_index`、`time` は
診断用です。非公開フィールドは編集可能なチェックポイント仕様ではありません。

TDDFTワーカー、埋め込みTDDFT、ネイティブE/H FDTD、MPIセッションにはこのAPIを追加していません。
ホストが対象部品を手動で結合する場合、**同じ受理済み結合時刻で場と物質の両方**を保存・復元します。
部品をまたぐ自動トランザクションはありません。単一状態の原子的な復元をMPI全体の
巻き戻しと読み替えません。`couple()`は引き続き自身で状態を生成・解放します。
可変dt・サブサイクル・ネイティブTDDFTの再入可能化は別作業です。基本はtransverseのA駆動です。
