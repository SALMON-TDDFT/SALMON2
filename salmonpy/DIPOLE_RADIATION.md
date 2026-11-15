# 分子双極子からFDTDへの一方向放射源

分子の物理的な電荷双極子pを実験室座標、e·bohr単位で渡します。
DipoleMapは双極子総量を保存して明示したマクロセルへ分配します。

P_i = sum_m weights[i,m] p_m / cell_volumes[i]

j_i = sum_m weights[i,m] (p_m(new)-p_m(old)) / (dt cell_volumes[i])

体積はFDTDセルのbohr³単位の値で、TDDFTの真空を含む計算箱の体積ではありません。
重みは非負で、分子ごとの和を1とします。自動規格化しません。
分子IDは重複のない文字列、場の点IDは非負整数で、指定した順序を保持します。
複数分子を同じセルへ重ねられます。座標回転と実際の分子数の反映は呼出側で明示してください。
均質な分子集団ならP=n*<p>ですが、数密度・配向平均を計算箱から推測しません。
時刻は原子単位、jはe/(bohr²·原子時間)です。SIへの暗黙変換はありません。

```python
from salmonpy.dipole_source import DipoleMap
mapping = DipoleMap(('molecule-1',), (10,20), [[.25],[.75]], [8.,8.])
j = mapping.current([[0.,0.,0.]], [[0.,0.,.001]], dt=.003)
```

入力はコピーされ、返す配列は読取専用です。一定の双極子は電流を作りません。
その静電場を初期化する機能ではありません。電気双極子近似で、分子の空間電荷分布、
磁気・高次多極子放射、局所場・自己場の帰還は扱いません。

## native E/Hへの接続

`drive_native_dipoles(field, mapping, dipoles, dt=...)`は初期状態のNativeEHFieldを借用します。
場グループの全ランクから同じ引数で呼び、他のグループからは呼びません。
双極子は(samples,molecules,3)、時刻は0からn*dtです。
読取専用のE/j履歴、それぞれの時刻、分子ID・点IDを返します。
Eはt=0を含み、jは区間平均を中点時刻で表します。
マクロセルからYee格子への分配とFDTD更新は既存native処理のままです。
Python側でもう一度Yee格子へ分配しないでください。

nativeの初期化条件により最初の区間電流はゼロが必要です。
すぐ変化する既存波形を使う場合は、呼出側で初期双極子を1回繰り返します。
これによるdtの遅延をAPI内に隠しません。

```python
import numpy as np
from salmonpy.dipole_source import drive_native_dipoles
# p: n*field.dtの電荷双極子、形状(samples,molecules,3)。
# field: 点ID・セル体積が一致する初期化済みNativeEHField。
delayed = np.concatenate((p[:1], p), axis=0)
result = drive_native_dipoles(field, mapping, delayed, dt=field.dt)
```

場はlen(delayed)-1ステップ以上の容量を持ち、点の順序・体積が一致する必要があります。
時間刻みも一致させ、補間は明示的に行います。全波形を最初の更新前に検査し、ランク間で照合します。
backend異常時に計算全体を巻き戻す契約ではありません。場のcloseは呼出側が担当します。
NativeEHFieldの明示的build record、ホストMPI context、境界から離れた真空近傍、
平面波入力という制約は引き継ぎます。分子の放射だけを見る場合はパルス振幅をゼロにします。
入射場を残す場合は出力に両方の寄与が含まれます。境界・CFL条件はFDTD側で扱います。

一方向接続なので放射したEをTDDFTへ戻しません。放射反作用、自己場と局所場の区別は未対応です。
実験的孤立系E-workerでは、応答の最後の3成分である全電荷双極子を使います。
周期系TDDFTの物質電流の符号規則を転用しないでください。FIELD_INPUT.mdも参照してください。

## サンプルと検証

```sh
python -m salmonpy.examples.dipole_radiation
```

wheel内の軽量例は電流生成と総量保存だけを確認します。
ソース内のtools/test_dipole_radiation_mpi.pyは明示したembedded library recordと新しい出力先で
native FDTDを実行します（合計2または3ランク）。tools/test_isolated_worker.pyが出力した
NPZを追加引数にすると分子の実測双極子も使えます。
検証記録はdocs/dipole-radiation-validation-2026-11-15.mdです。
短いC₂H₂波形は接続試験で、収束した分子発光スペクトルではありません。

## 細分化検証の現状

位置を固定した解析的な点双極子について、遅延電場の全成分と電流源の仕事を比較しました。
dt=.0015原子時間でdxを2→1→.75 bohrへ細かくすると、電場の相対誤差は46.5%→8.0%→4.3%、
源の仕事と連続体の放射エネルギーとの差は15.8%→2.5%→1.1%へ減りました。
閉曲面を通るエネルギー流の実測や、本番精度の収束保証ではありません。
時間刻みだけを半分にしても連続体解との一致は単調に改善しませんでした。
源セルのEは格子依存が強く、そのまま分子へ帰還する方法は未検証です。
再現ツールはソース内のtools/test_dipole_refinement_mpi.pyとtools/summarize_dipole_refinement.py、
条件と制限はdocs/radiation-refinement-validation-2026-11-15.mdにあります。
自己場差引きや真空の放射反作用モデルを自動的に有効化することはありません。

公開APIのisolated_responseサンプルが作った出力も、そのまま読み込めます。

```sh
python -m salmonpy.examples.dipole_radiation --archive ISOLATED_OUTPUT --cell-volume 8000
```

8000は使用するFDTDセルの実際の体積（bohr³）に置き換えてください。
単位・等間隔の時刻・分極電流生成と保存則を確認する例で、FDTD自体は起動しません。
native伝播までの実行は上記のMPI検証スクリプトに分けています。

## 分子へ渡す電場の明示選択

`select_molecular_field`は計算済みの電場`FieldValue`を選択する処理です。
双方向ソルバーは起動しません。方針を必ず指定します。

- `provided_drive`：利用者が選択済みの電場を保持し、provenanceの`meaning`に
  物理的な意味を記録します。参照場の併用は拒否します。
- `environment_scattered_only`：環境中の全電場から、対象分子自身だけによる
  真空参照場を引きます。外場・他分子からの場は残します。真空の放射反作用も
  除くため、完全な放射ダイナミクスや自然放出の計算とは呼びません。

両標本はE:au・bohr・point samplingとし、時刻・位置を厳密に一致させます。
参照のprovenanceには`source_ids=[molecule_id]`、`medium='vacuum'`、
`incident=False`が必要です。両provenanceの`footprint_id`、`sample_map_id`、
`own_source_history_id`も一致させます。これらは配置・観測・自己源履歴が同じという
利用者の宣言であり、別々のMaxwell計算の条件を検証済みという証明ではありません。
不足や不一致をゼロで補ったり、全電場へ暗黙に戻したりはしません。

戻り値`MolecularDrive`の`.field`と`.metadata`は不変です。`.record()`は入力・参照・
選択後のE、位置・時刻・方針・分子ID・出典を持つ独立したJSON互換データを返します。
補間・単位換算はしません。足し合わせた人工電場で動作を確認する例は次です。

```sh
python -m salmonpy.examples.molecular_field_selection
```

native場と参照場の因果的な更新、振動子による仕事の収支検証、孤立系TDDFTへの
実際の帰還は後続です。この選択処理だけで双方向結合が完成したとは扱いません。

### 開発用の帰還順序検証（公開TDDFTループではありません）

ソース側の`tools/test_molecular_feedback_mpi.py`は、独立した二つのnative真空場と
調和振動子で帰還を検証します。最初の区間を明示的に保持し、E(n)から次区間の双極子と
電流を求め、場のn+1ステップへ渡します。合計4ランク・各1スレッドです。
照射だけの別計算と比較して自己場の差し引きを確認します。この振動子の中点積分を、
孤立系TDDFTの端点E駆動と同じ更新法だと扱わないでください。

振動子内部の仕事の恒等式は満たしますが、電場を保持する結合には1次の時間誤差が
残ります。短いパルスの検証でdtを.003、.0015、.00075 auと下げると、仕事の差の絶対値は
3.60e-11、1.75e-11、8.63e-12 auとなりました。吸収仕事に比べて無視できず、
結合全体のエネルギー保存や実用的な分子帰還が確立したとは言えません。
合格した項目との区別はソース側の検証記録に残しています。
