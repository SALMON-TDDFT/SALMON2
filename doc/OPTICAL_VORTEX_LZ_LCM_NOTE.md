# Optical Vortex / Angular Momentum / LCM Note

対象 tree:
- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2`

対象:
- `single_scale_maxwell_tddft`
- optical vortex 入射
- 電子角運動量出力
- LCM 出力

## 1. 光渦の設定

主な入力変数:
- `yn_optical_vortex`
- `optical_vortex_charge`
- `optical_vortex_polarization`
- `optical_vortex_radius`
- `optical_vortex_center_x`
- `optical_vortex_center_y`
- `omega1`
- `tw1`
- `ae_shape1`

既定値:
- `yn_optical_vortex = 'n'`
- `optical_vortex_charge = 0`
- `optical_vortex_polarization = 'linear_x'`
- `optical_vortex_radius = -1d0`
- `optical_vortex_center_x = -1d30`
- `optical_vortex_center_y = -1d30`

制約:
- `yn_optical_vortex='y'` は `method_singlescale='3d'` 必須
- `optical_vortex_radius > 0`
- `omega1 > 0`
- `tw1 > 0`
- `optical_vortex_polarization` は `linear_x`, `linear_y`, `left_circular`, `right_circular`
- `ae_shape1` は `Acos2`, `Acos3`, `Acos4`, `Acos6`, `Acos8`

実装上の注意:
- 渦中心が未指定なら `x,y` ともセル中央が使われる
- 偏光の現在の規約は
  - `left_circular = (1, -i)/sqrt(2)`
  - `right_circular = (1, +i)/sqrt(2)`
- 位相は `omega1 * tt + 2*pi*phi_cep1 + charge * phi_xy`
- 半径外では入射場は 0

実務上のコツ:
- 円偏光と OAM の向きの対応は、名前より実際の `*_rt_micro.data` の `Lz_inc_*` で確認する
- `optical_vortex_center_x/y` は、角運動量と LCM の 2D マップの見え方に直結する
- 中心をずらすと星形や花弁状の非対称が出やすい

## 2. 角運動量関係パラメータ

主な入力変数:
- `yn_out_lz_rt`
- `out_lz_rt_step`
- `optical_vortex_center_x`
- `optical_vortex_center_y`
- `yn_spinorbit`

既定値:
- `yn_out_lz_rt = 'n'`
- `out_lz_rt_step = 100`

出力:
- `.../<sysname>_lz_xy/<sysname>_lz_xy_000010.data`
- `.../<sysname>_lz_total.data`

中身:
- 2D マップは `z` 方向積分済み
- total 側は
  - `local_lz_total`
  - `local_sz_total`
  - `local_jz_total`
  - `dlocal_jz_total_dt`

注意:
- `dlocal_jz_total_dt` は厳密微分ではなく有限差分
- `linear_x` の平面波では光の `Lz` は基本的に 0
- 電子角運動量の 2D マップは渦中心基準の座標で出る

実務上のコツ:
- `out_lz_rt_step` はまず `10` か `20` 程度にして挙動を見る
- 左右円偏光比較では、符号だけでなく絶対値も確認する
- 格子対称性が強い系では、見た目が円形でも 6 回対称や星型パターンが出る

## 3. LCM 関係パラメータ

主な入力変数:
- `yn_out_lcm_rt`
- `out_lcm_rt_step`

既定値:
- `yn_out_lcm_rt = 'n'`
- `out_lcm_rt_step = 100`

出力:
- `.../<sysname>_lcm_xy/<sysname>_lcm_xy_000010.data`

中身:
- `x`
- `y`
- `local_chern_marker_zint`

実装上の性格:
- occupied subspace を sharp occupation で作って Löwdin 直交化している
- 金属や smeared occupation 向きではない
- online RT 中の LCM は重い

現状の性能支配:
- `S1/S2` の GEMM
- `inv(S1), inv(S2)` の Löwdin 部分
- `ZT1/ZT2` はかなり改善済み

実務上のコツ:
- まず `out_lcm_rt_step` を粗くして、RT を止めない設定で様子を見る
- 本番サイズでは `yn_scalapack='y'` を優先する
- LCM が通るかの確認では、小さい系で correctness、大きい系で scaling を分ける

## 4. 計算設定のコツ

### 並列

LCM は両方効く:
- 空間並列: `ngrid_local` を減らす
- 軌道並列: `nocc_local` を減らす

優先:
- まず空間並列
- 次に軌道並列
- 大規模では両方必要

### ScaLAPACK

使う条件:
- build 時に ScaLAPACK 有効
- 実行時に `yn_scalapack='y'`

現状の方針:
- fallback の full 経路は残す
- `yn_scalapack='y'` のときだけ LCM Löwdin を reduced-memory ScaLAPACK 経路へ分岐
- `yn_scalapack='y'` と `yn_eigenexa='y'` が両方立っていても、現状の LCM は ScaLAPACK 優先

### 運用

まず見るべき出力:
- `*_rt_micro.data`
- `*_lz_total.data`
- `*_lz_xy/*.data`
- `*_lcm_xy/*.data`

切り分けの順:
1. 光渦の入射が `*_rt_micro.data` に出ているか
2. 角運動量 2D/total が出ているか
3. LCM が step ごとに完走するか
4. その後に時間チューニング

### パターン解釈

- hBN のように外形が円形でも、格子対称性で星型や花弁状は十分ありうる
- ただし artifact 切り分けのために
  - center
  - grid
  - `dt`
  - 偏光
  - charge
  の依存性を見るべき

## 5. 参照コード

- optical vortex:
  [optical_vortex_field.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/maxwell/optical_vortex_field.f90)
- input defaults / validation:
  [inputoutput.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/inputoutput.f90)
- electronic angular momentum:
  [rt_angular_momentum.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_angular_momentum.f90)
- LCM:
  [rt_local_chern_marker.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_local_chern_marker.f90)
- LCM 2D output:
  [main_tddft.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/main_tddft.f90)
