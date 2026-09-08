# DG Density Partition Weight Design

## Goal

`density_inv_weight_local = 1 / overlap_count` による単純平均をやめて、DC の buffer/supported-domain の考え方に沿った滑らかな partition-of-unity 重みで fragment 密度を合成する。

## Problem

現状の DG density 再構築では、同じ global grid 点に複数 fragment の寄与が重なると、その owner rank 上で `1/count` を掛けて平均している。

この処理は:

- fragment 中心と buffer 端を区別しない
- DC の「大きめの box で解いた上で中心を信頼する」という考え方を反映しない
- 各 fragment の寄与強度の違いを無視する
- 電荷保存を厳密に保証しない

その結果、再構築後の総電荷が `nelec` から数 % ずれる原因候補になっている。

## Chosen Approach

各 fragment に対して、buffer 端で滑らかに 0 へ落ちる `cosine taper` 型の未正規化重み `w_f(r)` を定義し、owner rank 上で各点ごとに

`W(r) = sum_f w_f(r)`

を作って、密度寄与を

`rho_f(r) * w_f(r) / W(r)`

で合成する。

## Why Cosine Taper

- `1/count` より DC 的
- `erf/erfc` より実装が軽い
- 端で 0 へ滑らかに落ちる
- 各方向の 1D 関数の積で 3D 重みを構成できる

## Weight Definition

各 fragment の local index `(ix, iy, iz)` に対して、各方向の 1D taper を作る。

- 中心領域では `1`
- buffer 領域では cosine で `1 -> 0`
- domain 外では `0`

3D 重みは

`w_f(ix,iy,iz) = wx(ix) * wy(iy) * wz(iz)`

とする。

buffer 幅は既存の fragment geometry から決まる量を使う。最初の実装では runtime の `nxyz_buffer` と fragment domain から構成し、追加の入力パラメータは導入しない。

## Normalization Strategy

owner rank 上で、同一点に寄与する全 fragment について未正規化重みの総和

`W(r) = sum_f w_f(r)`

を構築する。

実際の density 再構築では `density_inv_weight_local` の代わりに

`density_partition_weight_local = w_f(r) / W(r)`

を使う。

これにより各点で

`sum_f density_partition_weight_local(f,r) = 1`

となる。

## Data Model Changes

`s_dg_fragment_rt` に以下を追加する。

- fragment-local grid 上の未正規化 taper weight map
- owner rank 局所 grid 上の partition normalization sum
- 必要なら owner-local final normalized weight

既存の

- `density_weight_local`
- `density_inv_weight_local`

は削除対象または置換対象にする。

## Code Path Changes

### Map Build

`build_density_grid_owner_maps` で:

- 各 fragment 点の global owner map を作る
- 各 fragment 点の未正規化 taper weight を計算する
- owner rank 上で `W(r)` を集計する
- normalized partition weight を構成する

### Density Reconstruction

`calculate_density_from_fragments` で `rho_contrib` に掛ける重みを:

- 旧: `density_inv_weight_local(ixg,iyg,izg)`
- 新: `density_partition_weight_local(...)`

へ置き換える。

local/self/subgroup/recv の全足し込み経路で同じ重みを使う。

## Validation

見るべき指標:

- `total_charge` のずれが 3% からどこまで減るか
- `normalize` が不要に近づくか
- `density charge budget` で raw と weighted の差が妥当か
- Hartree / XC / reconstruct が壊れないか

## Risks

- taper が強すぎると中心領域まで削りすぎる
- taper が弱すぎると `1/count` と大差が出ない
- partition-of-unity 正規化を誤ると局所的な密度欠損を起こす

## Recommendation

最初は `cosine taper` を最小変更で入れ、総電荷ずれと密度形状を確認する。必要なら後で `erf/erfc` に差し替える。
