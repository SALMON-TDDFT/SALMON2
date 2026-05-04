# DG Mixed-Basis rho_mix Contract

## 目的
このメモは、DG mixed basis における rho_mix 構成の物理契約をコード実装に対応づけて固定する。

## 1. Mixed basis の定義

### 数式
raw basis を
- fragment basis: {phi_i}
- plane-wave basis: {p_mu}
とし、結合して
- B = [phi_1,...,phi_Nf,p_1,...,p_Np]
と書く。

重なり行列を S_raw とすると、実装は generalized EVP の Lowdin 変換で
- S_raw = U_s diag(lambda) U_s^H
- X = U_keep diag(1/sqrt(lambda_keep))
- H_ortho = X^H H_raw X
を作り、H_ortho の固有ベクトル V を解いて
- U = X V
を mixed-basis 変換として採用している。

よって mixed basis は
- chi = B U
で定義される。

### 対応コード
- S の固有分解と keep モード選別: src/rt/dg/rt_dg_plane_wave.f90
- X 構築: src/rt/dg/rt_dg_plane_wave.f90
- H_ortho 構築と固有分解: src/rt/dg/rt_dg_plane_wave.f90
- mixed_transform への保存: src/rt/dg/rt_dg_plane_wave.f90

## 2. Mixed basis の overlap

### 数式
上の構成なら
- S_mix = U^H S_raw U
は理論上 I に近い（切断誤差・数値誤差を除く）。

### 対応コード
- apply_overlap_operator_batch が S_raw を作用: src/rt/dg/rt_dg_fragment_ops.f90
- 今回追加した診断で S_mix = U^H S_raw U を直接計算し
  - ||S_mix - I||_F
  - diag min/max
  - offdiag norm
を出力: src/rt/dg/rt_dg_density_reconstruct.f90

## 3. coef_mix の意味

### 数式
raw 係数 c_raw から mixed 係数は
- c_mix = U^H (S_raw c_raw)
で定義される（S 内積での射影）。

逆変換は
- c_raw = U c_mix
である。

### 対応コード
- raw -> mixed: sync_mixed_coef_from_raw
  - overlap_all = S_raw * raw_all
  - mixed_all = U^H * overlap_all
  - src/rt/dg/rt_dg_fragment_ops.f90
- mixed -> raw: sync_raw_coef_from_mixed
  - raw_all = U * coef_mix
  - src/rt/dg/rt_dg_fragment_ops.f90
- 型コメント: coef_mix は canonical coefficients in orthonormal mixed basis
  - src/rt/dg/rt_dg_fragment_types.f90

## 4. rho_mix の物理定義（本件の判断）

### Case 判定
コード契約は「U により（切断後の）orthonormal mixed basis を作り、coef_mix はその係数」である。

したがって第一候補は Case A:
- rho_mix = C_mix f C_mix^H

である。

ただし成立条件は実測で
- ||U^H S_raw U - I||_F が十分小さい
で確認する。

今回追加した orthonormal_cc モードは、この条件を実行時に検証しながら同式で rho_mix を作る。

## 5. 実装上の mode 定義

- legacy:
  - 既存の rho_mix = C_mix f C_mix^H
- orthonormal_cc:
  - 数式は legacy と同じ
  - 追加で U^H S_raw U 診断を強制し、契約逸脱を警告
- metric_consistent:
  - 既存の局所 Gram-Schmidt 試験経路（比較用）

環境変数:
- SALMON_DG_RHO_MIX_MODE = legacy | orthonormal_cc | metric_consistent

## 6. 追加診断（rho_mix 構成前後）

- basis metric:
  - ||U^H S_raw U - I||_F
  - diag min/max
  - offdiag Frobenius norm
- coefficient norm (occupied):
  - ||c_mix||^2
  - c_mix^H S_mix c_mix
  - (U c_mix)^H S_raw (U c_mix)
- density diagnostics:
  - tr(rho_mix)
  - diag_min/diag_max
  - eig_min/eig_max (zheev)
  - purity = tr(rho^2)/tr(rho)
- block traces:
  - tr_ff, tr_fp, tr_pp, tr_total

出力実装:
- src/rt/dg/rt_dg_density_reconstruct.f90

## 7. 現時点の一行結論候補
現状の実装契約に従えば、最初に検証すべき命題は次。

- coef_mix is already expressed in an orthonormal mixed basis, and rho_mix should be built as C_mix f C_mix^H only if U^H S_raw U is numerically close to I.
