# Vh専用 dual rho 戦略 実装仕様書

## 0. 文書メタ情報
- 区分: 提案・設計案（未実装）
- 対象: SALMON v2.2.2（DC->DG連携経路）
- 目的: `rho_h` を Hartree (`Vh`) のみに使い、`Vxc` は従来どおり `rho` を使う dual rho 戦略を導入する

## 1. 背景と狙い
- DG-RT ではフラグメント境界近傍の密度不連続が Hartree 計算へ悪影響を与える可能性がある。
- 一方で `Vxc` は局所密度の物理量として、過剰な平滑化を避けるため `rho` を維持したい。
- `rho_h` は `rho` の派生量であり、`rho` 自体の規格化崩れや電子数ずれを補正するための経路ではない。
- したがって dual-rho は Hartree 安定化策であって、電子数補正策ではない。
- 本仕様では、Hartree 入力専用の補助密度 `rho_h` を導入し、以下を満たす。
  - `Vh = Vh[rho_h]`
  - `Vxc = Vxc[rho]`
  - DC seed 生成段階と DG-RT 段階の密度取り扱いポリシーを整合させる

## 2. 用語定義
- `rho`: 既存の全系電子密度（スピン和）
- `rho_s(ispin)`: 既存のスピン分解密度
- `rho_h`: Hartree 計算専用に構成した全系電子密度
- `rho_frag`: DGフラグメントごとに切り出した `rho`
- `rho_h_frag`: DGフラグメントごとに切り出した `rho_h`

## 3. スコープ
### 3.1 対象
- DC側変更: `yn_dc_for_dg='y'` のときのみ `rho_h` を選択的に構成
- DG側変更: `rho` 構成後に `rho_h` を作成し、`rho` と `rho_h` の両方をフラグメント分割してハミルトニアン行列要素を計算

### 3.2 非対象
- `yn_dc_for_dg='n'` の通常DC計算
- `Vxc` 入力を `rho_h` に切り替える変更
- HSE/RI専用カーネル最適化（本仕様では副作用を出さないことを優先）

## 4. 基本方針
1. データフローを明示的に二系統化する。
   - Hartree系: `rho_h` 系列
   - XC系: `rho` / `rho_s` 系列
2. `rho_h` は `rho` から生成する派生量とし、原データ `rho` は不変で保持する。
3. `rho_h` 生成は境界不連続抑制に限定し、総電荷保存を必須とする。
4. DCとDGで同一の `rho_h` 構成ロジック（または同等条件）を使う。
5. DG側では、各 RK stage において occupied-state renormalize 後の係数から `rho` を再構築し、その健全な `rho` を入力としてのみ `rho_h` を生成する。

## 5. 機能要件
## 5.1 新規サブルーチン要件（`rho -> rho_h`）
### 5.1.1 目的
- フラグメント境界由来の高周波成分を抑えた Hartree 用密度を生成する。

### 5.1.2 インタフェース案
```fortran
subroutine build_hartree_density_from_rho(mg, info, dg_frag, rho_in, rho_h_out, mode, preserve_charge)
```
- `rho_in`: 元密度（入力）
- `rho_h_out`: Hartree用密度（出力）
- `mode`: 平滑化モード（将来拡張用。初期は単一モード固定可）
- `preserve_charge`: 総電荷保存フラグ（初期値 `.true.`）

### 5.1.3 構成方式（初期実装方針）
- 初期実装では `rho_h` を実空間局所平均ではなく、全系密度 `rho` の Fourier 成分に対する高波数フィルタにより構成する。
- 基本手順は以下とする。
  1. `rho_h = rho` を初期値として複製する
  2. `rho_h` を全系 periodic 条件で Fourier 変換する
  3. 高波数成分 `|G| > Gcut` を減衰または除去する
  4. 逆 Fourier 変換して実空間の `rho_h` を得る
  5. 総電荷保存のため再正規化を行う
- `Vh` 用補助密度であるため、平滑化は全系 periodic 密度に対して一貫して適用し、フラグメントごとの独立平滑化は行わない。

### 5.1.4 必須要件
- 充電保存: $\int \rho_h(\mathbf{r}) d\mathbf{r} = \int \rho(\mathbf{r}) d\mathbf{r}$
- 有界性: `rho_h` が NaN/Inf を含まない
- 可逆不要: `rho_h -> rho` の逆変換は要求しない
- 決定性: 同一入力に対してMPI分割に依存しない結果（許容誤差内）
- 全系整合性: `rho_h` の生成は全系 periodic 境界条件の下で一意に定まり、フラグメント分割位置に依存しないこと
- 低波数保持: 低波数成分は極力保持し、除去対象は高波数成分に限定すること
- `rho_h` は `rho` の規格化不良を補正するための量ではない
- `rho` 再構築の電子数整合と occupied-state 正規化は、`rho_h` 生成前に満たされていることを前提とする
- 前提条件: 入力 `rho_in` は既に規格化・電子数整合が取れた健全な密度であること
- 禁止事項: `rho_h` を occupied-state renormalize の代替や電子数補正の吸収先として使わないこと

## 5.2 DC側要件
### 5.2.1 有効化条件
- `yn_dc_for_dg='y'` のときのみ `rho_h` 構成を実行する。

### 5.2.2 処理フロー
1. 従来どおり `rho` / `rho_s` を構成
2. `build_hartree_density_from_rho` で `rho_h` を生成
3. Hartree は `rho_h` から計算
4. XCは従来どおり `rho_s`（および `rho` 系）から計算
5. DG seed export に必要な密度系を出力（`rho` と `rho_h` の両方。互換モードあり）

### 5.2.3 互換性
- 既存の `yn_dc_for_dg='y'` ワークフローを壊さない。
- 出力フォーマット拡張時は「旧読込可能」またはバージョン識別子を必須とする。

## 5.3 DG側要件
### 5.3.1 有効化位置
- DG-RT 自己無撞着更新ループで、各 RK stage の occupied-state renormalize 完了後に `rho` を構成し、その直後に `rho_h` を生成する。

### 5.3.2 処理フロー
1. 各 RK stage で occupied-state renormalize を先に完了する
2. renormalize 後のフラグメント係数から `rho` / `rho_s` を再構成する
3. 健全性を満たした `rho` を入力として `build_hartree_density_from_rho` で `rho_h` を生成する
4. Hartree: `rho_h` を入力して `Vh` を計算
5. XC: `rho_s`（および `rho` 系）で `Vxc` を計算
6. `rho` と `rho_h` をそれぞれフラグメントへ分割
7. ハミルトニアン行列要素計算で
   - `Vh` 寄与は `rho_h_frag` 経路
   - `Vxc` 寄与は `rho_frag` / `rho_s_frag` 経路

### 5.3.3 実装注意
- `rho` と `rho_h` のメモリ領域を分離し、上書き事故を防止する。
- 既存の適応基底更新ロジックとの整合を保ち、更新判定に使う量は現行仕様を維持する。
- occupied-state renormalize 前の係数から `rho` を再構築しないこと。ここを外すと既知の 480 問題が再発し得るため、実装・レビュー時の必須確認項目とする。
- `rho_h` 導入の有無にかかわらず、電子数整合は `rho` 再構築側で満たす。`rho_h` 側で補正しない。

## 6. データ構造変更要件
- `s_dg_fragment_rt` または関連構造体に以下を追加（必要最小限）。
  - `type(s_scalar) :: rho_h`
  - `type(s_scalar), allocatable :: rho_h_frag(:)` または同等の分割バッファ
- DC->DG seed データにも `rho_h` チャネルを追加（互換フラグ付き）。

## 7. 入出力・フラグ仕様
## 7.1 新規フラグ（案）
- `yn_dual_rho_vh_only = 'y'|'n'`（デフォルト `'n'`）
- `dual_rho_h_mode = 'fft_lowpass'`（初期実装の既定）
- `dual_rho_h_gcut_ratio`（実数; `Gcut = dual_rho_h_gcut_ratio * Gnyquist`）
- `dual_rho_h_filter_shape = 'hard'|'soft'`（初期は `'soft'` 推奨）

## 7.2 有効条件
- DC側: `yn_dc_for_dg='y' .and. yn_dual_rho_vh_only='y'`
- DG側: `yn_dg_fragment_rt='y' .and. yn_dual_rho_vh_only='y'`

## 7.3 ログ
- ステップごとまたは所定間隔で以下を記録する。
  - `Q(rho)` と `Q(rho_h)` の差
  - `max|rho-rho_h|`
  - `Vh[rho]` との差分ノルム（デバッグ有効時）

## 8. 受け入れ基準
1. 保存則
   - $|Q(rho_h)-Q(rho)| < \varepsilon_Q$ を満たす
2. 数値安定性
   - no-kick で電流ドリフトが現行許容範囲内
   - NaN/発散を起こさない
3. 物理妥当性
   - 短時間スペクトル差分が基準ケースに対して許容範囲内
4. 後方互換
   - `yn_dual_rho_vh_only='n'` で既存結果と一致（丸め誤差レベル）

## 9. 実装フェーズ
1. Phase 1: 基盤
   - `rho_h` 生成サブルーチン追加
   - 共通ユーティリティ化（DC/DG共用）
2. Phase 2: DC適用
   - `yn_dc_for_dg='y'` 分岐で `rho_h` を導入
   - seed出力フォーマット拡張
3. Phase 3: DG適用
   - `rho` 再構成直後に `rho_h` 生成
   - `rho`/`rho_h` 二重分割と行列要素経路分離
4. Phase 4: 検証
   - no-kick / kick / 長時間安定性の回帰実行

## 10. リスクと緩和
- リスク: `rho_h` 平滑化により局所物理が過度に失われる
  - 緩和: `Vxc` を `rho` 維持、`rho_h` は Hartree 限定
- リスク: MPI分割依存差
  - 緩和: 全体和での正規化と再現性テストを必須化
- リスク: I/O互換破壊
  - 緩和: バージョンタグまたは互換読込分岐を実装

## 11. 実装完了の定義
- DCとDGの両経路で dual rho が同一ポリシーで有効化可能
- `Vh` のみ `rho_h`、`Vxc` は `rho` の分離がコード上で明示される
- 受け入れ基準を満たす回帰ログを添付してレビュー完了
