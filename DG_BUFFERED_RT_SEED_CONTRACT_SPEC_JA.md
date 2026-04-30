# DG Buffered RT Seed Contract Specification (JA)

## 1. 目的

本仕様書は、DC-LCFO から DG-RT へ受け渡す buffered seed 一式の契約を定義する。

主目的は次の 3 点である。

- RT 開始時の `coef`, `occ`, `esp` が同じ state set を参照することを保証する
- buffered basis を使う場合でも、「基底」と「状態係数」を混同しないことを保証する
- runtime/file の `nstate` 切り詰め時に silent mismatch を起こさないことを保証する

本仕様は、最近導入した以下の修正を前提とする。

- DC-LCFO が `eigenvalues_buffered.bin` を出力する
- RT が buffered `esp` を原則 strict 読み込みする
- RT が `nstate_keep = min(nstate_tot(runtime), nstate_tot(file))` を共通 cutoff として使う
- runtime cutoff を超える metadata 参照を strict では fatal とする

## 2. 用語

- fragment basis:
  各フラグメントで構築される局所基底 `|lambda_i>`
- buffered basis:
  fragment basis の定義はそのままに、実空間 support だけを拡張した基底表現
- state coefficient:
  ある state `|psi_n>` を fragment basis 上で展開した係数
- lineage:
  `coef`, `occ`, `esp` が同じ generalized eigenproblem から来ていること

## 3. 基本方針

### 3.1 基底と係数を分離して扱う

buffered RT で使う実空間基底は `basis_functions_buffered.bin` である。

一方、状態係数は「DC-LCFO で対角化された全系状態を fragment basis 上で表した係数」であり、基底そのものではない。

buffered basis モードを選ぶ場合、実空間 support 拡張に対応して、係数ファイルも `wavefunctions_buffered.bin` を使う。すなわち、次が整合したセットを構成する：

- `basis_functions_buffered.bin`: RT が使う buffered fragment basis
- `wavefunctions_buffered.bin`: buffered 基底に対応する state coefficients
- `occupations_buffered.bin`: 同じ state set の occupation
- `eigenvalues_buffered.bin`: 同じ state set の eigenvalue

一方、buffered basis を使わない（通常）モードでは：

- `basis_functions.bin`: fragment basis（通常）
- `wavefunctions.bin`: state coefficients
- `occupations.bin`: occupation（オプション）
- `eigenvalues.bin`: eigenvalue（オプション）

### 3.2 wavefunctions_buffered.bin の名称問題と将来の置換

`wavefunctions_buffered.bin` という名前は誤解を招きやすい。

#### 問題点

「buffered」という修飾語が「buffered 基底」を示唆するため、
「buffered 基底上での係数」だと誤解しやすい。

実質は「DC-LCFO で得た全系 state coefficients」であり、基底そのものではなく、
単に basis_functions_buffered.bin との組み合わせに対応する係数に過ぎない。

#### 将来の置換計画

より明瞭な名称への置き換えを計画している：

- **推奨候補1**: `lcfo_state_coefficients_on_buffered_basis.bin`
  - DC-LCFO で得た state coefficients であることを明示
  - buffered basis との対応を明示

- **推奨候補2**: `global_state_coefficients_on_buffered_basis.bin`
  - 全系（global）state coefficients であることを明示
  - buffered basis との対応を明示

#### 当面の運用

- RT が buffered basis モードを選ぶ場合、係数源は `wavefunctions_buffered.bin` である
- 移行方針は別仕様で定める

## 4. DC-LCFO 側の出力契約

DC-LCFO は、RT seed として少なくとも以下を同時生成しなければならない。

- `wavefunctions.bin`
- `occupations_buffered.bin`
- `eigenvalues_buffered.bin`
- `basis_functions_buffered.bin`

### 4.1 state set の同一性

上記 3 つの state 系列ファイル

- `wavefunctions.bin`
- `occupations_buffered.bin`
- `eigenvalues_buffered.bin`

は、同一の state set を参照しなければならない。

具体的には、次が一致しなければならない。

- `n_frag`
- `nspin`
- `nstate_frag`
- `nstate_tot`
- `n_mat`
- `n_basis`
- `index_basis`

`coef`, `occ`, `esp` のいずれかだけが別の state ordering を持つことを禁止する。

### 4.2 `eigenvalues_buffered.bin`

`eigenvalues_buffered.bin` は RT 側の `coef(:,i)` と同じ state labeling を持つこと。

つまり、`esp(i)` は `wavefunctions.bin` の state `i` に対応しなければならない。

`esp=0` のまま残る state を許してはならない。

## 5. RT 側の読込契約

### 5.1 strict を既定とする

buffered basis モードでは、RT は `eigenvalues_buffered.bin` の strict 読み込みを既定とする。

欠落、破損、metadata 不一致、または state cutoff と矛盾する参照が見つかった場合は fatal 停止する。

互換目的の例外として、次を許す。

- `SALMON_DG_ALLOW_MISSING_BUFFERED_ESP=1`
  欠落時に warning を出して fallback 継続
- `SALMON_DG_ENFORCE_STATE_TRUNCATION=0`
  旧来の緩い truncation 契約を許可

ただし、これらは debugging / migration 用であり、正規運用では使わない。

### 5.2 共通 cutoff

RT は

```text
nstate_keep = min(nstate_tot(runtime), nstate_tot(file))
```

を `coef`, `occ`, `esp` の共通 cutoff として使う。

このとき、metadata は cutoff 後の有効 state 領域を超えて参照してはならない。

strict mode では、`index_basis` 等が cutoff 超え state を参照した時点で fatal とする。

### 5.3 RT が使うべきファイル

buffered RT の正規入力は次とする。

#### 必須ファイル セット（buffered 基底モード）

1. `basis_functions_buffered.bin` — buffered fragment basis の定義
2. `wavefunctions_buffered.bin` — buffered 基底に対応する state coefficients **（正規source）**
3. `occupations_buffered.bin` — 上記 state set の occupation
4. `eigenvalues_buffered.bin` — 上記 state set の eigenvalue

#### 必須ファイル セット（非 buffered 基底モード）

1. `basis_functions.bin` — fragment basis の定義
2. `wavefunctions.bin` — state coefficients **（正規source）**
3. `occupations.bin` — 上記 state set の occupation（オプション）
4. `eigenvalues.bin` — 上記 state set の eigenvalue（オプション）

#### 重要: ファイル対応の混在禁止

- buffered 基底モード時に `wavefunctions.bin` を使用してはならない
- 非 buffered 基底モード時に `wavefunctions_buffered.bin` を使用してはならない
- `basis_functions_buffered.bin` と `wavefunctions.bin` の組み合わせは禁止

#### 理由と原則

buffered に切り替える対象は関連する 「**基底と周辺のすべての state metadata**」である。

- basis: `basis_functions.bin` → `basis_functions_buffered.bin`
- coefficient: `wavefunctions.bin` → `wavefunctions_buffered.bin`
- occupation: `occupations.bin` → `occupations_buffered.bin`
- eigenvalue: `eigenvalues.bin` → `eigenvalues_buffered.bin`

この一括対応により、

- DC-LCFO と RT の state correspondence を一意に保証できる
- state reordering や metadata mismatch を避けられる

#### 注意: wavefunctions_buffered.bin 名前の誤解

`wavefunctions_buffered.bin` という名前は「buffered 基底専用の係数」と誤解しやすいが、
実質は「DC-LCFO で得た state coefficients（buffered 基底表現）」である。
名称の明確化は将来のファイル置換時に行う予定である（セクション 3.2 参照）。

## 6. startup 変換と `esp`

startup Lowdin, occupied-subspace unitarity, stationary projection など、RT 起動時に state へ追加変換を入れる場合は注意が必要である。

### 6.1 許されること

- `coef` を変換する
- overlap/orthogonality を改善する
- occupied subspace 内で unitary 回転する

### 6.2 必須条件

`coef` の state ordering または state content を変えた場合、runtime `esp` を古い値のまま保持してはならない。

少なくとも次のいずれかを満たすこと。

- 変換後の `H,S` に対する Rayleigh quotient で `esp` を再定義する
- 変換後の generalized EVP を解き、`coef` と `esp` を同時更新する
- fixed-H debug mode では startup 変換を無効化する

## 7. fixed-H debug mode の推奨

frozen Hamiltonian 検証では、原因切り分けのため次を推奨する。

- startup 追加変換を無効化する
- `coef`, `occ`, `esp` の lineage が揃った state set をそのまま使う
- `Hc - eSc` 残差と Rayleigh-`esp` 差を起動直後に確認する

これにより

- DC→RT 受け渡しの不整合
- RT 起動変換による再ラベリング崩れ
- propagation 側の問題

を分離しやすくなる。

## 8. 非目標

本仕様は、Hartree 用 dual-rho 導入や `rho_h` 平滑化の方式そのものは定義しない。

また、本仕様は「DC 側の全系 LCFO state を RT 初期状態へそのまま持ち込むこと」が常に物理的に正しいとは主張しない。

本仕様が保証するのは、少なくとも RT が受け取る

- `coef`
- `occ`
- `esp`
- metadata

が同じ state set を指していることだけである。

## 9. 実装反映項目

本仕様に対応する必須実装項目は次である。

1. DC-LCFO が `eigenvalues_buffered.bin` を出力する
2. RT が buffered `esp` を strict に読み込む
3. RT が `nstate_keep` を `coef/occ/esp` 共通 cutoff として使う
4. RT が cutoff を超える metadata 参照を strict で fatal にする
5. buffered basis 使用時の canonical coefficient source を `wavefunctions.bin` とする

## 10. 運用メモ

旧データセットでは、次の理由で strict mode が停止することがある。

- `eigenvalues_buffered.bin` が存在しない
- runtime/file `nstate_tot` が一致しない
- buffered metadata が runtime cutoff を超える state index を参照している

これは仕様違反の早期検出であり、意図した防衛動作である。

旧データ互換が必要な場合のみ、一時的に legacy 環境変数を使う。
