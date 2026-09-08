# DG Fragment Cap Input Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Goal

DG-Fragment RT の行列次元 cap を入力ファイルから制御できるようにし、計算量とメモリ量を柔軟に調整できるようにする。

## Background

現状の `n_mat` 制限は環境変数 `SALMON_DG_NMAT_CAP` による全体一律 cap だけである。これは緊急回避には有効だが、入力ファイルに残らず再現性が低い。また、fragment ごとの占有軌道数に応じた cap を与えられないため、小さい fragment では過大な基底を保持しやすい。

## Design

新しい制御は `&dg_fragment` に追加する。既存の `&dc` にある `nstate_frag` は fragment basis 読み込み側の設定であり、RT の DG 行列次元制御とは責務を分ける。

追加する入力は以下の 3 つ。

- `dg_nmat_cap_mode`
  - `none`
  - `fixed`
  - `occ_multiple`
- `dg_nmat_cap_fixed`
  - 全体一律 cap
- `dg_nmat_cap_multiple`
  - 各 fragment の占有軌道数に掛ける倍率

適用規則は次の通り。

- `none`
  - cap を適用しない
- `fixed`
  - `n_mat(ispin) = min(n_mat(ispin), dg_nmat_cap_fixed)`
- `occ_multiple`
  - 各 fragment ごとに `cap_frag = min(n_basis(ifrag,ispin), floor(dg_nmat_cap_multiple * nocc_frag))`
  - その fragment で `index_basis` に残す global basis index を上位 `cap_frag` 個に制限する
  - 最終的な `n_mat(ispin)` は、生き残った global basis index の最大値から再圧縮して求める

## Occupation Definition

`nocc_frag` は「その fragment に局在する occupied states 数」とする。初期実装では、全 occupied state の fragment 重み

`w(ifrag, io, ispin) = sum_j |coef(j, io, ispin)|^2`

を用い、各 state を最大重みの fragment に 1 票で割り当てる方式を採る。これにより合計票数は occupied state 数と一致し、fragment ごとの cap が過度に膨らみにくい。

## Precedence

通常運用では namelist を使う。`SALMON_DG_NMAT_CAP` は後方互換と緊急デバッグ用 override として残し、設定されている場合は最優先で適用する。

優先順位:

1. `SALMON_DG_NMAT_CAP`
2. `dg_nmat_cap_mode`
3. cap なし

## Affected Areas

- `src/io/salmon_global.f90`
  - 新規入力変数追加
- `src/io/inputoutput.f90`
  - `&dg_fragment` namelist 追加
  - default 値設定
- `src/rt/dg/rt_dg_fragment.f90`
  - cap 適用ロジックを namelist + env override 対応に更新
- `src/rt/dg/rt_dg_fragment_soi.f90`
  - 同上

## Risks

- `occ_multiple` の占有数推定を雑にすると、必要な basis を切り過ぎる可能性がある
- fragment ごとの cap 適用後に `index_basis` と `n_mat` の整合性を必ず再圧縮する必要がある
- SOI/non-SOI の両経路で同じ規則を維持する必要がある

## Validation

- `fixed` 指定で従来の `SALMON_DG_NMAT_CAP` と同じ `n_mat_max` になること
- `occ_multiple` 指定で fragment ごとの有効 basis 数が想定どおり減ること
- build が通ること
- 起動ログに mode と適用後 `n_mat_max` を出して確認できること
