# DG Density Owner-Distributed Communication Design

**Date:** 2026-03-29

**Goal:** `calculate_density_from_fragments` の world 通信集中を解消し、`rank 0` に偏っている `recv_wait` / `unpack` を owner rank 側へ分散する。

## Background

現在の密度再構築では、fragment subgroup 内で `send_pack -> send_sum` を行った後、subgroup root だけが world 通信を担当している。その結果、owner grid を多く持つ root rank では `recv_wait` と `unpack` が集中し、実測では `rank 0` が `recv_wait ~ 1.36s` を抱えていた。一方、多くの他 root rank は `recv_wait=0` で `send_wait` を抱えており、通信パターンの偏在が主因と判断できる。

## Chosen Approach

subgroup 内 reduce は維持しつつ、reduce 後の world 通信担当を subgroup root 固定から subgroup 内の複数 rank に分散する。

- `density_subgroup_send_count` と `density_subgroup_send_slot_map` はそのまま利用する
- `send_pack -> send_sum` による subgroup 内集約も維持する
- ただし reduce 後の `send_sum` を world target ごとに subgroup 内担当 rank へ再配分し、その担当 rank が `comm_isend` / `comm_irecv` を実行する
- owner 側の受信後 unpack は既存 `density_recv_map` を使ってそのまま行う

## Why This Approach

この方式は、現行コードで最も複雑な前半の project 分散ロジックを大きく崩さずに、通信律速の本丸である root relay だけを切り出して改善できる。`send_pack` の構築や grid-to-owner slot 管理を作り直す必要がなく、world 側の受信分散だけを導入できるため、全面的な direct-send 化よりリスクが低い。

## High-Level Data Flow

1. subgroup 各 rank が block ごとの密度寄与を `send_pack` に積む
2. `comm_summation(..., dg_frag%icomm_frag)` で owner world rank ごとの payload を subgroup 内で集約する
3. world target `irank` ごとに subgroup 内担当 rank `handler_rank(irank)` を決める
4. `handler_rank(irank)` のみが `send_sum` から `rho_send(irank)` を取り出して world `isend` を出す
5. world owner rank 側も同じ規則で `irecv` を担当し、受信後に既存 `density_recv_map` で unpack する

## Handler Assignment

第一案では単純な静的規則を使う。

- `handler_id_frag = modulo(irank, dg_frag%isize_frag)`
- subgroup 内で `dg_frag%id_frag == handler_id_frag` の rank がその world target を担当する

この規則の利点は、sender/receiver の両方で独立に再計算でき、追加の通信や map 保存が不要なことにある。偏りが残る場合は将来、`density_recv_count` や `density_send_count` に基づく重み付き割当へ拡張する。

## Required Code Changes

### `src/rt/dg/rt_dg_density_reconstruct.f90`

- `rho_send` / `rho_recv` の確保条件を subgroup root 固定から subgroup handler 判定へ変更する
- `send_sum` から `rho_send` へコピーする箇所を handler rank 限定にする
- `comm_isend` / `comm_irecv` の発行条件を handler rank 限定にする
- trace に handler 観点の情報を追加し、担当分散が正しく働いているか見えるようにする

### `src/rt/dg/rt_dg_fragment_types.f90`

- 追加の persistent state が必要なら handler 用 map を保持する
- まずは静的規則を再計算する設計なので、必須変更は小さい見込み

### `src/rt/dg/rt_dg_fragment.f90`

- 必要に応じて handler map の deallocate/invalidate を追加する

## Error Handling

- `distribute_project=.false.` の経路は現状維持
- `icomm_frag == COMM_GROUP_NULL` または `isize_frag <= 1` の場合は現状の single-rank 振る舞いを維持
- handler 規則と送受信割当が合わない場合に備え、debug trace で `rank/id_frag/handler_id_frag/target_rank/npts` を出せるようにする

## Validation Plan

1. `make -C build -j2 salmon`
2. 同じケースで root 別 trace を再取得
3. `rank 0` の `recv_wait` が大きく低下し、複数 root に `recv_wait` / `unpack` が分散していることを確認
4. 総電荷規格化や Hartree 入口到達が変わらないことを確認

## Risks

- 静的 `mod` 割当だと owner 分布次第で偏りが残る可能性がある
- sender/receiver で handler 規則を完全一致させないと deadlock になる
- 既存 trace は root 前提の箇所があるため、ログ解釈を少し更新する必要がある

