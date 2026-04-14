# DG-RT Fragment ID Separation Spec

## 1. 目的

DG-RT で現在混在している fragment ID の意味を分離し、以下 2 系統を明示的に管理する。

- file fragment ID: DC-LCFO ファイル入出力で使う ID
- runtime spatial fragment ID: 実行時の空間所有と通信で使う ID

これにより、以下を同時に満たす。

- 基底/係数ファイルの整合性維持
- 実行時の空間 overlap 整合性維持
- Ne_coef_S, Ne_raw などの不変量回帰防止

## 2. 背景と現状課題

現状は single ifrag が以下の役割を兼任している。

- ファイルパス決定
- n_basis/index_basis メタデータ参照
- グリッド所有者判定
- halo 近傍 fragment 参照

このため、番号体系を runtime 側へ合わせる変更で file 側と不整合が発生しやすく、逆に file 側へ戻すと overlap 不整合が残る。

## 3. スコープ

### In scope

- DG-RT 初期化から time-step 1 までの fragment ID 解決経路
- 係数/基底/rg 読み込み時の ID 使用点
- density/hamiltonian/halo の owner 判定経路
- デバッグ出力の ID 表示整理

### Out of scope

- DC 側ファイルフォーマット変更
- 物理モデル式の変更
- DG-RT 以外の計算モードの仕様変更

## 4. 用語

- file_ifrag: DC-LCFO 生成物ディレクトリ番号に対応する ID (1..n_frag)
- spatial_ifrag: 実行時空間タイルに対応する ID (1..n_frag)
- root_rank(spatial_ifrag): spatial fragment root rank

## 5. 設計方針

single ifrag 依存をやめ、明示マッピングを導入する。

- spatial_to_file_ifrag(spatial_ifrag) -> file_ifrag
- file_to_spatial_ifrag(file_ifrag) -> spatial_ifrag

### 5.1 正規化ルール

- ファイルアクセス時は常に file_ifrag を使用
- 空間判定、owner 判定、halo 参照は常に spatial_ifrag を使用
- index_basis, n_basis は file メタデータとして読み込んだ直後に runtime 参照表へ再配置する

## 6. データ構造変更

対象: rt_dg_fragment 型

追加配列:

- ifrag_file_of_spatial(n_frag)
- ifrag_spatial_of_file(n_frag)
- id_array_spatial(n_frag)
- id_array_file(n_frag) (必要な場合のみ)

補助フラグ:

- mapping_valid

## 7. 主要処理フロー

### Phase A: マッピング確立

1. runtime 側で spatial_ifrag を既存ルールで確定
2. rg + basis ヘッダ overlap から spatial_to_file 候補を構築
3. 1 対 1 対応を検証
4. 双方向マップを確定

失敗時:

- fatal stop
- 不正ペア情報を全出力

### Phase B: 読み込み経路分離

1. 係数/基底/rg 読み込みループは file_ifrag で実施
2. 読み込み結果は spatial_ifrag 配列へ配置
3. n_basis/index_basis も同様に runtime 側配置を再構築

### Phase C: 実行時参照統一

1. density owner 判定: spatial_ifrag only
2. halo neighbor 参照: spatial_ifrag only
3. fragment root map: spatial_ifrag 基準

### Phase D: 診断強化

追加診断:

- assignment trace: spatial_ifrag, file_ifrag, overlap
- mapping table dump (rank0)
- metadata integrity checks

## 8. 実装タスク分解

### Task 1: 型と初期化

- rt_dg_fragment_types.f90 にマッピング配列追加
- allocate/deallocate 経路を更新

### Task 2: マッピング確立実装

- rt_dg_fragment.f90 初期化内に build_fragment_id_mapping() 追加
- overlap ベース照合 + 1:1 検証

### Task 3: file 読み込み API 分離

- ファイルパス生成は file_ifrag 専用関数化
- 既存 ifrag 直接書き込み箇所を置換

### Task 4: runtime 参照経路の分離

- owner 判定、halo、density reconstruction の ifrag 使用点を精査
- spatial_ifrag 固定に変更

### Task 5: 検証と回帰試験

- 1-step no-kick で baseline 比較
- trace on/off 双方で安定性確認

## 9. 受け入れ条件

必須条件:

- Ne_coef_S が baseline と一致
- Ne_raw が baseline と一致
- assignment trace で assigned_overlap == best_overlap
- ifrag ゼロ化問題の再発なし
- 既存 DG-RT 通常ケースでクラッシュしない

## 10. リスクと緩和策

リスク:

- index_basis 再配置の不整合
- owner 判定の hidden 連続 rank 仮定
- debug フラグ有効時のみ表面化する同期不整合

緩和策:

- マップ検証を fatal 条件にする
- owner 判定 API に ID 種別を明示
- trace と通常実行を分離して検証

## 11. 実行順序

推奨順:

1. Task 1
2. Task 2
3. Task 3
4. Task 4
5. Task 5

## 12. 初期実行コマンド

ビルド:

- make -j4

検証実行例:

- OMP_NUM_THREADS=1 SALMON_DG_ASSIGNMENT_TRACE=1 SALMON_DG_ELECTRON_PROBE=1 mpirun -np 8 build/salmon < inputfile_1step_pureDG_dt004_nokick_cmp

必要に応じて追加:

- SALMON_DG_DENSITY_TRACE=1
- SALMON_DG_BASIS_META_TRACE=1

## 13. 完了定義

以下が揃った時点で「2」の実装完了とする。

- マッピング分離コードが mainline 経路に統合
- 受け入れ条件を満たすログを保存
- 設計変更点を最小限の運用ドキュメントへ反映
