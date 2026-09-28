# ハイブリッド汎関数の共通実装名（2026-09-29）

HSE06、PBE0、PBEh(40)、PBEh(40)+rVV10で共有する厳密交換処理を
`exx_*` に整理した。計算式、分散方法、局所化・支持領域・ACE更新条件は変更していない。

| 旧モジュール／ファイル名 | 新モジュール／ファイル名 | 役割 |
|---|---|---|
| `hse_native` | `exx_native` | 実空間GS/RTから交換作用への接続 |
| `hse_ace` | `exx_ace` | ACE構築・作用 |
| `hse_wannier` | `exx_wannier` | MLWFと局所交換 |
| `hse_wannier_gauge` | `exx_wannier_gauge` | 局在化ゲージ |
| `hse_symmetry` | `exx_symmetry` | 対称性による軌道再構成 |
| `hse_spatial` | `exx_spatial` | 空間分割交換 |
| `hse_lcfo_rt` | `exx_lcfo_rt` | LCFO側の交換接続 |
| `hse_semilocal` | `hybrid_semilocal` | HSE/PBEh族の半局所成分 |
| `hse_ptcn`, `hse_ptcn_core` | `exx_ptcn`, `exx_ptcn_core` | PT-CN反復と実空間接続 |

対応する公開APIも `hse_refresh` → `exx_refresh`、`hse_ace_build` →
`exx_ace_build`、`lcfo_hse_refresh` → `lcfo_exx_refresh` などへ変更した。
共有型は `s_exx_ace`、`s_exx_wannier`、`s_exx_symmetry_map` とした。
既に役割が明確な `wannier_*`、`spatial_exx_*` は維持する。

`hse_exchange` と `hse_kernel_*` はHSE短距離カーネル専用なので変更しない。
`hybrid_semilocal` 内の `hse_semilocal_evaluate` もHSE固有評価なので維持する。
rVV10は非局所相関として独立した既存の `rvv10_*` モジュールに残る。

## 入力・保存データとの互換性

- namelist名・変数名・既存の互換入力は変更していない。
- `propagator='hse_ptcn'`、`'hse_taylor4'`、`'hse_taylor4_full'` は従来通り。
- ビルドオプション `USE_HSE`（configureのHSE有効化）も従来通り。
- 再開データ、スナップショットのファイル名・形式・識別情報は維持する。
- 共通局在化ログの `HSE_WANNIER` は `EXX_WANNIER`、LCFO共通交換ログの
  `LCFO HSE` は `LCFO EXX` に変更した。これらを読む外部スクリプトは更新が必要。
- Fortran内部APIの旧モジュール名を直接利用する外部コードは新名へ更新し、再ビルドする。

実行中のスペクトル計算用ビルドは変更せず、検証専用のビルドで確認する。

## 検証

検証専用ビルドはMPI/HSE/Libxc/ScaLAPACK有効のGNU Fortran 15。
MPI・HSE無効のCPU構成でもビルド成功を確認した。
改名に伴い、対象手続き37箇所へ明示的な `implicit none` を補った。
また422のHSE標準試験に計算完了fixtureを設定し、並列CTestで検証が先行する問題を修正した。

- MPI/HSE/Libxc/ScaLAPACK有効およびMPI/HSE無効のビルド：成功。
- 汎関数判定、PBE0/PBEh半局所項、HSE/PBEh/PBE0のsource ACE、
  fallback、PBE予備収束などと、PBEh+rVV10の実空間イオン運動・
  空間分割レーザー応答：25/25通過。
- PT-CN参照比較、HSE半局所項、MLWFゲージ勾配：3/3通過。
- 422（DC-HSE）、423/424（DC-PBE0→RT）の準備・計算・検証：
  並列CTestで9/9通過。422の検証は現在の既定仕様に合わせ、
  フラグメントSCFのcanonical交換とMLWF局所化なしを確認する。
- 改名前後24ソースの実行トークンを、識別子置換・文字列・コメント・
  `implicit none`追加を除いて比較し、一致を確認。
- 独立レビューと `git diff --check`：残る指摘なし。

ログは作業ディレクトリの `work/exx-rename-*.log` に保存。
アクセラレータ、他コンパイラ、長時間スペクトルの再計算は今回の改名検証には含まない。
