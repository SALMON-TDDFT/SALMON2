# ハイブリッド実装のCoding Rules対応

2026-09-29。`CODING_RULES.md`に対するレビューで見つかった項目を修正した。

## 修正範囲

- `exx_functional`へ、引数の汎関数名を分類するpure関数`is_hybrid(name)`と
  `is_global_hybrid(name)`を追加。GS・RT・入力検証・交換処理の重複判定を共用した。
  HF混合率、カーネル、データ形式は変えていない。
- 新しい判定関数のimportは利用手続き内の`use, only:`に限定した。
- PBE0追加で変更した132文字超の行を、判定共通化により短縮した。
  今回追加・変更したFortran行に132文字超はない。
- レビュー対象の`exx_functional`、`hse_native`、`hse_lcfo_rt`、`hse_semilocal`、
  `rvv10_distributed`とXC呼出し手続きに明示的な`implicit none`を補った。
- `yn_exx_dc_mlwf`の検証を既存の`yn_argument_check`へ統一した。
  この既存手続きは通常のFortran STOPを使うため、テストは終了コードだけに依存せず、
  不正入力メッセージと計算が未完了であることを確認する。
- `423_H_pbe0_dc_gs`と`424_H_pbe0_mesh_rt`を標準CTestへ登録。
  GS検証成功をRT準備のfixture依存にし、チェック済みGSデータをコピーする。
  CPU・HSE・MPI構成でのみ登録し、MPIなしの構成で存在しないテストに属性を付けない。

## 公式マニュアル

SALMON-DOCSの`source/input_keyword_list.rst`へ開発ブランチのnative hybrid入力を追加した。
ローカルブランチ`docs/native-hybrid-inputs`、コミット`be6aa8e`。
汎関数、半径、MLWF、PBE予備収束、ACE、rVV10と適用範囲を記載した。
MLWFの既定値は互換入力の解決後の値（10、200、1d-6）を記載し、計算例の指定値と区別した。
追加節はdocutilsで警告なしに構文検証済み。Sphinxサイト全体のビルドは未実施。

## 対象外

`hse_native`等の保存状態を計算コンテキスト型へ移す大規模な変更と、`hse_*`の一括改名は
行っていない。これは独立した設計・回帰検証が必要な構成改善として残る。
リポジトリ全体の既存コードを規約に合わせる一括整形も行っていない。

実行中のスペクトルは元のバイナリで継続し、検証には別ビルドを使用した。

## 検証結果

- GNU Fortran 15.2、Release、MPI/HSE/Libxc 5.2.3/ScaLAPACK有効の別ビルド：成功。
- 標準CTestの423/424（準備・計算・検証）：6/6通過。
- 汎関数判定・保存パラメータ・半局所項：3/3通過。
- HSE/PBEh/PBE0 source ACE、fallback、PBE予備収束：20/20通過。
- MPI/HSE/Libxc/ScaLAPACK無効のシリアルビルド：成功。
  既存101 C2H2 GSの計算は成功。CTestの検証起動は`python`実行名不在で失敗したが、
  同じ検証スクリプトを`python3`で同じ出力に対して実行し、固有値比較は通過した。
- HSE有効・MPI無効のCMake configure/generate：成功（この構成の全ビルド・実行は未実施）。
- 今回の追加・変更Fortran行の132文字検査、`git diff --check`：通過。
- OpenACC/CUDA、EigenExa、他コンパイラでは未検証。

独立レビューで指摘されたMPIなしのテスト登録条件を修正した。
