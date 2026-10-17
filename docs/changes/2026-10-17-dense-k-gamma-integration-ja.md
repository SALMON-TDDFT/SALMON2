# dense k と巨大 Gamma 点経路の統合

巨大 Gamma 点 snapshot e798eef4a と Coding Rules 修正 31c7c9eca に、
`hse-symmetry-fft-reuse` の b5fa41e91 までの履歴を merge。

- Gamma 点の共役ペア再利用、batch8/16、SR 分散 FFT、LCFO 再構築キャッシュを保持。
- dense k の対称性削減代表点での学習 ACE と、重み・履歴・因子の次元検査を採用。
- dense k の密度積 BLAS の一時的な単一スレッド化と復元を採用。
- SVD 品質診断を保持。時間計測の出力は dense k 側の形式へ統一し、
  時刻取得は SALMON の get_wtime を使用して OpenMP なしの構成を維持。
- 固定 LCFO 部分空間 RT の未使用ソースと説明を削除。実空間 RT を使用。
- dense 側 worktree の未 commit LCFO 2 ファイルは snapshot に既に含まれ、
  元 worktree は変更していない。

検証は `work/coding-rules-e798-20261017` に保存。
MPI/HSE/LibXC/ScaLAPACK、serial HSE、serial non-HSE のビルドが成功。
通常・シフト付き dense k oracle（strided/contiguous、並進あり）、
2/4 MPI 学習履歴、8/12 MPI 共役ペア probe、既存 422/423/424 の計9件を完走。
これは小規模の機能検証であり、富岳の性能測定や長時間物理計算ではない。
Fujitsu、GPU、EigenExa は今回未検証。過去 snapshot の provenance は変更しない。
