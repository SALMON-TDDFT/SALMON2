# 非DC有限温度HSEのRT

2026-10-06。GSで得たFermi分数占有を固定してRTに利用するCPUルートを追加。RT中の温度再熱平衡は行わない。電子密度・Hartree・半局所XCの既存更新は変更していない。

## 入力
`&system` の `temperature`（指定単位系のエネルギー）または `temperature_k`（K）をGSとRTに指定する。RTは通常通り `directory_read_data` にGSの `data_for_restart` を指定する。occupation.binの占有を既存処理で復元する。

適用範囲はHSE06、非DC、固定イオン、非偏極、完全一様kメッシュ、CPU、k点MPIのみ。Gammaも可能。nstateは電子数を収容する本数以上。`exx_mlwf_norm_fraction=0d0`、`exx_mlwf_radius=0d0`、`hse_sr_tolerance=0d0` を明示する。HSE本来のerfc遮蔽は維持。追加SR実空間打ち切りは未対応。

RTは `propagator='hse_taylor4'` と `yn_predictor_corrector='y'`。学習ACEの環境変数を指定しない。学習は既存の整数占有検証で拒否される。GPU、空間/軌道MPI、PT-CN、LCFO再構成、WF切り詰めは今回の適用外。yn_restartによるRT途中再開の既存制限は変更していない。

## 交換
重み付き源Q=psi sqrt(f/2)をcanonical frameで用いる。正の占有を捨てない。Wannier座標への変換は既存の交換実装を利用するが、非DC有限温度ではMLWF最適化・マスクはしない。ACEは全伝播軌道に対して構築し、Taylorの中点は二つの端点ACE作用を半分ずつ適用する既存処理を維持。各refreshで占有の有限性、0<=f<=2、k重み付き電子数を確認。

DCの温度依存占有と既存MLWF選択は変更していない。

## 検証
専用 work/si-k444-finite-temperature-20261006 にソース・ビルド・入力・ログ・結果を保存。
- GNU Fortran、MPI/OMP/HSE/Libxc/ScaLAPACKの専用ビルド成功。
- 既存入力で有限温度RTの拒否を再現してから変更。
- 独立なBloch和との分数占有交換作用比較: Gamma、複数k、整数占有の3試験通過。
- H2二分子、k=1x2x1、2MPI/1OMP、4軌道/4電子。kBT=.02 HaのGSが28SCFで収束。外場ゼロ/impulseの16step正常終了、有限、ノルム誤差9e-8。外場ゼロenergy変動3.89e-8 Ha。
- kBT=.1 Haの高温GS→impulse16stepも正常終了。数値はhot-validation.jsonを参照。

短時間経路検証であり、Si20fsスペクトル・レーザーの精度保証ではない。初回の構成不備(MPI無効になったbuild)と不正なmaxiter=0入力の記録も保存。正常検証はbuild-mpiと*-mpi/高温ケースのみ。
