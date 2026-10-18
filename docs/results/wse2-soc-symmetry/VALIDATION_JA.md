# HSE06 SOC NLCC・直交セル対称性の検証

エンジンの確定コミット: `7eca90d9eafc36b22c1d16dbef03d3f5e294c0f7`。
符号付き対角Cartesian操作、格子に一致する並進、CPU・固定核・非DC・k方向MPIに限定する。
一般の120度スピン回転や非直交primitive cellの対称性は今回の対象外。

## ローカル最終検証

- 単体・MPI・スカラー回帰・偏光部分群・非ハイブリッド互換性の15試験: PASS。
- 最終バイナリSHA256: `70697b0a88e5f304eb51dd17c01e80718744fb6cf412966eecd525ddaf9f1c3c`。
- OpenMP2、3×1×1 k、8³格子、Si2の直交開発セル。3つのGSは全て残差1e-8未満で収束。
- full群ON/OFFのGSエネルギー差2.0e-8 eV、x部分群/OFF差1.31e-6 eV。
- x/yのRTは各20ステップ、dt=.02 au、インパルス1e-4 au。OFF/ON差は最大6.09e-10 au、ACE/直接交換差は最大4.25e-11 au。
- 数値原本はlocal-results.json。実行入力・全ログは開発worktreeのvalidation/symmetry-reviewed-finalに保持。
- Si2は粗い人工的な開発セルであり、SiやWSe2の物性データではない。

```sh
# samples/hse_mlwf_reference 内
RUN_MPI_TESTS=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python3 -m unittest test_native_spin_nlcc test_native_exchange test_native_symmetry test_spinor_symmetry_field test_spinor_density_scope
# ソースルート内（新しい出力先を指定）
python3 samples/wse2_hse_soc_nlcc/validate_symmetry.py --salmon build-local/salmon --pseudo validation/Si_PBE19.vps --output validation/symmetry-reviewed-final --steps 20
```

複数kを同時に扱う既存ACEの小行列処理で、内部手続きの非再入可能な一時領域による行列/作用の破損を再現した。
構築・作用の小さいkループを逐次化し、BLASチームと主要Fock演算のスレッド並列を保つ。
修正前に失敗する多スレッドSPD行列・占有空間補間試験を追加した。

## 富岳の独立ビルド

専用ルート:
`/vol0001/hp260497/si-3d-weak-scaling/WSe2_HSE_SOC_NLCC_e1e1a89_20261018`

ディレクトリ名のe1e1a89は作成時のラベルであり、最終sourceは上記7eca90dに固定され、clean。
Fujitsu mpifrtpx/mpifccpx、tcsds-1.2.43、-Kfast/-Nalloc_assign、MPI/OpenMP/native HSE/ScaLAPACK。
別コピーのLibxc/FFTWをリンクし、configure/buildとも終了コード0。ELF64/AArch64を静的確認。

新しい `build/salmon` のSHA256:
`954a6e1e8670e9f3c17d96596734472e7d2487ee96d3003c949ebc8b96e3daed`

既存Si511バイナリの保護SHA:
`98bb3554ee4875e8eadf65918bec347648b1b7a41f8b211b259eee144aee8918`

ビルド前後と入力転送後に一致を確認。既存ジョブ52445014は変更していない。
ログインノードでSALMON/ctestを実行していない。富岳計算ノードでのGS/RT検証やWSe2の収束計算はまだ行っていない。

## WSe2準備の到達点

公式原本W_PBE19/Se_PBE19（W12/Se6価電子、j依存、NLCC非ゼロ）を取得し、48占有スピノルの単層W2Se4を準備。
14条件×4段階の入力、対応sym.dat、ジョブ、GS収束ゲート、シート応答比較、出典、SHAを保存。
基準候補は24×18×1、比較系列は12×8から60×36まで。実空間・真空・dt・インパルス・時間窓の試験も含む。
x偏光はnyy、y偏光はyny。x用GSを別に用意し、異なる代表k集合を跨ぐ再開を避ける。
全パッケージは富岳専用ルートのcases/に転送し、SHAと静的事前確認に合格。
ジョブ投入は行っていない。WSe2のk収束・光学精度・性能・必要計算時間を確認したとは主張しない。
