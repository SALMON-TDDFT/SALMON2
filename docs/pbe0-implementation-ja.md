# PBE0の実装と使い方

更新日：2026-09-29。対象ブランチ：`pbeh40-rvv10-water-md`。

## 汎関数と実装

`xc='pbe0'`で、25% Fock交換＋75% PBE交換＋100% PBE相関を使う。
全距離交換でありHSE06の短距離遮蔽は使わない。rVV10は加えない。
PBEh(40)とHSE06の既存の定義・既定値は変更していない。

交換係数は`src/xc/exx_functional.f90`の`exchange_fraction()`で管理する。
Fock交換と`salmon_xc.f90`の半局所項が同じ係数を参照する。
DC、実空間分散、MLWF、局所FFT、ACE、実空間RTは既存のPBEh経路を共用する。
LCFO保存データには汎関数名と係数、Coulomb cutoff等を記録し、PBEh(40)などの
異なる汎関数のデータをPBE0として再構成する入力は拒否する。

## ビルド

MPI、ネイティブHSE/PBEh機能（`USE_HSE=ON`）、Libxcが必要。
今回のビルドでは`USE_MPI=ON, USE_HSE=ON, USE_LIBXC=ON, USE_SCALAPACK=ON`を使用した。
LCFO分散全対角化を使う場合は`configure.py --enable-scalapack`を指定する。
アーキテクチャ設定によるHSE有効化とライブラリの準備は
[ビルドノート](hybrid-distribution-note-ja.md)および[HSEビルド](hse-build.md)を参照。

## DC-GS入力

既存の系・並列設定に対して以下を指定する。

```fortran
&functional
 xc='pbe0'
 pbeh_coulomb_radius=4d0
 exx_pre_scf_threshold=1d-4
 yn_exx_dc_mlwf='n'
 exx_mlwf_norm_fraction=1d0
 exx_mlwf_radius=0d0
/
&dc
 yn_dc_lcfo='y'
 lcfo_eigensolver='scalapack'
/
```

PBEで予備収束した後にPBE0へ切り替える。バッファ付き各フラグメントの交換は
局在化・支持切断なしで評価する。300 Kの電子占有はGSで共通化学ポテンシャルを
決めるための設定であり、レーザー下のRT温度制御ではない。

`pbeh_coulomb_radius`はPBE0でも同じ名前で使用する。入力の長さ単位に従い、
0はBorn–von Karmanセルの最短辺の半分を選ぶ。
**有限半径でのCoulomb切断はPBE0の定義ではなく数値近似**であり、物性計算には収束確認が必要。
以下のH2例の4 bohrは既存のPBEh(40)比較と揃えた値。

## 実空間RT入力

同じPBE0のDC-GSが出力した`data_dcdft`を使う。GSとRTでCoulomb cutoffを揃える。

```fortran
&calculation
 theory='tddft_response'
 yn_dc='n'
 yn_conventional_from_dcdft='y'
/
&functional
 xc='pbe0'
 pbeh_coulomb_radius=4d0
 yn_hse_wannier='y'
 exx_mlwf_interval=5
 exx_mlwf_maxiter=1000
 exx_mlwf_tolerance=1d-6
 exx_mlwf_radius=0d0
 exx_mlwf_norm_fraction=0.999d0
 exx_local_fft='auto'
 exx_ace_support='source'
/
```

波動関数は全体系の実空間メッシュで伝播する。上記の局所支持は交換・ACE構築に使う。
`source`方式のACEは固定イオンのnative RT用。移動イオンへの適用を検証した設定ではない。

## 32 H2の実行例

完全な入力を[samples/pbe0_h2](../samples/pbe0_h2)に置いた。
セル64×16×16 bohr、格子128×32×32、4×1×1フラグメント、32分子。
GSはMPI4、RTはMPI4・空間分割1×2×2。各プロセスのOpenMP/BLASスレッド数は1。

リポジトリ直下からの例（`build/salmon`は上記機能を有効にした実行ファイル）：

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p run-pbe0/gs run-pbe0/impulse run-pbe0/zero
cp samples/pbe0_h2/gs.inp run-pbe0/gs/inputfile
cp samples/pbe0_h2/impulse.inp run-pbe0/impulse/inputfile
cp samples/pbe0_h2/zero.inp run-pbe0/zero/inputfile
cp testsuites/pseudo/H_rps.dat run-pbe0/gs/
cp testsuites/pseudo/H_rps.dat run-pbe0/impulse/
cp testsuites/pseudo/H_rps.dat run-pbe0/zero/
(cd run-pbe0/gs && mpiexec -n 4 ../../build/salmon < inputfile > output 2>&1)
```

GSログの最終`DC #SCF`の`diff < 1d-10`、`end DC-LCFO complex`、`end SALMON`を確認してから：

```sh
cp -R run-pbe0/gs/data_dcdft run-pbe0/impulse/
cp -R run-pbe0/gs/data_dcdft run-pbe0/zero/
(cd run-pbe0/impulse && mpiexec -n 4 ../../build/salmon < inputfile > output 2>&1)
(cd run-pbe0/zero && mpiexec -n 4 ../../build/salmon < inputfile > output 2>&1)
```

二重起動や既存結果の上書きを避け、新しい実行ディレクトリを使う。

## 比較と検証状況

- 定義・保存パラメータ・Libxc半局所項の3テスト通過。
- PBE0/HSE06/PBEh(40)のGS→RT、MPI配置一致、source ACEとfallback、
  異なる汎関数データの拒否を含む8テスト通過。
- 32 H2 PBE0 GS：約106秒で収束。16ステップRTは正常終了、ACE33回すべて成功。
- PBE0の7000ステップ計算とHSE06の長時間計算は執筆時点で進行中。
  4汎関数の最終スペクトル・速度比較は未確定。

dt=.05 au、7000ステップは8.466 fs、公称Fourier間隔約0.489 eV。
今回のDC初期状態はRT Hamiltonianに対して完全には定常でないので、各汎関数で
`J_impulse - J_zero`を取り、三次窓を掛けて応答を比較する。
この差引きだけで平衡基底状態からの誘電関数になるとは主張しない。
共有計算機の同時実行数が変わるため、全経過時間だけから精密な速度比を決めない。

検証例（実行ファイル・MPIランチャーは環境に合わせて絶対パスを指定）：

```sh
export SALMON_TEST_EXE="$PWD/build/salmon"
export SALMON_TEST_MPIEXEC="$(command -v mpiexec)"
PYTHONPATH=developer_tests/653_functional python3 -m unittest test_source_ace test_source_ace_fallback
```
