# MIYABI-Gでの局所交換cuFFT検証

対象: SALMON2 `pbeh40-rvv10-water-md`、GPU転送再利用を含むコミット `ec7f565c`。
2026-09-29作成。以下は現地未実行の手順。最初の目標は1 GPUでのコンパイル・数値一致・再利用確認で、速度測定はその後に行う。

Si・8×8×8 k点のHSE06比較には、局所支持版ではなく新しい
`exx_kpoint_backend`を使う。[多k点のビルド・試験・PBSジョブ](../../testsuites/benchmark_si_kpoint_cufft/README.md)を参照。
下記の旧アーカイブ`ec7f565c`には多k点拡張は含まれないため、その試験では最新ブランチを使う。

## 1. ソースと環境

今回のソースアーカイブ `salmon-cufft-ec7f565c.tar.gz` をMIYABI-Gの作業領域へ転送し、展開する。
アーカイブにはコミット済みソースだけを含む。Macのバイナリ・実行中計算・結果ファイルは含まない。
GitHubから取得する場合は次のブランチを使用する。転送再利用の修正は `ec7f565c` に含まれる。

```sh
git clone --branch pbeh40-rvv10-water-md https://github.com/SALMON-TDDFT/SALMON2.git
cd SALMON2
git merge-base --is-ancestor ec7f565c HEAD
```

以下の展開手順はソースアーカイブを使う場合のみ実行する。

```sh
tar xzf salmon-cufft-ec7f565c.tar.gz
cp cufft-smoke.pbs salmon-cufft-ec7f565c/  # 添付のPBSを同じ場所へ転送してある場合
cd salmon-cufft-ec7f565c
```

MIYABI-GはGrace CPU＋Hopper GPUを搭載する。1ノードにGPUは1基。初回はMIGを使わず、`debug-g` の1ノードで試す。
公式のキュー表では `debug-g` は最大30分。以下は10分を要求する。

ログイン後、利用できるモジュールを確認する。

```sh
module avail nvidia
module load nvidia
module avail nv-hpcx
module avail fftw
```

`nvidia` と対応する `nv-hpcx`、`fftw` を選ぶ。以下の例はデフォルト名を用いるが、現地の `module avail` に合わせてバージョン付き名に置き換えてよい。選んだ名前はビルドと実行で統一する。

## 2. 最初は独立したcuFFT試験だけ実行

この段階で必要なのはNVHPC、FFTW、Python 3。SALMON全体、Libxc、ScaLAPACKのビルドは不要。
既存のPython試験がFortranの小さな検証プログラムをコンパイルし、GPUで実行する。

添付の `cufft-smoke.pbs` を使うか、リポジトリのルートに次の内容で作る。

```sh
#!/bin/bash
#PBS -N salmon-cufft
#PBS -q debug-g
#PBS -l select=1
#PBS -l walltime=00:10:00
#PBS -j oe

set -eo pipefail
cd "${PBS_O_WORKDIR:?}"
module purge
module load nvidia
module load fftw

# FFTWをpkg-configで発見できない場合は実際のインストール先を指定する。
# export FFTW_ROOT=/actual/path/to/fftw
# またはFFTW_FFLAGSとFFTW_LIBSを両方指定する。
# export FFTW_FFLAGS='-I/actual/path/to/fftw/include'
# export FFTW_LIBS='-L/actual/path/to/fftw/lib -lfftw3'

export NVFC=nvfortran
export CUFFT_TEST_FFLAGS='-gpu=cc90 -Minfo=accel'
export ACC_DEVICE_NUM=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# 初回だけ転送/カーネルの通知を残す。速度測定時にはunsetする。
export NVCOMPILER_ACC_NOTIFY=3

{
  date
  hostname
  uname -m
  module list
  nvfortran --version
  python3 --version
  nvidia-smi
  python3 testsuites/unit_pbeh_rvv10/test_cufft.py \
    --gpu -v CufftBackend.test_gpu_matches_fftw
} > "cufft-smoke.${PBS_JOBID}.log" 2>&1
```

自分の利用グループを指定して投入する。`YOUR_GROUP` は実際のグループ名に置き換える。

```sh
qsub -W group_list=YOUR_GROUP cufft-smoke.pbs
qstat
```

この試験は単一プロセスなので`mpiexec`は不要。GPUを使う実行は計算ノード上で行う。
`CUDA_VISIBLE_DEVICES`はジョブ環境が与える設定を維持する。

成功時はログに次のメッセージと unittest の `OK` が出る。

```text
PASS cuFFT/FFTW complex128 parity: one-shot, resident reuse, independent uploads, rebuilds and lifecycle
```

確認対象:

- 複素倍精度の交換作用がCPU FFTWと規格化誤差 `2e-12` 以下で一致。
- 非立方FFT・周期境界の支持点・ゼロ列・端数バッチ。
- 同じprepareでは再転送せず、source/filter/indicesの個別変更時に対応する転送回数だけ増える。
- 同じcapacityの端数バッチではcuFFTプランを再作成しない。
- geometry/capacity変更、明示release、二重release、finalization。

`--gpu`を付けたこの指定ではGPU試験をskipして成功扱いにしない。NVHPCやGPUが無い場合は失敗する。
コンパイルまたは実行の120秒タイムアウトが出た場合は、まずログを保存する。単体試験全体の所要時間にはコンパイルが含まれるため、交換処理の速度として扱わない。

## 3. 単体試験が通ったらSALMON全体をビルド

以下はNVHPC＋対応MPI環境で行う例。ログインノード上のビルド負荷・利用場所は現地手引きに従う。
Libxc/FFTW/ScaLAPACKは同じLinux Arm環境・互換なコンパイラ/MPI用を用いる。

```sh
module load nv-hpcx
mpifort --showme:command
mpicc --showme:command
```

前者が`nvfortran`、後者が`nvc`を使うことを確認する。GNU用MPI wrapperのまま進めない。

```sh
cmake -S . -B build-miyabi-cufft \
  -DSALMON_PLATFORM=generic \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS='-gpu=cc90' \
  -DOPENMP_FLAGS=-mp \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON \
  -DUSE_HSE=ON -DUSE_EXX_CUFFT=ON -DUSE_OPENACC=OFF
cmake --build build-miyabi-cufft -j4
```

`USE_OPENACC=OFF`は意図的。今回の局所交換だけにOpenACC/cuFFTを使う別バックエンドである。
`-gpu=mem:unified`や`mem:managed`は初回には加えず、実装した明示的な転送・常駐処理を検証する。

FFTW/Libxcの検出に失敗する場合は `-DFFTW_INSTALLDIR=/actual/path`、`-DLIBXC_INSTALLDIR=/actual/path` を追加する。CMakeログで実際の依存ライブラリを確認する。未検出の依存関係は自動取得・ビルドに進む場合があるので、その段階でネットワークエラーが出たら、計算ノードで繰り返さず依存ライブラリを事前準備する。

## 4. 全体計算でCPU/GPU比較

最初はΓ点・固定原子・1 MPIランク/1 GPUとする。既存の同じGS/RT入力と初期データを別の出力ディレクトリにコピーし、次のbackendだけを切り替える。

```fortran
&functional
 ! 既存の汎関数などの設定を維持して、以下を同じnamelistへ追加
 exx_local_fft='auto'
 exx_local_backend='cufft'  ! 比較側は'cpu'
 exx_gpu_batch_size=8
 exx_mlwf_norm_fraction=0.999d0
 exx_mlwf_radius=0d0
/
```

これはnamelistの差分例であり、単独で実行できる入力ではない。`num_kgrid=1,1,1`、`nproc_k=1`とする。CPU側も同じ支持切り詰めを使い、GPU効果と近似誤差を混ぜない。

まず短いGS/RTでエネルギー・電流・最終状態を比較。その後、これまでのH₂セルのimpulse 16ステップで壁時計時間・RT1ステップ時間・CPU RSS・GPUメモリを比較する。単体試験の合格はSCF/RT全体の合格を意味しない。

- 同じGPU対応実行ファイルで`cpu`/`cufft`を切り替える。
- 初期化時間とRT時間を分ける。バッチ数8を基準に1/4/16も比較する。
- `EXX_ADAPTIVE`のlocal pairsが0なら、GPU局所FFTの速度試験にはならない。
- 性能測定では`NVCOMPILER_ACC_NOTIFY`を解除し、同じCPUスレッド数・配置で複数回測る。
- 1ノードに複数GPUはない。複数GPUは2/4ノード、基本1 MPIランク/ノードから進める。MPI配置は現地のNVHPC/HPC-X手引きに合わせる。
- CPU RSSとGPUメモリは分けて記録。GH200のCPU–GPU接続だけから高速化率を予測しない。

現在の常駐期間は交換作用1回の内部。RTステップ間には保持しない。CPU側MPIのため、targetと結果はバッチ毎に転送する。

## 5. 最初の実行後に確認するログ

`cufft-smoke.<jobid>.log` とPBS出力を保存する。失敗した場合もコンパイル診断を省略しない。
全体ビルドへ進んだ場合はCMakeの設定出力、ビルド出力、`CMakeCache.txt`、入力・実行ログも残す。

## 参照

- [MIYABI公式システム仕様](https://www.cc.u-tokyo.ac.jp/en/supercomputer/miyabi/system.php)
- [公式キュー制限](https://www.cc.u-tokyo.ac.jp/supercomputer/miyabi/service/job.php)
- [公式Miyabi実践講習資料：PBS・モジュール環境](https://www.cc.u-tokyo.ac.jp/events/lectures/244/20250425-3.pdf)
- [NVIDIA OpenACCガイド：GPU指定・実行診断](https://docs.nvidia.com/hpc-sdk/compilers/openacc-gs/)
