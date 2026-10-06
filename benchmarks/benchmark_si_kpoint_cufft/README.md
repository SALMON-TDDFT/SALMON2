# MIYABI-G: Si HSE06 8×8×8 k点 CPU/GPU比較

この試験は従来型Siセル（8原子、32電子、占有16状態）、格子定数10.26 bohr、
実空間12³点、HSE06の**GS固定SCF反復**を測定する。粗い実空間格子なので物性値の
収束試験ではない。Nkは各軸の点数で、`--nk 8`は合計512 k点。

## 実装と適用範囲

`&functional` の `exx_kpoint_backend='cpu'/'cufft'` で切り替える。
Γ点局所支持版の `exx_local_backend` とは別である。ここではMLWF切り詰めを使わない。

- このベンチマークはHSE06を使用。共通エンジンはPBE0/PBEhも扱い、完全な一様kメッシュ（総数2以上）と直交セル、k点並列に対応。
- GPU：MPI受信密度タイルの並べ替え、k格子FFT、同じ離散交換カーネルの乗算、逆FFT、送信用並べ替え。
- CPU：密度行列と作用のBLAS、MPI通信、ACE、その他のHamiltonian。
- k格子FFTの逆変換は未規格化。既存の作用BLASにある`-1/Nk`を維持する。
- 定数配列・FFTプラン・GPUバッファはkernelオブジェクトの寿命中再利用し、同じデータは再転送しない。
- 各MPIタイルの入出力は転送が必要。GPU-aware MPIは不要。
- `hse_block_rows`が行タイル数。`exx_gpu_batch_size`はこの経路には使わない。
- 分数占有・空状態を含むWannier経路、DC、空間分割、縮約k点、Wannierスナップショットには未対応。

手元にはNVHPC/NVIDIA GPUがなく、実GPUコンパイル・実行・性能は未検証。

## 1. ビルド

[MIYABI-G手順](../../docs/inputs/miyabi-g-cufft-ja.md)のNVHPC＋対応HPC-X環境を使う。
FFTW、Libxc、ScaLAPACK/BLASも同じ環境に揃える。

```sh
module load nvidia nv-hpcx fftw
mpifort --showme:command
mpicc --showme:command
cmake -S . -B build-miyabi-cufft \
  -DSALMON_PLATFORM=generic \
  -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS='-gpu=cc90' \
  -DOPENMP_FLAGS=-mp -DUSE_MPI=ON -DUSE_SCALAPACK=ON \
  -DUSE_HSE=ON -DUSE_EXX_CUFFT=ON -DUSE_OPENACC=OFF
cmake --build build-miyabi-cufft -j4
```

モジュール名/バージョンと依存ライブラリの場所は現地環境に合わせる。
`USE_OPENACC=OFF`はSALMON全体の既存OpenACCルートを無効にする指定。
今回の交換用ソースには別途OpenACCが付く。GPU用ビルド1本でCPU/GPUを切り替える。

## 2. 計算ノード上のGPU数値試験

既存の `cufft-smoke.pbs` と同じ1ノードPBS環境で、実行コマンドを以下に置き換える。
FFTWは`FFTW_ROOT`または`FFTW_FFLAGS`/`FFTW_LIBS`、あるいはpkg-configで指定する。

```sh
export CUFFT_TEST_FFLAGS='-gpu=cc90'
python3 experiments/kpoint_cufft/test_gpu.py --gpu -v
```

奇数3³ k点、非連続slot、端数、kernel/座標/slot変更、プラン再利用を
独立FFTW参照と比較する。成功メッセージは`PASS k-cuFFT/FFTW parity`。
GPUが無い場合は明示的に失敗する。

続いてMPI接続全体を試す（1 GPUなので最初は1ランク）。
この試験はFFTWに加えてBLASのリンク指定が必要。NVHPC付属BLASを使う例：

```sh
export BLAS_LIBS='-lblas'
# module show fftw で確認した実際の場所を指定する。
export FFTW_ROOT=/actual/path/to/fftw
python3 experiments/kpoint_cufft/run.py --gpu --ranks 1
```

FFTW_ROOTの代わりにFFTW_FFLAGS/FFTW_LIBSを指定してもよい。
BLAS_LIBSは実際にリンクできるライブラリに合わせる。2/4ランクのGPU検証は
対応するノードを確保し、各ノードに1ランクずつ配置する。たとえば2ノードでは：

```sh
python3 experiments/kpoint_cufft/run.py --gpu --ranks 2 \
  --mpi-args='--map-by ppr:1:node --bind-to core'
```

## 3. まず2³で全体確認、次に8³

リポジトリのルートから投入する。`YOUR_GROUP`を自分の利用グループへ置き換える。
PBS内のmodule行はビルド時と同じバージョンに揃える。

```sh
qsub -W group_list=YOUR_GROUP benchmarks/benchmark_si_kpoint_cufft/miyabi.pbs
```

既定は1ノード、1 MPI×8 OpenMP、2³ k点、3 SCF反復、CPU/GPU各1回。
GPU試験と小さいSi比較が通ったら、本命の8³を実行する。
初回所要時間は未測定なので、次はshort-g・1時間の例とする。

```sh
qsub -W group_list=YOUR_GROUP -q short-g -l walltime=01:00:00 \
  -v NK=8,SCF_STEPS=3 benchmarks/benchmark_si_kpoint_cufft/miyabi.pbs
```

1時間はジョブ全体の上限。足りなければ部分ログを確認して設定を調整する。
1回が十分短ければ`run.py`に`--repeats 3`を追加して中央値を取る。
反復ごとにCPU/GPUの実行順を反転する。

2ノード・2 GPUへ進める例（最初の1ノード検証後）：

```sh
qsub -W group_list=YOUR_GROUP -q short-g \
  -l select=2 -l walltime=01:00:00 -v NK=8,SCF_STEPS=3,MPI_RANKS=2 \
  benchmarks/benchmark_si_kpoint_cufft/miyabi.pbs
```

MIYABI-Gは1ノード1 GPU。PBS例はHPC-X/Open MPIの
`--map-by ppr:1:node:PE=8 --bind-to core`を使う。現地MPIの対応オプションを確認する。

## 4. 保存結果と判定

出力先`si-k<NK>-<jobid>/`には以下を保存する。

- `settings.json`：格子、k点、反復数、プロセス/スレッド数、実行ファイル。
- `cpu-1/`と`cufft-1/`：実際の入力、擬ポテンシャル、SALMON出力。
- `summary.json`：プロセス起動込み壁時計時間、SALMON全体時間、SCF最大時間、SCF反復当たり時間、反復エネルギー。
- `comparison.json`：CPU/GPUの反復エネルギー最大差と全体壁時計速度比。
- `host-memory.rank-*.txt`：各ランクのGNU time最大RSS（LinuxではkB）。MPI launcherだけのRSSではない。
- `gpu-memory.rank-*.csv`：1秒間隔のnode-wide GPU使用量（MiB）。短いピークを取り逃がし得るため厳密な割当最大値とは呼ばない。

両者のSCF反復数が指定値と一致すること、GPU側に`EXX_KPOINT_BACKEND=cufft`が出ることを確認する。
既定のエネルギー一致条件は各反復で`1e-6 eV`以下。物理精度の保証値ではない。
不一致や異常終了は成功として集計しない。反復数を増やす前にログを調べる。

`HSE_PROFILE`のCPU側はdensity/comm/FFT/kernel/packing/actionに分かれる。
GPU側のFFT欄は転送・並べ替え・FFT・kernelをまとめた区間で、kernel欄は0になる。
rank0の診断なので、純粋なcuFFT時間や全ランク最大値とは解釈しない。
初回のプラン作成・定数転送は別途全体時間に含まれる。

速度比較は同じ格子・反復数・MPI/OpenMP・BLASスレッド数で行う。
CPU版とGPU版は同じ初期化設定から別ディレクトリで開始する。
この初回試験はGSのみであり、RTの1ステップ時間ではない。
