# MIYABI cuFFTのAlltoallクラッシュへの対応

報告：MPI 1 rank / OpenMP 8 threads、NVHPC/HPC-Xでreturncode 139。`exx_k_exchange`のGPU状態通知→`exx_native`の内部コールバック→`communication`のMPI_Alltoall。gdbのC側PMPI_AlltoallではcommunicatorポインタがNULL。

## 調査と変更

同じコールバックによる初期ハンドシェイクは、GPU prepareより前にも実行される。従って単に「MPI 1ではコミュニケータを作らない」とは断定できない。Fortranの整数ハンドル0とC側のNULLポインタも区別する必要がある。

疑わしい経路は、カーネル内の内部手続きが、呼出し側の内部コールバックを呼び、さらに外側の`info%icomm_k`を参照する部分。GPUによるメモリ破損やMPIライブラリの不整合など、他の原因はまだ排除できていない。

修正ではcommunicatorをカーネルとコールバックの明示引数にし、SALMON側のコールバックをモジュール手続きへ移した。GPU prepare/apply後の状態通知はカーネル本体から直接呼び、二重の内部手続き参照を除いた。1ランクでも状態確認と通信は維持し、失敗したGPU処理を成功扱いしない。

## 検証

GNU MPIビルド成功。CPU上の模擬GPUでMPI 1/2/3/4、OMP 8のprepare/apply・失敗注入・端数・再利用試験に合格。複製コミュニケータを明示的に渡し、通信先の一致も検証した。CPU交換と参照解のMPI 1〜4比較も合格。

NVHPC/HPC-X/cuFFTでのSIGSEGVそのものは手元で再現・修正確認できていない。コンパイラの不具合と確定したわけではない。

## MIYABIでの再確認

コード更新後、既存のGPUビルド設定で再ビルドする。個別に適用する場合は[この問題だけのパッチ](patches/miyabi-k-cufft-communicator.patch)をリポジトリ直下で`git apply --check`してから適用する。パッチには同時開発中のメモリ削減変更を含めていない。

計算ノード上で、既存のNVHPC/HPC-X、FFTW、BLAS設定を使う：

```sh
python3 experiments/kpoint_cufft/run.py --gpu --ranks 1 --threads 8
```

試験が通った後、報告時と同じSi入力・MPI 1・OMP 8で再実行する。単体試験だけで実計算側の問題が解消したと判断しない。再発した場合はprepare前後とコールバック入口でcommunicator値が保たれるかを比較し、デバイス処理後の破損とMPI ABI/ライブラリ混在も切り分ける。
