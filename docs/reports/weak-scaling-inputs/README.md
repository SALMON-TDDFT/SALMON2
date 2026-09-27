# Diamond弱スケーリング：入力と再実行

[統合ノートへ](../../../DEVELOPMENT_NOTES.md)

| ディレクトリ | 原子数 | 配置 | MPI | GS nstate | RT nstate |
|---|---:|---|---:|---:|---:|
| [c32](c32) | 32 | 4×1×1 | 4 | 128 | 64 |
| [c64](c64) | 64 | 8×1×1 | 8 | 256 | 128 |
| [c128](c128) | 128 | 16×1×1 | 16 | 512 | 256 |

各ディレクトリに`gs-inputfile`、`rt-inputfile`、`atom.dat`を保存。擬ポテンシャルはリポジトリの[testsuites/pseudo/C_rps.dat](../../../testsuites/pseudo/C_rps.dat)を使用します。GSバイナリデータは含めていません。新規GSからのRTは可能ですが、初期状態を再生成するため、公開測定の保存済みGSとバイト単位で同じ結果になるとは限りません。

## ビルド条件

測定バイナリの実装は`88385397`。MPI/HSE/ScaLAPACKを有効化、GNU Fortran15、Release/O3、外部BLAS、AArch64のループベクトル化を無効化（`-fexternal-blas -fno-tree-loop-vectorize`）。ScaLAPACK/OpenBLASのリンク先は環境に合わせて指定してください。[ビルド仕様](../../hse-build.md)を参照。全構成にこのコンパイラ固有フラグを適用するものではありません。

## C32の例

リポジトリ直下で、ビルド済み実行ファイルを`build/salmon`とした例です。未使用の`benchmark-c32`ディレクトリで実行します。MPIの配置オプションは使用するMPI/計算機に合わせてください。

```sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=FALSE

# GSにRT用の設定を持ち込まない。
unset SALMON_LCFO_RT SALMON_LCFO_RT_MLWF SALMON_LCFO_RT_DIRECT_WF
unset SALMON_LCFO_RT_RADIUS SALMON_LCFO_RT_ACE_INTERVAL SALMON_LCFO_RT_U_INTERVAL
unset SALMON_LCFO_RT_CONTINUITY SALMON_LCFO_RT_FFT_BATCH

mkdir -p benchmark-c32/gs benchmark-c32/rt
cp docs/reports/weak-scaling-inputs/c32/gs-inputfile benchmark-c32/gs/inputfile
cp docs/reports/weak-scaling-inputs/c32/rt-inputfile benchmark-c32/rt/inputfile
cp docs/reports/weak-scaling-inputs/c32/atom.dat benchmark-c32/gs/
cp docs/reports/weak-scaling-inputs/c32/atom.dat benchmark-c32/rt/
cp testsuites/pseudo/C_rps.dat benchmark-c32/gs/
cp testsuites/pseudo/C_rps.dat benchmark-c32/rt/
(cd benchmark-c32/gs && mpirun -np 4 ../../build/salmon < inputfile > run.log 2>&1)
```

GSの`end SALMON`と最後の`DC #SCF ... diff`が1e-7未満であることを確認してから進めます。今回のC32は1147反復、9.9556585e-8でした。

```sh
ln -s ../gs/data_dcdft benchmark-c32/rt/data_dcdft
export SALMON_LCFO_RT=1
export SALMON_LCFO_RT_MLWF=1
export SALMON_LCFO_RT_DIRECT_WF=1
export SALMON_LCFO_RT_RADIUS=6
export SALMON_LCFO_RT_ACE_INTERVAL=1
export SALMON_LCFO_RT_U_INTERVAL=1
export SALMON_LCFO_RT_FFT_BATCH=1
(cd benchmark-c32/rt && mpirun -np 4 ../../build/salmon < inputfile > run.log 2>&1)
```

C64/C128はそれぞれ対応する入力とMPI8/16を使います。すべて逐次実行し、同時に別SALMONジョブやビルドを走らせません。測定時はOpenMPIの`--map-by slot --bind-to none --nooversubscribe`を使用。CPU bindingは無効でした。

## 比較するもの

`rt iterations`の最大時間を使用し、弱効率を100×T32/TNと定義。GS・初期局在化は別集計。電流16時刻、有限値、最終密度積分4×原子数、Gram誤差、ACE構築35回、Taylor呼び出し32回、初期＋終点Gram検査17回を確認します。OS負荷・CPU配置の違いを記録してください。

[測定JSON](../diamond64-mlwf-support/weak-32-64-128-results.json)と[入力・バイナリのハッシュ](manifest.json)を同梱。後者のローカルパスは公開用の識別子で、実在パスではありません。
