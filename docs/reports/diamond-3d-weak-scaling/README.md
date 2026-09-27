# Diamond：3次元弱スケーリング入力

既存の1次元Diamond系列を3次元に展開した **未実行の入力セット**。
[開発ノート](../../../DEVELOPMENT_NOTES.md) / [富岳ビルド](../../hse-platforms.md#fugaku)

## ケース一覧

1フラグメント＝Diamond通常立方晶セル8原子、1 MPI/fragment、軌道MPI1。

| ケース | 原子数 | 電子数 | MPI（GS/RT共通） | 全格子 | セル辺長 bohr | GS / RT nstate |
|---|---:|---:|---:|---|---:|---:|
| [4×4×4](4x4x4) | 512 | 2048 | 64 | 64³ | 26.88 | 2048 / 1024 |
| [6×6×6](6x6x6) | 1728 | 6912 | 216 | 96³ | 40.32 | 6912 / 3456 |
| [8×8×8](8x8x8) | 4096 | 16384 | 512 | 128³ | 53.76 | 16384 / 8192 |
| [10×10×10](10x10x10) | 8000 | 32000 | 1000 | 160³ | 67.20 | 32000 / 16000 |

各ケースに `gs/inputfile`、`rt/inputfile`、`atom.dat`。共通の `C_rps.dat` はリポジトリの `testsuites/pseudo/C_rps.dat` と同一。実行時の相対パスを使うため、このディレクトリ構成を維持してください。生成条件・SHA256は `manifest.json`。

## 固定条件

- 格子定数6.72 bohr、コア16³、格子間隔0.42 bohr、Gammaのみ。
- **バッファーは8,8,8格子**。フラグメント周期セルは32³格子・辺長13.44 bohr・64原子。従来1Dの8,0,0と異なるので、1Dの実行時間との直接比較はしない。
- `nstate_frag=256`（局所128占有＋128空状態）。全体GSは4状態/原子、RTは2占有状態/原子。DC密度由来のGS→LCFOの差は従来どおり受け入れる。
- HSE06、SCF閾値1e-7、最大1800反復、従来のmixrate0.01等を継承。収束はこの新しい3D入力では未検証。
- `lcfo_eigensolver='chefsi'`、filter degree60、最大200cycle、残差許容1e-7。ScaLAPACK有効ビルドが必要。`yn_scalapack='n'`は各fragmentの通常SALMON対角化の設定であり、LCFO CheFSIやRTの自動分散線形代数を無効にしない。
- RTは直接WF係数Taylor4、dt0.02 a.u.、16steps、x方向impulse1e-4。MLWF局所積分R6 bohr、ACE1/U1、FFT batch1、実測FFT計画OFF。`rt-env.sh`で設定する（namelistだけでは実験的LCFO RTにならない）。
- エネルギー出力間隔40、最終密度出力16。16stepsは速度比較の短時間試験であり、収束した誘電関数を得る長さではない。

## 現実装のメモリ制約：大きい入力を直ちに投入しない

`src/xc/lcfo_rt_wannier.f90` の初期局在化で `raw_local(no,no,6,1)` と `raw(no,no,6,1)` が各空間rankに確保される。complex(8)を16byteとして **この2配列だけ** の容量は以下。U、ACE、WF、基底、rootへの係数集約などは別途必要。

| ケース | 占有軌道数 | 上記2配列だけのGiB/rank |
|---|---:|---:|
| 4³ | 1024 | 0.188 |
| 6³ | 3456 | 2.136 |
| 8³ | 8192 | 12.000 |
| 10³ | 16000 | 45.776 |

10³はこの下限だけで約45.8 GiB/rankなので、32 GiB/node構成では1rank/nodeでも収まらない。8³も4rank/nodeではこの2配列だけで48 GiB/node。MPI数を増やすだけではこの複製配列は減らない。**入力の生成完了は現在のコードで全サイズが実行可能なことを意味しない。** 特に8³/10³の実測前には初期MLWFと残る密行列のメモリ分散・削減が必要。4³から開始し、最大メモリを記録する。GS/CheFSIのメモリと初期局在化時間も別に確認する。

## 実行例：4×4×4

計算資源の確保はサイトの手順に従う。以下の `mpiexec -n 64` はランチャーの例で、Slurm環境ならサイト指定の `srun -n 64` などへ置き換える。ログインノードで数値計算を走らせない。MPI/ノード・OMPスレッド・CPU配置は全ケースで固定し、メモリ要件を満たす設定を選ぶ。ここではジョブスクリプトや課金グループを仮定しない。

```sh
# このREADMEのあるディレクトリから。実行ファイルには絶対パスを設定。
export SALMON_EXE=/absolute/path/to/SALMON2/build/salmon
# OMP_NUM_THREADSやランク配置は確保した資源に合わせて全ケース同一に設定。
export OMP_DYNAMIC=FALSE
. ./gs-env.sh
(cd 4x4x4/gs && mpiexec -n 64 "$SALMON_EXE" < inputfile > run.log 2>&1)
```

`end SALMON`、最終DC-SCF残差<1e-7、CheFSI対角化完了と `data_dcdft` 出力を確認してからRTへ進む。未収束/失敗GSをそのまま再利用しない。

```sh
(cd 4x4x4/rt && ln -s ../gs/data_dcdft data_dcdft)
. ./rt-env.sh
(cd 4x4x4/rt && mpiexec -n 64 "$SALMON_EXE" < inputfile > run.log 2>&1)
```

他のケースはディレクトリとMPI数を表の値へ変更する。GS/RTはそれぞれ対応するサイズの初期状態を使用し、古い出力を混ぜない。同じ計算機・バイナリ・環境設定で逐次測定する。

## 記録と静的検証

RTループの `rt iterations` 最大時間を取り、弱効率は **100×T(4³)/T(N³)** とする。GS・初期MLWF・出力時間は別集計。実時間、最大メモリ、MPI/ノード、OMP数、配置、ACE/U/FFT設定、コンパイラ、バイナリ版を保存する。ACE35構築・Taylor32呼び出し・Gram17検査と有限値を確認。最終密度積分の目標は電子数の列。サイズ間の電流差には有限サイズと異なるGSの影響があるため、実装誤差とは呼ばない。

```sh
python3 validate.py
```

原子数・一意性・セル内配置・8原子/core・格子/MPI/状態数・参照パス・ハッシュの静的チェック済み。SALMONの実入力読込、SCF収束、3D RTの安定性、富岳の実行と速度測定は未実施。`generate.py` は入力を再生成するので、手編集後に実行すると上書きされる。
