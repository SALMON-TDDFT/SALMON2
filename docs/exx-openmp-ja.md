# 空間分割交換計算のOpenMP対応

2026-09-29。Γ点の空間分割交換が使うFFTW pencil経路と、FFT前後の対密度生成・カーネル乗算・作用加算をOpenMP化した。

## 実装

- 独立した1次元FFTを連続した区間へ分割し、スレッドごとに専用のFFTW計画で処理する。計画の作成・破棄は直列。スレッド数変更時には対応するキャッシュを再構築する。
- 1回の変換を通して同じOMPチームを維持し、転置用の詰め替えも分担する。MPI通信は主スレッドだけが行う（既存のMPI_THREAD_FUNNELEDを維持）。
- 小さいFFTは同期負担が勝ったため、従来の逐次処理を維持する。FFT担当数はOMP上限と仕事量で決まり、目安は1担当あたり65,536複素成分以上。これは性能上の粒度で、物理近似や系のメモリ上限ではない。
- FFT前後の点ごとの処理は262,144成分以上でOMP化し、担当スレッド数をバッチ列数以下に抑える。
- バッチ幅は従来の4を維持する。端数の未使用列をゼロで埋め、対数が変化するたびに別のFFT作業配列を保持することを避ける。
- 新たなnamelistやFFTWスレッドライブラリは不要。`OMP_NUM_THREADS`で上限を設定する。OMPなしのコンパイルでは逐次経路になる。

16対・64対バッチ、毎FFTのチーム生成、待機方式変更も計測したが、Si小セルで総時間が悪化する構成は採用しなかった。途中の測定も[記録](reports/exx-openmp/)に保存している。

## 速度

MPI4、128×64×64格子、4チャネルのFFT往復を10回。計画作成は時間から除く。最終版の同一実行内でOMP1/2/4を順に測定した。

| OMP | FFT往復全体 ms | 局所FFT ms | MPI ms | 詰め替え等 ms |
|---:|---:|---:|---:|---:|
| 1 | 14.096 | 5.129 | 2.167 | 8.014 |
| 2 | 10.603 | 2.650 | 2.579 | 5.930 |
| 4 | 8.841 | 1.671 | 2.524 | 4.728 |

OMP1→4で往復全体は約1.59倍、局所FFT部分は約3.07倍。各列はランク最大値なので、内訳の和は全体時間と一致しない。この結果はFFT単体であり、ハイブリッドRT全体の高速化率ではない。

Si 4×1×1（32原子、MPI4）の局所FFTは小さく、自動選択では逐次FFTが使われる。HSE06の同一GS・16ステップ比較は次の通り。各条件1回で、環境によるばらつきを含む。

| OMP | 変更前 秒/step | 最終構成 秒/step |
|---:|---:|---:|
| 1 | 1.240 | 1.236 |
| 2 | 1.139 | 1.168 |
| 4 | 1.154 | 1.192 |

最終バイナリの接続確認では、OMP4でHSE06が1.063秒/step・224.23 MiB/ランク、PBE0が1.113秒/step・223.88 MiB/ランクだった。試行間の幅があるため、このSiセルで高速化したとは結論しない。変更前に対する電流差はHSE06で8.63e-17以下、PBE0で4.10e-20以下、出力されたエネルギー差は両者0だった。

## 検証

- GNU MPI/ScaLAPACK版、MPIなし版のビルド。
- FFTのOMPなしコンパイルとMPI1/2/4、OMP1→3→4→1の切替、大小格子、端数バッチ、両スペクトル配置、逐次結果との一致。
- 大格子でHSE/Coulomb交換作用のOMP1/4一致（許容2e-12）、端数対象、対象軌道0列。
- source ACE 7試験、対省略6試験。RTはOMP4、試験用GSは従来のOMP1で用意。
- 共有FFTを使うrVV10のMPI2/4・OMP1/2回帰。
- 独立した読み取りレビュー。富岳・NVHPCでの今回の変更は未検証。

試験用H2のHSE GSをOMP4にした際、1e-10収束条件に届かなかった。変更前の実行ファイルでも同じ未収束を再現したため、RT比較のGS条件は元のOMP1に固定した。閾値は緩めていない。

```sh
python3 testsuites/unit_fftw_pencils_omp/run.py --build /path/to/build
python3 testsuites/unit_fftw_pencils_omp/run.py --build /path/to/build --no-openmp
python3 testsuites/unit_fftw_pencils_omp/run.py --build /path/to/build --action
python3 testsuites/unit_fftw_pencils_omp/run.py --build /path/to/build --benchmark
```

ベンチマークの入力・実行ファイルSHA256、ランク別RSS、数値差は[最終接続確認](reports/exx-openmp/si-laser-omp-verified.json)に記録。FFT実測は[ログ](reports/exx-openmp/fftw-omp-benchmark-final.log)を参照。
