# ACE vendor threading implementation plan

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

**Goal:** Extend the approved single-k / multiple-k threading policy to Fujitsu BLAS/LAPACK and NVIDIA NVPL without changing the production benchmark binary.

**Architecture:** Keep optional control callbacks. Carry an explicit three-integer scope state (active, prior BLAS count, prior LAPACK count), with request zero restoring that state. Enter and leave scopes on the calling thread and each k worker. OpenBLAS changes its process-wide setting only outside the parallel region; NVPL changes BLAS and LAPACK thread-local settings independently; Fujitsu changes the OpenMP task setting. Unknown providers retain the serial-k fallback.

**Tech stack:** Fortran, OpenMP, CMake link probes, OpenBLAS / Fujitsu SSL2BLAMP / NVPL.

1. Add a failing OpenMP scope test; test distinct worker settings, nested scope restoration, and independent NVPL defaults (including zero).
2. Adapt the callback contract and k-worker scopes; extend vendor adapters and configuration detection. Fujitsu selection requires its compiler plus SSL2BLAMP link flags. NVPL requires both local-control symbols. No new namelist input.
3. Run unit probes with OpenMP on/off, OpenBLAS real library, NVPL API test doubles and Fujitsu OpenMP semantics under GNU. Build MPI and non-MPI; run related MPI regressions.
4. Document commands, supported link choices and limitations. Vendor ABI/performance verification remains an on-machine task; test doubles do not establish vendor correctness.

Sources: NVIDIA NVPL BLAS and LAPACK service APIs; RIKEN Fujitsu BLAS/LAPACK and SSL II manuals. The sequential-library fallback and deferred CTest/manual registration policy are preserved. No push in this task.

## 利用方法と検証範囲

既存のCMakeビルドを再configureすると、選択したリンク構成から制御APIを検出する。追加のnamelist指定は不要。configure出力の `EXX thread control:` を確認する。k点数はMPI分割後のランク内の値であり、全体のk点数ではない。通信callbackを使うACE構築ではk点OMPを有効にしない。

- 富岳: 既存 `platforms/fugaku.cmake` の `-Kopenmp -Nfjomplib` と `-SSL2BLAMP` を維持する。富士通コンパイラとSSL2BLAMPの指定を確認してOpenMP task制御を選択する。OpenMP無効ビルドでは新たな並列化を行わない。
- MIYABI-G: CPU BLAS/LAPACKの両方をNVPLへリンクする。NVHPCの `-Mnvpl` 等はインストール済みSDKの構成に従い、CMakeのリンク設定にも含める。両方の `*_set_num_threads_local` がリンクできた場合のみ有効。NVPLの既定値0とBLAS/LAPACKの異なる設定も保存・復元する。cuFFTの有効化だけではNVPL使用もACEのGPU実行も意味しない。
- OpenBLAS: 既存の実行時制御を維持する。プロセス全体に影響するため、別のアプリケーションスレッドが同時にBLASを使わない通常のACE呼出しを前提とする。
- その他/未検出: 既存のライブラリ設定を保持し、k点ループは逐次実行する。MKL/ArmPLの専用制御は今回追加していない。

`python3 developer_tests/651_hybrid_exchange/threads/run.py` はGNU FortranとOpenBLASで実行する手元検証。OpenMP有効/無効、OpenBLAS制御有効/無効、富士通向けOpenMP task設定の入れ子と復元、NVPL API模擬実装による各ワーカーの設定・独立した既定値・ACE作用・異常終了時の復元を確認する。NVPL模擬試験およびGNUのOpenMP試験は、富岳/MIYABIの実ライブラリのABIや性能の検証ではない。実機で同一入力をOMP1/複数スレッドで比較する必要がある。進行中の5 fs計測バイナリには変更を適用しない。

参照:
- https://www.r-ccs.riken.jp/fugaku/docs/manual/en/lang/math/j2ul-2575-01enz0.pdf
- https://docs.nvidia.com/nvpl/latest/blas/api/service.html
- https://docs.nvidia.com/nvpl/latest/lapack/api/service.html

## 手元検証結果

2026-09-30: 専用試験6構成、GNU MPI/ScaLAPACKビルド、GNU非MPIビルド、CTest 438/439とfixtureを含む11試験が成功。`work/exx-block3d/ace-vendor-{unit,build,serial,ctest}.log` に結果を保存。読取りレビューで重大な問題なし。実ベンダーライブラリと性能は未検証。

## ACE適用への拡張（承認済み）

構築と同じoptional thread_controlを適用にも渡す。通信callbackなし・ランク内nk>1・OMP複数スレッド・制御対応時にk点OMP、それ以外は逐次kとBLAS内部並列。各workerがoverlap(no,nt)を一度確保し、処理後に解放して設定を復元する。MPIを含むcallbackは並列領域の外に維持。任意ターゲット数、ゼロ交換、通信呼出し回数、OMP有無、制御有無、設定復元を既存単体試験へ追加し、MPI/nonMPIビルドと既存回帰で確認する。

適用側の追加検証: 専用試験6構成（矩形target、格子点ゼロ、無効target、構築/適用それぞれのNVPL worker設定復元）、GNU MPI/非MPIビルド、関連CTest11件に成功。`work/exx-block3d/ace-apply-{unit,build,serial,ctest}.log`。読取りレビューで重大な問題なし。実行速度は未計測。

## 分散ゲージへの適用（承認済み）

`gauge_tiles_refresh` の全体を一つのスレッド制御スコープとする。内部処理を `refresh_core` に置き、早期returnがあっても外側で必ず復元する。これにより通常のpzgemm、分散SVD/固有値計算、および極分解の代替反復内のpzgemmを同じ設定で実行する。反復やMPIを呼ぶ処理そのものをOMPループにしない。既存並列領域内またはOMP無効時は設定変更なし。通常・異常経路のスコープ実行/復元をMPI試験で確認し、分散ゲージの既存数値比較をOMP1/2で実施する。周辺の明示オブジェクトリンク式の試験には新しいモジュール依存を加える。

分散ゲージ検証: MPI2/4×OMP1/2、各ScaLAPACK有効/無効の8構成に成功。OpenBLAS実APIで設定変更と復元を確認。MPI/nonMPIビルド、関連CTest11件に成功。`work/exx-block3d/gauge-threads-*.log`。読取りレビューで重大な問題なし。速度とベンダー実機は未検証。

## 分散ACE構築と逆変換（承認済み）

`orbital_ace_build` に外側wrapperを設け、元の内部packerと早期returnを保ちながら、分散metricのpzheevと複製metricのzheevを同じスレッド設定で実行し、必ず復元する。分散/複製、packed/dense、sparse training、異常入力、空ランクを既存試験で検証する。

逆変換はexx_nativeの非分散gaugeのMATMUL前後と、ブロック逆変換のownerループ前後に限定して設定する。MATMULがBLASを使うかはコンパイラ設定による（手元GNUビルドは-fexternal-blas）。ブロック経路のuse_blocked_inverse=.false.は維持。通常のgauge_tiles_rotate/orbital_rotateは列通信と配列更新で、BLAS制御だけでは高速化されない。多k点交換、LCFO GS対角化、PTCN等は今回の対象外。性能未計測、進行中のバイナリは維持。

分散ACE追加検証: ScaLAPACK有効MPI1/2/4×OMP1/2×軌道数2/7、無効MPI2/4×OMP1/2×軌道数2/7が成功。通常/packed/sparse training、異常入力、空ランク、実OpenBLAS設定復元を含む。MPI/nonMPIビルドと関連CTest11件、条件数回帰も成功。`work/exx-block3d/orbital-threads-*.log`。無効のブロック逆変換はビルド確認のみ。読取りレビューで重大な問題なし。

## Wannier正変換・逆変換の局所更新（承認済み）

orbital_rotate / gauge_tiles_rotate のoutput更新を明示的な格子点×出力軌道ループにし、OMP collapse(2)で並列化する。列のbroadcast/reductionは並列領域の外に置き、入力列順序と各出力要素の加算順序を維持する。追加の全メッシュ配列は不要。複素U、非一様weights、正変換/随伴変換、軌道を持たないrank、OMP1/2/4とOMP無効を検証する。列ごとのOMP領域起動コストを伴うため、速度改善は実測まで断定しない。

正逆変換検証: 複素U・非一様weights・随伴変換でMPI1/2/4、OMP1/2/3/4の出力が完全一致。OMP指示無効でも参照値と一致。ScaLAPACK有効/無効、既存ゲージ更新試験、MPI/nonMPIビルド、関連CTest11件成功。`work/exx-block3d/rotation-omp-*.log`。読取りレビューで重大な問題なし。
