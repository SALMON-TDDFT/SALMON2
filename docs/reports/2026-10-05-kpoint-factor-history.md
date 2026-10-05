# k点別学習ACEの入力整合とMac検証

2026-10-05。k点別履歴がGamma専用SR条件を必須にしていたため通常の多k点入力と両立しなかった。Gammaは従来のSR距離閾値1e-3、支持率1、source ACEを維持。多k点はSR距離閾値0、支持率0、occupied ACEを要求する。HSE遮蔽核は維持。多k点のSR距離打切り・局所FFTは未実装。

静的HSE06、整数占有、全均一kメッシュ、対称性削減なし、k-only MPI、Taylor4 predictor/correctorの制約を維持。Gammaの既存空間MPI経路も維持。SALMON_FACTOR_HISTORY=learnedで予測、strictで厳密比較。

## 再現条件と結果

Si8原子、10.26bohr立方セル、32電子16軌道、12³格子、8MPI×1OMP。GS閾値1e-8。RT横方向z impulse1e-4、dt0.16au、256step（約0.991fs）、3履歴、厳密補正8step、warmup80。OMP_NUM_THREADS=1、OPENBLAS_NUM_THREADS=1。区間時間はGFORTRAN_UNBUFFERED_ALL=yでstdoutのstep終了時刻を記録し、end80→end256の176stepを比較（厳密補正込み）。短い単発試験であり長時間・収束スペクトルの精度保証ではない。

| kメッシュ | GS回数/残差 | 厳密秒step | 学習秒step | 学習後速度比 | 学習後Jz L2誤差 |
|---|---|---:|---:|---:|---:|
|2³|69 / 8.4002938e-9|0.122154|0.035562|3.435|1.377%|
|4³|88 / 7.3539855e-9|1.230729|0.333793|3.687|1.061%|

全8/64k点で教師更新、81stepから154step予測。両方式正常終了、電流256行/energy257行、有限・時刻一致。2³時刻再測定は元計算と電流差0。全区間電流L2誤差は2³1.236%、4³0.776%。最大energy差1.838e-8 / 1.279e-8Ha。ただし予測ACE期待値energyと厳密瞬間交換汎函数energyを区別する。

試験実行ファイルは採用LAPACK版オブジェクト集合に修正exx_nativeをリンク。SHA256 2bac28ff4fd8c66b8c841e206c686b3014521440620e40ed5bdcfb3d8c6cd327。公開版はエラーメッセージのみ132文字以内に短縮。LCFO-only ae6a4e11変更は試験バイナリに含まない。

証拠はローカルwork/si-k222-learned-mpi8-v5-20261005、work/si-k222-learned-mpi8-timed-20261005、work/si-k444-learned-mpi8-20261005のplan/inputfile/run.py/analyze.py/output.log/timing.json/comparison.json。先行入力失敗v1〜v4を保持。機械依存バイナリ・計算データはGitに含めない。
