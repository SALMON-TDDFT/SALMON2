# WSe₂用スピノルpair FFT拡張

ソース037c301。既存HSE06＋SOC＋NLCC・SU(2)対称性展開を維持し、全4スピン密度ブロックに(r,r′) pair分散を追加した。起動時SALMON_EXX_FFT_ORBITS=pairで有効、既定は従来経路。

非等方kメッシュ、64bit pair番号、両FFT配置、不均等/空MPIランク、NaN・経路不一致・過大タイルの拒否を検査。17回帰試験合格。独立レビューの整数変換上限指摘を修正済み。

同一最終バイナリSHAはextension-record.json。小さい矩形Si₂のSOC/NLCC fixtureでGS3条件の残差<1e-8、x/y各20stepの対称性ON/OFFとACE/直接交換が合格。同一保存済みGSからpair ON/OFFを比較し、全RT出力の最大差はx=1.6822e-16、y=1.9022e-15。原本はresults.json、same-gs-pair-legacy.jsonとworktree validation/spinor-pair-final-20261019/。

この版は実空間点の対称性によるFFT削減を実装せず、全pairをタイル分配する。Si用pair FFTの性能結果は転用しない。WSe₂の数値収束・富岳実行時RSS・高速化率は未確認。GPU、q削減、学習ACEは今回追加していない。

富岳の独立配置先は既存WSe₂領域のpairfft-spinor-20261019/。既存build/salmonとSi用binaryは保持する。ジョブ取消・再投入はしない。

富岳の別ビルドは終了コード0。新binary SHA: eb0329041a9ddd121c9004fe6f8a1aaf53bb6097febfdf06e69fa8bddb303998。2026-10-19の明示承認を受け、ジョブ52458963・52458965・52458967の起動先を専用wrapperで切替えた。2026-10-19T13:55:49.061387+09:00の検証時は全件QUE。旧binary、wrapper、input、preflightは保持。入力・資源・ジョブIDは不変、取消・再投入なし。実使用ログ・全rank終了・RSSは計算開始後の確認事項。

[富岳の切替検証原本](fugaku-records/deployment-verified.json)とビルド原本を回収し、転送前後のアーカイブSHA一致を確認済み。
