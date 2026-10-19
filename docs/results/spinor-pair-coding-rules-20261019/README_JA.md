# スピノルpair FFTのCoding Rules対応（2026-10-19）

今回のpair FFT拡張をCODING_RULES.mdに合わせて修正した。

- 全変更手続き・内部手続き・追加インターフェースに明示的implicit none。依存は局所use/onlyへ限定。
- 実験用環境変数選択を廃止し、&functional の yn_exx_spinor_pair_fft（既定n）へ移行。
  default/read/broadcast/lowercase/variables.log/yn_argument_check/互換性検証を追加。
  yはCPU・全q・既存固定イオンSOC HSE06条件の全pair経路。nは従来経路。
- 標準testsuites/442・443・444を追加。CPU/HSE、MPI条件で登録。
  単体のSU(2)/ACE/基準交換比較、両FFT配置、不均等/空ランク、集団エラー拒否、
  SOCなしでの不正な有効化を確認。準備・実行・検証のfixtureとラベルを設定。
- 初回CTestで準備処理がprobeを削除した問題を検出。実行ファイルとmodule出力を
  ケースディレクトリ外の専用ビルド領域へ移し、全項目を再実行して修正を確認。
- 入力仕様はdocs/inputs/hse-feature-selection.mdとSALMON-DOCSの
  source/input_keyword_list.rstへ反映。SALMON-DOCS commit 20f15a0b24ea5f0d746f5cc53d3418f641870178。
  追加したhybrid feature表をdocutilsで警告なしに構文検証（サイト全体の生成は未実施）。

GNU Fortran・MPI/OpenMP/HSE有効のout-of-sourceビルド成功。CTest 15/15、既存回帰17/17合格。
単体MPIは4ランクと8ランク（k点を持たないランクを含む）。追加変更Fortran行は132文字以下。
独立レビューで入力処理とテスト登録を再確認し、未解決の必須指摘なし。

同一バイナリSHA256: 206fe1398df56c6cbf6b937b9639edaf2553ead6cd87adf4cfc2ef18b1abceda。
小さな矩形Si2のSOC/NLCC開発fixtureで、namelistをyにしてGS3条件が残差<1e-8で収束。
x/y各20stepのRT、対称性ON/OFF、ACE/直接交換が合格（results.json）。
偏光方向の対称性はx=nyy、y=yny。
同じGSからのRTで旧環境変数だけでは既定nのまま、大文字Yはyへ正規化されることを確認。
既定経路とpair経路の最大差1.6822e-16（selector-results.json）。

442/443は交換カーネル単体試験であり、全体のGS/RT試験とは区別する。
WSe2物理収束・性能の証明ではない。今回の修正版についてFujitsu再ビルド、
GPU/OpenACC、q削減、学習ACE、MPI無効構成の再ビルドは実施していない。
再起動形式や物理カーネルは変更していない。

富岳の配置済み037c301と既存バイナリ・ジョブ入力・起動先は今回変更していない。
富岳の037c301は旧環境変数方式であり、本修正版のnamelist方式と混同しないこと。
将来の富岳反映では独立ビルドと対応入力を一緒に更新し、既存バイナリを保持する。
