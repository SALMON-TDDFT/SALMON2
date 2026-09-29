# EXX作業領域の寿命整理

承認済みの優先順1（コピー・同時保持削減）、2（LCFO再構成後の解放）を調査する。

1. native Γ点refreshでtarget_workをlocalへmove_allocして再利用する。
2. unchangedの場合はtarget_workへ戻し、既存キャッシュは維持する。
3. changedの場合のみ旧cached_sourceを解放する。MLWF更新後にaction_workをwへ受け渡す。
4. DCではwをcached_actionへ、通常RTではaction_workへ移す。作用側は2つの作業配列を独立に確保する。
5. 新たな通信・再計算・入力変更は導入しない。全k点経路は変更しない。
6. LCFO再構成の配列寿命を確認し、既に局所変数として解放されているものは変更しない。
7. MPI RT回帰と128 H2のimpulse 16ステップを各1回実行し、採用済みACE分散の保存結果と比較する。

## 結果

1の作業領域再利用を実装し、回帰13試験合格。128 H2の1回比較でRSS約11%削減、電流・エネルギー一致。2のLCFO再構成後の配列は既に解放済みで追加変更なし。詳しくは[測定ノート](../exx-workspace-memory-ja.md)。軌道ブロック化と支持球内保存は次段階。
