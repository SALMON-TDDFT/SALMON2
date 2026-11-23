# 別プロセス型高速PR版HSEの必要性と性能モード（2026-11-23）

## 判断と対応範囲

この経路は、呼び出し元のMPI・Fortranのプロセス共有状態・作業領域を子プロセスへ分離する明示的な入口として維持する。内部MPIや縮約k点の高速HSEには、検証済みの呼び出し元MPI経路を優先する。別プロセス型は従来から全k点の範囲で、縮約入力は変更前ワーカーもOPEN 1で初期化できない。今回は物理対応範囲を拡張せず、この拒否が変わらないことを確認した。

## 実装と検証

性能モードをMODE_PERFORMANCE／SALMONPY_HSE_PERFORMANCE 1で更新前に合意する。rt_stepの既存checked引数へ接続し、新しいプロセス共有変数や数値更新式は追加しない。通常モード・初期化・SHA・容量/形状/ステップ番号・明示的保存復元の検査は維持する。未対応ワーカーは拒否し、黙って別経路へ下げない。

同じ最新のGit-only準備済み計算核で、変更前プログラムと変更後プログラムをビルドした。ソース・全アプリケーションオブジェクトをワークスペースの禁止ソルバーチェッカーで監査し、CMakeの検査も維持した。既存GSを試験用の読取リンクで再利用し、GS計算・GSコピーは行わない。製品のコピー/目録処理は変更せず、試験時だけ同一ハッシュを照合するリンクへ置換した。

全k点の学習付きTaylor法とPT-CN、各32ステップで、変更前・通常・性能モードの電流とエネルギーが完全一致した。複数状態を交互に更新して兄弟状態が変わらないこと、逆符号の駆動が独立対照と一致することも確認した。8ステップの時間発展後状態を一時保存し、元ワーカー終了後に通常/性能の新しいワーカーへ直接復元して、9〜32ステップの結果が一致した。保存形式/入力・容量は同じで、復元先のモードを維持する。未対応プログラムの性能要求拒否と、失敗した性能ワーカーを再起動しないことも確認した。

ソース/wheelで同じ結果。両スイートは710件中709成功・任意スキップ1。小さい合成軌道での接続試験であり、物理収束や速度向上、施設MPI/多ノード/他環境の証拠ではない。記録はcalc2026/salmonpy-hse-process-performance-20261123。ローカル保存を確認し、クラウド同期は未確認。

```json
{
  "schema": "salmonpy.serial-fast-pr-hse-performance/1",
  "decision": "retain explicit isolated serial CPU route; prefer host-MPI for internal parallelism/reduced symmetry",
  "supported_profiles": [
    "taylor-full",
    "ptcn-full"
  ],
  "existing_unsupported_profiles": [
    "taylor-reduced",
    "ptcn-reduced"
  ],
  "steps": 32,
  "source_and_wheel_current_energy_vs_checked_previous": "array_equal",
  "independent_states": true,
  "worker_exit_direct_restore": "array_equal; receiving mode retained",
  "old_worker_rejects_performance": true,
  "failed_worker_no_automatic_respawn": true,
  "GS_collected": false,
  "tests_each": {
    "total": 710,
    "passed": 709,
    "skipped": 1
  },
  "worker_sha256": "c3866966035ee7711eb2baeefe6fad7a5b95f1651ada51d169a90e61d8bed347",
  "scope": "macOS arm64 synthetic full-k evolution, not physical convergence or speed benchmark",
  "wheel_sha256": "44b3d1cde8f6bea26a981a1b0d837028846a44ad1abe2cfe8348ea4e12336544"
}
```
