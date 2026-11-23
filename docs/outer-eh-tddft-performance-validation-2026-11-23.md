# E/H・LD・TDDFT結合の性能モード検証（2026-11-23）

既存の合成GSを読取リンクで再利用し、新規GS計算・GSコピーは行わない。試験時だけコピーと目録関数を置換し、製品側のコピー／検査契約は変更していない。電子・電磁場の更新と明示的な保護更新は製品コードを使用した。

{
  "schema": "salmonpy.outer-eh-tddft-performance/1",
  "steps": 96,
  "ranks": [
    3,
    6
  ],
  "field_ranks": [
    1,
    2
  ],
  "ranks_per_tddft_group": [
    1,
    2
  ],
  "tddft_points": 2,
  "source_and_wheel_cases": 12,
  "failed_step": 65,
  "accepted_step_after_restore": 64,
  "comparison": "array_equal: A/E/j and LD-cell field observation",
  "GS_collected": false,
  "scope": "PZ synthetic existing orbitals; actual E/H and LD and TDDFT kernels, not physical convergence or speed benchmark",
  "LD_cell_max_abs_E": 0.00117801626960232,
  "cross_layout_max_abs_error": 4.336808689942018e-19,
  "wheel_tests": {
    "total": 707,
    "passed": 706,
    "skipped": 1
  },
  "wheel_sha256": "61ebd7d7642e31e076b695d335a4c00cb65c67059f092fd6def088df5eff3b9f",
  "library_sha256": "f6a964f52b013f337c85566f28ee15381b97fc65cc109f042fa867219f8fcb2f"
}

通常の連続実行、通常／性能モードでの途中例外・明示的な保護復元・再試行を比較。65番目の更新中に実際の電子状態を進めた後で例外を返し、場・LD・電子状態を64ステップ終了時へ戻した。元の初期状態から履歴を再生する経路は使用していない。LD領域内の電場が非ゼロであることも確認した。

記録: calc2026/salmonpy-outer-eh-tddft-performance-20261123。波形・要約・ログだけを保存した。クラウド同期完了は未確認。高速HSE・分子の外側スケジュール・他環境・物理的収束・速度向上の証拠にはしない。
