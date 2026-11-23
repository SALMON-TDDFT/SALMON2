# E/H・LD・高速PR版HSEの外側性能モード（2026-11-23）

学習付きTaylor法／PT-CNと全k点／縮約k点の4条件で、FDTD内のLDと2点のHSEを結合した。既存の監査済みライブラリをSHA照合して使用した。新しいGS計算・GSコピー・電子状態のディスク保存は行っていない。検証時だけコピー・目録・GSパス解決を、指定済みGS配下の読取リンクへ置換し、ファイル内容・擬ポテンシャル・対称性・外場・計算・通信・退避復元の検査は製品コードを使った。製品の資産経路契約は変更しない。

通常モードの連続実行、通常／性能モードでの途中失敗・明示的な保護復元・再試行を比較した。65番目の更新で実際のHSE状態を進めた後に例外を返し、64ステップ終了時へ場・LD・電子／交換資源を直接復元して再試行した。初期状態から履歴を再生する経路は使っていない。各条件内ではA/E/J、LD領域の電場、HSEの電流・エネルギー出力が完全一致し、エネルギー出力は初期値＋96ステップの97行だった。

```json
{
  "schema": "salmonpy.outer-fast-pr-hse-performance/1",
  "profiles": [
    "taylor-full",
    "taylor-reduced",
    "ptcn-full",
    "ptcn-reduced"
  ],
  "source_cases": 24,
  "wheel_cases": 24,
  "steps": 96,
  "ranks": [
    3,
    6
  ],
  "field_ranks": [
    1,
    2
  ],
  "ranks_per_HSE_group": [
    1,
    2
  ],
  "material_points": 2,
  "failed_update": 65,
  "restored_accepted_step": 64,
  "comparison": "array_equal: A/E/j, LD-cell fields, HSE current and energy files",
  "energy_rows": 97,
  "full_reduced_max_abs_errors": {
    "taylor-3": 1.328147661294743e-18,
    "taylor-6": 1.0842021724855044e-18,
    "ptcn-3": 3.5101045334218206e-18,
    "ptcn-6": 5.7462715141731735e-18
  },
  "cross_layout_max_abs_errors": {
    "taylor-full": 1.0028870095490916e-18,
    "taylor-reduced": 7.589415207398531e-19,
    "ptcn-full": 4.4994390158148434e-18,
    "ptcn-reduced": 3.6591823321385775e-18
  },
  "comparison_rtol": 2e-10,
  "comparison_atol": 2e-12,
  "GS_collected": false,
  "tests": {
    "total": 707,
    "passed": 706,
    "skipped": 1
  },
  "library_sha256": "4d7e4c98e87b631790f20e0584fe34d707cf770f8526bbaddb6e46eb6c1b7f88",
  "wheel_sha256": "e1eb0dfdb1b8dac7f35cd112a13126291f86c90911d29128989047bef80321e1",
  "scope": "synthetic initial orbitals; not physical convergence, disk/process-exit restart, multinode execution or speed benchmark"
}
```

記録: calc2026/salmonpy-outer-fast-hse-performance-20261123。記録対象は波形、エネルギー、ログ、検証コード、wheelで、GSは含まない。ローカル保存を確認し、OneDriveのクラウド同期完了は未確認。小さな合成初期軌道による接続試験であり、物理的収束、分子の外側スケジュール、ディスク／プロセス終了後の再開、多ノードや速度向上を保証しない。
