# A/E/H受け渡し・分子外側結合の性能モード（2026-11-23）

HSE検証は57dfb441でcommit/push済み。その後、A/E/H用の値の受け渡しと分子の外側結合へ性能モードを接続した。数値演算は既存Fortranカーネルを使い、物理的な更新順は変えていない。A/E/Hの自作応答は合成試験用の計算である。

A/E/H: 実E/H＋LDと2つの自作物質群。Hを取得する前の合意、単位・位置・ステップ番号の契約を維持した。65番目の更新中の例外後、場と自作状態だけでなく受け渡しアダプターのステップ番号・失敗状態も64ステップ終了時へ戻して再試行した。

分子: 既存C2H2・PZのGSを読取リンクで再利用し、新しいGS計算やコピーを行っていない。試験時だけコピー関数をリンクへ置換し、実ワーカー通信と本計算・自己真空参照のFDTDは製品コードを使用した。内部ワーカーの構築入口を使う検証であり、新たな公開ランタイム認証記録を発行したものではない。真空／LD環境、FDTD各群1／2ランクで比較した。33番目の場の更新の検証を拒否した後、両FDTDを32ステップ、分子を31ステップ終了時へ直接復元し、同じワーカーで再試行した。分子更新回数は失敗した1回を含む64回で、初期状態から履歴を再計算していない。

性能モードでも、分子の試行・受理に必要な退避、自己真空参照の差し引き、予測／受理電場の一致条件、明示的な検証を維持する。プロセス停止やMPI障害に対する一般的な復帰保証ではない。明示的ディスク復元で生成する分子候補にも、復元先ワーカーのモードを渡すよう修正したが、今回の実計算はメモリ内復元の検証である。

```json
{
  "schema": "salmonpy.packet-molecular-outer-performance/1",
  "packet": {
    "source_cases": 6,
    "wheel_cases": 6,
    "ranks": [
      3,
      6
    ],
    "ports": [
      "A",
      "E",
      "H"
    ],
    "steps": 96,
    "failed_update": 65,
    "restored_accepted_step": 64
  },
  "molecule": {
    "source_cases": 12,
    "wheel_cases": 12,
    "ranks": [
      4,
      6
    ],
    "environments": [
      "vacuum",
      "ld"
    ],
    "field_steps": 64,
    "molecular_steps": 63,
    "rejected_field_update": 33,
    "restored_field_step": 32,
    "restored_molecular_step": 31,
    "worker_updates_after_retry": 64,
    "prefix_replay": false,
    "scope": "internal existing C2H2 worker validation, not a newly issued public runtime receipt"
  },
  "comparison": "array_equal: transferred fields/current, molecular dipole/electron number/energy, accepted field output streams",
  "tests_each": {
    "total": 709,
    "passed": 708,
    "skipped": 1
  },
  "GS_collected": false,
  "limitations": [
    "no physical convergence or speed claim",
    "no disk/process-exit resume campaign",
    "no process/MPI failure recovery"
  ],
  "wheel_sha256": "4fac18af3cab63c0154d59343859b8ca10075d215f128454a68fc39632f44708"
}
```

記録: calc2026/salmonpy-packet-molecular-performance-20261123。波形・要約・ログ・検証コード・wheelのみ保存し、GSや電子状態のディスク保存は含まない。ローカル保存を確認したがクラウド同期完了は未確認。物理的収束、速度向上、多ノードの証拠ではない。
