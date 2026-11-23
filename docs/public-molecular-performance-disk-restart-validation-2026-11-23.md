# 公開ランタイムと分子性能モードの一括ディスク再開（2026-11-23）

## 検証記録

既存の監査済み実行ファイルをSHA照合し、既存のC2H2/PZ参照に対する独立状態、通常／性能モード、明示的な履歴再生、同じワーカー内の直接復元、ワーカー終了後の直接ディスク復元を測定した。検証ツールは全条件を確認してからrevision5の記録を確定する。その記録をload_isolated_runtimeで読み、公開の初期状態ファクトリーと両モードの公開入口を実際に使った。ソース／wheelで生成した記録はバイト一致する。公開APIの検証記録であり、インターネット上のReleaseへ配布したという意味ではない。

## 別MPIプロセスでの再開

公開molecular_feedbackサンプルを、真空／LD環境と4／6ランクで実行した。各条件は通常の連続実行、性能モードの連続実行、性能モードでの32ステップ保存、MPI終了後の性能モードでの再開の4つ。ソース16・wheel16回のローカルMPI起動を行った。

本計算FDTDと自己真空参照FDTDを32、分子を31ステップ終了時に保存した。復元中は場と分子の時間発展メソッドを呼ぶと試験を失敗させる仕組みを入れ、全ランクで呼び出し0回だった。Fortran内の実番号を照合し、分子ワーカーの実測更新回数は復元直後0回、64ステップまでの続行後は新しい32回だけだった。復元先の性能モードも保持した。

全フレームの双極子・電子数・エネルギー・場・電流・エネルギー収支観測と、両FDTDの受理済み観測ファイルが連続実行と完全一致した。ソースとwheelも一致した。解析用履歴は再計算に使用していない。

## 保存と範囲

記録はcalc2026/salmonpy-molecular-performance-disk-restart-20261123。公開APIに渡せる実測記録はruntime/report.json。時間発展後の一括チェックポイントは各条件のcheckpointに残した。新規GS計算・GSコピーは行っていない。試験時だけGSコピーを読取リンクへ置換し、製品のコピー処理は変更しない。ローカル保存を確認し、クラウド同期完了は未確認。

同じ実行ファイル・MPI分割でのMac上の短時間C2H2接続試験であり、物理的収束・速度向上・他環境・多ノードの証拠ではない。両スイートは709件中708成功・任意スキップ1。重いCIの手動起動なし。

```json
{
  "schema": "salmonpy.public-molecular-performance-disk-restart/1",
  "runtime_record": "runtime/report.json",
  "executable_sha256": "15c9cb71230289c9b3743208dfdd7a5c2b69e2e24ff2b08e6d3ff4b5b08e5f9e",
  "execution_modes": [
    "checked",
    "performance"
  ],
  "public_reference_error": 0.0,
  "public_energy_error": 0.0,
  "source_jobs": 16,
  "wheel_jobs": 16,
  "conditions": [
    {
      "kind": "source",
      "environment": "vacuum",
      "ranks": 4,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "source",
      "environment": "vacuum",
      "ranks": 6,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "source",
      "environment": "ld",
      "ranks": 4,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "source",
      "environment": "ld",
      "ranks": 6,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "wheel",
      "environment": "vacuum",
      "ranks": 4,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "wheel",
      "environment": "vacuum",
      "ranks": 6,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "wheel",
      "environment": "ld",
      "ranks": 4,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    },
    {
      "kind": "wheel",
      "environment": "ld",
      "ranks": 6,
      "frames_exact": true,
      "field_steps": [
        32,
        32
      ],
      "molecular_step": 31,
      "restore_updates": 0,
      "continuation_updates": 32,
      "mode": "performance",
      "forbidden_restore_calls": 0
    }
  ],
  "field_steps": 64,
  "molecular_steps": 63,
  "GS_calculated": false,
  "GS_copied": false,
  "excited_state_checkpoints_retained": true,
  "tests_each": {
    "total": 709,
    "passed": 708,
    "skipped": 1
  },
  "scope": "macOS arm64 single-machine fixed-ion PZ C2H2 short-run integration; same executable and MPI layout; no physical convergence or speed claim",
  "wheel_sha256": "1c17ab5b86dd348db95baf0566d89a82487c243c2748aa894d04f5f8eae6c001"
}
```
