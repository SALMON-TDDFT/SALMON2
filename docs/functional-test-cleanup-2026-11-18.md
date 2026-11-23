# Functional test cleanup / 機能単位のテスト整理

The previous pass reduced637→626 methods. This pass reduces626→620 and87→83
unit-test files. Counts describe test methods, not independent product features.
LD analytic/owned-state/input contracts now live in one file, and SBE reference/
owned-state contracts in one file. Independent test classes preserve failure names.

前段の637→626件に続き、今回は626→620件、87→83ファイルへ整理しました。
数値は製品の機能数ではなく検証メソッド数です。LD・SBEは機能別にまとめています。

| Deleted old check / 削除した旧検査 | Retained check / 残した検査 |
|---|---|
| LD batch/step equality and independence | Native stateless/batch exact parity + public snapshot independent branches |
| LD partial-batch overflow rollback duplicate | Persistent LD full-snapshot rollback + public input/failed-native-call transaction |
| SBE original independence test | Shared-owner independence, same-drive equality, immutable settings + nonfinite rollback |
| TDDFT initial replay mismatch duplicate | Candidate-cleanup failure with accepted current/time/PID preservation; intermediate mismatch separate |
| GS model-source existence only | Read actual packaged input fields and verify capture/activate content |
| Separate parser-type packaging test | Same assertion inside GS input-ownership check |

Unique SBE immutable-settings checks were transferred before deleting the old test.
Duplicate TDDFT startup-failure injection was removed from the intermediate-mismatch
case because the shared-bundle test checks both creation/restore failure and leases.
Parser ownership native round-trip tools/Fortran drivers are retained; these are
actual behavior checks, distinct from the removed unit packaging smoke test.

LD analytic Drude/damping/dielectric tests, SBE archived SALMON reference and all8
correction combinations, snapshot corruption, native finalizer failures, output
retention, checkpoints and MPI integration tests remain. No fixed target count
was used to justify removing a unique error path or physical comparison.

Before/after the selected suite runs49→43 tests. Python line tracing of
ld/sbe/native-handle/snapshot/tddft code records the identical581 executable lines.
This is a limited guard against accidental lost paths, not branch coverage or
physical accuracy evidence. Exact file/line sets are retained in the diagnostic
record. No production solver code, native binary, wheel format or CI was changed.

代表49→43件の実行を比較し、対象実装の実行行581行が一致しました。これは限定的な補助検査で、
分岐網羅や物理精度の証明ではありません。独自の数値参照・異常系は残しています。

Canonical diagnostic record: calc2026/salmonpy-test-consolidation-20261118/functional-cleanup.
Only logs/JSON are collected, no GS data. Local OneDrive write only; sync unverified.
