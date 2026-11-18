# Component regression tests / 部品テストのまとめ実行

Use an existing validated library/build and existing GS bundle. The runner does
not build, download, compute SCF, or start GitHub Actions. Each case has its own
MPI process lifetime; combining these in one Python process would violate native
attachment restrictions. LD/SBE now share their point-MPI test implementation.

既存の検証済みライブラリ・ビルド記録・GSを明示して使います。ビルド、取得、SCF、
GitHub Actionsは実行しません。native接続の寿命制約を守るため、ケースごとにMPIを起動します。
LD/SBEは準備・通信・時系列比較・結果保存の共通部分を一つにまとめました。

```sh
python tools/run_component_tests.py NEW_OUTPUT --suite handles --library LIBRARY
python tools/run_component_tests.py NEW_OUTPUT --suite memory-gs --build BUILD --fixture EXISTING_GS
python tools/run_component_tests.py NEW_OUTPUT --suite handles memory-gs --library LIBRARY --build BUILD --fixture EXISTING_GS
```

Use an environment with NumPy/mpi4py and the package available. `--python` and
`--mpiexec` select executables. `--fixture` can repeat for different XC bundles.
`--dry-run` prints commands without creating output or launching processes.
Tests run sequentially;
OMP/OpenBLAS/MKL threads are1 in child processes. No parent settings are modified.

NumPy/mpi4pyとパッケージを利用できる環境で実行します。必要なら`--python`と
`--mpiexec`を指定します。XCごとに`--fixture`を繰り返せます。`--dry-run`ではコマンドを
表示するだけでファイルやプロセスを作りません。ケースは順次実行し、子プロセスのOMP/BLASは1スレッドです。

| Suite | Cases / ケース | Preserved checks / 保持した検査 |
|---|---|---|
| handles | LD/SBE × MPI2/3 | used-model broadcast, bitwise serial parity; LD finalizers closed; SBE shared native model / 使用済みモデル転送・直列との完全一致・LD解放・SBE共有 |
| memory-gs | each fixture × MPI1/orbital2/spatial2 per material group | file-route parity, independent states, source-close lifetime, interleaved load, direct restore, unavailable import force / 従来経路一致・独立性・GS終了後の寿命・読み込み切替・直接復元・未計算力拒否 |

`NEW_OUTPUT/summary.json` indexes case status, command and elapsed time. Individual
logs and original result records are retained. Failure is recorded, other
selected cases continue, and the overall command exits nonzero if any case fails.
An existing output directory is never overwritten. Large GS files are not copied
into collected reports; the old file-route control still makes its normal local
working copy inside the selected scratch output. Do not collect that working GS
bundle as results; retain only diagnostics/JSON.

`summary.json`に状態・コマンド・時間をまとめ、個別ログと結果JSONも残します。
失敗の後も他の選択ケースを実行し、1件でも失敗すると全体は失敗終了します。
既存出力先は上書きしません。従来経路の比較側は通常どおり作業用GSコピーを作るため、
結果の収集時には作業用GS本体を含めず、診断とJSONだけを保存します。

The redundant LD/SBE-specific scripts have been removed. Use the common entry
with `--model ld` or `--model sbe`.
旧LD/SBE個別スクリプトは削除しました。共通入口でモデルを指定します。

```sh
mpiexec -n 3 python tools/test_material_handles_mpi.py LIBRARY NEW_JSON --model sbe
```

Unit tests with different error cases and numerical/native tests with different
oracles remain separate. The runner consolidates orchestration, not coverage.
単体テストの異常系や、異なる比較対象を持つ数値テストは残しています。検査範囲を減らす整理ではありません。

## Unit-test cleanup / 単体テストの整理

The original637 methods are now620. Six repeated import-side-effect tests are
covered in the public-API import test (including the missing-MPI environment).
Five name-existence-only tests are removed because the public inventory or real
operation checks cover them. Numerical checks and individual failure scenarios
remain. GS link tests are in test_native_gs.py; standard/HSE shared ownership
classes are in test_shared_models.py. Unit test files decrease from89 to83.

637件から620件へ整理しました。重複するimport検査6本は公開API検査へ集約し、
公開一覧や実動作で確認済みの存在確認5本を削除しました。数値・異常系の検査は維持します。
GSリンク検査はGSテストへ、標準TDDFT/HSEの共有所有権検査は同じファイルへまとめ、
単体テストファイルも89本から83本へ減らしました。

The second cleanup groups LD analytic/ownership/input checks in test_ld.py and
SBE reference/ownership checks in test_sbe.py. Superseded tests were deleted only
after mapping their unique assertions to retained checks; this is not simply
merging method names. [Detailed mapping](functional-test-cleanup-2026-11-18.md).

第2段階ではLDとSBEを機能単位のファイルに整理し、包含関係を確認した旧検査を削除しました。
[削除した検査と残した検査の対応](functional-test-cleanup-2026-11-18.md)を記録しています。
