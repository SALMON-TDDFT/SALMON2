# Bundled SALMON solver sources / 同梱ソルバーソース

The repository and wheel contain complete standard and fast-PR solver snapshots
under `sources/`, with original licenses and acknowledgements. The only supported
HSE implementation is the fast PR profile. HSE development prioritizes symmetry
and this optimized profile; see [current symmetry scope](../salmonpy/FAST_HSE_SYMMETRY.en.md). `overrides/` holds completed source
files for experimental stages. All copies are bound to inventories in
`source-trees/`; there are no compatibility patch files or git-apply steps.

```sh
python tools/prepare_salmon_runtime.py NEW_SOURCE
python -m salmonpy.runtime build NEW_BUILD --jobs 2 --cmake-arg=-DBLA_VENDOR=OpenBLAS
python tools/prepare_hse_rt_state.py NEW_HSE_SOURCE --explicit-routing
```

Choose explicit compiler/library flags using [the runtime guide](../salmonpy/RUNTIME.en.md).
Preparation does not download, compile, change an external checkout or overwrite
a destination. Build output must be outside this checkout. The historical
UPSTREAM OUTPUT form is accepted but the repository-owned source is used.
The HSE reference-only profile is for numerical comparison, not an alternative
supported implementation. Standard USE_HSE=ON is rejected.

Update source files and inventory digests together. Compare shared LD/SBE kernels
against native/, run source/object prohibition audits and relevant regressions.
The upstream commit and legacy patch_sha256 preserve provenance and existing
binary receipt compatibility; current files are bound by source_tree_sha256 and
replacements_sha256. Numerical changes require renewed capability validation.
No source is compiled on import, and full solver compilation is not part of
ordinary wheel installation. Prefer compatible prebuilt runtimes for users.

## 日本語

HSEは高速PR版を本筋とし、対称性を中心に開発・検証します。対応範囲は[対称性ガイド](../salmonpy/FAST_HSE_SYMMETRY.md)を参照してください。

本体ソースをリポジトリとwheelに同梱し、SALMONpyと同じ版で更新・検証します。
固定した本体を検証して外部の新規ディレクトリへコピーし、段階別の完成形を配置します。
外部Gitチェックアウトと互換パッチ適用は不要です。HSEは高速PR版だけを保守します。
旧patch_sha256は由来情報であり、現在のファイルの識別には別のハッシュを使います。
ライセンス・元コミットを保持し、import時にはコンパイルしません。

## NLCC consistency patch (CPU, unpolarized NC)

The current bundled source aligns Libxc's NLCC energy weighting and GGA
core-density gradients, and adds the explicit force from the sampled core-density
interpolant. Atoms are the outer OpenMP work items. The core force is reduced
over spatial ranks only, before force symmetry, without k/orbital multiplicity.
Invalid scalar/table values and arithmetic overflow are rejected collectively.

Validation tools are `tools/run_nlcc_native_tests.py` with the prescribed-density
programs in `tests/runtime/`, and scalar-only Si probes in
`tools/measure_nlcc_force_work.py`, `tools/test_nlcc_mpi_scf.py` and
`tools/test_nlcc_no_core_regression.py`. Each newly linked test executable audits
all application objects. Solver probes discard temporary eigen/k/info outputs;
only requested scalar reports are exported.

Piecewise-linear radial knots and the existing hard core cutoff need not be
differentiable; tests distinguish interior derivatives from boundary crossings.
Spin-polarized Libxc-NLCC remains unsupported. Meta-GGA force validation is not
claimed, and an OpenACC NLCC-force path explicitly fails as unvalidated.
Existing binary bundles are not changed by this source update. The Calculator
enables NLCC forces only when a matching validation receipt is loaded with
`load_validated_build` or `load_runtime`. The receipt binds the executable and
patch hashes to native checks and a whole-application regression; relocation
requires repeating the latter. Source/object audit alone does not authorize an
old binary's NLCC forces. This capability does not certify physical convergence.
