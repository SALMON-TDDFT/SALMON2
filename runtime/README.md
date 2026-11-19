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

## Source unification: shared SBE, FDTD and LD stage

Current explicit-routing fast-PR preparations reuse the standard SBE type/velocity
kernels/wrappers and the FDTD/LD/EH files through complete-file inventories.
Native/solver kernels are byte-checked before approval; no third source copy or
patch is added. Reference-only snapshots remain unchanged.

The same source builds with HSE ON/OFF. OFF needs no HSE FFTW specification or
linkage. `build_fast_hse_mpi.py` now emits two distinct records:

| Record | Use |
|---|---|
| build.json | Host library for HSE material groups and, when field_abi is declared, NativeEHField field groups |
| executable-build.json | Audited executable for SalmonExecutable.from_mpi_build and batch MPI FDTD |

field_abi protocol1 declares salmonpy.eh-centered-current/1 and is checked against
the compiled checkpoint protocol during explicit build. Legacy embedded field
records remain accepted. Existing prebuilt runtimes do not acquire the capability
from a Python upgrade. Same-build field/material groups still run in separate
processes; this adds no concurrent same-process field/electronic guarantee.

The default standard GS/TDDFT preparation remains. Current counts include15
standard-only paths and71 differing paths. The five standard alternative-HSE
files are not migration targets:only fast-PR HSE is supported. Shared bridge/output
code and equivalent audit guards should be organized, not counted as missing
physics features. Retained GS/input/standard TDDFT initialization still needs porting.

## 日本語

HSEは高速PR版を本筋とし、対称性を中心に開発・検証します。対応範囲は[対称性ガイド](../salmonpy/FAST_HSE_SYMMETRY.md)を参照してください。

本体ソースをリポジトリとwheelに同梱し、SALMONpyと同じ版で更新・検証します。
固定した本体を検証して外部の新規ディレクトリへコピーし、段階別の完成形を配置します。
外部Gitチェックアウトと互換パッチ適用は不要です。HSEは高速PR版だけを保守します。
旧patch_sha256は由来情報であり、現在のファイルの識別には別のハッシュを使います。
ライセンス・元コミットを保持し、import時にはコンパイルしません。

### SBE・FDTD・LD共有化の段階

現在のexplicit-routing高速PR版は、通常版・wheelと同じSBEカーネル、およびFDTD／LD／
E/Hの既存ファイルを使います。完成形の一覧から参照し、第三のソース複製やパッチは
追加しません。参照専用ソースは保持し、import時のコンパイルも増やしません。

`build_fast_hse_mpi.py` が出す `build.json` はHSEと、field_abi宣言がある場合の
NativeEHFieldに使います。別の `executable-build.json` はMPIバッチFDTDの
SalmonExecutable.from_mpi_buildへ渡します。ライブラリと実行ファイルは別のハッシュです。
同じビルドを場・物質の別MPI群で使えますが、同一プロセスでの同時利用を保証するもの
ではありません。旧embedded記録も使え、配布済みバイナリは自動更新しません。

同じソースのHSE有効／無効ビルドを検証しました。無効側はHSE用FFTW指定・リンクが
不要です。標準GS／TDDFTの既定入口はまだ切り替えません。固有15パスのうち標準側の
代替HSE5ファイルは移植対象外で、HSEは高速PR版に限定します。残るGS／入力／通常TDDFT
の初期化・所有権と、通信・出力・監査コードの整理を進めます。

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

### Shared initialization input and output infrastructure

Fast preparation also compiles the existing standard `io/model_inputs.f90` and
`io/step_output.f90` through the shared IO CMake list. The duplicate output helper
in `salmonpy/native_adapters` is removed. MPI/serial host builders reuse the solver
object; when linking an older prepared build without that object, they compile
this same maintained source explicitly and include it in the object audit.

The ordinary input record is a value-only prerequisite for retained GS/TDDFT
migration. It does **not** replace fast HSE's larger `hse_model_inputs` snapshot,
and does not make ordinary embedded GS/TDDFT available in the fast build yet.
No public lifecycle or numerical update order changes in this stage.

### Fast model-owned parser

Fast inputoutput now provides the standard s_input_parser workspace contract.
Explicit model registration reads through a local isolated parser, and captures
its full fast settings separately as before. Existing batch reads retain named
scratch files for the MD reader. Root-owned scratch units stay local; coordinate
flags and diagnostics are broadcast so all ranks see the same parser metadata.
The complete fast default block, unit initialization, KY resolution and input
validation are preserved (only trailing whitespace normalized in moved defaults).
Do not replace the complete HSE model snapshot with the common parser value view.

Validation entry: tools/test_parser_ownership.py BUILD OUTPUT --fast-hse
--reference-build AUDITED_OLD_BUILD; it compares every variables.log entry with
fresh legacy parses from an unchanged library as well as owned/interleaved reads.
No new GS/SCF calculation is needed for this input-only test.

### Internal retained fast GS import

The fast GS context can import compatible existing wavefunction restarts without
entering SCF/warmup/output loops. gs_export_seed makes value-only detached copies
through the standard gs_seed module now shared in the fast CMake list. The common
seed does not export HSE plans/owners or complete fast model settings. Internal
GS import retains target-HSE phase, context-owned communicators and cache cleanup;
its probe omits unconverged energy/force/residual fields. Legacy SCF iteration and
output body is unchanged. This is not yet a Python NativeGSBackend fast profile.
The no-SCF extension to tools/test_gs_lifecycle.py reuses the existing native probe
under a compile flag, retaining older builds' legacy probe compilation.

### Python fast GS handles

A compiled sp_gs_protocol probe now advertises gs_abi protocol1/import mode and
its three arrays in the fast host build record. NativeGSBackend.load reuses the
existing GS Python API, preserves material IDs and returns density/occupations/
orbitals. Per-import symmetry and density-symmetry scratch are switched/restored;
sym.dat assets are carried for reduced meshes. Native monotonic tokens invalidate
closed owners, and matching calls are collective within the host subgroup.

Fast imported metadata has None for energy/residual/tolerance/convergence; native
placeholder values are not physical zeros. SCF creation, computed eigenvalues/
forces, portable export and GS→TDDFT initialization remain standard-profile
features, not fast GS features yet. The retained_gs sample has an import-only
mode rather than a new sample file. See NATIVE_GS guides for the Python example.
