# Optional full SALMON runtime / 外部SALMON実行環境

SALMONpy's wheel contains LD/SBE kernels only. DFT, TDDFT, FDTD and band workflows
need a compatible **external** full executable. The full SALMON source is not
part of this repository's working tree. It is prepared outside the checkout.

`source.json` pins upstream SALMON commit
`3322d0561c0cad5180679136b3eddbe8b909487c` and the SHA256 of
`patches/salmonpy-runtime.patch`. The patch preserves the already validated worker,
GS status/force, band-path, LD/SBE extraction and PZHEEVD changes. It is not a
request to merge these changes into SALMON upstream.

## Prepare and build

Use paths appropriate for your system. UPSTREAM and RUNTIME must be distinct;
RUNTIME must not exist and must be outside UPSTREAM.

```sh
git clone https://github.com/SALMON-TDDFT/SALMON2.git UPSTREAM
python tools/prepare_salmon_runtime.py UPSTREAM RUNTIME
cmake -S RUNTIME -B RUNTIME_BUILD -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DBLA_VENDOR=OpenBLAS -DCMAKE_BUILD_TYPE=Release
cmake --build RUNTIME_BUILD --target salmon --parallel 2
python tools/check_forbidden_eigensolvers.py --object-tree RUNTIME_BUILD/src
```

Preparation archives only the pinned commit (ignoring uncommitted files), checks
and applies the patch, checks that the shared kernels match `native/`, and runs
the forbidden-eigensolver source audit. It never changes the upstream checkout.
It records provenance in the prepared tree. No private Git history is required:
only the explicitly obtained upstream checkout and files in this repository.

Compiler flags, OpenMP/MPI and Libxc must be configured for the intended runtime.
On the tested Apple Silicon GNU Fortran build use LP64 OpenBLAS for full SALMON,
not Accelerate. Libxc functionals require a compatible `USE_LIBXC=ON` build.
The SHA256 used by `SalmonExecutable` must be recorded from the audited executable.
A fresh hash alone does not establish compatibility. Existing validated external
executables remain usable; source cleanup does not change their interface.

The CI workflow prepares this same source and compares ordinary TDDFT execution
with its interactive worker. For local native adapter checks set
`SALMONPY_RUNTIME_SOURCE` to RUNTIME; the tests obtain the Si pseudopotential from
that pinned source. `tools/test_band_path_native.py RUNTIME` tests the prepared
band implementation. SALMONpy does not bundle the general SALMON sample library.

## 日本語

本体の全ソース・機種別設定・一般サンプルはSALMONpyには同梱しません。
固定した公開SALMON版に互換パッチを適用し、別の場所で本体を構築します。
上記の準備コマンドは元チェックアウトを変更せず、新規保存先だけに出力します。
必要な本体変更を削除したわけではなく、パッチと由来情報へ分離しています。
wheelのみを使うLD/SBE・構造・解析では、この準備は不要です。

## Maintenance

When updating a shared kernel, regenerate the corresponding patch hunk and
record its new hash in source.json; preparation rejects mismatched copies.
Updating the upstream commit requires explicitly revalidating all native adapters.
Keep upstream LICENSE, LICENSE.THIRD-PARTY and acknowledgements in the prepared
runtime; do not substitute SALMONpy's narrower third-party notice for them.
