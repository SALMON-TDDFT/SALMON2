# HSE build portability and Nk=4^3 tests

## Status

The implementation uses the existing CPU toolchain selection. The Fugaku and
Linux configurations are source-audited, not certified by runs on those hosts:
this development session has no Fujitsu compiler or Linux runtime available.
Apple Silicon / GNU Fortran 15 is the executed platform. Passing its tests does
not replace tests with Fujitsu, Intel or a Linux GNU toolchain.

## Fugaku

The recommended entry point follows the [official SALMON build instructions](https://salmon-tddft.jp/webmanual/current/html/install_and_run.html#build-and-install):

```sh
mkdir build
cd build
python3 ../configure.py --arch=fujitsu-a64fx-ea --enable-scalapack
make -j 8
```

The executable is `build/salmon`. To install it elsewhere, add
`--prefix=/absolute/path/to/install` to configure.py, then run `make install`;
the executable is installed under that prefix's `bin/`. Use CMake 3.14 or later
and a fresh build directory. The existing architecture name resolves through
the compatibility alias to the updated Fugaku toolchain. MPI and HSE are enabled;
`--enable-scalapack` states the distributed-linear-algebra requirement explicitly.
No HSE-specific library flags or separate dependency builds are required when
compatible libraries or source downloads are available. `--enable-libxc` is not
required for native HSE: its C-library dependency is handled automatically.

The configure.py command construction and short architecture-name resolution
are tested locally; this does not execute the Fujitsu compiler.

Direct CMake is also supported. On a Fugaku login host with `mpifrtpx` and `mpifccpx` on PATH and no explicit
compiler/toolchain override, ordinary CMake selects the target toolchain:

```sh
cmake -S . -B build
cmake --build build -j 8
```

The automatic selection defaults to Release, MPI, native HSE and vendor
ScaLAPACK. CMake prints `SALMON: selecting Fugaku MPI/HSE toolchain` and records
`platforms/fugaku.cmake` in `CMakeCache.txt`. `USE_MPI=OFF` also disables the
ScaLAPACK default; explicit `USE_SCALAPACK=OFF` is respected. The old
`platforms/fujitsu-a64fx-ea.cmake` remains a compatibility alias.

If your shell sets CC/FC or a package manager specifies compilers, automatic
selection deliberately leaves those choices alone. For a Fugaku cross build,
use a fresh build directory and explicitly select the toolchain:

```sh
cmake -S . -B build-fugaku \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/platforms/fugaku.cmake"
cmake --build build-fugaku -j 8
```

`SALMON_PLATFORM=generic` disables automatic detection.
`SALMON_PLATFORM=fugaku` requires the two compiler wrappers and rejects
conflicting compiler overrides. An explicit toolchain always takes precedence.
Do not reuse a build directory configured for another compiler/architecture.

The toolchain selects Fujitsu OpenMP and SSL II / ScaLAPACK
(`-Kopenmp -Nfjomplib`, `-SCALAPACK -SSL2BLAMP`). Vendor library selection follows
the [RIKEN/RIST usage guide](https://www.r-ccs.riken.jp/fugaku/docs/workshop/2025/en/seminar_for_fugaku_users_beginner_course_en_202511.pdf).
Dependencies receive the selected compilers and absolute target toolchain.
HSE library probes compile/link and do not execute target programs on the login
host. Libxc native-host optimization is disabled; FFTW itself is serial, with
independent per-worker plans when SALMON uses OpenMP.

Compatible installed FFTW and Libxc are used when detected. Otherwise CMake's
build downloads hash-pinned sources and builds them automatically with the target
toolchain. First-time fallback builds require network access; for offline sites,
provide compatible installed libraries through CMAKE_PREFIX_PATH or the documented
LIBXC_INSTALLDIR / FFTW_INSTALLDIR hints. No manual dependency build is needed when
the downloads are available. Do not run MPI/CTest target executables on the login
host: submit execution in the site's compute allocation.

The new platform-selection tests use fake compiler paths to check configuration
logic only. Fresh generic MPI/HSE and HSE-off builds are tested on Apple Silicon.
Actual Fujitsu compilation, linking against SSL II/ScaLAPACK, and 3D numerical
execution on Fugaku remain to be verified by the user; these are not certified by
local tests.

## Linux

For GNU C/Fortran plus installed MPI, ordinary CMake with `USE_MPI=ON` selects
compatible installed libraries or builds missing dependencies. For Intel oneAPI,
use `platforms/intel-oneapi.cmake` to retain its MPI/MKL/OpenMP settings. No
Apple-only OpenBLAS selection is applied on Linux. Linux x86_64 and AArch64
compiler/runtime tests remain required before claiming support verified there.

## Numerical dependency regression

During the Nk=4^3 test, the auto-built Netlib 3.12.1 path produced wrong ZHEEV
eigenvectors despite correct eigenvalues and INFO=0. A saved 16x16 SCF matrix
reproduced a normalized residual of 0.155; recompiling ZLARF1L without loop
vectorization reduced it to 7.1e-16. OpenBLAS gave 5.6e-16. An independent generated
Hermitian matrix reproduces the failure, without storing the SCF checkpoint.

The build now disables loop vectorization only for the Netlib fallback with
GNU Fortran 15.x on arm64/aarch64. This is a conservative architecture/compiler
guard based on the observed Mac failure; it is not a claim that the same failure
has been reproduced on Linux. The installed vendor-library and SALMON build flags
are unchanged. `hse_lapack_eigenvectors` checks the actual selected library's
residual and orthonormality before the GS fixture starts.

## Si integration tests

Cases 420/421 use a fresh HSE06 Si8 GS, 12^3 real-space grid, full shifted Nk=4^3,
16 occupied states, and a z-polarized impulse of 1e-4 au. GS convergence below
1e-8 is mandatory. Successful GS verification is the fixture required by RT
preparation, which also checks the producer's checkpoint and convergence log.

The 64-step Taylor4+ACE trajectory (dt=.16 au, 10.24 au total) is a CI regression.
It checks finite current/energy, bounded current, post-kick energy drift, and
the generated dielectric/conductivity relation. It does not resolve optical
peaks or establish k-point convergence.

```sh
OMP_NUM_THREADS=1 ctest --test-dir build \
  -R '(420_bulk_Si_hse_gs|421_bulk_Si_hse_rt)' --output-on-failure
```

Selecting only case 421 pulls in the LAPACK and GS fixtures automatically.
Python 3 is needed for verification, but not for building SALMON.

## Executed validation (2026-09-25/26)

Apple Silicon, GNU Fortran 15, AppleClang 17, Open MPI 5.0.9; Release builds.

| Dependency path | MPI ranks × OMP threads | GS | 64-step impulse |
| --- | --- | --- | --- |
| Automatically built Netlib 3.12.1 (workaround), Libxc 5.2.3, FFTW 3.3.10 | 4 × 1 | 110.47 s; residual 3.5022627e-9 | 170.04 s; verification passed |
| Installed OpenBLAS, Libxc 7, FFTW | 4 × 2 | 59.57 s; residual 3.4765297e-9 | 88.46 s; verification passed |

Both GS calculations converged at reported iteration 73. The installed-library
run selected only `421_bulk_Si_hse_rt`; CTest automatically executed all seven
LAPACK/GS/RT checks, with zero failures (150.37 s total). It used
`OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1`; SALMON reported two OpenMP threads.
The bundled run verified GS and RT in separate CTest invocations after correcting
a verifier that counted the trailing k-grid coefficient row as a k point.
These timings involve different libraries and thread counts, so they are not
an isolated OpenMP speedup measurement.

The GS/RT checks cover convergence, 64 k points with normalized weights, finite
outputs, the default HSE screening and Taylor propagator, electron count 32,
nonzero bounded impulse current, post-kick energy drift, and the dielectric/
conductivity identity. The original failing Netlib eigenvector path was reproduced
before applying the workaround. No new Fugaku or Linux execution is claimed.

### frtpxがhse_lcfo_rt.f90でSIGSEGVになる場合（2026-09-27調査中）

利用者のtcsds-1.2.43環境で、直列ビルドでも`hse_lcfo_rt.f90`の翻訳中に
`Compilation abnormally ended due to SIGSEGV`、続いて`flist: Invalid format`が発生。
`-Kfast`を`-O0`へ置換してもstatus=11で再現した。最適化レベル低下による
回避は未成立。現在の情報では、原因となる構文・処理段階は特定していない。

既存のbuildを使い、ソース直下で次を実行する：

```sh
python3 tools/diagnose_frtpx_lcfo.py --build build --compiler mpifrtpx
```

診断はコンパイルのみ。元コードのO0、OpenMP指定を外した対照、すべての手続きの
実行文をstub化した対照、および各手続きの実行文だけを復元した10ケースを逐次実行する。
宣言と手続きのインターフェースは保持する。全生成物と既存moduleのコピーは新規の
一時ディレクトリに隔離し、本番object/moduleの置換・リンク・数値計算は行わない。
各ケースのstatus、コンパイラ版、コマンド、ログを保存する。非ゼロstatusが構文エラーか
SIGSEGVかはログで区別する。一手続きだけで再現しない場合は手続き間の組合せや
module宣言/import処理も候補に残るため、単一の原因を断定しない。

ローカルGNU Fortran 15で13ケースすべての構文・コンパイルを確認した。
富士通コンパイラでの診断結果・回避策は未確認。診断用stubは実行してはならない。
