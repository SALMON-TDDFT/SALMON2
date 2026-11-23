# Build and distribute SALMONpy wheels

Users with a matching wheel should start with [README.en.md](README.en.md).
This guide is for building the package and validating distributions.
LD/SBE libraries are included; the separate DFT/TDDFT/FDTD executable is not.
The corresponding solver sources are included for explicit builds described in [RUNTIME.en.md](RUNTIME.en.md).

## Tested wheel configurations

The 2026-11-23 [wheel CI](https://github.com/otobetomohito/SALMONpy/actions/runs/45127091962)
passed at `35a7f6c5`, merged into SALMONpy main as `82194785`.
All 20 CPython/OS wheels were built, repaired, installed and tested from one checked source archive.

| Platform | CPU | CPython | Scope |
|---|---|---|---|
| macOS 15 | arm64 | 3.10–3.14 | LD/SBE wheel and Python interface tests passed |
| macOS 15 | x86_64 | 3.10–3.14 | Same |
| Linux manylinux_2_28 | x86_64 | 3.10–3.14 | Same |
| Windows 2022 | AMD64 | 3.10–3.14 | Same |

The CI macOS deployment target is 15.0. Inspect the actual wheel tags and bundled dependency
requirements; successful CI is not a claim of support for every OS release or machine.
Linux arm64, Windows arm64, musl, free-threaded Python and PyPy are not in this matrix.
Wheel results do not validate full electronic/FDTD solvers, physical observables or site MPI.
The [Linux source/host-contract CI](https://github.com/otobetomohito/SALMONpy/actions/runs/45127091916)
separately built the MPI library from tracked files and passed C/Fortran query and parser checks without GS calculations.

## Build from a source checkout

Use a Fortran/C compiler, CMake and Python with scikit-build-core.
The lightweight package links Accelerate on macOS and LP64 OpenBLAS on Linux/Windows.
The separately built GNU Fortran full solver uses its own documented BLAS configuration;
do not infer its compatibility from the smaller LD/SBE libraries.
Windows source builds need a matching compiler/binutils/LP64 OpenBLAS toolchain; CI uses MSYS2 UCRT64.

```sh
python -m pip install build
python tools/check_forbidden_eigensolvers.py --source-root native
python -m build
```

This builds an sdist and then its wheel. Set compiler paths through the build environment.
CMake retains configure, incremental source and application-object eigensolver checks;
do not bypass them. Source audits inspect only the current build source, not historical snapshots.
The compiler's object inspection tool must be available (`nm` in the GNU configuration).

## Bundle binary dependencies before redistribution

A locally built wheel may still reference libraries installed on the build machine.
Use delocate on macOS, auditwheel on Linux or delvewheel on Windows, then run:

```sh
python tools/repair_ld_wheel.py dist/WHEEL_FILENAME.whl wheelhouse
python -m pip install --only-binary=:all: wheelhouse/WHEEL_FILENAME.whl
```

Replace the filename with an actual wheel. Keep compiler/BLAS dependency licenses with the package;
see `packaging/ld/licenses` and LICENSE.THIRD-PARTY. The repair tool checks ctypes-loaded dependencies,
refreshes installed-library SHA256 records after relocation, and updates wheel RECORD hashes.
CMake also writes checksums after its install-time binary changes: source installation does not
require repair merely to make hashes match. Repair is still needed for portable redistribution.
Hashes check integrity, not publisher authenticity.

## Verify the installed product outside the checkout

```sh
python -m salmonpy.quickstart
python "SOURCE_DIRECTORY/tools/test_installed_ld.py" "SOURCE_DIRECTORY" --suite full
```

Replace `SOURCE_DIRECTORY` with the source directory’s absolute path. Run from outside it using the target installed environment.
The runner rejects importing the source checkout. `--suite core` checks public APIs, LD/SBE,
manual/example availability and the first-use entry. `--suite full` runs the complete suite;
optional dependency checks are skipped only when the feature is not installed.
Install task-specific extras when validating them. Restart Python after replacing compiled libraries.
A missing numerical symbol such as `sp_convergence_check` means the library is not from the matching
source revision; rebuild explicitly or install a matching wheel. Do not add a Python numerical fallback
or silently skip required kernel tests.

The installed examples, manuals and task catalog are part of the wheel. These build/test helper
commands require the source distribution or checkout. None runs during ordinary import.

## Automatic and manual CI

`ld-wheels.yml` builds/checks one sdist. On PR/manual runs all four OS jobs use the same archive:
CPython3.12 runs the full suite with ASE/plot extras once per OS; all remaining interpreters run core.
A separate NumPy-only environment also runs core. Pushes retain one Linux sdist-to-wheel full suite.
The unrepaired source-built wheel is checked before the matrix. A failed archive blocks downstream jobs.

The full-suite invocation count for PR/manual runs changed from 25 to 4; this is not an estimate
of Actions minutes saved. All 20 builds and per-wheel basic checks remain.
New automatic runs replace stale runs for the same workflow/branch or PR; explicit manual runs are not canceled.
Push and PR checks remain separate because their Git trees can differ.

`native-contract.yml` separately checks inventories, a tracked-only export, audited MPI build,
and GS-free C/Fortran query/parser contracts. Keep this required for source, inventory and packaging changes;
such commits must not use `[skip ci]`. This verifies a different product from LD/SBE wheels.
Real-material observable reproduction and convergence campaigns are not required automatic CI.

`runtime-bundles.yml` is manual: `scope=worker` is one Linux serial worker check;
`scope=bundle` performs broader integration and distribution. NLCC validation is also manual.
Do not launch physical campaigns for documentation/configuration edits. Input preparation never dispatches CI.

For configuration-only checks use actionlint and
`python tools/test_ci_configuration.py` (Python3.12+, PyYAML). Guide text is synchronized by
`python tools/sync_manual_guides.py --check`.

## Availability and historical evidence

The package is a development release and has not been published on PyPI.
Actions artifacts and tagged releases are distinct; access follows repository permissions and artifact retention.
Do not advertise a download URL or runtime capability without checking the actual artifact/version.
The fixed 2026-11-14 Si runtime is historical and needs worker-protocol revalidation for the current Python package;
see [SI_INTRO.en.md](SI_INTRO.en.md). Its earlier numerical results do not certify later binaries.

Earlier platform/build/Libxc evidence remains in the Git history and dated records under docs/.
Their machine-specific compiler paths and old patch commands are not current installation instructions.
Current solver sources live in `runtime/sources/solver`; shared-kernel paths are recorded in `runtime/source.json`.
Fast HSE has additional FFTW dependencies; see [runtime/README.md](../runtime/README.md).
