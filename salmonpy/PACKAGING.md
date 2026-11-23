# LD and SBE wheels

Full-runtime setup: prefer a compatible prebuilt binary, with an explicit source-build fallback. See [runtime guide](RUNTIME.en.md). Branch CI artifacts and public releases are distinct.

> Solver sources are bundled under `runtime/sources` (inside `_runtime_source`
> in wheels). Prepare a new external build directory using the current runtime
> guide. No upstream checkout or compatibility-patch application is required.
> Historical validation commands below describe earlier layouts.

The package uses scikit-build-core and a separate CMake project. LD and velocity-gauge SBE calculation routines, C ABI adapters and Python APIs are installed.
The SALMON executable and MPI are not required. SBE links system Accelerate on
macOS and bundled OpenBLAS on Linux/Windows. Source builds on Linux/Windows
require an LP64 OpenBLAS development installation. Fortran/C libraries are built
for each platform, with the same Python interface.

The table below is verified for 0.1.0a4, including the Python TDDFT adapter.
[CI run 43363121675](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43363121675)
passed all four platform jobs and the source-distribution build at `e837ebac`.
Each of 20 installed wheels (CPython 3.10–3.14) passed 34 tests (680 total).
Windows OpenBLAS discovery uses CMAKE_PREFIX_PATH pointing to UCRT64.

The separately built TDDFT worker executable is not included in these wheels.
Its solver regression passed on Linux CI and local macOS arm64; Windows and
Intel macOS TDDFT called directly from Python execution are not claimed by wheel test results.

## Target matrix

| Platform | Architecture | CPython | Status |
| --- | --- | --- | --- |
| macOS | arm64 | 3.10–3.14 | CI passed (macOS 15) |
| macOS | x86_64 | 3.10–3.14 | CI passed |
| Linux manylinux_2_28 | x86_64 | 3.10–3.14 | CI passed |
| Windows | AMD64 | 3.10–3.14 | CI passed |

The CI macOS deployment target is 15.0; dependent runtimes can impose additional
requirements, which the repair/test stage must validate. Linux arm64, Windows
arm64, musl, free-threaded Python and PyPy are not initial release targets.

## Build and install

With a Fortran compiler, CMake and Python available:

```sh
python -m pip install build
python -m build
```

This builds an sdist and then builds the wheel from that sdist. For compiler
selection set FC before configuring. Source builds require the compiler's `nm`
tool for mandatory application-object inspection (Windows CI uses MSYS2 UCRT64
GNU tools). Existing SALMON build configurations are unaffected.

A raw wheel may still depend on compiler libraries on the build computer. Before
redistribution, repair it using the appropriate tool installed in the build
environment: delocate on macOS, auditwheel on Linux, delvewheel on Windows.

```sh
python tools/repair_ld_wheel.py dist/<wheel-name>.whl wheelhouse
python -m pip install wheelhouse/<repaired-wheel-name>.whl
```

The repair step bundles dependencies, refreshes the LD checksum after relocation,
and updates the wheel RECORD hashes. The checksum is an integrity check, not a
signature. Windows repair explicitly analyzes the ctypes-loaded DLL. Compiler
runtime license notices must accompany redistributed dependencies according to
their licenses; review release artifacts before public distribution.

## CI and testing

`.github/workflows/ld-wheels.yml` builds a checked source distribution once.
On PR/manual runs, all four platform jobs build, repair and install CPython
3.10–3.14 wheels from that same archive. Every interpreter runs the public API,
LD/SBE and AI-guide tests (`--suite core`); CPython3.12 runs the full suite with
ASE/plot extras once per OS. A separate CPython3.12 NumPy-only environment runs
core tests. The runner rejects source-checkout imports and empty selections.
Ordinary pushes keep the Linux sdist-to-wheel check and full suite. CI stores
artifacts; it does not publish to PyPI.

Local verification on 2026-11-13: gfortran 15, CPython 3.14, macOS arm64; built
from sdist, repaired with delocate, installed in a separate venv, and all 13 API
and numerical tests passed outside the source checkout. Wheel-record repair
checks cover .dll/.so/.dylib names. Windows AMD64, Linux x86_64 and both macOS architectures subsequently passed
all CPython 3.10–3.14 builds and installed-wheel tests in run 43354726041.

The wheel includes AI_USAGE.md, AGENTS.md, INPUT_CONTRACT.md, the manual search,
environment diagnosis and examples. Tests and build tools remain in the source distribution. Restart the interpreter after rebuilding
solver code. No upload/public release has been performed.

Verified CI: https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43354726041
Code commit: 887391751b6918c7edc6ac819e49fc07dc82202a.
20 wheels, each with 13 passing tests; source-distribution rebuild also passed.

Version 0.1.0a2 adds non-LD Python workflow tools. Its final code commit
045f34faef532ceb106ea00274319a70d27eca8d passed the same 20-wheel matrix,
27 tests per wheel (540 checks):
https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43357199430
Full FDTD executable integration is a separate local test, not part of the wheel CI.


## Local Mac full GS-to-TDDFT validation (2026-11-13)

The full GNU Fortran SALMON executable must use LP64 OpenBLAS on the tested
Apple Silicon Mac. Accelerate passed the LD/SBE and synthetic TDDFT checks but
its ZDOTC complex-return interface crashed in a minimal GNU Fortran probe and
in full GS initialization. The same probe returned 30+0i with OpenBLAS.
Do not infer full-GS compatibility from the smaller wheel calculation routines.

The following paths are the Apple Silicon validation calling program’s paths, not portable defaults.
For another calling program, locate its compiler and LP64 OpenBLAS installation explicitly.
For a new full application build on that calling program, use CMake LAPACK discovery with
`-DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas -DBLA_VENDOR=OpenBLAS
-DBLA_SIZEOF_INTEGER=4`, and leave LAPACK_VENDOR_FLAGS unset (no Accelerate).
The verified compiler is `/opt/homebrew/bin/gfortran`; use GNU C for the OpenMP
build. Run both source and complete application-object audits before execution.
Keep `OPENBLAS_NUM_THREADS=1` when running OpenMP SALMON to avoid nested threading.

With this configuration a real Si GS converged in 37 iterations, and the
1000-step TDDFT and SBE Python-interface versus Fortran comparisons passed. See the local-run
section of `docs/tddft-validation-2026-11-13.md` for numerical evidence.

## Optional Libxc for the full TDDFT executable (Mac validation)

The external SALMON worker now accepts `pz`, `tbmbj`, `libxc_pz`, and `libxc_pbe`.
The two Libxc choices require `USE_LIBXC=ON`. They are not bundled into the
LD/SBE wheel, and this does not enable arbitrary Libxc functionals in the worker.
The installed Mac Libxc lacked `libxcf90`/`xc_f90_lib_m.mod`; the tested compatible
installation is Libxc **5.2.3**, built separately with the same GNU Fortran compiler.
Do not combine that f90 module with a different Libxc library installation.

After obtaining and extracting the official Libxc 5.2.3 source archive:

The following Python driver avoids shell-specific line continuation. Run it in
an explicitly selected build workspace with the extracted source at `libxc-5.2.3`.
Set `compiler` to the actual compiler path if it is not available on PATH.
This is a portable invocation pattern, not verification of full SALMON on every OS.

```python
from pathlib import Path
import shutil
import subprocess

compiler = shutil.which('gfortran')
if compiler is None:
    raise RuntimeError('Select an installed GNU Fortran compiler')
source = Path('libxc-5.2.3').resolve()
build = Path('libxc-build').resolve()
prefix = Path('libxc-install').resolve()
subprocess.run(['cmake', '-S', str(source), '-B', str(build),
    '-DCMAKE_POLICY_VERSION_MINIMUM=3.5',
    f'-DCMAKE_INSTALL_PREFIX={prefix}', f'-DCMAKE_Fortran_COMPILER={compiler}',
    '-DENABLE_FORTRAN=ON', '-DENABLE_XHOST=OFF', '-DBUILD_TESTING=OFF',
    '-DCMAKE_BUILD_TYPE=Release'], check=True)
subprocess.run(['cmake', '--build', str(build), '--config', 'Release', '--parallel', '4'], check=True)
subprocess.run(['cmake', '--install', str(build), '--config', 'Release'], check=True)
```

Add `-DUSE_LIBXC=ON` and `-DLIBXC_INSTALLDIR=<absolute prefix selected above>` to the
full OpenBLAS-backed SALMON configuration above. Check configuration reports
`Found Libxc version 5.2.3`, then audit source and all application objects before
running. The source archive SHA256 used here was
`851a45aee9ddaafea49f684fe3e3ab2fbd79f1c1289b8c1cea216330b27e887b`.
See `docs/libxc-validation-2026-11-13.md` for actual tests and limitations.

solver TB-mBJ does not need Libxc. Import-time compilation is not used.
The optional full-runtime wheel and capability-selection flow in
`docs/plans/2026-11-13-runtime-distribution-proposal.md` are proposals, not
currently available install commands or APIs.

## Current CI scope and independent distribution

The manual `runtime-bundles.yml` workflow owns both runtime distribution and
ordinary/interactive TDDFT comparison. Choose `scope=worker` for one Linux serial
build and the ordinary/predictor-corrector worker checks; choose `scope=bundle`
for relocation/distribution and broader integration. The duplicate standalone
`tddft-worker.yml` workflow is removed. These comparisons use the same audited
current source,not an unavailable historical baseline. They establish adapter
agreement,not historical solver regression. Worker diagnostics upload only logs,
not generated GS fixture directories.

Wheel extras are installed only for each OS's CPython3.12 full suite; other
versions still exercise compiled calculation routines and public/sample syntax
without optional extras. Required Python dependencies remain NumPy only. These
are configured checks; historical runs do not certify the changed configuration.

SALMONpy is intended for a separate repository, not a merge into the SALMON main
repository. solver source provenance and full-executable compatibility must be
preserved when moving it. The three LD/SBE calculation routines live in native/. Full solver sources and verified completed replacements are bundled under runtime/ and copied into a new external build directory.

solver checksums are generated after CMake installs the shared library, including
its RPATH/install-name changes. Direct source installation therefore does not
require preparing wheels with their required shared libraries for checksum validity. Portable redistribution still needs
dependency repair; that step refreshes the checksum and RECORD again.
The sdist CI checks and installs the unrepaired source-built wheel as well.

## Matching calculation routines for source tests

The numerical suite needs calculation routines from this source revision, including the
convergence ABI. An older installed `.so`/`.dylib` can load LD/SBE successfully
while lacking `sp_convergence_check`. Build a fresh wheel explicitly and run the
source test files against it from outside the checkout, as in the installed-wheel
checks above. Do not add a pure-Python numerical fallback, compile on import or
silently skip those numerical tests. A missing ABI is an environment/version
mismatch, not evidence that the convergence algorithm failed. Full fast-HSE
runtime builds additionally require FFTW; see runtime/README.md in the checkout.

## Required source-integrity CI

Changes to `runtime/sources/**`, `runtime/source-trees/**`, `packaging/**` or `pyproject.toml` must not use `[skip ci]`. `native-contract.yml` checks source inventories with pytest, exports tracked Git files only, prepares/builds the MPI library from that export, and runs GS-free C/Fortran query and explicit-parser contracts. The query itself is short; preparing its full runtime still costs a build. Ordinary pushes also retain the Linux sdist/wheel job. OS wheel matrices remain on PR/manual triggers. Full runtime campaigns, TDDFT solver comparisons and NLCC physics checks are manual-only; real-material observable reproduction is not required CI. This respects the Actions budget without skipping source-integrity checks. No workflow is manually dispatched by input preparation or imports. CI configuration is not evidence that a run passed; attach the actual run result when available.


## CI consolidation (2026-11-23)

The PR/manual wheel path previously invoked the full suite25 times:20 platform/
Python wheels,four minimal environments and one source-built wheel. It now runs
four full suites,with core checks on every other wheel and minimal environment.
This is an84% reduction in *full-suite invocations*,not an84% prediction of
Actions minutes. All20 wheel builds remain,and the unrepaired Linux wheel still
has its installed-library hash and core tests checked before matrix jobs begin.
A failed source archive blocks downstream platform builds.

`native-contract.yml` remains distinct: it builds the MPI library from a Git-only
export and checks C/Fortran query/parser contracts without GS or physical solver
runs. This does not duplicate the standalone LD/SBE wheel product. Its triggers
also cover source metadata,shared kernels,build/audit helpers and ABI fixtures.
NLCC and runtime physics campaigns remain manual; incompatible serial/MPI/Libxc
build settings are not treated as interchangeable artifacts.

Push and PR checks remain separate because the branch head and proposed merge
can differ. Newer automatic checks replace obsolete runs for the same workflow/
branch or PR; explicit manual runs use unique groups and are not canceled by
this policy. No existing remote job is modified by local validation.

Use `python tools/test_ci_configuration.py` (Python3.12+,PyYAML) and actionlint
for local configuration checks. Build an sdist and install its wheel outside the
checkout,then use `tools/test_installed_ld.py SOURCE --suite core` or `--suite full`.
Do not dispatch expensive physics jobs merely to validate workflow edits. The
new GitHub configuration must still pass a real run before claiming CI success.
