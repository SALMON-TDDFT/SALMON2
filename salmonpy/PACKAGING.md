# LD wheels

The package uses scikit-build-core and a separate CMake project. Only the LD
kernel, C ABI adapter and Python API are installed; the SALMON executable, MPI,
BLAS and ground-state functionality are not required. Native libraries are built
for each platform, with the same Python interface.

## Target matrix

| Platform | Architecture | CPython | Status |
| --- | --- | --- | --- |
| macOS | arm64 | 3.10–3.14 | Local 3.14/macOS 26 wheel verified; CI configured |
| macOS | x86_64 | 3.10–3.14 | CI configured; execution pending |
| Linux manylinux_2_28 | x86_64 | 3.10–3.14 | CI configured; execution pending |
| Windows | AMD64 | 3.10–3.14 | CI configured; execution pending |

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

A raw wheel may still depend on compiler libraries on the build host. Before
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

`.github/workflows/ld-wheels.yml` uses cibuildwheel to build, repair and install
wheels, then runs all numerical and input-lifecycle tests against the installed
package. The runner rejects imports from the source checkout. CI stores wheels
and the sdist as artifacts; it does not publish to PyPI.

Local verification on 2026-11-13: gfortran 15, CPython 3.14, macOS arm64; built
from sdist, repaired with delocate, installed in a separate venv, and all 13 API
and numerical tests passed outside the source checkout. Wheel-record repair
checks cover .dll/.so/.dylib names. Other OS results require actual CI execution;
configuration alone is not evidence of support.

The wheel includes AI_USAGE.md and INPUT_CONTRACT.md. Development examples and
tests live in the source distribution. Restart the interpreter after rebuilding
native code. No upload/public release has been performed.
