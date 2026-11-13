# Full solver runtime

LD/SBE wheels already include their compiled library. DFT, TDDFT, FDTD and band
calculations use a **linked executable plus its runtime libraries**. An object
file alone would still require a linker and dependencies on the user's machine.
Importing SALMONpy never downloads or builds anything.

## Prefer an existing or prebuilt runtime

First establish the user's memory limit, CPU budget and desired workflow. Check
OS/architecture and required features; LD wheel compatibility does not establish
full-runtime compatibility. Reuse a previously audited executable via
`SalmonExecutable(path, recorded_sha256)`, or verify a managed runtime:

```python
from salmonpy.runtime import load_runtime
engine = load_runtime('my-runtime', required=['tddft', 'libxc'])
# Pass executable=engine to the existing DFT/TDDFT/band loaders.
# For FDTD: engine.run(calculation, 'new-output', threads=2).
```

For public releases, use the trusted publisher's catalog URL:

```sh
python -m salmonpy.runtime status --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime install my-runtime --catalog CATALOG_URL --require tddft libxc
python -m salmonpy.runtime verify my-runtime --require tddft libxc
```

`CATALOG_URL` is an explicit placeholder, not a currently published URL. The
bundled catalog is currently empty. Public hosting awaits the independent
SALMONpy repository; do not invent a download URL or claim an asset is published.
Multiple matches require `--artifact ID`. Destinations must be new. A missing
match raises `RuntimeUnavailable` and prints the explicit build entry point;
it never silently drops Libxc, changes architecture or starts compilation.
Downloads use HTTPS and SHA256. Archive paths and all installed files are checked.
Catalogs must come from a trusted publisher: their hashes are not signatures.

## Current branch artifacts

The **Full runtime bundles** Actions workflow builds macOS arm64, macOS x86_64
and Linux x86_64 profiles and uploads a bundle only after relocation and native
TDDFT tests succeed. A configured runner is not evidence of a passing run. These
profiles use a single MPI rank and PZ/TB-mBJ; they do not include Libxc. Windows
full-runtime bundles are not yet provided (Windows LD/SBE wheels are separate).
The locally tested macOS arm64 Libxc bundle is a separate artifact/profile.

Download the matching successful run's `salmonpy-runtime-...` artifact from the
current branch's Actions page. Extract the **outer Actions ZIP**, retaining the
inner runtime ZIP and `catalog.json`. Then:

```sh
python -m salmonpy.runtime status --catalog downloaded/catalog.json --require tddft
python -m salmonpy.runtime install my-runtime --catalog downloaded/catalog.json --archive downloaded/RUNTIME_ZIP --require tddft
```

`RUNTIME_ZIP` means the actual ZIP filename inside that artifact. CI catalogs
have no public asset URL; installation requires this local archive. No compiler
or linker is needed. OS versions older than the recorded build baseline are
rejected conservatively. Native calculations still need compatible inputs and
pseudopotentials, which are not part of the executable bundle.

## Explicit fallback build

Only offer a build when no appropriate prebuilt/existing runtime is available.
Explain required tools and confirm the user's CPU/memory budget. Obtain the
pinned upstream source, then explicitly build in a new directory:

```sh
git clone https://github.com/SALMON-TDDFT/SALMON2.git UPSTREAM
python -m salmonpy.runtime build UPSTREAM NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

These are GNU compiler options, not universal platform flags. Supply the actual
GNU C compiler path on macOS; `/usr/bin/gcc` is normally Apple Clang. Install
CMake >=3.22, Git, GNU Fortran/C and LP64 OpenBLAS beforehand. The command never installs
system dependencies. Windows native full-solver support is unvalidated; use a
supported Linux environment/WSL or a separately validated site build.
`--jobs` controls build concurrency, not simulation threads or a memory limit.
For Libxc supply `--cmake-arg=-DUSE_LIBXC=ON` and
`--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX` for a compatible Fortran-enabled Libxc.

The wheel includes the pinned patch, shared-kernel references and audit scripts,
so no SALMONpy checkout is needed. Source preparation, compilation and object
checks are logged in `NEW_BUILD/build.log`; `build.json` records status, exact
commands and the resulting executable SHA256. Failures preserve the log and
source; retry in a new directory. Success means **built and object-audited**,
not independently validated scientific accuracy. Run the relevant native
regression before production and use the recorded hash:

```python
import json
from pathlib import Path
from salmonpy.fdtd import SalmonExecutable
record = json.loads(Path('NEW_BUILD/build.json').read_text())
if record['status'] != 'built-and-audited':
    raise RuntimeError('runtime build did not complete')
engine = SalmonExecutable(record['executable'], record['sha256'])
```

MPI runtimes must match the cluster's MPI stack and remain explicit site builds.
The current Python adapters' single-rank scope is unchanged. This acquisition
layer does not alter the DFT/TDDFT update order, state isolation or A→J interface.
