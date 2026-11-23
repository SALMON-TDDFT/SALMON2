# Prepare the calculation programs for DFT, TDDFT and FDTD

[日本語](RUNTIME.md)

The first LD/SBE test needs only a matching SALMONpy wheel:
`python -m salmonpy.quickstart`. DFT, TDDFT, FDTD and bands need an additional
executable or shared library compatible with the selected SALMONpy interface.

## Choose by how the program will be called

| Use | Required program | Guide |
|---|---|---|
| LD, SBE, pulses, geometry and analysis | SALMONpy wheel and task-specific extras | [First steps](README.en.md) |
| DFT/TDDFT/FDTD/bands in a separate process | Compatible executable | This page |
| Electronic and field calculations within one MPI job | Shared library built for the caller's MPI | [Coupling](COUPLED.en.md) |
| Fast-HSE ground states and propagation | Compatible fast-HSE calculation library | [Ground states](NATIVE_GS.en.md), [HSE](SHARED_HSE.en.md) |

Prefer a compatible binary distribution. Build explicitly when one is unavailable.
Importing or checking inputs never downloads, compiles, initializes MPI or starts a calculation.
Wheel compatibility does not establish compatibility of the additional calculation programs.

## Use an existing distribution

Check the OS, CPU, required features and compatibility with the installed Python package.
A standard SALMON executable need not contain SALMONpy's step/state interfaces.
Verify an installed distribution with:

```python
from salmonpy.runtime import load_runtime
engine = load_runtime('my-runtime', required=['tddft'])
```

Pass `engine` to the corresponding DFT/TDDFT/band `executable=` argument.
For Libxc-PZ/PBE, require `['tddft', 'libxc']`; do not replace it with a non-Libxc build.
Hashes and feature labels alone do not certify every input or the worker protocol.

Install from a trusted local catalog and its matching ZIP:

```sh
python -m salmonpy.runtime status --catalog downloaded/catalog.json --require tddft
python -m salmonpy.runtime install my-runtime --catalog downloaded/catalog.json --archive downloaded/RUNTIME_ZIP --require tddft
python -m salmonpy.runtime verify my-runtime --require tddft
```

Replace `RUNTIME_ZIP` with the actual calculation-program ZIP. The destination must be new.
The installer checks hashes, file inventory and ZIP paths. Hashes are not publisher signatures;
obtain the catalog from a trusted distributor. Public HTTPS catalogs are also supported.
For private assets, download through an authorized account before installing locally.

For Actions artifacts, unpack the outer ZIP and use the inner runtime ZIP plus `catalog.json`.
The dedicated `prepare` in the [Si starter](SI_INTRO.en.md) instead takes the outer ZIP.
For the Si starter, use a matching ZIP from its packaged distribution inventory. See [Si starter](SI_INTRO.en.md) for platform availability and obtaining the files.

## For developers preparing distributions

The manual **Runtime bundles and TDDFT integration** workflow, `scope=bundle`, builds
Linux x86_64, Apple Silicon and Intel Mac configurations and checks relocation before uploading.
`scope=worker` is a shorter Linux-only check and does not produce distributions.
Configured jobs are not evidence of a passing run at a particular revision.
The ordinary wheel workflow checks the separate LD/SBE product. This full-runtime workflow
does not currently provide a Windows distribution.

## Build when a compatible binary is unavailable

Sources, inventories and audit tools are included in wheels. No upstream checkout is required.
Prepare the compilers, CMake, BLAS/LAPACK, and any required MPI/Libxc in the build environment.
The following example uses GNU Fortran/C and LP64 OpenBLAS; adapt compiler paths and library locations.

```sh
python -m salmonpy.runtime build NEW_BUILD --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=gfortran \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DBLA_VENDOR=OpenBLAS \
  --cmake-arg=-DOPENMP_FLAGS=-fopenmp \
  --cmake-arg=-DCMAKE_BUILD_TYPE=Release \
  --cmake-arg='-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -ffree-line-length-none -fallow-argument-mismatch'
```

`--jobs` controls build concurrency, not calculation concurrency or a memory limit.
On macOS `/usr/bin/gcc` is usually Apple Clang; specify the actual GNU C compiler for this OpenMP example.
This guide does not establish a validated Windows full-solver build.

- Libxc: prepare its Fortran interface and add `--cmake-arg=-DUSE_LIBXC=ON` and `--cmake-arg=-DLIBXC_INSTALLDIR=LIBXC_PREFIX`.
- MPI: use the site's matching `mpifort` and `mpicc`, adding `--cmake-arg=-DUSE_MPI=ON`. Building MPI does not make a single-rank Python entry internally parallel.
- Fast HSE: add `--cmake-arg=-DUSE_HSE=ON`; the default is OFF. FFTW is required. The build can download a pinned FFTW if no installed copy is usable. Offline builds must provision it first and pass `--cmake-arg=-DFFTW_INSTALLDIR=FFTW_PREFIX`. See [details](../runtime/README.md).

Successful builds record `built-and-audited`, executable location and SHA256 in `build.json`.
Failures retain `build.log` and staged sources. Do not overwrite validated binaries.
A successful build/audit is not a material-convergence or all-features validation.

```python
import json
from pathlib import Path
from salmonpy.fdtd import SalmonExecutable
record = json.loads(Path('NEW_BUILD/build.json').read_text())
if record['status'] != 'built-and-audited':
    raise RuntimeError('runtime build did not complete')
engine = SalmonExecutable(record['executable'], record['sha256'])
```

The single maintained tree is `runtime/sources/solver`. State-ownership changes are already
integrated; do not apply historical extraction patches. Preserve the source and complete
application-object eigensolver audits during builds. Reusing an unchanged validated program
requires its recorded SHA256, not a new startup source/object audit.

## Compute nodes without internet access

1. Identify the target OS, CPU, Python, MPI and site modules.
2. Stage wheels, dependencies, calculation programs, inputs, structures, potentials and initial states on an authorized connected machine. Record provenance, licensing and hashes; keep credentials out of calculation inputs.
3. Transfer while preserving relative paths. Install local wheels with `pip --no-index --find-links`, and calculation programs with local catalog and ZIP files.
4. Run a small example within allocated compute resources. Successful management-node preparation or checking never automatically starts the calculation.

Resolve missing files before calculation; there is no implicit download or build.
Supply pseudopotentials as local files. GS/2 stores their references, SHA256 and provenance;
explicit acquisition is available through `salmonpy.pseudopotentials`; see [downloads](PSEUDOPOTENTIALS.en.md).
See [preparation-only operations](PREPARATION.en.md) and [structure files](STRUCTURE.en.md).

## Forces with NLCC potentials

NLCC forces require validation evidence tied to the binary and supported conditions.
Use `load_runtime(..., required=['nlcc-forces'])`, or
`load_validated_build(build_json, validation_json)` for a validated local build.
See [ASE Calculator](ASE_CALCULATOR.en.md) for verification and packaging instructions.
Implementation checks do not replace grid, k-point, SCF and force convergence studies.

## Reuse verification within one job

For many material points, explicitly reuse the initial complete verification:

```python
import secrets
from salmonpy.runtime import prepare_runtime_job, load_runtime
validation = prepare_runtime_job(
    'my-runtime', 'new-job-validation.json',
    job_id='CURRENT_JOB_ID', job_token=secrets.token_urlsafe(32),
    required=['tddft', 'single-rank'],
)
engine = load_runtime('my-runtime', required=['tddft', 'single-rank'], validation=validation)
print(engine.runtime_validation_report)
```

Neither call executes SALMON, initializes MPI, downloads or builds. Prepare once on the job leader
and distribute within that job. Never log or save the secret in inputs/checkpoints; create a new one
for the next job. Omitting `validation` uses complete verification on every load.
The metadata-based fast path is POSIX-specific; unsupported/unusable cases warn and use full verification.
Every reused load still hashes the executable. It can skip rehashing ancillary files and scanning the tree,
but per-file metadata checks remain. This is not a security guarantee against malicious changes or
shared filesystems reporting unreliable change metadata.

For process-isolated TDDFT, pass `runtime_validation=validation` to both provider and backend.
Custom implementations need not participate. Shared-library loading is a separate operation.
The `mpi_runtime_receipt` example does not initialize MPI or run a calculation without explicit `--run`.
See [memory and measurement](MEMORY.en.md) and [MPI layout](DISTRIBUTED.en.md).

## State ownership and restart

| Goal | Current guide |
|---|---|
| Independent TDDFT states sharing one model | [SHARED_TDDFT.en.md](SHARED_TDDFT.en.md) |
| Different models used sequentially within one MPI group | [TDDFT_REGISTRY.en.md](TDDFT_REGISTRY.en.md) |
| Fast HSE and symmetry | [SHARED_HSE.en.md](SHARED_HSE.en.md), [FAST_HSE_SYMMETRY.en.md](FAST_HSE_SYMMETRY.en.md) |
| Ground states and interrupted SCF | [NATIVE_GS.en.md](NATIVE_GS.en.md) |
| Restart and failed-restore handling | [CHECKPOINTING.en.md](CHECKPOINTING.en.md) |
| External program calls | [API.en.md](API.en.md), [C_ABI.en.md](C_ABI.en.md) |

Sequential state switching is not concurrent threaded execution. Check each guide for supported
functionals, MPI partitions, formats and recoverable failures.
