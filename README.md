# SALMONpy

Use Python to calculate material response, electronic dynamics and light propagation,
or to prepare structures, pulses and analysis for an existing calculation.

- [English: start here](salmonpy/README.en.md) · [日本語：はじめに](salmonpy/README.md)
- [Python examples](salmonpy/examples/README.en.md) · [サンプル集](salmonpy/examples/README.md)
- [Installing the calculation programs](salmonpy/RUNTIME.en.md) · [計算プログラムの準備](salmonpy/RUNTIME.md)
- [Task catalog for AI assistants](salmonpy/manual.json)

The installation wheel (a Python distribution package) includes the compiled LD/SBE
libraries. Install a matching wheel to try these calculations without a Fortran compiler.
DFT, TDDFT, bands and FDTD additionally need a compatible SALMON executable or shared
library: use a prebuilt distribution when available, or explicitly build the bundled sources.
Importing SALMONpy never downloads or compiles software.

For coupled calculations, Python can call FDTD and electronic calculations one step
at a time in the same processes, using MPI groups to share the work. Start with the
basic examples before the [coupling guide](salmonpy/COUPLED.en.md).

This is the independent SALMONpy repository. It includes the solver sources used by
its Python interfaces; an upstream SALMON checkout is not required.

## Repository layout

| Directory | Purpose |
|---|---|
| `salmonpy/` | Python APIs, C ABI wrappers, examples, bilingual guides and tests |
| `native/` | Fortran calculation routines shared with the bundled solver |
| `packaging/` | Wheel build/repair configuration and runtime license notices |
| `runtime/` | Solver sources, file checksums and build configuration |
| `tests/native/` | Numerical reference routines and C/Fortran interface checks |
| `tools/` | SALMONpy packaging, runtime preparation and integration checks |
| `docs/` | Design and validation history; older source paths describe past layouts |

## Developer checks

```sh
python -m salmonpy.build_native
python -m unittest discover -s salmonpy/tests -v
python tools/test_ld_packaging.py
python tools/check_forbidden_eigensolvers.py --source-root native
```

Use `python -m build` for source/wheel distribution; see
[packaging](salmonpy/PACKAGING.md) for compilers, BLAS and preparing wheels with their required shared libraries.
Do not commit build products, ground states or visualization datasets.

SALMONpy is derived from SALMON (Apache-2.0). See [LICENSE](LICENSE),
[NOTICE](NOTICE) and [third-party notices](LICENSE.THIRD-PARTY).

## Migration provenance

Moved from `SALMON-TDDFT/SALMON2:salmonpy-field-components` on 2026-11-14 at
`8006261b9e48541b591e496a813b3f06ca5fc1b1`, retaining the full 5,282-commit history.
The Si starter milestones include `78e06940` and `811bbe1c`; see [NOTICE](NOTICE).
A shallow clone hides earlier history. `runtime/source.json` deliberately retains
the pinned upstream SALMON source URL; it is not a distribution link.

The calculation sources are maintained in `runtime/sources/solver` for standard exchange-correlation functionals and fast PR HSE. See [source layout and build instructions](runtime/README.md).
