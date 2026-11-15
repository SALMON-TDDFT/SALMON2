# SALMONpy

Full-runtime setup: prefer a compatible prebuilt binary, with an explicit source-build fallback. See [runtime guide](salmonpy/RUNTIME.en.md). Branch CI artifacts and public releases are distinct.

Python components for material response, electronic dynamics, electromagnetic
simulation and visualization. This is the independent SALMONpy repository; it includes pinned solver sources maintained with the Python integration.

- [English guide](salmonpy/README.en.md) · [日本語ガイド](salmonpy/README.md)
- [Python examples](salmonpy/examples/README.en.md) · [サンプル集](salmonpy/examples/README.md)
- [Optional full SALMON runtime](runtime/README.md)
- [AI task catalog](salmonpy/manual.json)

The wheel compiles LD/SBE native kernels and also carries the full solver sources for explicit builds. DFT, TDDFT, bands and FDTD use
a separately built compatible SALMON executable. No compilation occurs on import.

## Repository layout

| Directory | Purpose |
|---|---|
| `salmonpy/` | Python APIs, C ABI wrappers, examples, bilingual guides and tests |
| `native/` | Three shared LD/SBE Fortran kernel files |
| `packaging/` | Wheel build/repair configuration and runtime license notices |
| `runtime/` | Complete solver sources, hash inventories and stage-specific replacement files |
| `tests/native/` | Independent numerical reference routines and native checks |
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
[packaging](salmonpy/PACKAGING.md) for compilers, BLAS and wheel repair.
Do not commit build products, ground states or visualization datasets.

SALMONpy is derived from SALMON (Apache-2.0). See [LICENSE](LICENSE),
[NOTICE](NOTICE) and [third-party notices](LICENSE.THIRD-PARTY).

## Migration provenance

Moved from `SALMON-TDDFT/SALMON2:salmonpy-field-components` on 2026-11-14 at
`8006261b9e48541b591e496a813b3f06ca5fc1b1`, retaining the full 5,282-commit history.
The Si starter milestones include `78e06940` and `811bbe1c`; see [NOTICE](NOTICE).
A shallow clone hides earlier history. `runtime/source.json` deliberately retains
the pinned upstream SALMON source URL; it is not a distribution link.
