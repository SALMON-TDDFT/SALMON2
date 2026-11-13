# SALMONpy

Python components for material response, electronic dynamics, electromagnetic
simulation and visualization. This repository is being prepared for independent
SALMONpy distribution; it is not a SALMON application source checkout.

- [English guide](salmonpy/README.en.md) · [日本語ガイド](salmonpy/README.md)
- [Python examples](salmonpy/examples/README.en.md) · [サンプル集](salmonpy/examples/README.md)
- [Optional full SALMON runtime](runtime/README.md)
- [AI task catalog](salmonpy/manual.json)

The wheel bundles only the LD/SBE native kernels. DFT, TDDFT, bands and FDTD use
a separately built compatible SALMON executable. No compilation occurs on import.

## Repository layout

| Directory | Purpose |
|---|---|
| `salmonpy/` | Python APIs, C ABI wrappers, examples, bilingual guides and tests |
| `native/` | Three shared LD/SBE Fortran kernel files |
| `packaging/` | Wheel build/repair configuration and runtime license notices |
| `runtime/` | Pinned upstream revision and compatibility patch for optional full SALMON |
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
