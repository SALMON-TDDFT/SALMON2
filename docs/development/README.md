# SALMONpy development

Run build/test commands from the repository root. For using the package, start with the [user guide](../../salmonpy/README.en.md).

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
[packaging](../../salmonpy/PACKAGING.md) for compilers, BLAS and preparing wheels with their required shared libraries.
Do not commit build products, ground states or visualization datasets.

## Migration provenance

Moved from `SALMON-TDDFT/SALMON2:salmonpy-field-components` on 2026-11-14 at
`8006261b9e48541b591e496a813b3f06ca5fc1b1`, retaining the full 5,282-commit history.
The Si starter milestones include `78e06940` and `811bbe1c`; see [NOTICE](../../NOTICE).
A shallow clone hides earlier history. `runtime/source.json` deliberately retains
the pinned upstream SALMON source URL; it is not a distribution link.

The calculation sources are maintained in `runtime/sources/solver` for standard exchange-correlation functionals and fast PR HSE. See [source layout and build instructions](../../runtime/README.md).

Implementation notes: [日本語](hse-symmetry-notes.ja.md) / [English](hse-symmetry-notes.en.md).

Additional implementation and distribution notes (historical context):

- GS model inputs: [日本語](gs-model-input-notes.ja.md) / [English](gs-model-input-notes.en.md)
- Si starter distributions: [日本語](si-intro-distribution-notes.ja.md) / [English](si-intro-distribution-notes.en.md)

- Process/output implementation notes: [日本語](output-process-notes.ja.md).
