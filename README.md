# SALMONpy

Use Python to calculate material response, electronic dynamics and light propagation,
or to prepare structures, pulses and analysis for an existing calculation.

- [English: start here](salmonpy/README.en.md) · [日本語：はじめに](salmonpy/README.md)
- [Python examples](salmonpy/examples/README.en.md) · [サンプル集](salmonpy/examples/README.md)
- [Installing the calculation programs](salmonpy/RUNTIME.en.md) · [計算プログラムの準備](salmonpy/RUNTIME.md)
- [Task catalog for AI assistants](salmonpy/manual.json)

Already using SALMONpy? Choose [a feature to add](salmonpy/examples/README.en.md#add-to-existing) / [既存計算への追加](salmonpy/examples/README.md#add-to-existing).

This is a development release and is not published on PyPI. Obtain a wheel
matching your OS, CPU and Python version; see the language guides above.

With that wheel installed, finish the first short LD/SBE calculations with:

```sh
python -m salmonpy.quickstart
```

The command needs no additional input and prints `SALMONpy first calculations: PASS`
with the results directory. See the language guides above for wheel installation.
For AI assistants: `python -m salmonpy.manual` prints a capability summary and first-test instructions;
search a task afterward instead of reading the entire JSON catalog.

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

For building distributions, source layout and implementation notes, see the
[developer guide](docs/development/README.md).

SALMONpy is derived from SALMON (Apache-2.0). See [LICENSE](LICENSE),
[NOTICE](NOTICE) and [third-party notices](LICENSE.THIRD-PARTY).
