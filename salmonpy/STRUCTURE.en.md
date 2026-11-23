# Structures and maximum grid spacing

[日本語](STRUCTURE.md)

`salmonpy.structure` provides `StructureInput`, `GridSelection` and
`grid_from_spacing`. ASE is optional, required only when reading/copying structures.
No network, solver execution, repetition, cell standardization or coordinate
wrapping is added. Use separately prepared supercell files.

```python
from salmonpy.structure import StructureInput, grid_from_spacing
structure = StructureInput.from_file("supercell.xyz", format="extxyz", index=0)
atoms = structure.atoms  # independent copy, Å; modifications do not affect source
selection = grid_from_spacing(atoms.cell.array, spacing=0.4, unit="angstrom")
print(selection.as_dict())
print(structure.record())
```

The numeric width is an example, not a converged material recommendation.
`StructureInput.from_atoms(atoms)` is the in-memory route and also copies geometry.
`from_file` requires an explicit `cif` or `extxyz` format and nonnegative integer
frame index (default 0, not ASE's default last frame). File bytes are hashed and
parsed from the same snapshot; filename `@` has no indexing meaning. ASE's MPI
broadcast is disabled. CIF symmetry expansion is part of parsing, so asymmetric
unit site count need not equal returned atom count. The helper adds no supercell
expansion. Use extxyz for explicitly enumerated cell/positions. Readers may raise
ASE/format or OS exceptions; no fallback format or guessed units is attempted.

`record()` returns a detached summary with cell in Å, periodicity, atom count and,
for files, basename, format, index and SHA256 of the entire file. It is not a full
geometry archive; retain the referenced file. `atoms` preserves the parsed cell,
positions, species and ASE settings in a fresh copy. The loader does not assert
DFT eligibility: existing charged/magnetic/nonperiodic DFT restrictions still apply.

## Grid contract

Supply a scalar maximum spacing or exactly three widths and explicit `angstrom`
or `bohr` units. Cell vectors are always supplied in Å. SALMON's Bohr conversion
is 0.52917721067 Å. For each axis, start with ceil(L/h), enforce minimum (default
8) and round upward to an explicitly specified multiple (default 1). Actual width
is L/N and never exceeds the requested width. The input cell is unchanged.
`GridSelection` is returned only by the helper and contains immutable tuples;
`as_dict()` gives requested/actual widths in Å, counts, minimum, multiples and total
points. A 5×6×7 Å cell at maximum 0.4 Å gives (13,15,18).

Initial support is positive axis-aligned orthogonal cells. Nonorthogonal/rotated
cells are rejected for width resolution; no hidden rotation is performed. The
standalone helper receives only a cell, so the caller must check periodicity.
DFT's ASE entry checks three periodic axes. Counts are limited to 2**20 per axis;
this is an input bound, not a memory-budget check. Explicit multiples are user/
calculation implementation constraints, never inferred from rank count. For a calculation implementation needing special
alignment, resolve explicitly and pass `selection.counts` as `grid`.

## Existing DFT / TDDFT entry points

Both `salmonpy.ase.dft_from_atoms` and `tddft_from_atoms` now accept
`grid_spacing=...` plus `grid_spacing_unit=...` instead of `grid=(...)`.
Specifying both, neither, or an unrelated spacing unit is rejected. Existing
explicit-count calls keep their behavior, including their existing geometry scope.
`states`, `kmesh`, pseudopotentials and XC selection remain separate; no physical
convergence or runtime capability is inferred. The DFT model stores adopted counts;
use the grid helper before creation when you need its detailed report. These new
keywords do not extend `SalmonCalculator`; use the helper's counts there.

## Offline material-ID example

```sh
python -m salmonpy.examples.prepare_structure supercell.xyz new-inputs --format extxyz --material-id film --spacing 0.4 --unit angstrom --kmesh 2 2 2
python -m salmonpy.preparation new-inputs
```

The example registers a structure/grid resolver with the common material registry,
stages the supplied file, and preserves material ID, source hash, selected frame,
cell, adopted grid and explicit kmesh in resolved.json. Change settings per material
in Python using the same resolver/registry and preparation APIs. The sample uses
one selected file per invocation; it is not a DFT calculation implementation registration and does
not choose potentials, run GS or certify complete scientific input. It requires
ASE already available offline. Outputs are outside the input bundle.

Mixed/partial CIF or ASE occupancy is rejected at structure loading and the DFT boundary; supply an explicit ordered structure instead of silently choosing the dominant species or filling vacancies. Nested ASE info metadata are also copied independently.
