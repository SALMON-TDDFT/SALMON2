# CIF to SALMON structure converter

`cif2salmon.py` converts one periodic crystal structure in a crystallographic
information file (CIF) into a SALMON input fragment. It expands crystallographic
symmetry and writes all atoms in the unit cell as an `&atomic_red_coor` or
`&atomic_coor` namelist.

The converter is an optional utility. SALMON itself does not depend on Python or
Gemmi.

## Requirements

- Python 3.8 or later
- [Gemmi](https://gemmi.readthedocs.io/)

Install Gemmi into the Python environment used to run the converter:

```sh
python3 -m pip install gemmi
```

## Usage

Write a structure fragment to standard output:

```sh
python3 utility/cif2salmon/cif2salmon.py structure.cif
```

Write it to a new file:

```sh
python3 utility/cif2salmon/cif2salmon.py structure.cif -o structure_salmon.inp
```

An existing output file is not overwritten unless `--force` is specified.

If a CIF contains multiple structure data blocks, select one by the name after
`data_` (both `sample` and `data_sample` are accepted):

```sh
python3 utility/cif2salmon/cif2salmon.py structures.cif --block sample
```

To write only the atomic-coordinate namelist and species-index comments:

```sh
python3 utility/cif2salmon/cif2salmon.py structure.cif --format atoms
```

Reduced coordinates are written by default. Use `--coordinates cartesian` to
write Cartesian coordinates in angstrom:

```sh
python3 utility/cif2salmon/cif2salmon.py structure.cif \
  --coordinates cartesian
```

The default duplicate-position tolerance is `1e-5` angstrom. It can be changed
with `--symprec ANGSTROM`.

## Output and conversion rules

The default output contains the following structure information:

- `&units` with `unit_system = 'A_eV_fs'`;
- `&system` with the cell, `yn_periodic`, `nelem`, and `natom`;
- comments showing each SALMON species index, atomic number, and atom count;
- a comment showing the total atom count as `natom`;
- the selected atomic-coordinate namelist with symmetry-expanded positions.

An orthorhombic (rectangular) cell is written using `al`. Other cells are written using
`al_vec1`, `al_vec2`, and `al_vec3`. The CIF conventional cell is preserved; the
converter does not reduce it to a primitive cell.

`--coordinates reduced` writes fractional positions in `&atomic_red_coor` and
is the default. `--coordinates cartesian` writes positions in angstrom in
`&atomic_coor`. Both forms use the same symmetry-expanded unit cell.

The output is a structure fragment, not a complete SALMON input. In particular,
the converter cannot infer pseudopotential files, valence-electron counts,
orbital counts, real-space grids, or k-point grids. The species index in the
fifth column of `&atomic_red_coor` or `&atomic_coor` must correspond to the same
index in the user-provided `&pseudo` namelist.

## Validation and limitations

The converter deliberately stops instead of guessing when a CIF contains a
structure that cannot be represented unambiguously by SALMON atomic-coordinate
input. It rejects:

- explicit partial or unknown occupancies;
- disorder groups and mixed-site structures;
- atom sites whose chemical element cannot be determined;
- missing or invalid unit-cell data;
- missing space-group metadata and symmetry operations;
- inconsistent space-group information;
- overlapping sites with different labels or elements.

If `_atom_site_occupancy` is absent, full occupancy is assumed. Fractional atom
coordinates are preferred, but Cartesian atom coordinates are also accepted and
converted using the CIF unit cell. Coordinates are wrapped into `[0, 1)` after
symmetry expansion.

The converter accepts ordinary CIF 1.1 files supported by Gemmi, including
gzip-compressed `.cif.gz` files. CIF 2.0-specific syntax is not guaranteed.
