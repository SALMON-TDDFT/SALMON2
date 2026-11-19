#!/usr/bin/env python3
#
# Copyright 2026 SALMON developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Convert a crystallographic information file to a SALMON input fragment."""

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import sys
from typing import Any, List, Optional, Sequence, Tuple

try:
    import gemmi
except ImportError:
    gemmi = None


FRACTIONAL_TAGS = (
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
)
CARTESIAN_TAGS = (
    "_atom_site_Cartn_x",
    "_atom_site_Cartn_y",
    "_atom_site_Cartn_z",
)
CELL_TAGS = (
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
)


class ConversionError(Exception):
    """An expected error caused by an unsupported or invalid CIF structure."""


@dataclass(frozen=True)
class Atom:
    """One atom in the expanded unit cell."""

    symbol: str
    atomic_number: int
    label: str
    fractional: Tuple[float, float, float]


@dataclass(frozen=True)
class SalmonStructure:
    """Validated information needed for a SALMON structure fragment."""

    block_name: str
    vectors: Tuple[Tuple[float, float, float], ...]
    atoms: Tuple[Atom, ...]
    species: Tuple[Tuple[str, int], ...]


def _require_gemmi() -> None:
    if gemmi is None:
        raise ConversionError(
            "the optional dependency 'gemmi' is required; install it with "
            "'python3 -m pip install gemmi'"
        )


def _column(block: Any, tag: str) -> List[str]:
    return list(block.find_values(tag))


def _has_complete_columns(block: Any, tags: Sequence[str]) -> bool:
    sizes = [len(_column(block, tag)) for tag in tags]
    return bool(sizes[0]) and len(set(sizes)) == 1


def _looks_like_structure(block: Any) -> bool:
    has_cell = all(bool(block.find_value(tag)) for tag in CELL_TAGS)
    has_coordinates = _has_complete_columns(block, FRACTIONAL_TAGS)
    has_coordinates = has_coordinates or _has_complete_columns(block, CARTESIAN_TAGS)
    return has_cell and has_coordinates


def _select_block(document: Any, requested_name: Optional[str]) -> Any:
    if requested_name:
        name = requested_name[5:] if requested_name.lower().startswith("data_") else requested_name
        for block in document:
            if block.name.lower() == name.lower():
                if not _looks_like_structure(block):
                    raise ConversionError(
                        "CIF data block '{}' does not contain a complete periodic structure".format(
                            block.name
                        )
                    )
                return block
        available = ", ".join(block.name for block in document) or "(none)"
        raise ConversionError(
            "CIF data block '{}' was not found; available blocks: {}".format(
                requested_name, available
            )
        )

    candidates = [block for block in document if _looks_like_structure(block)]
    if not candidates:
        raise ConversionError("the CIF contains no data block with a complete periodic structure")
    if len(candidates) > 1:
        names = ", ".join(block.name for block in candidates)
        raise ConversionError(
            "the CIF contains multiple structure blocks ({}); select one with --block".format(
                names
            )
        )
    return candidates[0]


def _parse_number(raw_value: str, description: str) -> float:
    if gemmi.cif.is_null(raw_value):
        raise ConversionError("{} is unknown or missing".format(description))
    value = gemmi.cif.as_number(raw_value)
    if not math.isfinite(value):
        raise ConversionError("{} is not a finite number".format(description))
    return value


def _validate_cell(structure: Any) -> None:
    cell = structure.cell
    lengths = (cell.a, cell.b, cell.c)
    angles = (cell.alpha, cell.beta, cell.gamma)
    if any(not math.isfinite(value) or value <= 0.0 for value in lengths):
        raise ConversionError("cell lengths must be finite and positive")
    if any(not math.isfinite(value) or not 0.0 < value < 180.0 for value in angles):
        raise ConversionError("cell angles must be finite and between 0 and 180 degrees")
    if not math.isfinite(cell.volume) or cell.volume <= 1.0e-12:
        raise ConversionError("the unit-cell volume is zero or invalid")


def _set_cartesian_coordinates(block: Any, structure: Any) -> None:
    if _has_complete_columns(block, FRACTIONAL_TAGS):
        return
    if not _has_complete_columns(block, CARTESIAN_TAGS):
        raise ConversionError("atom sites require complete fractional or Cartesian coordinates")

    columns = [_column(block, tag) for tag in CARTESIAN_TAGS]
    if len(columns[0]) != len(structure.sites):
        raise ConversionError("the Cartesian coordinate table does not match the atom-site table")
    for index, site in enumerate(structure.sites):
        xyz = tuple(
            _parse_number(
                columns[axis][index],
                "Cartesian coordinate for site {}".format(index + 1),
            )
            for axis in range(3)
        )
        site.fract = structure.cell.fractionalize(gemmi.Position(*xyz))


def _validate_sites(block: Any, structure: Any) -> None:
    if not structure.sites:
        raise ConversionError("the selected CIF data block contains no atom sites")

    occupancies = _column(block, "_atom_site_occupancy")
    if occupancies and len(occupancies) != len(structure.sites):
        raise ConversionError("the occupancy column does not match the atom-site table")

    for index, site in enumerate(structure.sites):
        description = site.label or "site {}".format(index + 1)
        if occupancies:
            occupancy = _parse_number(
                occupancies[index], "occupancy of '{}'".format(description)
            )
        else:
            occupancy = site.occ
        if abs(occupancy - 1.0) > 1.0e-8:
            raise ConversionError(
                "partial occupancy is unsupported: '{}' has occupancy {}".format(
                    description, occupancy
                )
            )
        if site.disorder_group:
            raise ConversionError(
                "disordered atom sites are unsupported: '{}' belongs to disorder group {}".format(
                    description, site.disorder_group
                )
            )
        if site.element.atomic_number <= 0:
            raise ConversionError(
                "the chemical element of '{}' cannot be determined".format(description)
            )
        coordinates = (site.fract.x, site.fract.y, site.fract.z)
        if any(not math.isfinite(value) for value in coordinates):
            raise ConversionError(
                "fractional coordinates of '{}' are missing or invalid".format(description)
            )


def _prepare_spacegroup(structure: Any) -> None:
    has_metadata = bool(
        structure.symops
        or structure.spacegroup_hall
        or structure.spacegroup_hm
        or structure.spacegroup_number
    )
    if not has_metadata:
        raise ConversionError(
            "space-group metadata or explicit symmetry operations are required; "
            "the converter will not silently assume P1"
        )

    structure.determine_and_set_spacegroup("S.H21N")
    inconsistency = structure.check_spacegroup().strip()
    if inconsistency:
        raise ConversionError("inconsistent space-group information: {}".format(inconsistency))


def _cell_vectors(cell: Any) -> Tuple[Tuple[float, float, float], ...]:
    matrix = cell.orth.mat.tolist()
    return tuple(
        tuple(float(matrix[row][column]) for row in range(3)) for column in range(3)
    )


def _normalize_fractional(value: float) -> float:
    value = value - math.floor(value)
    if value < 1.0e-12 or 1.0 - value < 1.0e-12:
        return 0.0
    return value


def _periodic_distance(
    first: Sequence[float],
    second: Sequence[float],
    vectors: Sequence[Sequence[float]],
) -> float:
    delta = [first[axis] - second[axis] for axis in range(3)]
    delta = [value - round(value) for value in delta]
    cartesian = [
        sum(vectors[axis][component] * delta[axis] for axis in range(3))
        for component in range(3)
    ]
    return math.sqrt(sum(value * value for value in cartesian))


def _expanded_atoms(
    structure: Any,
    vectors: Sequence[Sequence[float]],
    symmetry_tolerance: float,
) -> Tuple[Atom, ...]:
    atoms: List[Atom] = []
    for site in structure.get_all_unit_cell_sites():
        fractional = tuple(
            _normalize_fractional(value)
            for value in (site.fract.x, site.fract.y, site.fract.z)
        )
        atom = Atom(
            symbol=site.element.name,
            atomic_number=site.element.atomic_number,
            label=site.label,
            fractional=fractional,
        )
        duplicate = None
        for previous in atoms:
            if _periodic_distance(fractional, previous.fractional, vectors) <= symmetry_tolerance:
                duplicate = previous
                break
        if duplicate is None:
            atoms.append(atom)
        elif duplicate.atomic_number != atom.atomic_number:
            raise ConversionError(
                "different elements occupy the same position: '{}' and '{}'".format(
                    duplicate.label, atom.label
                )
            )
        elif duplicate.label != atom.label:
            raise ConversionError(
                "distinct atom sites overlap after symmetry expansion: '{}' and '{}'".format(
                    duplicate.label, atom.label
                )
            )
        # Equal labels are symmetry images of the same source site and are deduplicated.
    return tuple(atoms)


def read_cif_structure(
    filename: Path,
    block_name: Optional[str] = None,
    symmetry_tolerance: float = 1.0e-5,
) -> SalmonStructure:
    """Read, expand, and validate one periodic structure from a CIF."""

    _require_gemmi()
    if symmetry_tolerance <= 0.0 or not math.isfinite(symmetry_tolerance):
        raise ConversionError("--symprec must be a finite positive distance in angstrom")
    try:
        document = gemmi.cif.read(str(filename))
    except (OSError, RuntimeError, ValueError) as error:
        raise ConversionError("cannot read '{}': {}".format(filename, error)) from error

    block = _select_block(document, block_name)
    try:
        structure = gemmi.make_small_structure_from_block(block)
    except (RuntimeError, ValueError) as error:
        raise ConversionError(
            "cannot interpret CIF data block '{}': {}".format(block.name, error)
        ) from error

    _validate_cell(structure)
    _set_cartesian_coordinates(block, structure)
    _validate_sites(block, structure)
    _prepare_spacegroup(structure)
    vectors = _cell_vectors(structure.cell)
    atoms = _expanded_atoms(structure, vectors, symmetry_tolerance)
    if not atoms:
        raise ConversionError("symmetry expansion produced no atoms")

    species: List[Tuple[str, int]] = []
    seen_atomic_numbers = set()
    for atom in atoms:
        if atom.atomic_number not in seen_atomic_numbers:
            species.append((atom.symbol, atom.atomic_number))
            seen_atomic_numbers.add(atom.atomic_number)

    return SalmonStructure(
        block_name=block.name,
        vectors=vectors,
        atoms=atoms,
        species=tuple(species),
    )


def _format_number(value: float) -> str:
    if abs(value) < 5.0e-13:
        value = 0.0
    return "{:.12f}".format(value)


def _orthorhombic_lengths(
    vectors: Sequence[Sequence[float]],
) -> Optional[Tuple[float, float, float]]:
    off_diagonal = (
        vectors[0][1],
        vectors[0][2],
        vectors[1][0],
        vectors[1][2],
        vectors[2][0],
        vectors[2][1],
    )
    if all(abs(value) <= 1.0e-10 for value in off_diagonal):
        return tuple(math.sqrt(sum(value * value for value in vector)) for vector in vectors)
    return None


def _fractional_to_cartesian(
    fractional: Sequence[float],
    vectors: Sequence[Sequence[float]],
) -> Tuple[float, float, float]:
    return tuple(
        sum(vectors[axis][component] * fractional[axis] for axis in range(3))
        for component in range(3)
    )


def render_salmon(
    structure: SalmonStructure,
    source_name: str,
    atoms_only: bool = False,
    coordinate_system: str = "reduced",
) -> str:
    """Render a structure as a pasteable SALMON namelist fragment."""

    if coordinate_system not in ("reduced", "cartesian"):
        raise ConversionError(
            "coordinate system must be either 'reduced' or 'cartesian'"
        )

    species_index = {
        atomic_number: index
        for index, (_, atomic_number) in enumerate(structure.species, start=1)
    }
    species_count = {
        atomic_number: sum(
            atom.atomic_number == atomic_number for atom in structure.atoms
        )
        for _, atomic_number in structure.species
    }
    lines = [
        "! Generated by utility/cif2salmon/cif2salmon.py from {}".format(
            Path(source_name).name
        ),
        "! CIF data block: {}".format(structure.block_name),
        "! This is a structure fragment, not a complete SALMON input.",
    ]

    if not atoms_only:
        lines.extend(["", "&units", "  unit_system = 'A_eV_fs'", "/", "", "&system"])
        lines.append("  yn_periodic = 'y'")
        lengths = _orthorhombic_lengths(structure.vectors)
        if lengths is not None:
            lines.append("  al = {}".format(", ".join(_format_number(x) for x in lengths)))
        else:
            for index, vector in enumerate(structure.vectors, start=1):
                lines.append(
                    "  al_vec{} = {}".format(
                        index, ", ".join(_format_number(x) for x in vector)
                    )
                )
        lines.extend(
            [
                "  nelem = {}".format(len(structure.species)),
                "  natom = {}".format(len(structure.atoms)),
                "/",
            ]
        )

    lines.append("")
    for index, (symbol, atomic_number) in enumerate(structure.species, start=1):
        lines.append(
            "! species index: {} = {} (Z={}), count = {}".format(
                index, symbol, atomic_number, species_count[atomic_number]
            )
        )
    lines.append("! total atom count: natom = {}".format(len(structure.atoms)))
    if coordinate_system == "reduced":
        namelist_name = "atomic_red_coor"
    else:
        namelist_name = "atomic_coor"
        lines.append("! Cartesian atomic coordinates are in angstrom.")
    lines.extend(["", "&{}".format(namelist_name)])
    for atom in structure.atoms:
        if coordinate_system == "reduced":
            coordinates = atom.fractional
        else:
            coordinates = _fractional_to_cartesian(
                atom.fractional, structure.vectors
            )
        lines.append(
            "  '{}'  {}  {}  {}  {:d}".format(
                atom.symbol,
                *(_format_number(value) for value in coordinates),
                species_index[atom.atomic_number]
            )
        )
    lines.extend(["/", ""])
    return "\n".join(lines)


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert one periodic CIF structure into a SALMON atomic-coordinate "
            "namelist and its companion lattice information."
        )
    )
    parser.add_argument("cif_file", type=Path, help="input .cif or .cif.gz file")
    parser.add_argument("-o", "--output", type=Path, help="output file (default: stdout)")
    parser.add_argument(
        "--force", action="store_true", help="allow --output to overwrite an existing file"
    )
    parser.add_argument(
        "--block", help="CIF data-block name; required when multiple structures are present"
    )
    parser.add_argument(
        "--format",
        choices=("fragment", "atoms"),
        default="fragment",
        help="write a structure fragment or only the atomic namelist (default: fragment)",
    )
    parser.add_argument(
        "--coordinates",
        choices=("reduced", "cartesian"),
        default="reduced",
        help=(
            "write &atomic_red_coor or angstrom-valued &atomic_coor "
            "(default: reduced)"
        ),
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=1.0e-5,
        metavar="ANGSTROM",
        help="periodic duplicate tolerance in angstrom (default: 1e-5)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _argument_parser().parse_args(argv)
    try:
        structure = read_cif_structure(args.cif_file, args.block, args.symprec)
        output = render_salmon(
            structure,
            str(args.cif_file),
            atoms_only=args.format == "atoms",
            coordinate_system=args.coordinates,
        )
        if args.output is None:
            sys.stdout.write(output)
        else:
            mode = "w" if args.force else "x"
            try:
                with args.output.open(mode, encoding="utf-8") as stream:
                    stream.write(output)
            except FileExistsError as error:
                raise ConversionError(
                    "output file '{}' already exists; use --force to replace it".format(
                        args.output
                    )
                ) from error
            except OSError as error:
                raise ConversionError(
                    "cannot write output file '{}': {}".format(args.output, error)
                ) from error
    except ConversionError as error:
        print("cif2salmon: error: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
