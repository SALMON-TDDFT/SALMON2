#!/usr/bin/env python3
"""Align a periodic structure so crystal symmetries preserve fragment boxes."""

import argparse
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import os
from pathlib import Path
import shlex
import sys
import tempfile


_DIMENSION = 3
_MAX_POINT_OPERATION_ORDER = 12
_IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def _fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, float):
        return Fraction(str(value))
    return Fraction(value)


def _validate_vector(values, name):
    values = tuple(values)
    if len(values) != _DIMENSION:
        raise ValueError(f"{name} must have exactly {_DIMENSION} entries")
    return values


def _validate_mesh(mesh):
    mesh = _validate_vector(mesh, "mesh")
    if any(not isinstance(size, int) or isinstance(size, bool) or size <= 0 for size in mesh):
        raise ValueError("mesh dimensions must be positive integers")
    return mesh


def _validate_integer_matrix(matrix, name):
    rows = tuple(tuple(row) for row in matrix)
    if len(rows) != _DIMENSION or any(len(row) != _DIMENSION for row in rows):
        raise ValueError(f"{name} must be a 3x3 integer matrix")
    if any(
        not isinstance(value, int) or isinstance(value, bool)
        for row in rows
        for value in row
    ):
        raise ValueError(f"{name} must be a 3x3 integer matrix")
    return rows


def _determinant_3x3(matrix):
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _matrix_product(left, right):
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(_DIMENSION))
            for column in range(_DIMENSION)
        )
        for row in range(_DIMENSION)
    )


def _has_bounded_finite_order(matrix):
    power = _IDENTITY
    for _ in range(_MAX_POINT_OPERATION_ORDER):
        power = _matrix_product(power, matrix)
        if power == _IDENTITY:
            return True
    return False


@dataclass(frozen=True)
class SymOp:
    """A crystal symmetry operation x' = W x + tau in fractional coordinates."""

    rotation: tuple
    translation: tuple

    def __post_init__(self):
        rows = _validate_integer_matrix(self.rotation, "rotation")
        if abs(_determinant_3x3(rows)) != 1:
            raise ValueError("rotation determinant must be +/-1")
        if not _has_bounded_finite_order(rows):
            raise ValueError(
                "rotation must have finite order no greater than "
                f"{_MAX_POINT_OPERATION_ORDER}"
            )
        translation = tuple(
            _fraction(value) % 1
            for value in _validate_vector(self.translation, "translation")
        )
        object.__setattr__(self, "rotation", rows)
        object.__setattr__(self, "translation", translation)


@dataclass(frozen=True)
class PeriodicAtom:
    """An atom labeled by species at an exact fractional position."""

    species: str
    position: tuple
    label: str = ""
    species_index: str = ""
    trailing: tuple = ()

    def __post_init__(self):
        if not isinstance(self.species, str) or not self.species:
            raise ValueError("atom species must be a non-empty string")
        position = tuple(
            _fraction(value) % 1
            for value in _validate_vector(self.position, "atom position")
        )
        object.__setattr__(self, "position", position)
        object.__setattr__(self, "trailing", tuple(self.trailing))


@dataclass(frozen=True)
class AlignmentResult:
    """Validated periodic translation and its aligned geometry."""

    translation: tuple
    atoms: tuple
    operations: tuple
    grid_maps: tuple
    fragment_maps: tuple
    self_mapped_fragments: int


@dataclass(frozen=True)
class IntegerGridMap:
    """An integer affine permutation of a finite periodic grid."""

    mesh: tuple
    coefficients: tuple
    shift: tuple

    def __post_init__(self):
        mesh = _validate_mesh(self.mesh)
        coefficients = _validate_integer_matrix(self.coefficients, "coefficients")
        shift = _validate_vector(self.shift, "shift")
        if any(
            not isinstance(value, int) or isinstance(value, bool) for value in shift
        ):
            raise ValueError("shift must contain exactly 3 integers")
        shift = tuple(shift[axis] % mesh[axis] for axis in range(_DIMENSION))
        object.__setattr__(self, "mesh", mesh)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "shift", shift)

        images = {
            self.map_index(index)
            for index in itertools.product(*(range(size) for size in mesh))
        }
        if len(images) != mesh[0] * mesh[1] * mesh[2]:
            raise ValueError(
                "grid coefficients and shift do not define a bijection of the "
                "periodic discrete mesh"
            )

    def map_index(self, index):
        index = _validate_vector(index, "grid index")
        if any(
            not isinstance(value, int)
            or isinstance(value, bool)
            or value < 0
            or value >= self.mesh[axis]
            for axis, value in enumerate(index)
        ):
            raise ValueError("grid index is outside the periodic mesh")
        return tuple(
            (
                sum(self.coefficients[alpha][beta] * index[beta] for beta in range(_DIMENSION))
                + self.shift[alpha]
            )
            % self.mesh[alpha]
            for alpha in range(_DIMENSION)
        )


def integer_grid_map(operation, mesh):
    """Return the exact integer grid permutation induced by ``operation``.

    Fractional coordinates obey x_beta = i_beta / N_beta, hence

      j_alpha = sum_beta W(alpha,beta) i_beta N_alpha/N_beta
                + N_alpha tau_alpha  (mod N_alpha).

    Every scaled coefficient and shift must be integral, and the resulting
    finite map must be a bijection.
    """

    if not isinstance(operation, SymOp):
        raise TypeError("operation must be a SymOp")
    mesh = _validate_mesh(mesh)
    coefficients = []
    for alpha in range(_DIMENSION):
        row = []
        for beta in range(_DIMENSION):
            coefficient = Fraction(
                operation.rotation[alpha][beta] * mesh[alpha], mesh[beta]
            )
            if coefficient.denominator != 1:
                raise ValueError(
                    "symmetry induces a non-integral grid coefficient "
                    f"at ({alpha},{beta}): {coefficient}"
                )
            row.append(coefficient.numerator)
        coefficients.append(tuple(row))

    shift = []
    for alpha in range(_DIMENSION):
        component = mesh[alpha] * operation.translation[alpha]
        if component.denominator != 1:
            raise ValueError(
                "symmetry induces a non-integral grid shift "
                f"on axis {alpha}: {component}"
            )
        shift.append(component.numerator % mesh[alpha])

    return IntegerGridMap(mesh, tuple(coefficients), tuple(shift))


def inversion_center_diagnostic(grid_map):
    """Return inversion-center coordinates measured in grid spacings.

    The center is q/2 and may lie at a half-grid position. It is diagnostic
    information, never an integrality condition for accepting the map.
    """

    if not isinstance(grid_map, IntegerGridMap):
        raise TypeError("grid_map must be an IntegerGridMap")
    negative_identity = ((-1, 0, 0), (0, -1, 0), (0, 0, -1))
    if grid_map.coefficients != negative_identity:
        raise ValueError("grid map is not an inversion")
    return tuple(Fraction(component, 2) for component in grid_map.shift)


def fragment_target_enumeration(grid_map, fragment_shape):
    """Enumerate target fragments reached from every source fragment."""

    if not isinstance(grid_map, IntegerGridMap):
        raise TypeError("grid_map must be an IntegerGridMap")
    fragment_shape = _validate_mesh(fragment_shape)
    if any(
        grid_map.mesh[axis] % fragment_shape[axis] != 0
        for axis in range(_DIMENSION)
    ):
        raise ValueError("fragment dimensions must divide the periodic mesh")

    fragment_counts = tuple(
        grid_map.mesh[axis] // fragment_shape[axis]
        for axis in range(_DIMENSION)
    )
    targets = {
        source: set()
        for source in itertools.product(*(range(count) for count in fragment_counts))
    }
    for index in itertools.product(*(range(size) for size in grid_map.mesh)):
        source = tuple(
            index[axis] // fragment_shape[axis] for axis in range(_DIMENSION)
        )
        image = grid_map.map_index(index)
        target = tuple(
            image[axis] // fragment_shape[axis] for axis in range(_DIMENSION)
        )
        targets[source].add(target)
    return {source: frozenset(values) for source, values in targets.items()}


def periodic_same_species_atom_bijection(operation, atoms):
    """Return the atom permutation induced by a fractional symmetry operation."""

    if not isinstance(operation, SymOp):
        raise TypeError("operation must be a SymOp")
    atoms = tuple(atoms)
    if any(not isinstance(atom, PeriodicAtom) for atom in atoms):
        raise TypeError("atoms must contain PeriodicAtom values")

    unused = set(range(len(atoms)))
    permutation = []
    for atom in atoms:
        transformed = tuple(
            (
                sum(
                    operation.rotation[alpha][beta] * atom.position[beta]
                    for beta in range(_DIMENSION)
                )
                + operation.translation[alpha]
            )
            % 1
            for alpha in range(_DIMENSION)
        )
        matches = [
            index
            for index in sorted(unused)
            if atoms[index].species == atom.species
            and atoms[index].position == transformed
        ]
        if not matches:
            raise ValueError(
                "symmetry does not define a periodic same-species atom bijection"
            )
        target = matches[0]
        permutation.append(target)
        unused.remove(target)
    if unused:
        raise ValueError("symmetry does not define a periodic same-species atom bijection")
    return tuple(permutation)


def _signed_axis_sources(operation):
    """Return (input axis, sign) for each output axis."""

    sources = []
    used = set()
    for alpha, row in enumerate(operation.rotation):
        nonzero = [(beta, value) for beta, value in enumerate(row) if value]
        if len(nonzero) != 1 or abs(nonzero[0][1]) != 1:
            raise ValueError(
                "operation is not a signed-axis permutation "
                f"on output axis {alpha}"
            )
        beta, sign = nonzero[0]
        if beta in used:
            raise ValueError("operation is not a signed-axis permutation")
        used.add(beta)
        sources.append((beta, sign))
    return tuple(sources)


def validate_lattice_metric(operation, cell):
    """Require exact preservation of an orthogonal-cell lattice metric."""

    cell = tuple(_fraction(value) for value in _validate_vector(cell, "cell"))
    if any(length <= 0 for length in cell):
        raise ValueError("cell lengths must be positive")
    metric = tuple(cell[axis] ** 2 for axis in range(_DIMENSION))
    for beta in range(_DIMENSION):
        for gamma in range(_DIMENSION):
            transformed = sum(
                operation.rotation[alpha][beta]
                * metric[alpha]
                * operation.rotation[alpha][gamma]
                for alpha in range(_DIMENSION)
            )
            expected = metric[beta] if beta == gamma else Fraction(0)
            if transformed != expected:
                raise ValueError(
                    "operation does not preserve the orthogonal lattice metric"
                )


def _validate_buffer(buffer):
    buffer = _validate_vector(buffer, "buffer")
    if any(
        not isinstance(width, int) or isinstance(width, bool) or width < 0
        for width in buffer
    ):
        raise ValueError("buffer widths must be nonnegative integers")
    return buffer


def validate_buffer_geometry(operation, buffer):
    """Require a signed-axis operation to map complete buffer boxes."""

    buffer = _validate_buffer(buffer)
    for alpha, (beta, _sign) in enumerate(_signed_axis_sources(operation)):
        if buffer[alpha] != buffer[beta]:
            raise ValueError(
                "buffer width on output axis "
                f"{alpha} must equal input axis {beta}"
            )


def _compose_operations(left, right):
    rotation = _matrix_product(left.rotation, right.rotation)
    translation = tuple(
        (
            sum(
                left.rotation[alpha][beta] * right.translation[beta]
                for beta in range(_DIMENSION)
            )
            + left.translation[alpha]
        )
        % 1
        for alpha in range(_DIMENSION)
    )
    return SymOp(rotation, translation)


def validate_symmetry_group(operations):
    """Validate identity uniqueness, inverses, and affine group closure."""

    operations = tuple(operations)
    if not operations:
        raise ValueError("symmetry operation set is empty")
    identity = SymOp(_IDENTITY, (0, 0, 0))
    identity_count = sum(operation == identity for operation in operations)
    if identity_count != 1:
        raise ValueError("symmetry identity must occur exactly once")
    operation_set = set(operations)
    if len(operation_set) != len(operations):
        raise ValueError("symmetry operation set contains duplicates")
    for left in operations:
        for right in operations:
            if _compose_operations(left, right) not in operation_set:
                raise ValueError("symmetry operation set fails group closure")
        if not any(
            _compose_operations(left, right) == identity
            and _compose_operations(right, left) == identity
            for right in operations
        ):
            raise ValueError("symmetry operation has no inverse in the group")


def parse_symmetry_file(path):
    """Parse ordered three-row affine operations from ``sym.dat``."""

    rows = []
    for line_number, line in enumerate(Path(path).read_text(encoding="ascii").splitlines(), 1):
        content = line.split("#", 1)[0].strip()
        if not content:
            continue
        fields = content.split()
        if len(fields) != 4:
            raise ValueError(
                f"symmetry row {line_number} must contain W(row,1:3), tau(row)"
            )
        try:
            rows.append((tuple(int(value) for value in fields[:3]), _fraction(fields[3])))
        except ValueError as error:
            raise ValueError(f"invalid symmetry row {line_number}: {error}") from error
    if not rows or len(rows) % _DIMENSION:
        raise ValueError("symmetry file must contain exactly three rows per operation")
    operations = []
    for start in range(0, len(rows), _DIMENSION):
        operations.append(
            SymOp(
                tuple(rows[start + axis][0] for axis in range(_DIMENSION)),
                tuple(rows[start + axis][1] for axis in range(_DIMENSION)),
            )
        )
    return tuple(operations)


def parse_atom_file(path, cell):
    """Parse SALMON Cartesian atom rows into exact fractional coordinates."""

    cell = tuple(_fraction(value) for value in _validate_vector(cell, "cell"))
    if any(length <= 0 for length in cell):
        raise ValueError("cell lengths must be positive")
    atoms = []
    for line_number, line in enumerate(Path(path).read_text(encoding="ascii").splitlines(), 1):
        content = line.split("#", 1)[0].strip()
        if not content:
            continue
        try:
            fields = shlex.split(content, posix=False)
        except ValueError as error:
            raise ValueError(f"invalid atom row {line_number}: {error}") from error
        if len(fields) < 5:
            raise ValueError(
                f"atom row {line_number} must contain label, x, y, z, species index"
            )
        label = fields[0]
        species = label.strip("'\"")
        try:
            cartesian = tuple(_fraction(fields[axis + 1]) for axis in range(_DIMENSION))
        except ValueError as error:
            raise ValueError(f"invalid Cartesian atom row {line_number}: {error}") from error
        atoms.append(
            PeriodicAtom(
                species,
                tuple(cartesian[axis] / cell[axis] for axis in range(_DIMENSION)),
                label,
                fields[4],
                tuple(fields[5:]),
            )
        )
    if not atoms:
        raise ValueError("atom file contains no atom rows")
    return tuple(atoms)


def _translate_operation(operation, translation):
    return SymOp(
        operation.rotation,
        tuple(
            (
                operation.translation[alpha]
                + translation[alpha]
                - sum(
                    operation.rotation[alpha][beta] * translation[beta]
                    for beta in range(_DIMENSION)
                )
            )
            % 1
            for alpha in range(_DIMENSION)
        ),
    )


def _translate_atoms(atoms, translation):
    return tuple(
        PeriodicAtom(
            atom.species,
            tuple((atom.position[axis] + translation[axis]) % 1 for axis in range(_DIMENSION)),
            atom.label,
            atom.species_index,
            atom.trailing,
        )
        for atom in atoms
    )


def _candidate_operation_data(operation, mesh):
    """Precompute exact integer affine data without constructing a grid map."""

    sources = _signed_axis_sources(operation)
    base_shift = []
    coefficients = []
    for alpha, (beta, sign) in enumerate(sources):
        coefficient = Fraction(sign * mesh[alpha], mesh[beta])
        if coefficient.denominator != 1 or abs(coefficient.numerator) != 1:
            raise ValueError(
                "signed-axis map must preserve grid spacing across mapped axes"
            )
        shift = mesh[alpha] * operation.translation[alpha]
        if shift.denominator != 1:
            raise ValueError(
                f"symmetry induces a non-integral grid shift on axis {alpha}: {shift}"
            )
        coefficients.append(coefficient.numerator)
        base_shift.append(shift.numerator)
    return sources, tuple(coefficients), tuple(base_shift)


def _candidate_grid_shifts(operation_data, numerators, mesh):
    """Apply half-grid translation using integer congruences only."""

    sources, coefficients, base_shift = operation_data
    shifts = []
    for alpha, (beta, _sign) in enumerate(sources):
        doubled = (
            2 * base_shift[alpha]
            + numerators[alpha]
            - coefficients[alpha] * numerators[beta]
        )
        if doubled % 2:
            return None
        shifts.append((doubled // 2) % mesh[alpha])
    return tuple(shifts)


def _fast_complete_fragment_map(operation_data, shifts, mesh, fragment_counts):
    """Map fragment intervals from integer shifts, without any grid enumeration."""

    shape = tuple(mesh[axis] // fragment_counts[axis] for axis in range(_DIMENSION))
    sources, _coefficients, _base_shift = operation_data
    mapping = {}
    for source in itertools.product(*(range(count) for count in fragment_counts)):
        target = []
        for alpha, (beta, sign) in enumerate(sources):
            if shape[alpha] != shape[beta]:
                raise ValueError("mapped fragment axes must have equal grid widths")
            start = source[beta] * shape[beta]
            if sign > 0:
                image_start = (start + shifts[alpha]) % mesh[alpha]
            else:
                image_start = (
                    shifts[alpha] - (start + shape[beta] - 1)
                ) % mesh[alpha]
            if image_start % shape[alpha]:
                return None
            target.append(image_start // shape[alpha])
        mapping[source] = tuple(target)
    return mapping


def _wrapped_component(value):
    value %= 1
    return value - 1 if value > Fraction(1, 2) else value


def find_fragment_compatible_translation(
    atoms, operations, cell, mesh, fragment_counts, buffer
):
    """Search deterministic half-grid translations and fully validate the winner."""

    atoms = tuple(atoms)
    operations = tuple(operations)
    cell = tuple(_fraction(value) for value in _validate_vector(cell, "cell"))
    mesh = _validate_mesh(mesh)
    fragment_counts = _validate_mesh(fragment_counts)
    buffer = _validate_buffer(buffer)
    validate_symmetry_group(operations)
    for operation in operations:
        validate_lattice_metric(operation, cell)
        validate_buffer_geometry(operation, buffer)
        periodic_same_species_atom_bijection(operation, atoms)

    denominators = tuple(2 * size for size in mesh)
    operation_data = tuple(
        _candidate_operation_data(operation, mesh) for operation in operations
    )
    best = None
    for numerators in itertools.product(*(range(value) for value in denominators)):
        maps = []
        for data in operation_data:
            shifts = _candidate_grid_shifts(data, numerators, mesh)
            if shifts is None:
                break
            fragment_map = _fast_complete_fragment_map(
                data, shifts, mesh, fragment_counts
            )
            if fragment_map is None:
                break
            maps.append(fragment_map)
        if len(maps) != len(operations):
            continue
        self_count = sum(
            source == target for fragment_map in maps for source, target in fragment_map.items()
        )
        translation = tuple(
            Fraction(numerators[axis], denominators[axis])
            for axis in range(_DIMENSION)
        )
        wrapped = tuple(_wrapped_component(value) for value in translation)
        norm = sum((wrapped[axis] * cell[axis]) ** 2 for axis in range(_DIMENSION))
        rank = (-self_count, norm, translation)
        if best is None or rank < best[0]:
            best = (rank, translation, tuple(maps), self_count)
    if best is None:
        raise ValueError("no half-grid translation preserves complete fragment boxes")

    _rank, translation, fast_maps, self_count = best
    aligned_operations = tuple(
        _translate_operation(operation, translation) for operation in operations
    )
    aligned_atoms = _translate_atoms(atoms, translation)
    validate_symmetry_group(aligned_operations)
    grid_maps = tuple(integer_grid_map(operation, mesh) for operation in aligned_operations)
    shape = tuple(mesh[axis] // fragment_counts[axis] for axis in range(_DIMENSION))
    full_maps = []
    for operation, grid_map, fast_map in zip(aligned_operations, grid_maps, fast_maps):
        validate_lattice_metric(operation, cell)
        validate_buffer_geometry(operation, buffer)
        periodic_same_species_atom_bijection(operation, aligned_atoms)
        targets = fragment_target_enumeration(grid_map, shape)
        if any(len(values) != 1 for values in targets.values()):
            raise ValueError("final grid validation found a split fragment box")
        full_map = {source: next(iter(values)) for source, values in targets.items()}
        if full_map != fast_map:
            raise ValueError("fast and full fragment maps disagree")
        full_maps.append(tuple(sorted(full_map.items())))
    return AlignmentResult(
        translation,
        aligned_atoms,
        aligned_operations,
        grid_maps,
        tuple(full_maps),
        self_count,
    )


def _format_fraction(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _format_decimal(value, digits=24):
    value = _fraction(value)
    denominator = value.denominator
    twos = 0
    fives = 0
    while denominator % 2 == 0:
        denominator //= 2
        twos += 1
    while denominator % 5 == 0:
        denominator //= 5
        fives += 1
    if denominator == 1:
        places = max(twos, fives)
        scaled = abs(value.numerator) * 2 ** (places - twos) * 5 ** (places - fives)
        divisor = 10 ** places
        integer, fractional = divmod(scaled, divisor)
        sign = "-" if value < 0 else ""
        if places == 0:
            return f"{sign}{integer}.0"
        return f"{sign}{integer}.{fractional:0{places}d}"
    with localcontext() as context:
        context.prec = digits + 8
        decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        return format(decimal_value, f".{digits}f")


def _render_atoms(atoms, cell):
    lines = []
    for atom in atoms:
        cartesian = tuple(atom.position[axis] * cell[axis] for axis in range(_DIMENSION))
        fields = [
            f"  {atom.label or repr(atom.species)}",
            *(_format_decimal(value) for value in cartesian),
            atom.species_index,
            *atom.trailing,
        ]
        lines.append(" ".join(fields).rstrip())
    return "\n".join(lines) + "\n"


def _render_symmetry(operations, translation, buffer):
    translation_text = " ".join(_format_fraction(value) for value in translation)
    buffer_text = " ".join(str(value) for value in buffer)
    lines = [
        "# Generated by align_periodic_structure_to_fragments.py",
        f"# translation={translation_text}",
        f"# buffer={buffer_text}",
    ]
    for operation in operations:
        for alpha in range(_DIMENSION):
            lines.append(
                " ".join(
                    [*(f"{value:2d}" for value in operation.rotation[alpha]),
                     _format_decimal(operation.translation[alpha])]
                )
            )
    return "\n".join(lines) + "\n"


def _periodic_difference(left, right):
    difference = abs((left - right) % 1)
    return min(difference, 1 - difference)


def _validate_output_paths(input_atom, input_sym, output_atom, output_sym, force):
    inputs = {Path(input_atom).resolve(), Path(input_sym).resolve()}
    outputs = (Path(output_atom).resolve(), Path(output_sym).resolve())
    if outputs[0] == outputs[1]:
        raise ValueError("atom and symmetry outputs must be different dedicated files")
    if any(path in inputs for path in outputs):
        raise ValueError("outputs must be dedicated files and must not alias inputs")
    existing = [path for path in outputs if path.exists()]
    if existing and not force:
        raise FileExistsError(f"output already exists: {existing[0]}")


def _write_validated_outputs(
    output_atom, output_sym, atom_text, sym_text, cell,
    expected_atoms, expected_operations, tolerance, force
):
    output_atom = Path(output_atom)
    output_sym = Path(output_sym)
    if output_atom.resolve() == output_sym.resolve():
        raise ValueError("atom and symmetry outputs must be different files")
    existing = [path for path in (output_atom, output_sym) if path.exists()]
    if existing and not force:
        raise FileExistsError(f"output already exists: {existing[0]}")
    for path in (output_atom, output_sym):
        path.parent.mkdir(parents=True, exist_ok=True)
    temporary_paths = []
    destinations = (output_atom, output_sym)
    originals = {
        destination: destination.read_bytes() if destination.exists() else None
        for destination in destinations
    }
    published = []
    try:
        for destination, content in ((output_atom, atom_text), (output_sym, sym_text)):
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="ascii", dir=destination.parent,
                prefix=f".{destination.name}.", suffix=".tmp", delete=False
            ) as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
                temporary_paths.append(Path(stream.name))
        reread_atoms = parse_atom_file(temporary_paths[0], cell)
        reread_operations = parse_symmetry_file(temporary_paths[1])
        if len(reread_atoms) != len(expected_atoms):
            raise ValueError("atom output re-read changed the atom count")
        for actual, expected in zip(reread_atoms, expected_atoms):
            if (
                actual.species != expected.species
                or actual.label != expected.label
                or actual.species_index != expected.species_index
                or actual.trailing != expected.trailing
                or any(
                    _periodic_difference(actual.position[axis], expected.position[axis])
                    > tolerance
                    for axis in range(_DIMENSION)
                )
            ):
                raise ValueError("atom output failed coordinate or metadata re-validation")
        if len(reread_operations) != len(expected_operations):
            raise ValueError("symmetry output re-read changed the operation count")
        for actual, expected in zip(reread_operations, expected_operations):
            if actual.rotation != expected.rotation or any(
                _periodic_difference(actual.translation[axis], expected.translation[axis])
                > tolerance
                for axis in range(_DIMENSION)
            ):
                raise ValueError("symmetry output failed affine-operation re-validation")
        try:
            for destination in destinations:
                os.replace(temporary_paths[0], destination)
                temporary_paths.pop(0)
                published.append(destination)
        except OSError:
            for destination in published:
                try:
                    destination.unlink()
                except FileNotFoundError:
                    pass
            for destination, content in originals.items():
                if content is not None:
                    destination.write_bytes(content)
            raise
    finally:
        for path in temporary_paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass


def _build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Translate a periodic SALMON structure so supplied crystal "
            "symmetries map regular fragment and buffer boxes completely."
        )
    )
    parser.add_argument("--input-atom", required=True, type=Path)
    parser.add_argument("--input-sym", required=True, type=Path)
    parser.add_argument("--output-atom", required=True, type=Path)
    parser.add_argument("--output-sym", required=True, type=Path)
    parser.add_argument("--cell", required=True, nargs=3, metavar=("LX", "LY", "LZ"))
    parser.add_argument("--mesh", required=True, nargs=3, type=int, metavar=("NX", "NY", "NZ"))
    parser.add_argument(
        "--fragments", required=True, nargs=3, type=int, metavar=("FX", "FY", "FZ")
    )
    parser.add_argument(
        "--buffer", required=True, nargs=3, type=int, metavar=("BX", "BY", "BZ"),
        help="nonnegative regular per-fragment buffer widths in grid points",
    )
    parser.add_argument(
        "--tolerance", default="1e-12",
        help="maximum periodic coordinate error accepted when outputs are re-read",
    )
    parser.add_argument("--force", action="store_true", help="replace existing output files")
    return parser


def main(argv=None):
    arguments = _build_parser().parse_args(argv)
    cell = tuple(_fraction(value) for value in arguments.cell)
    mesh = tuple(arguments.mesh)
    fragments = tuple(arguments.fragments)
    buffer = _validate_buffer(tuple(arguments.buffer))
    tolerance = _fraction(arguments.tolerance)
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    _validate_output_paths(
        arguments.input_atom, arguments.input_sym,
        arguments.output_atom, arguments.output_sym, arguments.force
    )
    if any(mesh[axis] % fragments[axis] for axis in range(_DIMENSION)):
        raise ValueError("fragment counts must divide the periodic mesh")
    atoms = parse_atom_file(arguments.input_atom, cell)
    operations = parse_symmetry_file(arguments.input_sym)
    result = find_fragment_compatible_translation(
        atoms, operations, cell, mesh, fragments, buffer
    )
    atom_text = _render_atoms(result.atoms, cell)
    sym_text = _render_symmetry(result.operations, result.translation, buffer)
    _write_validated_outputs(
        arguments.output_atom, arguments.output_sym, atom_text, sym_text,
        cell, result.atoms, result.operations, tolerance, arguments.force
    )
    translation_text = " ".join(_format_fraction(value) for value in result.translation)
    print(f"translation={translation_text}")
    print("buffer=" + " ".join(str(value) for value in buffer))
    for index, (operation, grid_map, fragment_map) in enumerate(
        zip(result.operations, result.grid_maps, result.fragment_maps), 1
    ):
        shifts = " ".join(_format_fraction(value) for value in operation.translation)
        print(f"operation={index} aligned_shift={shifts} grid_shift={grid_map.shift}")
        print(f"operation={index} fragment_map={dict(fragment_map)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
