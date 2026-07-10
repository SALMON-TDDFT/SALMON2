#!/usr/bin/env python3
"""Exact discrete geometry helpers for SAWF fragment alignment.

This module intentionally contains no file parsing or command-line interface.
It models symmetry in fractional coordinates and validates whether that
symmetry induces a permutation of a finite periodic real-space grid.
"""

from dataclasses import dataclass
from fractions import Fraction
import itertools


_DIMENSION = 3


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


@dataclass(frozen=True)
class SymOp:
    """A crystal symmetry operation x' = W x + tau in fractional coordinates."""

    rotation: tuple
    translation: tuple

    def __post_init__(self):
        rows = tuple(tuple(row) for row in self.rotation)
        if len(rows) != _DIMENSION or any(len(row) != _DIMENSION for row in rows):
            raise ValueError("rotation must be a 3x3 integer matrix")
        if any(
            not isinstance(value, int) or isinstance(value, bool)
            for row in rows
            for value in row
        ):
            raise ValueError("rotation must be a 3x3 integer matrix")
        translation = tuple(
            _fraction(value) for value in _validate_vector(self.translation, "translation")
        )
        object.__setattr__(self, "rotation", rows)
        object.__setattr__(self, "translation", translation)


@dataclass(frozen=True)
class PeriodicAtom:
    """An atom labeled by species at an exact fractional position."""

    species: str
    position: tuple

    def __post_init__(self):
        if not isinstance(self.species, str) or not self.species:
            raise ValueError("atom species must be a non-empty string")
        position = tuple(
            _fraction(value) % 1
            for value in _validate_vector(self.position, "atom position")
        )
        object.__setattr__(self, "position", position)


@dataclass(frozen=True)
class IntegerGridMap:
    """An integer affine permutation of a finite periodic grid."""

    mesh: tuple
    coefficients: tuple
    shift: tuple

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

    grid_map = IntegerGridMap(mesh, tuple(coefficients), tuple(shift))
    images = {
        grid_map.map_index(index)
        for index in itertools.product(*(range(size) for size in mesh))
    }
    point_count = mesh[0] * mesh[1] * mesh[2]
    if len(images) != point_count:
        raise ValueError(
            "symmetry does not induce a bijection of the periodic discrete mesh"
        )
    return grid_map


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
