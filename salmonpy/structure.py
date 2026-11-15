"""Optional ASE structure input and explicit grid-width resolution; no execution."""
from dataclasses import dataclass
from copy import deepcopy
import hashlib
import io
import math
from pathlib import Path
import numpy as np

from .errors import InputError, InputTypeError

__all__ = ['StructureInput', 'GridSelection', 'grid_from_spacing']
_BOHR_ANGSTROM = 0.52917721067  # SALMON input conversion, also used by salmonpy.ase.


def _validate_occupancy(atoms):
    occupancy = atoms.info.get('occupancy')
    if occupancy is None:
        return
    if not isinstance(occupancy, dict):
        raise InputError('unsupported occupancy metadata; provide an ordered structure')
    for site in occupancy.values():
        if not isinstance(site, dict) or len(site) != 1:
            raise InputError('mixed occupancy is unsupported; provide an explicit ordered structure')
        value = next(iter(site.values()))
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)) or value != 1:
            raise InputError('partial occupancy is unsupported; provide an explicit ordered structure')


def _copy_atoms(atoms):
    result = atoms.copy()
    result.info = deepcopy(atoms.info)
    return result


@dataclass(frozen=True, init=False)
class StructureInput:
    """Detached ASE geometry. Accessing atoms returns a new copy each time."""
    _atoms: object
    _metadata: tuple

    def __init__(self):
        raise InputError('use StructureInput.from_atoms or from_file')

    @classmethod
    def from_atoms(cls, atoms):
        from ase import Atoms
        if not isinstance(atoms, Atoms):
            raise InputTypeError('expected ASE Atoms')
        if (not len(atoms) or not np.isfinite(atoms.positions).all()
                or not np.isfinite(atoms.cell.array).all()):
            raise InputError('nonempty finite structure required')
        _validate_occupancy(atoms)
        value = object.__new__(cls)
        object.__setattr__(value, '_atoms', _copy_atoms(atoms))
        object.__setattr__(value, '_metadata', (('source', 'ASE Atoms'),))
        return value

    @classmethod
    def from_file(cls, path, *, format, index=0):
        """Read one local CIF/extxyz frame from hashed bytes, without ASE MPI.

        No repetition/wrapping/standardization. CIF symmetry expansion is part
        of ASE parsing. File hash identifies all bytes, not just the chosen frame.
        """
        if format not in ('cif', 'extxyz'):
            raise InputError('explicit format must be cif or extxyz')
        if type(index) is not int or index < 0:
            raise InputError('index must be a nonnegative frame integer')
        path = Path(path)
        data = path.read_bytes()
        from ase.io import read
        stream = io.BytesIO(data) if format == 'cif' else io.StringIO(data.decode('utf-8'))
        atoms = read(stream, format=format, index=index, parallel=False)
        result = cls.from_atoms(atoms)
        object.__setattr__(result, '_metadata', (
            ('source', 'local file'), ('filename', path.name), ('format', format),
            ('index', index), ('sha256', hashlib.sha256(data).hexdigest())))
        return result

    @property
    def atoms(self):
        return _copy_atoms(self._atoms)

    def record(self):
        """Detached summary/provenance; not a full geometry or restart archive."""
        return {**dict(self._metadata), 'atom_count': len(self._atoms),
                'cell_angstrom': self._atoms.cell.array.tolist(),
                'pbc': self._atoms.pbc.tolist()}


@dataclass(frozen=True, init=False)
class GridSelection:
    counts: tuple
    requested_spacing_angstrom: tuple
    actual_spacing_angstrom: tuple
    minimum: int
    multiples: tuple

    def __init__(self):
        raise InputError('use grid_from_spacing to create a GridSelection')

    def as_dict(self):
        return {'counts': list(self.counts),
                'requested_spacing_angstrom': list(self.requested_spacing_angstrom),
                'actual_spacing_angstrom': list(self.actual_spacing_angstrom),
                'minimum': self.minimum, 'multiples': list(self.multiples),
                'total_points': math.prod(self.counts)}


def grid_from_spacing(cell_angstrom, *, spacing, unit, minimum=8, multiples=(1, 1, 1)):
    """Ceiling grid for positive axis-aligned cells, with explicit alignment.

    Width means a maximum. No MPI layout is inferred; counts are bounded by the
    existing DFT API's 2**20 per-axis limit, not a memory or convergence guarantee.
    Caller must separately require periodic boundaries for electronic DFT use.
    """
    if unit not in ('angstrom', 'bohr'):
        raise InputError('grid spacing unit must be angstrom or bohr')
    try:
        cell = np.array(cell_angstrom, dtype=float, copy=True)
        widths = np.asarray(spacing)
        if widths.dtype.kind not in 'iuf' or widths.shape not in ((), (3,)):
            raise InputError('grid spacing must be numeric, not boolean')
        widths = np.broadcast_to(widths, (3,)).astype(float, copy=True)
    except (TypeError, ValueError) as exc:
        raise InputError('expected a 3x3 cell and scalar or three numeric spacings') from exc
    if (cell.shape != (3, 3) or not np.isfinite(cell).all()
            or np.any(np.diag(cell) <= 0)
            or np.any(cell != np.diag(np.diag(cell)))):
        raise InputError('width resolution requires positive axis-aligned orthogonal cell')
    if unit == 'bohr':
        widths *= _BOHR_ANGSTROM
    if not np.isfinite(widths).all() or np.any(widths <= 0):
        raise InputError('grid spacing must be finite and positive')
    if type(minimum) is not int or not 1 <= minimum <= 2**20:
        raise InputError('minimum must be an integer in [1, 2**20]')
    if (not isinstance(multiples, (list, tuple)) or len(multiples) != 3
            or any(type(n) is not int or not 1 <= n <= 2**20 for n in multiples)):
        raise InputError('multiples must be three positive integers <= 2**20')
    counts = []
    for length, width, multiple in zip(np.diag(cell), widths, multiples):
        ratio = float(length) / float(width)
        if not math.isfinite(ratio) or ratio > 2**20:
            raise InputError('requested grid exceeds per-axis DFT limit')
        count = max(minimum, math.ceil(ratio))
        count = ((count + multiple - 1) // multiple) * multiple
        # Ensure the stated maximum even at a floating-point division boundary.
        if length / count > width:
            count += multiple
        if count > 2**20:
            raise InputError('aligned grid exceeds per-axis DFT limit')
        counts.append(count)
    result = object.__new__(GridSelection)
    values = {'counts': tuple(counts),
              'requested_spacing_angstrom': tuple(float(v) for v in widths),
              'actual_spacing_angstrom': tuple(float(v/n) for v, n in zip(np.diag(cell), counts)),
              'minimum': minimum, 'multiples': tuple(multiples)}
    for key, value in values.items():
        object.__setattr__(result, key, value)
    return result
