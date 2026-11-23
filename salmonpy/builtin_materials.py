"""Explicit LD/SBE/DFT configuration adapters; registration never initializes native code."""
from dataclasses import dataclass
import hashlib
from pathlib import Path

from .config import MaterialCatalog, MaterialRegistry
from .errors import InputError
from .preparation import _relative

__all__ = ['builtin_material_registry', 'material_from_catalog',
           'ConfiguredMaterial', 'ConfiguredState', 'ConfiguredGroundState']
_VERSION = '1'
_LD_KEYS = {'units', 'plasma_frequency', 'strength', 'resonance', 'damping'}
_SBE_KEYS = {'units', 'directory', 'sysname', 'bands', 'kpoints', 'electrons', 'lattice'}
_SBE_OPTIONS = {'correction_order', 'nonlocal_correction', 'active_bands'}


def _keys(parameters, required, optional=()):
    if type(parameters) is not dict or not required <= parameters.keys() or parameters.keys() - required - set(optional):
        raise InputError('missing or unknown built-in material parameters')


def _ld(parameters):
    from .ld import LorentzDrude
    _keys(parameters, _LD_KEYS, {"epsilon_infinity","valid_frequency_range"})
    if parameters['units'] != 'atomic_gaussian':
        raise InputError('LD units must explicitly be atomic_gaussian')
    model = LorentzDrude(**{k: v for k, v in parameters.items() if k != 'units'})
    c = model.config
    return model, {'units': 'atomic_gaussian', 'plasma_frequency': c.plasma_frequency,
                   'strength': list(c.strength), 'resonance': list(c.resonance),
                   'damping': list(c.damping), **({'epsilon_infinity':c.epsilon_infinity} if c.epsilon_infinity!=1.0 else {}),
                   **({'valid_frequency_range':list(c.valid_frequency_range)} if c.valid_frequency_range is not None else {})}


def _hash(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _sbe(parameters, base_directory):
    from .sbe import SBEInitialState
    _keys(parameters, _SBE_KEYS, _SBE_OPTIONS)
    if parameters['units'] != 'atomic':
        raise InputError('SBE units must explicitly be atomic (lattice in bohr)')
    directory, stem = parameters['directory'], parameters['sysname']
    if directory != '.':
        _relative(directory)
    if len(_relative(stem).parts) != 1:
        raise InputError('SBE sysname must be a filename stem')
    if type(parameters.get('nonlocal_correction', False)) is not bool:
        raise InputError('nonlocal_correction must be boolean')
    base = Path(base_directory).resolve()
    root = (base/directory).resolve()
    if not root.is_relative_to(base):
        raise InputError('SBE directory escapes base_directory')
    paths = {}
    for suffix in ('_k.data', '_eigen.data', '_tm.data'):
        path = root/(stem+suffix)
        if not path.resolve().is_relative_to(base):
            raise InputError('SBE file escapes base_directory')
        name = (Path(directory)/(stem+suffix)).as_posix()
        paths[name] = path
    before = {name: _hash(path) for name, path in paths.items()}
    kwargs = {k: v for k, v in parameters.items() if k not in ('units', 'directory', 'sysname')}
    model = SBEInitialState.from_salmon(root, stem, **kwargs)
    if before != {name: _hash(path) for name, path in paths.items()}:
        raise InputError('SBE files changed while reading')
    resolved = {**parameters, 'lattice': model.lattice.tolist(),
                'active_bands': model.eigen.shape[0],
                'correction_order': model.correction_order,
                'nonlocal_correction': model.nonlocal_correction,
                'file_sha256': before}
    return model, resolved


def _local_file(name, base_directory):
    _relative(name)
    base = Path(base_directory).resolve()
    path = (base/name).resolve()
    if not path.is_relative_to(base) or not path.is_file():
        raise InputError(f'expected a local file inside base_directory: {name!r}')
    return path


def _dft(parameters, base_directory):
    from .ase import dft_from_atoms
    from .dft import Pseudopotential
    from .structure import StructureInput, grid_from_spacing
    _keys(parameters, {'structure', 'pseudopotentials', 'xc', 'states', 'kmesh', 'grid'}, {'electrons','cval'})
    structure_ref = parameters['structure']
    _keys(structure_ref, {'path', 'format'}, {'index'})
    refs = parameters['pseudopotentials']
    if type(refs) is not list or not refs:
        raise InputError('pseudopotentials must be a nonempty list')
    paths = {structure_ref['path']: _local_file(structure_ref['path'], base_directory)}
    for ref in refs:
        _keys(ref, {'path', 'atomic_number', 'local_channel'})
        paths[ref['path']] = _local_file(ref['path'], base_directory)
    before = {name: _hash(path) for name, path in paths.items()}
    structure = StructureInput.from_file(paths[structure_ref['path']],
                                        format=structure_ref['format'], index=structure_ref.get('index', 0))
    atoms = structure.atoms
    pseudos = [Pseudopotential.from_file(paths[ref['path']], atomic_number=ref['atomic_number'],
                                       local_channel=ref['local_channel']) for ref in refs]
    for pseudo in pseudos:
        if pseudo.suffix == '.upf':
            import xml.etree.ElementTree as ET
            from ase.data import atomic_numbers
            element = ET.fromstring(pseudo.data).find('PP_HEADER').get('element', '').strip().capitalize()
            if atomic_numbers.get(element) != pseudo.atomic_number:
                raise InputError('UPF header element does not match declared atomic_number')
    grid = parameters['grid']
    if type(grid) is not dict or set(grid) not in ({'counts'}, {'spacing'}):
        raise InputError('grid must specify exactly counts or spacing')
    selection = None
    if 'spacing' in grid:
        _keys(grid['spacing'], {'value', 'unit'})
        selection = grid_from_spacing(atoms.cell.array, spacing=grid['spacing']['value'],
                                      unit=grid['spacing']['unit'])
        counts = selection.counts
    else:
        counts = grid['counts']
    if type(parameters['xc']) is not str:
        raise InputError('xc must be explicit string')
    model = dft_from_atoms(atoms, pseudopotentials=pseudos, states=parameters['states'],
                           grid=counts, kmesh=parameters['kmesh'], xc=parameters['xc'],
                           electrons=parameters.get('electrons'),cval=parameters.get('cval'))
    if before != {name: _hash(path) for name, path in paths.items()}:
        raise InputError('DFT files changed while reading')
    resolved = {**parameters, 'structure': {**structure_ref, 'index': structure_ref.get('index', 0)},
                'electrons': model.electrons, 'states': model.states, 'kmesh': list(model.kmesh),
                'xc': model.xc, 'file_sha256': before, 'structure_record': structure.record(),
                'grid_resolution': selection.as_dict() if selection is not None else
                {'mode': 'explicit_counts', 'counts': list(model.grid)}}
    if model.xc=='tbmbj':resolved['cval']=model.cval
    return model, resolved


def builtin_material_registry():
    """Return a new explicit registry containing LD/SBE/DFT parameter resolvers.

    Custom resolvers may be added normally; material_from_catalog only dispatches
    the built-ins. No solver execution, library loading or compilation here.
    """
    registry = MaterialRegistry()
    registry.register('salmonpy.ld', version=_VERSION, models=['ld'],
                      resolve_parameters=lambda model, parameters, **kw: _ld(parameters)[1])
    registry.register('salmonpy.sbe', version=_VERSION, models=['sbe'],
                      resolve_parameters=lambda model, parameters, *, base_directory:
                      _sbe(parameters, base_directory)[1])
    registry.register('salmonpy.dft', version=_VERSION, models=['dft'],
                      resolve_parameters=lambda model, parameters, *, base_directory:
                      _dft(parameters, base_directory)[1])
    return registry


@dataclass(frozen=True, init=False)
class ConfiguredState:
    """Identity-carrying proxy; numerical calls retain the original model contract."""
    material_id: str
    definition_digest: str
    _state: object
    initial_state_id: str | None = None
    initial_state_digest: str | None = None

    def __init__(self):
        raise InputError('use ConfiguredMaterial.new_state')

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return getattr(self._state, name)

    def __enter__(self):
        self._state.__enter__()
        return self

    def __exit__(self, *args):
        return self._state.__exit__(*args)


def _make(cls, **values):
    result = object.__new__(cls)
    for key, value in values.items():
        object.__setattr__(result, key, value)
    return result


@dataclass(frozen=True, init=False)
class ConfiguredMaterial:
    """Resolved model and ID; obtain through material_from_catalog or the GS bridge.

    model exposes the original LD/SBE/DFT/TDDFT object. new_state delegates arguments and
    returns existing independent time-evolution states with immutable identity fields.
    No initial_states recipe, point/run ID, MPI or output/restart propagation is
    selected implicitly. LD remains E-driven; SBE remains A/c-driven.
    """
    material_id: str
    definition_digest: str
    model: object
    initial_state_id: str | None = None
    initial_state_digest: str | None = None

    def __init__(self):
        raise InputError('use material_from_catalog')

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return getattr(self.model, name)

    def new_state(self, **kwargs):
        if not callable(getattr(self.model, 'new_state', None)):
            raise InputError('this model has no time-evolution state; run DFT and prepare TDDFT first')
        return _make(ConfiguredState, material_id=self.material_id,
                     definition_digest=self.definition_digest,
                     _state=self.model.new_state(**kwargs),
                     initial_state_id=self.initial_state_id, initial_state_digest=self.initial_state_digest)

    def record(self):
        result = {'material_id': self.material_id, 'definition_digest': self.definition_digest}
        if self.initial_state_id is not None:
            result.update(initial_state_id=self.initial_state_id, initial_state_digest=self.initial_state_digest)
        return result

    def run(self, *args, **kwargs):
        from .dft import DFT
        if not isinstance(self.model, DFT):
            raise InputError('run is available only for a configured DFT model')
        if isinstance(kwargs.get('initial_guess'), ConfiguredGroundState):
            kwargs['initial_guess'] = kwargs['initial_guess'].result
        result = self.model.run(*args, **kwargs)
        return _make(ConfiguredGroundState, material_id=self.material_id,
                     definition_digest=self.definition_digest, result=result)

    def close(self):
        close = getattr(self.model, 'close', None)
        if close is not None:
            close()

    def __enter__(self):
        enter = getattr(self.model, '__enter__', None)
        if enter is not None:
            enter()
        return self

    def __exit__(self, *args):
        self.close()


@dataclass(frozen=True, init=False)
class ConfiguredGroundState:
    """Identity-bearing result of ConfiguredMaterial.run; no disk ID persistence."""
    material_id: str
    definition_digest: str
    result: object

    def __init__(self):
        raise InputError('use ConfiguredMaterial.run')

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return getattr(self.result, name)

    def record(self):
        return {'material_id': self.material_id, 'definition_digest': self.definition_digest}

    def tddft_initial_state(self, **kwargs):
        initial = self.result.tddft_initial_state(**kwargs)
        return _make(ConfiguredMaterial, material_id=self.material_id,
                     definition_digest=self.definition_digest, model=initial)


def material_from_catalog(catalog, material_id, *, base_directory):
    """Build an LD/SBE/DFT model after revalidating resolved parameters/files.

    This factory does not execute calculations. DFT uses run; evolution uses
    new_state. OS/file errors propagate.
    Initial-state source/prepare recipes in the catalog are not consumed; these
    models use their existing zero-history LD / closed-shell SBE initial states.
    """
    if not isinstance(catalog, MaterialCatalog) or not isinstance(material_id, str):
        raise InputError('expected MaterialCatalog and string material_id')
    definitions = catalog.as_dict()['materials']
    if material_id not in definitions:
        raise InputError(f'unknown material_id: {material_id!r}')
    definition = definitions[material_id]
    if definition['backend_version'] != _VERSION:
        raise InputError('unsupported built-in adapter version')
    parameters = definition['parameters']
    pair = (definition['backend'], definition['model'])
    if pair == ('salmonpy.ld', 'ld'):
        model, resolved = _ld(parameters)
    elif pair == ('salmonpy.dft', 'dft'):
        model, resolved = _dft({k: v for k, v in parameters.items()
                                if k not in ('file_sha256', 'structure_record', 'grid_resolution')}, base_directory)
    elif pair == ('salmonpy.sbe', 'sbe'):
        model, resolved = _sbe({k: v for k, v in parameters.items() if k != 'file_sha256'}, base_directory)
    else:
        raise InputError('material_from_catalog supports only registered built-in LD/SBE/DFT adapters')
    if resolved != parameters:
        raise InputError('resolved material parameters or source hashes changed; resolve a new catalog explicitly')
    return _make(ConfiguredMaterial, material_id=material_id,
                 definition_digest=catalog.definition_digests[material_id], model=model)
