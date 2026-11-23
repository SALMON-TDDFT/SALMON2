"""Explicit catalog selection of a pinned SALMONpy GS/1 for TDDFT."""
import hashlib
import json
import math
from pathlib import Path
import re
import shlex

from .builtin_materials import ConfiguredMaterial, material_from_catalog, _hash, _keys, _make
from .config import MaterialCatalog, _copy_json
from .dft import DFT, GroundStateResult
from .errors import InputError, InputTypeError
from .fdtd import SalmonExecutable
from .preparation import _relative

__all__ = ['tddft_from_catalog']


def _manifest(path, expected):
    with path.open('rb') as stream:
        content = stream.read(16 * 1024 * 1024 + 1)
    if len(content) > 16 * 1024 * 1024:
        raise InputError('GS manifest exceeds 16 MiB')
    if hashlib.sha256(content).hexdigest() != expected:
        raise InputError('GS manifest hash mismatch')
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise InputError('duplicate GS manifest key')
            result[key] = value
        return result
    try:
        record = _copy_json(json.loads(content, object_pairs_hook=pairs))
    except (ValueError, RecursionError) as exc:
        raise InputError('invalid GS manifest') from exc
    if type(record) is not dict:
        raise InputError('GS manifest must be an object')
    if record.get('schema') != 'salmonpy.gs/1' or record.get('status') != 'complete':
        raise InputError('expected a completed GS/1 manifest')
    files = record.get('files')
    if type(files) is not dict or any(
        type(value) is not str or not re.fullmatch('[0-9a-f]{64}', value)
        for value in files.values()
    ):
        raise InputError('GS files must be an object of SHA256 digests')
    energy = record.get('energy_hartree')
    if type(energy) not in (int, float) or not math.isfinite(energy):
        raise InputError('GS energy_hartree must be finite')
    iterations = record.get('iterations')
    if type(iterations) is not int or iterations < 1:
        raise InputError('GS iterations must be a positive integer')
    return record


def _check_input(model, root):
    """Compare actual GS scientific input, not only manifest declarations."""
    import f90nml
    from .bands import _input_groups
    expected = _input_groups(model.input_text())
    actual = _input_groups((root/'native/gs.inp').read_text())
    if set(actual) != set(expected):
        raise InputError('GS input groups mismatch with configured DFT')
    ignored = {'control', 'scf', 'atomic_red_coor'}
    try:
        for group in set(expected) - ignored:
            wanted = f90nml.reads('&'+group+expected[group]+'/')
            found = f90nml.reads('&'+group+actual[group]+'/')
            if wanted != found:
                raise InputError(f'GS input mismatch in {group}')
        def rows(body):
            parsed = []
            for line in body.splitlines():
                tokens = shlex.split(line.split('!', 1)[0])
                if not tokens:
                    continue
                if len(tokens) != 5:
                    raise InputError('invalid GS input atomic row')
                parsed.append((int(tokens[4]), *(float(x.replace('D', 'e').replace('d', 'e')) for x in tokens[1:4])))
            return parsed
        if rows(actual['atomic_red_coor']) != rows(expected['atomic_red_coor']):
            raise InputError('GS input atomic positions/species mismatch')
    except (ValueError, KeyError, IndexError) as exc:
        raise InputError(f'GS input does not match configured DFT: {exc}') from exc
    for index, pseudo in enumerate(model.pseudopotentials):
        if _hash(root/'native'/model._pseudo_name(index)) != pseudo.sha256:
            raise InputError('GS pseudo bytes mismatch with configured material')


def tddft_from_catalog(catalog, initial_state_id, *, base_directory, executable,
                       max_steps=1000, threads=1, timeout=300.):
    """Select a local pinned GS/1, verify it against the DFT material, load TDDFT.

    Source must be {format: salmonpy.gs/1, path: relative-directory,
    manifest_sha256: lowercase-SHA256}. Requires the same executable bytes used
    for GS creation. No native worker starts until new_state; no downloading or
    SCF is implicit. Exact geometry/order/discretization matching is deliberate.
    """
    from .sbe import _integer
    from .ld import _positive
    if not isinstance(catalog, MaterialCatalog) or not isinstance(initial_state_id, str):
        raise InputError('expected MaterialCatalog and string initial_state_id')
    if not isinstance(executable, SalmonExecutable):
        raise InputTypeError('expected explicit SalmonExecutable')
    max_steps = _integer(max_steps, 'max_steps', 1, 2**31-2)
    threads = _integer(threads, 'threads', 1, 2**31-1)
    timeout = _positive(timeout, 'timeout')
    entries = catalog.as_dict()['initial_states']
    if initial_state_id not in entries:
        raise InputError(f'unknown initial_state_id: {initial_state_id!r}')
    entry = entries[initial_state_id]
    if set(entry) != {'material_id', 'source'}:
        raise InputError('TDDFT selection requires an explicit saved GS source; prepare is unsupported')
    source = entry['source']
    _keys(source, {'format', 'path', 'manifest_sha256'})
    if source['format'] != 'salmonpy.gs/1':
        raise InputError('only saved SALMONpy GS/1 sources are supported')
    digest = source['manifest_sha256']
    if type(digest) is not str or not re.fullmatch('[0-9a-f]{64}', digest):
        raise InputError('manifest_sha256 must be lowercase SHA256')
    _relative(source['path'])
    base = Path(base_directory).resolve()
    candidate = base/source['path']
    root = candidate.resolve()
    if candidate.is_symlink() or not root.is_relative_to(base) or not root.is_dir():
        raise InputError('GS source must be a local directory inside base_directory')
    path = root/'ground-state.json'
    if path.is_symlink():
        raise InputError('GS manifest symlinks are unsupported')
    record = _manifest(path, digest)
    if record.get('executable_sha256') != executable.sha256 or _hash(executable.path) != executable.sha256:
        raise InputError('GS executable identity mismatch; native restart conversion is not implicit')
    material = material_from_catalog(catalog, entry['material_id'], base_directory=base)
    if not isinstance(material.model, DFT):
        raise InputError('saved GS selection requires a configured DFT material')
    model = material.model
    wanted = {**model._restart_signature(), 'fractional_positions': model.fractional_positions.tolist()}
    for name, value in wanted.items():
        if record.get(name) != value:
            raise InputError(f'GS material mismatch: {name}')
    gs = GroundStateResult.open(root)
    _check_input(model, root)
    initial = gs.tddft_initial_state(executable=executable, max_steps=max_steps, threads=threads, timeout=timeout)
    try:
        snapshot = dict(initial.manifest)
        expected = {name.removeprefix('native/'): value for name, value in record['files'].items()}
        if any(not name.startswith('native/') for name in record['files']):
            raise InputError('unsupported GS file inventory')
        if set(snapshot) != set(expected) | {'rt.inp'} or any(snapshot.get(name) != value for name, value in expected.items()):
            raise InputError('copied GS snapshot mismatch')
        if _hash(path) != digest:
            raise InputError('GS manifest changed during loading')
        return _make(ConfiguredMaterial, material_id=material.material_id,
                     definition_digest=material.definition_digest, model=initial,
                     initial_state_id=initial_state_id, initial_state_digest=digest)
    except BaseException:
        initial.close()
        raise
