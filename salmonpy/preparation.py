"""Explicit local input staging and offline transfer verification.

No acquisition, installation, solver execution or backend discovery. Trusted
configuration resolvers run only during prepare_inputs, not verify_inputs.
"""
import argparse
import ast
import tokenize
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
import unicodedata

from .config import _copy_json, resolve_materials
from .errors import InputError

__all__ = ['prepare_inputs', 'prepare_fdtd', 'prepare_dft', 'prepare_tddft', 'verify_inputs', 'preflight_inputs', 'check_job']
_SCHEMA = 'salmonpy.input-bundle/1'
_GENERATED = {'materials.json', 'resolved.json', 'manifest.json'}
_DIGEST = re.compile(r'[0-9a-f]{64}\Z')


def _relative(name):
    if not isinstance(name, str) or not name or '\\' in name:
        raise InputError('asset name must be a portable relative POSIX path')
    path = PurePosixPath(name)
    if path.is_absolute() or path.as_posix() != name:
        raise InputError(f'unsafe asset name: {name!r}')
    for part in name.split('/'):
        stem = part.split('.')[0].upper()
        if (part in ('', '.', '..') or part.endswith((' ', '.'))
                or any(c in '<>:"|?*' or ord(c) < 32 for c in part)
                or stem in {'CON', 'PRN', 'AUX', 'NUL'}
                or re.fullmatch(r'(COM|LPT)[1-9]', stem)):
            raise InputError(f'unsafe asset name: {name!r}')
    return path


def _names(names, *, assets=False):
    """Reject cross-platform filename and parent-directory collisions."""
    canonical = {}
    files = set(names)
    for name in names:
        path = _relative(name)
        if assets and path.parts[0].casefold() in _GENERATED:
            raise InputError(f'reserved asset name: {name}')
        for i in range(1, len(path.parts) + 1):
            prefix = '/'.join(path.parts[:i])
            key = unicodedata.normalize('NFC', prefix).casefold()
            if key in canonical and canonical[key] != prefix:
                raise InputError(f'case-colliding asset path: {name}')
            canonical[key] = prefix
            if i < len(path.parts) and prefix in files:
                raise InputError(f'file/directory collision: {name}')


def _hash(path):
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _record(path, provenance):
    digest, size = _hash(path)
    return {'sha256': digest, 'size': size, 'provenance': provenance}


def _write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, allow_nan=False,
                               sort_keys=True, indent=2) + '\n', encoding='utf-8')


def _read(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise InputError(f'duplicate JSON key: {key!r}')
            result[key] = value
        return result
    try:
        if path.stat().st_size > 16 * 1024 * 1024:
            raise InputError('input bundle manifest exceeds 16 MiB')
        return _copy_json(json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=pairs))
    except (ValueError, RecursionError) as exc:
        raise InputError(f'invalid input bundle manifest: {exc}') from exc


def _publish_stage(stage,dest):
    """Reserve a new destination; never delete a pre-existing directory."""
    dest.mkdir(exist_ok=False)
    try:
        for child in stage.iterdir():child.rename(dest/child.name)
    except BaseException:
        shutil.rmtree(dest)
        raise


def _prepare_native(directory,kind,input_name,input_data,assets,metadata,*,expected_hashes=None):
    """Common input-only publication for native payloads and explicit files."""
    dest=Path(directory).expanduser().absolute()
    if dest.exists() or dest.is_symlink():raise FileExistsError(dest)
    manifest=f'salmonpy-prepared-{kind}.json'
    if set(assets)&{input_name,manifest}:raise InputError('asset collides with a prepared input filename')
    _names([*assets,input_name,manifest])
    payloads={input_name:input_data,**assets}
    dest.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.salmonpy-{kind}-',dir=dest.parent) as temporary:
        stage=Path(temporary);entries={}
        for name,data in payloads.items():
            target=stage/name;target.parent.mkdir(parents=True,exist_ok=True)
            if isinstance(data,bytes):target.write_bytes(data)
            else:shutil.copyfile(data,target)
            digest,size=_hash(target)
            if expected_hashes is not None and name in expected_hashes and digest!=expected_hashes[name]:
                raise InputError('initial guess changed while staging: '+name)
            entries[name]={'sha256':digest,'size':size}
        _write(stage/manifest,{'schema':f'salmonpy.prepared-{kind}/1','status':'prepared-not-run',
                              'input':input_name,**metadata,'files':entries})
        _publish_stage(stage,dest)
    return dest


def prepare_fdtd(calculation,directory,*,mpi=None):
    """Write native FDTD inputs only; return a NEW directory, never execute."""
    from .fdtd import FDTD,FDTDMPI
    from .errors import InputTypeError
    if not isinstance(calculation,FDTD):raise InputTypeError('calculation must be FDTD')
    if mpi is not None and not isinstance(mpi,FDTDMPI):raise InputTypeError('mpi must be FDTDMPI')
    dest=Path(directory).expanduser().absolute()
    if dest.exists() or dest.is_symlink():raise FileExistsError(dest)
    text=calculation.input_text() if mpi is None else mpi.input_text(calculation)
    assets={name:value.encode('utf-8') for name,value in calculation.input_assets().items()}
    return _prepare_native(directory,'fdtd','inputfile',text.encode('utf-8'),assets,
                           {'ranks':1 if mpi is None else mpi.ranks,
                            'mpi_grid':None if mpi is None else list(mpi.grid)})


def prepare_dft(calculation,directory,*,scf=None,save_restart=True,initial_guess=None):
    """Write a fresh DFT input and immutable PP/symmetry assets, without SCF.

    Optional explicit initial_guess stages the authenticated wavefunction seed.
    No runtime is loaded or started; its required checksum is recorded.
    """
    from .dft import DFT,SCFSettings,_restart_output_text,_initial_guess_source,_initial_guess_text,_GUESS_FILES
    from .errors import InputTypeError
    if not isinstance(calculation,DFT):raise InputTypeError('calculation must be DFT')
    scf=SCFSettings() if scf is None else scf
    text=_restart_output_text(calculation.input_text(scf),save_restart)
    metadata={'xc':calculation.xc,'save_restart':save_restart}
    if calculation.xc=='hse06':metadata['required_runtime_profile']='fast-pr'
    assets=dict(calculation.input_assets);expected={}
    if initial_guess is not None:
        source,record=_initial_guess_source(calculation,initial_guess,directory)
        text=_initial_guess_text(text)
        metadata['required_executable_sha256']=record['executable_sha256']
        metadata['initial_guess']={'mode':'wavefunction','source_directory':str(source.directory),
                                   'manifest_sha256':_hash(source.directory/'ground-state.json')[0]}
        for name in _GUESS_FILES:
            key='initial_guess/'+name;relative='native/data_for_restart/'+name
            assets[key]=source.directory/relative;expected[key]=record['files'][relative]
    return _prepare_native(directory,'dft','gs.inp',text.encode('utf-8'),assets,metadata,expected_hashes=expected)


def prepare_tddft(input_file,directory,*,files=None):
    """Stage supplied RT input bytes and explicitly selected local files only.

    No executable/worker/GS loader, discovery or namelist rewriting. Restart
    paths and dependencies remain user-managed; an offline management node
    need not have the compute node's GS directory. Copies stream to disk.
    """
    source=Path(input_file).expanduser()
    if not source.is_file():raise InputError('TDDFT input_file must be a local regular file')
    files={} if files is None else files
    if type(files) is not dict:raise InputError('files must map relative names to explicit local paths')
    assets={}
    for name,path in files.items():
        path=Path(path).expanduser()
        if not path.is_file():raise InputError(f'missing local TDDFT asset: {name}')
        assets[name]=path
    return _prepare_native(directory,'tddft','rt.inp',source,assets,
                           {'restart_policy':'external-user-managed','validation_scope':'explicit-files-only'})


def prepare_inputs(data, directory, *, registry, base_directory, files, provenance=None):
    """Prepare a NEW local directory and return its absolute Path.

    ``files`` maps desired relative paths to explicitly selected local source
    paths (relative to base_directory unless absolute). No directory recursion,
    symlinks, URL fetching or inferred file dependencies. ``provenance`` optionally
    maps the same paths to JSON objects; omitted provenance is marked unknown.
    Model resolvers receive the staged directory as base_directory and must check
    their references there. This is not an automatic dependency/physics validator.
    """
    dest = Path(directory).absolute()
    if dest.exists() or dest.is_symlink():
        raise FileExistsError(dest)
    if type(files) is not dict:
        raise InputError('files: expected relative-name to source-path mapping')
    _names(files, assets=True)
    try:
        recipe = _copy_json(data)
        metadata = _copy_json({} if provenance is None else provenance)
    except RecursionError as exc:
        raise InputError('input nesting is too deep') from exc
    if type(metadata) is not dict or set(metadata) - set(files):
        raise InputError('provenance: expected metadata only for declared assets')
    if any(type(value) is not dict for value in metadata.values()):
        raise InputError('provenance: each entry must be a JSON object')
    base = Path(base_directory).resolve()
    sources = {}
    missing = []
    for name, source in files.items():
        try:
            path = base / source
        except TypeError as exc:
            raise InputError(f'files[{name!r}]: expected local source path') from exc
        if path.is_symlink() or not path.is_file():
            missing.append(f'{name}: {source}')
        else:
            sources[name] = path
    if missing:
        raise InputError('missing or non-regular input files: ' + '; '.join(missing))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.salmonpy-input-', dir=dest.parent) as temporary:
        stage = Path(temporary)
        entries = {}
        for name, source in sources.items():
            target = stage/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            entries[name] = _record(target, metadata.get(name, {'source': 'unknown'}))
        catalog = resolve_materials(recipe, registry=registry, base_directory=stage)
        _write(stage/'materials.json', recipe)
        _write(stage/'resolved.json', catalog.as_dict())
        for name in ('materials.json', 'resolved.json'):
            entries[name] = _record(stage/name, {'source': 'salmonpy.preparation'})
        _write(stage/'manifest.json', {'schema': _SCHEMA, 'files': entries})
        verify_inputs(stage)
        # mkdir reserves the new name; never overwrite even if another job won a race.
        _publish_stage(stage,dest)
    return dest


def verify_inputs(directory, *, manifest_sha256=None):
    """Check inventory, sizes and SHA256 without importing/calling any backend.

    Returns a detached manifest dictionary. Optional externally recorded manifest
    hash anchors the transferred manifest. Without it this is an internal
    consistency check, not authentication. Does not certify solver compatibility,
    semantic completeness, physical accuracy or package dependencies.
    """
    root = Path(directory).absolute()
    if root.is_symlink() or not root.is_dir():
        raise InputError('input bundle must be a real directory')
    manifest = root/'manifest.json'
    if manifest.is_symlink() or not manifest.is_file():
        raise InputError('missing or non-regular manifest.json')
    if manifest_sha256 is not None:
        if not isinstance(manifest_sha256, str) or not _DIGEST.fullmatch(manifest_sha256):
            raise InputError('manifest_sha256: expected lowercase SHA256')
        if _hash(manifest)[0] != manifest_sha256:
            raise InputError('manifest SHA256 mismatch')
    data = _read(manifest)
    if type(data) is not dict or set(data) != {'schema', 'files'} or data['schema'] != _SCHEMA:
        raise InputError('unsupported input bundle manifest')
    entries = data['files']
    if type(entries) is not dict or not {'materials.json', 'resolved.json'} <= entries.keys():
        raise InputError('manifest must inventory materials.json and resolved.json')
    if 'manifest.json' in entries:
        raise InputError('manifest cannot inventory itself')
    _names([*entries, 'manifest.json'])
    directories = set()
    for name, entry in entries.items():
        if (type(entry) is not dict or set(entry) != {'sha256', 'size', 'provenance'}
                or type(entry['size']) is not int or entry['size'] < 0
                or type(entry['sha256']) is not str or not _DIGEST.fullmatch(entry['sha256'])
                or type(entry['provenance']) is not dict):
            raise InputError(f'invalid manifest entry: {name}')
        directories.update(str(parent) for parent in PurePosixPath(name).parents if str(parent) != '.')
    found = set()
    for path in root.rglob('*'):
        name = path.relative_to(root).as_posix()
        if path.is_symlink():
            raise InputError(f'symlink in input bundle: {name}')
        if path.is_dir():
            if name not in directories:
                raise InputError(f'unexpected directory: {name}')
        elif path.is_file():
            found.add(name)
        else:
            raise InputError(f'non-regular bundle member: {name}')
    expected = set(entries) | {'manifest.json'}
    if found != expected:
        raise InputError(f'input inventory mismatch; missing={sorted(expected-found)}, extra={sorted(found-expected)}')
    for name, entry in entries.items():
        actual, size = _hash(root/name)
        if actual != entry['sha256'] or size != entry['size']:
            raise InputError(f'input checksum/size mismatch: {name}')
    return data


def preflight_inputs(directory, *, manifest_sha256=None, registry=None):
    """Return a JSON-ready input-only report, without running a solver.

    With an explicit trusted registry, re-resolve the staged material recipe
    and compare it with the prepared definitions. No registry means file
    integrity only. This does not validate native namelist semantics, runtime
    availability, MPI layout, GS compatibility, convergence or memory capacity.
    """
    root = Path(directory).absolute()
    report = {'schema': 'salmonpy.input-preflight/1', 'verified': False,
              'directory': str(root), 'status': 'failed', 'files': 0,
              'scope': 'file integrity only' if registry is None else
                       'file integrity and registered material parameters',
              'errors': [], 'solver_started': False,
              'limitations': ['Native parser/runtime and MPI layout are not checked.',
                              'GS compatibility, physical convergence and memory capacity are not checked.']}
    phase = 'integrity'
    try:
        manifest = verify_inputs(root, manifest_sha256=manifest_sha256)
        report['files'] = len(manifest['files'])
        report['manifest_sha256'] = _hash(root/'manifest.json')[0]
        if registry is not None:
            phase = 'materials'
            catalog = resolve_materials(_read(root/'materials.json'),
                                        registry=registry, base_directory=root)
            if catalog.as_dict() != _read(root/'resolved.json'):
                raise InputError('material definitions differ from prepared resolved.json; prepare a new bundle')
        report['verified'] = True
        report['status'] = 'passed'
    except Exception as exc:
        report['errors'].append({'phase': phase, 'type': type(exc).__name__,
                                 'message': str(exc)})
    return report


def check_job(code, inputs, *, manifest_sha256=None):
    """Statically check a Python source and local JSON/bundle; never run either.

    Does not import the supplied code, call its main(), follow its imports,
    execute plugins, or infer dynamic file references. A passing report covers
    syntax/integrity only, not code-input agreement or a runnable calculation.
    """
    code, inputs = Path(code), Path(inputs)
    report = {'schema': 'salmonpy.job-static-check/1', 'status': 'failed',
              'code': str(code.absolute()), 'inputs': str(inputs.absolute()),
              'code_executed': False, 'solver_started': False, 'errors': [],
              'checks': {}, 'not_checked': [
                  'Agreement between Python code and input parameters or paths.',
                  'Imported modules, runtime/MPI availability and native namelist semantics.',
                  'Physical convergence, GS compatibility and memory capacity.']}
    try:
        with tokenize.open(code) as stream:
            text = stream.read()
        tree = ast.parse(text, filename=str(code))
        compile(tree, str(code), 'exec', dont_inherit=True)
        report['checks']['python'] = {'status': 'passed', 'sha256': _hash(code)[0],
                                     'scope': 'syntax only; no imports or execution'}
    except Exception as exc:
        report['errors'].append({'phase': 'python', 'type': type(exc).__name__,
                                 'message': str(exc), 'line': getattr(exc, 'lineno', None)})
    if inputs.is_dir():
        result = preflight_inputs(inputs, manifest_sha256=manifest_sha256)
        report['checks']['inputs'] = result
        report['errors'].extend(result['errors'])
    else:
        try:
            if manifest_sha256 is not None:
                raise InputError('manifest SHA256 applies to an input bundle directory only')
            _read(inputs)
            report['checks']['inputs'] = {'status': 'passed', 'sha256': _hash(inputs)[0],
                                         'scope': 'JSON syntax and finite values only; no scientific schema validation'}
        except Exception as exc:
            report['errors'].append({'phase': 'inputs', 'type': type(exc).__name__,
                                     'message': str(exc)})
    if not report['errors']:
        report['status'] = 'passed'
    return report


def _main():
    parser = argparse.ArgumentParser(description='Check staged inputs before submission or after offline transfer; no solver launch.')
    parser.add_argument('directory')
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--code', type=Path, help='Statically check this Python file with the JSON or bundle input; never execute it')
    parser.add_argument('--resolve-builtins', action='store_true',
                        help='Recheck built-in material parameters with local files; no solver execution')
    parser.add_argument('--report', type=Path, help='New JSON report outside the immutable input bundle')
    args = parser.parse_args()
    if args.report is not None:
        root = Path(args.directory).resolve()
        if args.report.resolve().is_relative_to(root):
            parser.error('--report must be outside the input bundle')
    registry = None
    if args.code is not None and args.resolve_builtins:
        parser.error('--code does not execute material resolvers; do not combine with --resolve-builtins')
    if args.resolve_builtins:
        from .builtin_materials import builtin_material_registry
        registry = builtin_material_registry()
    report = (check_job(args.code, args.directory, manifest_sha256=args.manifest_sha256)
              if args.code is not None else
              preflight_inputs(args.directory, manifest_sha256=args.manifest_sha256, registry=registry))
    text = json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + '\n'
    if args.report is not None:
        with args.report.open('x', encoding='utf-8') as stream:
            stream.write(text)
    print(text, end='')
    raise SystemExit(0 if report['status'] == 'passed' else 2)


if __name__ == '__main__':
    _main()
