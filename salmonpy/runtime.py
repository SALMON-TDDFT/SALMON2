"""Explicit acquisition of a compatible full solver; import is offline.

Catalogs are a trust boundary: use publisher-provided catalogs, not arbitrary
third-party URLs. Hashes establish integrity, not scientific correctness.
"""
from .errors import InputError, RuntimeFailure
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from urllib.parse import urlparse
from urllib.request import urlopen
import zipfile

INTERFACE = 'salmonpy-runtime/2'
_MAX_ARCHIVE = 1024 * 1024 * 1024


class RuntimeUnavailable(RuntimeFailure):
    """No compatible published binary; build is never started implicitly."""


def host_platform():
    system = platform.system().lower()
    machine = platform.machine().lower()
    machine = {'amd64': 'x86_64', 'aarch64': 'arm64'}.get(machine, machine)
    if system == 'linux':
        libc, version = platform.libc_ver()
    elif system == 'darwin':
        libc, version = '', platform.mac_ver()[0]
    else:
        libc, version = '', platform.version()
    return dict(system=system, machine=machine, libc=libc, version=version)


def _version(value):
    if not re.fullmatch(r'\d+(\.\d+)*', value):
        raise InputError('missing or invalid platform version')
    return tuple(int(p) for p in value.split('.'))


def _compatible(target, host):
    if any(target.get(k) != host.get(k) for k in ('system', 'machine', 'libc')):
        return False
    try:
        have, need = _version(host['version']), _version(target['version'])
        size = max(len(have), len(need))
        return have + (0,)*(size-len(have)) >= need + (0,)*(size-len(need))
    except (KeyError, ValueError, TypeError):
        return False


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _https(url):
    parsed = urlparse(url)
    if parsed.scheme != 'https' or not parsed.netloc or parsed.username or parsed.password:
        raise InputError('runtime downloads require an HTTPS URL without credentials')


def _download(url, target):
    _https(url)
    with urlopen(url, timeout=60) as response, Path(target).open('wb') as out:
        _https(response.url)
        size = 0
        while block := response.read(1024*1024):
            size += len(block)
            if size > _MAX_ARCHIVE:
                raise InputError('runtime archive exceeds size limit')
            out.write(block)


def _catalog(catalog=None):
    if catalog is None:
        catalog = Path(__file__).with_name('runtime_catalog.json')
    if isinstance(catalog, dict):
        data = catalog
    elif str(catalog).startswith('https://'):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'catalog.json'
            _download(str(catalog), path)
            data = json.loads(path.read_text(encoding='utf-8'))
    else:
        data = json.loads(Path(catalog).read_text(encoding='utf-8'))
    if data.get('schema') != 'salmonpy.runtimes/1' or not isinstance(data.get('artifacts'), list):
        raise InputError('unsupported runtime catalog')
    return data


def inspect_runtime(catalog=None, *, required=(), host=None):
    """List matching binaries. Default bundled catalog is read offline."""
    host = host_platform() if host is None else host
    available = [e for e in _catalog(catalog)['artifacts']
                 if e.get('interface') == INTERFACE
                 and _compatible(e.get('platform', {}), host)
                 and set(required).issubset(e.get('features', []))]
    return {'host': host, 'required': list(required), 'available': available,
            'fallback': 'python -m salmonpy.runtime build --help',
            'note': 'No automatic compilation. An empty list means no matching published asset in this catalog.'}


def _relative(name):
    if not isinstance(name, str) or not name or '\\' in name or ':' in name:
        raise InputError('unsafe runtime path')
    path = PurePosixPath(name)
    if path.is_absolute() or any(p in ('', '.', '..') for p in name.split('/')):
        raise InputError('unsafe runtime path')
    return path


def _manifest(directory, *, required=(), host=None):
    directory = Path(directory).expanduser().resolve()
    meta = json.loads((directory/'runtime.json').read_text(encoding='utf-8'))
    if meta.get('schema') != 'salmonpy.runtime/1' or meta.get('interface') != INTERFACE:
        raise InputError('unsupported runtime interface')
    if not _compatible(meta.get('platform', {}), host_platform() if host is None else host):
        raise RuntimeUnavailable('runtime platform is incompatible with this host')
    if not set(required).issubset(meta.get('features', [])):
        raise RuntimeUnavailable('runtime lacks required features')
    files = meta.get('files')
    if not isinstance(files, dict) or meta.get('executable') not in files:
        raise InputError('runtime manifest lacks executable inventory')
    for name, expected in files.items():
        path = directory.joinpath(*_relative(name).parts)
        if any(p.is_symlink() for p in [path, *path.parents] if p != directory.parent):
            raise InputError('runtime symlinks are not supported')
        if not path.is_file() or _hash(path) != expected:
            raise InputError('runtime file checksum mismatch: '+name)
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob('*') if p.is_file()}
    if actual != set(files) | {'runtime.json'}:
        raise InputError('runtime contains unlisted files')
    return meta


def load_runtime(directory, *, required=(), host=None):
    """Verify all bundle files and return the existing SalmonExecutable adapter."""
    from .fdtd import SalmonExecutable
    meta = _manifest(directory, required=required, host=host)
    path = Path(directory).expanduser().resolve()/meta['executable']
    if not os.access(path, os.X_OK):
        raise InputError('runtime executable is not executable')
    engine = SalmonExecutable(path, meta['files'][meta['executable']])
    claims = meta.get('validated_capabilities', {})
    if type(claims) is not dict: raise InputError('invalid runtime capability inventory')
    record = claims.get('nlcc-forces')
    feature = 'nlcc-forces' in meta.get('features', [])
    if feature != (record is not None): raise InputError('NLCC feature requires a validation record')
    if record is not None:
        if type(record) is not str: raise InputError('invalid NLCC validation path')
        relative = _relative(record)
        if relative.as_posix() not in meta['files']: raise InputError('unverified NLCC validation record')
        receipt = json.loads((Path(directory).expanduser().resolve()/relative).read_text())
        from ._runtime_capabilities import validate_nlcc_receipt
        scope = validate_nlcc_receipt(receipt, executable_sha256=engine.sha256,
            source_patch_sha256=meta.get('provenance', {}).get('patch_sha256'),
            source_tree_sha256=meta.get('provenance', {}).get('source_tree_sha256'))
        object.__setattr__(engine, '_nlcc_force_xc', scope)
    return engine


def load_validated_build(build_record, validation_record):
    """Load an audited serial build plus binary-bound NLCC validation.

    Records must come from trusted explicit build/validation commands. This does
    not run a solver or promise material convergence. Ordinary executables and
    old runtime bundles keep an empty NLCC capability.
    """
    from .fdtd import SalmonExecutable
    from ._runtime_capabilities import validate_nlcc_receipt
    path = Path(build_record).expanduser().resolve()
    record = json.loads(path.read_text())
    exe = Path(record.get('executable','')).resolve()
    if (record.get('status') != 'built-and-audited' or exe != path.parent/'build/salmon'
            or record.get('build_configuration',{}).get('mpi_enabled') is not False
            or not exe.is_file() or _hash(exe) != record.get('sha256')):
        raise InputError('expected unchanged audited serial runtime build')
    receipt = json.loads(Path(validation_record).read_text())
    if receipt.get('build_executable_sha256') != record['sha256']:
        raise InputError('NLCC evidence belongs to another build')
    scope = validate_nlcc_receipt(receipt, executable_sha256=record['sha256'],
                                source_patch_sha256=record.get('source',{}).get('patch_sha256'),
                                source_tree_sha256=record.get('source',{}).get('source_tree_sha256'))
    engine = SalmonExecutable(exe,record['sha256'])
    object.__setattr__(engine,'_nlcc_force_xc',scope)
    return engine


def install_runtime(directory, *, catalog=None, required=(), artifact=None, host=None, archive=None):
    """Download a matching immutable bundle into a NEW directory; never build."""
    dest = Path(directory).expanduser().resolve()
    if dest.exists():
        raise FileExistsError(dest)
    matches = inspect_runtime(catalog, required=required, host=host)['available']
    if artifact is not None:
        matches = [e for e in matches if e['id'] == artifact]
    if not matches:
        raise RuntimeUnavailable('No compatible prebuilt runtime. Use an audited existing runtime, '
                                 'or explicitly run: python -m salmonpy.runtime build --help')
    if len(matches) != 1:
        raise InputError('Multiple compatible runtimes; select an explicit artifact ID')
    entry = matches[0]
    if archive is None:
        if not entry.get('url'):
            raise RuntimeUnavailable('This CI artifact has no public download URL. Download the CI artifact and pass --archive PATH.')
        _https(entry['url'])
    local_archive = Path(archive).expanduser().resolve() if archive is not None else None
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.salmonpy-runtime-', dir=dest.parent) as tmp:
        tmp = Path(tmp)
        archive, stage = tmp/'download.zip', tmp/'runtime'
        if local_archive is None:
            _download(entry['url'], archive)
        else:
            shutil.copyfile(local_archive, archive)
        if _hash(archive) != entry['sha256']:
            raise InputError('runtime archive checksum mismatch')
        stage.mkdir()
        with zipfile.ZipFile(archive) as z:
            names = set()
            total = 0
            for info in z.infolist():
                name = info.filename
                _relative(name)
                mode = info.external_attr >> 16
                if name in names or info.is_dir() or stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in (0, stat.S_IFREG)):
                    raise InputError('unsupported or duplicate runtime archive member')
                names.add(name)
                total += info.file_size
                if total > _MAX_ARCHIVE:
                    raise InputError('expanded runtime exceeds size limit')
                target = stage.joinpath(*PurePosixPath(name).parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as source, target.open('wb') as out:
                    shutil.copyfileobj(source, out)
            meta = _manifest(stage, host=host)
            if any(meta.get(k) != entry.get(k) for k in ('interface', 'platform', 'features')):
                raise InputError('runtime catalog and manifest disagree')
            (stage/meta['executable']).chmod(0o755)
        # Reserve destination without replacing an existing runtime, even in a race.
        dest.mkdir(exist_ok=False)
        try:
            for child in stage.iterdir():
                shutil.move(str(child), dest/child.name)
            return load_runtime(dest, required=required, host=host)
        except BaseException:
            shutil.rmtree(dest)
            raise


def build_runtime(upstream=None, directory=None, *, jobs, cmake_args=()):
    """Explicit local build. Retain logs/source on failure; no dependency install.

    Returns the audited executable path and hash, not a portable distribution or
    a claim that all scientific workflows have been validated on this machine.
    """
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise InputError('jobs must be a positive integer chosen for the memory budget')
    if any(not s.startswith('-D') or '=' not in s for s in cmake_args):
        raise InputError('CMake options must be explicit -DNAME=VALUE definitions')
    resources = Path(__file__).with_name('_runtime_source')
    if not resources.exists():
        resources = Path(__file__).resolve().parents[1]
    prepare = resources/'tools/prepare_salmon_runtime.py'
    checker = resources/'tools/check_forbidden_eigensolvers.py'
    if not prepare.is_file() or not checker.is_file():
        raise RuntimeFailure('runtime build resources missing; reinstall a complete SALMONpy wheel')
    if directory is None:
        directory, upstream = upstream, None
    if directory is None:
        raise InputError('a new build directory is required')
    work = Path(directory).expanduser().resolve()
    if work.is_relative_to(resources) or (upstream is not None and work.is_relative_to(Path(upstream).expanduser().resolve())):
        raise InputError('build directory must be outside source checkouts')
    work.mkdir(parents=True, exist_ok=False)
    commands = [
        [sys.executable, str(prepare), str(work/'source')],
        [sys.executable, str(checker), '--source-root', str(work/'source/src')],
        ['cmake', '-S', str(work/'source'), '-B', str(work/'build'), '-DUSE_HSE=OFF', *cmake_args,
         '-DCMAKE_REQUIRE_FIND_PACKAGE_LAPACK=TRUE', '-DCMAKE_REQUIRE_FIND_PACKAGE_LIBXCF90=TRUE'],
        ['cmake', '--build', str(work/'build'), '--target', 'salmon', '--parallel', str(jobs)],
        [sys.executable, str(checker), '--object-tree', str(work/'build')],
    ]
    record = {'schema': 'salmonpy.runtime-build/1', 'status': 'started', 'commands': commands,
              'platform': host_platform(), 'validation': 'source/object audit only; run workflow-specific tests before production'}
    record_path = work/'build.json'
    record_path.write_text(json.dumps(record, indent=2)+'\n')
    try:
        with (work/'build.log').open('w') as log:
            for command in commands:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        executable = work/'build'/('salmon.exe' if os.name == 'nt' else 'salmon')
        config=(work/'build/config.h').read_text()
        record['build_configuration']={'hse_enabled':bool(re.search(r'^\s*#define\s+USE_HSE(?:\s|$)',config,re.M)),
                                       'mpi_enabled':bool(re.search(r'^\s*#define\s+USE_MPI(?:\s|$)',config,re.M)),
                                       'config_sha256':_hash(work/'build/config.h')}
        record.update(status='built-and-audited', executable=str(executable), sha256=_hash(executable))
        record['source'] = json.loads((work/'source/salmonpy-runtime-source.json').read_text())
        return record
    except BaseException:
        record['status'] = 'failed'
        raise
    finally:
        record_path.write_text(json.dumps(record, indent=2)+'\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('status', 'install'):
        p = sub.add_parser(name)
        p.add_argument('--catalog', help='trusted local catalog JSON or HTTPS catalog URL')
        p.add_argument('--require', nargs='*', default=[], help='required feature names, e.g. tddft libxc')
        if name == 'install':
            p.add_argument('directory'); p.add_argument('--artifact')
            p.add_argument('--archive', help='local ZIP downloaded from CI; verified against the catalog')
    p = sub.add_parser('build', help='explicit source build; requires compilers, libraries and CMake; uses bundled sources')
    p.add_argument('upstream', help='new output directory, or legacy upstream argument'); p.add_argument('directory', nargs='?')
    p.add_argument('--jobs', required=True, type=int)
    p.add_argument('--cmake-arg', action='append', default=[], help='repeat --cmake-arg=-DNAME=VALUE')
    p = sub.add_parser('verify')
    p.add_argument('directory'); p.add_argument('--require', nargs='*', default=[])
    args = parser.parse_args(argv)
    try:
        if args.command == 'status':
            result = inspect_runtime(args.catalog, required=args.require)
        elif args.command == 'install':
            exe = install_runtime(args.directory, catalog=args.catalog, required=args.require, artifact=args.artifact, archive=args.archive)
            result = {'executable': str(exe.path), 'sha256': exe.sha256}
        elif args.command == 'verify':
            exe = load_runtime(args.directory, required=args.require)
            result = {'executable': str(exe.path), 'sha256': exe.sha256}
        else:
            result = build_runtime(args.upstream, args.directory, jobs=args.jobs, cmake_args=args.cmake_arg)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, zipfile.BadZipFile) as exc:
        parser.exit(2, str(exc)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
