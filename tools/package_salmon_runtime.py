"""Create a relocatable single-rank runtime bundle (macOS/Linux).

Run native integration against the relocated result BEFORE publishing the ZIP.
The base URL must refer to the intended SALMONpy release, never inferred upstream.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from salmonpy.runtime import INTERFACE, host_platform, load_runtime


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def linux_dependencies(binary):
    output = subprocess.check_output(['ldd', str(binary)], text=True)
    if 'not found' in output:
        raise RuntimeError(output)
    # glibc and its loader are supplied by the target OS, all other DSOs bundled.
    system = re.compile(r'^(ld-linux[^/]*|lib(c|m|pthread|dl|rt|resolv|util)\.so\..*)$')
    result = []
    for line in output.splitlines():
        match = re.search(r'=> (/\S+)', line) or re.match(r'\s*(/\S+)', line)
        if match:
            path = Path(match.group(1))
            if not system.fullmatch(path.name): result.append(path)
    return result


def package(build, source, destination, *, base_url=None, license_dirs=()):
    root = Path(__file__).resolve().parents[1]
    build, source, dest = (Path(p).resolve() for p in (build, source, destination))
    if base_url is not None and not base_url.startswith('https://'):
        raise ValueError('base URL must be the HTTPS release asset directory')
    if sys.platform not in ('darwin', 'linux'):
        raise RuntimeError('portable full-runtime packaging is currently implemented for macOS/Linux only')
    cache = (build/'CMakeCache.txt').read_text()
    def enabled(name): return bool(re.search(r'^'+name+r':\w+=(ON|TRUE|1|YES)$', cache, re.M | re.I))
    if enabled('USE_MPI'):
        raise ValueError('MPI runtimes must be built for the cluster; this bundle profile is single rank')
    provenance = json.loads((source/'salmonpy-runtime-source.json').read_text())
    expected = json.loads((root/'runtime/source.json').read_text())
    if provenance != expected:
        raise ValueError('prepared source does not match this SALMONpy runtime revision')
    record_path = build.parent/'build.json'
    if not record_path.is_file():
        raise ValueError('package only fresh outputs of python -m salmonpy.runtime build (build.json required)')
    record = json.loads(record_path.read_text())
    if record.get('build_configuration',{}).get('mpi_enabled') is not False:
        raise ValueError('portable bundle requires a completed serial capability record; rebuild explicitly')
    if (record.get('status') != 'built-and-audited' or record.get('source') != expected
            or Path(record.get('executable', '')).resolve() != build/'salmon'
            or record.get('sha256') != digest(build/'salmon')
            or source != build.parent/'source'):
        raise ValueError('build provenance or executable checksum does not match')
    platform = host_platform()
    if sys.platform == 'darwin':
        arches = subprocess.check_output(['lipo', '-archs', str(build/'salmon')], text=True).split()
        if arches != [platform['machine']]:
            raise ValueError('runtime binary architecture differs from the builder host')
    else:
        header = (build/'salmon').read_bytes()[:20]
        if header[:4] != b'\x7fELF' or header[5] != 1:
            raise ValueError('expected a little-endian ELF runtime')
        machine = int.from_bytes(header[18:20], 'little')
        if {'x86_64': 62, 'arm64': 183}.get(platform['machine']) != machine:
            raise ValueError('runtime binary architecture differs from the builder host')
    subprocess.run([sys.executable, str(root/'tools/check_forbidden_eigensolvers.py'),
                    '--object-tree', str(build/'src')], check=True)
    dest.mkdir(parents=True, exist_ok=False)
    bundle = dest/'runtime'
    (bundle/'bin').mkdir(parents=True)
    executable = bundle/'bin/salmon'
    shutil.copy2(build/'salmon', executable)
    (bundle/'lib').mkdir()
    if sys.platform == 'darwin':
        from delocate.delocating import delocate_path
        delocate_path(str(bundle), str(bundle/'lib'), sanitize_rpaths=True)
    else:
        copies = {}
        for path in linux_dependencies(executable):
            if path.name in copies and path.resolve() != copies[path.name].resolve():
                raise ValueError('dependency basename collision: '+path.name)
            copies[path.name] = path
            shutil.copy2(path, bundle/'lib'/path.name)
        subprocess.run(['patchelf', '--set-rpath', '$ORIGIN/../lib', str(executable)], check=True)
        for path in (bundle/'lib').iterdir():
            subprocess.run(['patchelf', '--set-rpath', '$ORIGIN', str(path)], check=True)
        if any(not p.resolve().is_relative_to(bundle) for p in linux_dependencies(executable)):
            raise RuntimeError('unbundled runtime dependencies remain')
    license_root = bundle/'licenses'
    license_root.mkdir()
    for name in ('LICENSE', 'LICENSE.THIRD-PARTY', 'ACKNOWLEDGEMENTS'):
        shutil.copy2(source/name, license_root/name)
    shutil.copytree(root/'packaging/ld/licenses', license_root/'compiler-blas')
    for i, path in enumerate(license_dirs):
        shutil.copytree(path, license_root/f'additional-{i}')
    if enabled('USE_LIBXC') and not license_dirs:
        raise ValueError('Libxc build requires --license-dir with the Libxc license')
    features = ['dft', 'tddft', 'fdtd', 'bands', 'pz', 'tbmbj', 'single-rank', 'material-id-independent']
    if enabled('USE_LIBXC'): features += ['libxc', 'libxc-pz', 'libxc-pbe']
    name = 'salmonpy-runtime-'+platform['system']+'-'+platform['machine']
    meta = {'schema': 'salmonpy.runtime/1', 'interface': INTERFACE,
            'platform': platform, 'features': features, 'executable': 'bin/salmon',
            'provenance': provenance, 'validation': 'object-audited; native relocation tests required before publication',
            'files': {p.relative_to(bundle).as_posix(): digest(p) for p in sorted(bundle.rglob('*')) if p.is_file()}}
    (bundle/'runtime.json').write_text(json.dumps(meta, indent=2)+'\n')
    load_runtime(bundle)
    archive = dest/(name+'.zip')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for path in sorted(bundle.rglob('*')):
            if path.is_file(): z.write(path, path.relative_to(bundle).as_posix())
    entry = {'id': name, 'interface': INTERFACE, 'platform': platform, 'features': features,
             'url': (base_url.rstrip('/')+'/'+archive.name) if base_url else None, 'sha256': digest(archive)}
    (dest/'catalog.json').write_text(json.dumps({'schema': 'salmonpy.runtimes/1', 'artifacts': [entry]}, indent=2)+'\n')
    print(json.dumps({'bundle': str(bundle), 'archive': str(archive), 'catalog': str(dest/'catalog.json')}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build'); p.add_argument('source'); p.add_argument('destination')
    p.add_argument('--base-url', help='omit for branch CI artifacts without a public release URL')
    p.add_argument('--license-dir', action='append', default=[])
    a = p.parse_args()
    package(a.build, a.source, a.destination, base_url=a.base_url, license_dirs=a.license_dir)
