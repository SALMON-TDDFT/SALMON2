"""Fixed Si starter: verify a prebuilt runtime, create GS, and compute charge j.

python -m salmonpy.examples.si_intro prepare NEW_DIR --runtime-artifact ACTIONS_ZIP
python -m salmonpy.examples.si_intro run NEW_DIR --threads 1 --memory-budget-gib 4

No download or compilation occurs. The fixed small case verifies reproducibility,
not spatial/k-point/time convergence or production optical accuracy.
"""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import stat
import tempfile
import time
import zipfile

import numpy as np
from salmonpy.crystal import respond
from salmonpy.dft import DFT, Pseudopotential, SCFSettings
from salmonpy.fields import AcosPulse
from salmonpy.runtime import install_runtime, load_runtime

_ASSET_NAMES = ('input.json', 'expected.json', 'provenance.json', 'Si_rps.dat', 'LICENSE.txt')
_REQUIRED = ('dft', 'tddft', 'pz', 'single-rank')
_MAX_ARTIFACT = 512*1024*1024


def _assets():
    directory = Path(__file__).with_name('si_intro_data')
    assets = {name: (directory/name).read_bytes() for name in _ASSET_NAMES}
    provenance = json.loads(assets['provenance.json'])
    if provenance.get('schema') != 'salmonpy.si-intro-provenance/1':
        raise ValueError('unsupported Si starter provenance schema')
    expected_names = set(_ASSET_NAMES)-{'provenance.json'}
    if set(provenance.get('files', {})) != expected_names:
        raise ValueError('incomplete starter provenance file inventory')
    for name in expected_names:
        if hashlib.sha256(assets[name]).hexdigest() != provenance['files'][name]:
            raise ValueError('packaged starter checksum mismatch: '+name)
    return assets


def _digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def _model(configuration, directory):
    if configuration.get('schema') != 'salmonpy.si-intro/1':
        raise ValueError('unsupported Si starter input schema')
    pseudo = configuration['pseudo']
    if pseudo['filename'] != 'Si_rps.dat':
        raise ValueError('the fixed Si starter requires Si_rps.dat')
    model = DFT(**configuration['dft'], pseudopotentials=[Pseudopotential.from_file(
        Path(directory)/pseudo['filename'], atomic_number=pseudo['atomic_number'],
        local_channel=pseudo['local_channel'])])
    return model, SCFSettings(**configuration['scf'])


def _clock():
    try:
        value = time.perf_counter()
        return float(value) if np.isfinite(value) else None
    except Exception:
        return None


def _elapsed(start):
    end = _clock()
    if start is None or end is None or end < start:
        return {'status': 'unavailable', 'seconds': None}
    return {'status': 'measured', 'seconds': end-start}


def compare_result(expected, energy, times, current):
    """Compare all endpoint currents, their signs/times, and the total GS energy."""
    tolerances = {}
    for key in ('energy_atol', 'current_atol', 'current_rtol'):
        value = expected[key]
        if isinstance(value, bool) or not np.isfinite(value) or value < 0:
            raise ValueError('reference tolerances must be finite and nonnegative')
        tolerances[key] = float(value)
    target_time = np.asarray(expected['time'], dtype=float)
    target_current = np.asarray(expected['current'], dtype=float)
    target_energy = float(expected['energy_hartree'])
    if (target_time.ndim != 1 or len(target_time) < 2 or np.any(np.diff(target_time) <= 0)
            or target_current.shape != (len(target_time), 3)
            or not np.isfinite(target_time).all() or not np.isfinite(target_current).all()
            or not np.isfinite(target_energy)):
        raise ValueError('invalid expected energy/time/current arrays')
    times = np.asarray(times, dtype=float); current = np.asarray(current, dtype=float)
    matching_shape = times.shape == target_time.shape and current.shape == target_current.shape
    finite = bool(np.isfinite(energy) and np.isfinite(times).all() and np.isfinite(current).all())
    # Time acceptance is roundoff-only and cannot hide a shifted step.
    dt = float(np.min(np.diff(target_time)))
    time_ok = matching_shape and finite and bool(np.allclose(times, target_time, rtol=0, atol=dt*1e-10))
    energy_error = abs(float(energy)-target_energy) if np.isfinite(energy) else None
    current_error = float(np.max(np.abs(current-target_current))) if matching_shape and finite else None
    energy_ok = energy_error is not None and energy_error <= tolerances['energy_atol']
    current_ok = matching_shape and finite and bool(np.allclose(
        current, target_current, atol=tolerances['current_atol'], rtol=tolerances['current_rtol']))
    return {'passed': bool(time_ok and energy_ok and current_ok),
            'energy_passed': bool(energy_ok), 'time_passed': bool(time_ok),
            'current_passed': bool(current_ok), 'energy_absolute_error_hartree': energy_error,
            'current_max_absolute_error_au': current_error, 'tolerances': tolerances,
            'current_convention': 'physical charge current density j = -SALMON electronic matter current J',
            'scope': 'fixed-input numerical reproducibility; not physical convergence'}


def prepare(directory, runtime_artifact):
    """Verify a known Actions ZIP and prepare a new, self-contained input folder."""
    destination = Path(directory).expanduser().resolve()
    if destination.exists():
        raise FileExistsError(destination)
    archive = Path(runtime_artifact).expanduser().resolve()
    assets = _assets()
    provenance = json.loads(assets['provenance.json'])
    sha = _digest(archive)
    matches = [entry for entry in provenance['runtime_artifacts'] if entry['outer_sha256'] == sha]
    if len(matches) != 1:
        raise ValueError('runtime artifact checksum is not a pinned Si starter artifact; see SI_INTRO guide')
    entry = matches[0]
    inner = entry['archive_member']
    if Path(inner).name != inner or not inner.endswith('.zip') or '\\' in inner:
        raise ValueError('invalid pinned runtime archive member')
    with zipfile.ZipFile(archive) as zipped:
        members = zipped.infolist()
        if len(members) != 2 or {m.filename for m in members} != {'catalog.json', inner}:
            raise ValueError('unexpected outer artifact members')
        if sum(m.file_size for m in members) > _MAX_ARTIFACT:
            raise ValueError('expanded runtime artifact exceeds starter size limit')
        if any(m.is_dir() or stat.S_ISLNK(m.external_attr >> 16) for m in members):
            raise ValueError('invalid outer artifact members')
        # Never extract supplied paths. Only the two known names are copied.
        with tempfile.TemporaryDirectory(prefix='salmonpy-si-artifact-') as temporary:
            temporary = Path(temporary)
            for name in ('catalog.json', inner):
                with zipped.open(name) as source, (temporary/name).open('wb') as target:
                    shutil.copyfileobj(source, target)
            destination.mkdir(parents=True, exist_ok=False)
            record = {'schema': 'salmonpy.si-intro-preparation/1', 'status': 'started',
                      'runtime_artifact_sha256': sha, 'runtime_artifact': entry,
                      'asset_sha256': {name: hashlib.sha256(data).hexdigest() for name, data in assets.items()}}
            try:
                engine = install_runtime(destination/'runtime', catalog=temporary/'catalog.json',
                                         archive=temporary/inner, required=_REQUIRED)
                for name, data in assets.items():
                    (destination/name).write_bytes(data)
                model, scf = _model(json.loads(assets['input.json']), destination)
                (destination/'gs.inp').write_text(model.input_text(scf), encoding='utf-8')
                record.update(status='complete', executable_sha256=engine.sha256)
            except BaseException as exc:
                record.update(status='failed', error=str(exc))
                raise
            finally:
                _write_json(destination/'preparation.json', record)
    return destination


@contextmanager
def _thread_environment(threads):
    names = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
    previous = {name: os.environ.get(name) for name in names}
    try:
        for name in names:
            os.environ[name] = str(threads if name == 'OMP_NUM_THREADS' else 1)
        yield
    finally:
        for name, value in previous.items():
            if value is None: os.environ.pop(name, None)
            else: os.environ[name] = value


def run(directory, *, threads=1, memory_budget_gib):
    """Run exactly the packaged starter inputs; retain failures and comparisons."""
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 4:
        raise ValueError('threads must be an integer from 1 to 4 within the user CPU budget')
    if isinstance(memory_budget_gib, bool) or not np.isfinite(memory_budget_gib) or memory_budget_gib <= 0:
        raise ValueError('positive finite memory_budget_gib required')
    root = Path(directory).expanduser().resolve()
    assets = _assets()
    for name, data in assets.items():
        path = root/name
        if path.is_symlink() or not path.is_file() or path.read_bytes() != data:
            raise ValueError('prepared starter asset differs from the packaged reference: '+name)
    configuration = json.loads(assets['input.json'])
    expected = json.loads(assets['expected.json'])
    model, scf = _model(configuration, root)
    if (root/'gs.inp').read_text(encoding='utf-8') != model.input_text(scf):
        raise ValueError('prepared gs.inp differs from the fixed starter input')
    # Validate expected arrays and tolerance schema before any calculation.
    compare_result(expected, expected['energy_hartree'], expected['time'], expected['current'])
    preparation = json.loads((root/'preparation.json').read_text(encoding='utf-8'))
    if preparation.get('schema') != 'salmonpy.si-intro-preparation/1' or preparation.get('status') != 'complete':
        raise ValueError('starter preparation is incomplete')
    for name, data in assets.items():
        if preparation.get('asset_sha256', {}).get(name) != hashlib.sha256(data).hexdigest():
            raise ValueError('prepared asset inventory differs: '+name)
    pinned = json.loads(assets['provenance.json'])['runtime_artifacts']
    if preparation.get('runtime_artifact') not in pinned:
        raise ValueError('unrecognized prepared runtime artifact')
    engine = load_runtime(root/'runtime', required=_REQUIRED)
    if engine.sha256 != preparation.get('executable_sha256'):
        raise ValueError('runtime executable changed since preparation')
    output = root/'result'
    output.mkdir(exist_ok=False)
    report = {'schema': 'salmonpy.si-intro-result/1', 'status': 'started',
              'platform': platform.platform(), 'python': platform.python_version(),
              'executable_sha256': engine.sha256,
              'setup': {'status': 'passed', 'runtime_artifact': preparation['runtime_artifact']},
              'calculation': {'status': 'not_started'}, 'validation': {'status': 'not_performed'},
              'resources': {'memory_budget_gib': float(memory_budget_gib), 'memory_enforced': False,
                            'threads': threads, 'simultaneous_workers': 1,
                            'peak_memory': {'status': 'not_measured', 'bytes': None,
                                            'scope': 'whole job including child processes'}},
              'scope': 'small fixed Si starter; not grid/k-point/time-converged optics'}
    started = _clock()
    try:
        report['calculation']['status'] = 'running'
        with _thread_environment(threads):
            gs = model.run(engine, output/'gs', scf=scf, threads=threads, timeout=600.)
            rt = configuration['rt']
            drive = AcosPulse(**rt['pulse']).sample(np.arange(rt['steps']+1)*rt['dt'])
            with gs.tddft_initial_state(executable=engine, max_steps=rt['steps'], threads=threads, timeout=600.) as initial:
                response = respond(initial, drive, output_directory=output/'rt')
        np.savez(output/'current.npz', time=response.time, A=drive.values, current=response.current)
        np.savetxt(output/'current.csv', np.column_stack((response.time, drive.values, response.current)),
                   delimiter=',', header='time_au,Ac_x_au,Ac_y_au,Ac_z_au,j_x_au,j_y_au,j_z_au', comments='')
        report['calculation'] = {'status': 'passed', 'energy_hartree': gs.energy_hartree,
                                 'scf_iterations': gs.iterations, 'time_samples': len(response.time),
                                 'current_convention': response.field_convention}
        comparison = compare_result(expected, gs.energy_hartree, response.time, response.current)
        report['validation'] = {'status': 'passed' if comparison['passed'] else 'failed', **comparison}
        report['status'] = 'passed' if comparison['passed'] else 'validation_failed'
        return report
    except BaseException as exc:
        report['status'] = 'calculation_failed'
        report['calculation'].update(status='failed', error=str(exc))
        raise
    finally:
        report['resources']['wall_time'] = _elapsed(started)
        _write_json(output/'report.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prepare_parser = sub.add_parser('prepare', help='verify and unpack one pinned prebuilt Actions artifact')
    prepare_parser.add_argument('directory')
    prepare_parser.add_argument('--runtime-artifact', required=True)
    run_parser = sub.add_parser('run', help='run fixed GS and A-driven current; no build or download')
    run_parser.add_argument('directory')
    run_parser.add_argument('--threads', type=int, default=1, choices=range(1, 5))
    run_parser.add_argument('--memory-budget-gib', type=float, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        print(prepare(args.directory, args.runtime_artifact))
        return 0
    result = run(args.directory, threads=args.threads, memory_budget_gib=args.memory_budget_gib)
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
