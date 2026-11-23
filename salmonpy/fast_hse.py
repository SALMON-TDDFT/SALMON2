"""Process-isolated serial fast-PR HSE material response.

Each live material point owns one explicitly built native worker. There is no
compilation, download, MPI import or MPI initialization here. Native checkpoints
require the same worker, initial bundle and configuration; callers authenticate
checkpoint files before loading them. Native text output remains in work_root.
"""
from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
import hashlib
import json
import os
import queue
import re
import subprocess
import tempfile
import threading

import numpy as np

from .errors import InputError, InputTypeError, RuntimeFailure
from .sbe import _integer
from .tddft import _copy_bundle, _digest
from .embedded_tddft import _manifest
from .bands import _input_groups, _scalar
from .material_response import _response_timing

__all__ = ['FastHSEMaterialProvider', 'FastHSEResponseBackend']


def _vector(value, name):
    try:
        raw = np.asarray(value)
        if np.iscomplexobj(raw):
            raise InputError(name + ' must be real')
        result = np.array(raw, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise InputError(name + ' must be a finite real three-vector') from exc
    if result.shape != (3,) or not np.isfinite(result).all():
        raise InputError(name + ' must be a finite real three-vector')
    return result


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise InputError(name + ' must be positive and finite')
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(name + ' must be positive and finite') from exc
    if not np.isfinite(result) or result <= 0:
        raise InputError(name + ' must be positive and finite')
    return result


def _relative(name):
    value = Path(name)
    if (not name or '\x00' in name or '\n' in name or '\r' in name
            or value.is_absolute() or PureWindowsPath(name).is_absolute()
            or '..' in value.parts or '\\' in name):
        raise InputError('bundle input references must be relative paths')
    return value


def _input_text(text):
    if '!SALMONPY_END_INPUT' in text or any(len(line) > 256 for line in text.splitlines()):
        raise InputError('unsupported input marker or line length')
    groups = _input_groups(text)
    if _scalar(groups.get('units',''),'unit_system',default='au') != 'au':
        raise InputError('fast HSE worker input must use atomic units (unit_system=au)')
    em = groups.get('emfield','')
    if _scalar(em,'trans_longi',default='tr') != 'tr':
        raise InputError('fast HSE material requires transverse A drive')
    shape = _scalar(em,'ae_shape1',default='none')
    second = _scalar(em,'ae_shape2',default='none')
    if shape not in ('none','impulse') or second != 'none':
        raise InputError('fast HSE input must disable internal pulses; A is supplied externally')
    if shape == 'impulse':
        try:
            impulse = float(_scalar(em,'e_impulse',default='.01').replace('d','e'))
        except ValueError as exc:
            raise InputError('explicit zero impulse required for external A drive') from exc
        if impulse != 0:
            raise InputError('explicit zero impulse required for external A drive')
    if _scalar(groups.get('functional', ''), 'xc') != 'hse06':
        raise InputError('fast HSE material requires xc=hse06')
    if _scalar(groups.get('calculation', ''), 'theory') != 'tddft_response':
        raise InputError('fast HSE material requires a tddft_response input')
    propagation = _scalar(groups.get('propagation', ''), 'propagator')
    if propagation not in ('hse_taylor4', 'hse_taylor4_full', 'hse_ptcn'):
        raise InputError('an explicit fast-PR HSE propagator is required')
    # Input references form part of the bundle identity; external relative/absolute
    # assets cannot silently escape it. Quoted Fortran paths are supported.
    token = re.compile(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*")
    clean = token.sub(lambda m: '' if m[0].startswith('!') else m[0], text)
    refs = re.compile(r'\b(file_pseudo|file_atom_coor|file_atom_red_coor|file_kw|directory_read_data|base_directory)'
                      r'\s*(?:\([^)]*\))?\s*=\s*(\'(?:\'\'|[^\'])*\'|"(?:""|[^"])*")', re.I)
    for match in refs.finditer(clean):
        _relative(match[2][1:-1].replace("''", "'").replace('""', '"'))
    return groups, propagation


@dataclass(frozen=True)
class _Initial:
    directory: Path
    input_file: str
    manifest: tuple
    max_steps: int
    propagator: str


class FastHSEMaterialProvider:
    """Identify a fixed GS/RT bundle without converting or downloading assets.

    Keep the source bundle unchanged until states are created. Every state gets
    its own verified copy. max_steps is the worker's fixed propagation capacity.
    """
    def __init__(self, directory, input_file='rt.inp', *, max_steps=1000):
        root = Path(directory).expanduser().resolve()
        inp = _relative(os.fspath(input_file))
        if not root.is_dir():
            raise InputError('existing initial bundle directory required')
        manifest = _manifest(root)
        for name in ('info.bin', 'occupation.bin', 'wfn.bin'):
            if not (root / 'data_for_restart' / name).is_file():
                raise InputError('missing GS record: ' + name)
        if not (root / inp).is_file():
            raise InputError('missing fast HSE input file')
        _, propagator = _input_text((root / inp).read_text())
        steps = _integer(max_steps, 'max_steps', 1, 2**31 - 2)
        self.initial = _Initial(root, str(inp), manifest, steps, propagator)
        self.identity = 'salmonpy.fast-hse-initial/1/' + hashlib.sha256(json.dumps(
            [str(inp), manifest, steps], sort_keys=True).encode()).hexdigest()

    def prepare(self, context):
        if _manifest(self.initial.directory) != self.initial.manifest:
            raise InputError('initial fast HSE bundle changed after preparation')
        return self.initial


class FastHSEResponseBackend:
    """One serial CPU native worker per material point, one rank per solver group.

    build_record must come from the explicit fast-HSE worker builder. timeout is
    the maximum seconds for initialization or any protocol response. threads
    configures each child only; the host process's thread settings are unchanged.
    """
    execution = 'leader_only'
    reduction = 'leader'

    def __init__(self, build_record, work_root, *, threads=1, timeout=300., execution_mode="checked"):
        from .ld import _execution_mode_name
        self._execution_mode=_execution_mode_name(execution_mode)
        if os.name != 'posix':
            raise InputError('fast HSE worker requires POSIX /dev/stdin (macOS or Linux)')
        record_path = Path(build_record).expanduser().resolve()
        record = json.loads(record_path.read_text())
        if record.get('schema') != 'salmonpy.fast-hse-worker/1' or record.get('status') != 'built-and-audited':
            raise InputError('audited fast HSE worker build record required')
        self.executable = Path(record['executable']).expanduser().resolve()
        self.sha256 = record['sha256']
        self.origin_sha256 = record.get('origin_sha256')
        if not isinstance(self.origin_sha256, str) or re.fullmatch('[0-9a-f]{64}', self.origin_sha256) is None:
            raise InputError('fast HSE worker requires its audited origin hash')
        if not self.executable.is_file() or _digest(self.executable) != self.sha256:
            raise InputError('fast HSE worker checksum mismatch')
        self.work_root = Path(work_root).expanduser().resolve()
        if not self.work_root.is_dir():
            raise InputError('existing work root required')
        self.threads = _integer(threads, 'threads', 1, 2**31 - 1)
        self.timeout = _positive(timeout, 'timeout')
        self.identity = 'salmonpy.fast-hse-worker/1/' + hashlib.sha256(json.dumps(
            [self.sha256, self.origin_sha256, self.threads, self.execution_mode], sort_keys=True).encode()).hexdigest()

    @property
    def execution_mode(self):return self._execution_mode

    def create_state(self, initial, point_id, dt, context):
        if context is not None and context.solver_comm.size != 1:
            raise InputError('fast HSE requires one rank per material solver group; no internal MPI')
        if not isinstance(initial, _Initial):
            raise InputTypeError('prepared FastHSEMaterialProvider state required')
        dt = _positive(dt, 'dt')
        if _manifest(initial.directory) != initial.manifest:
            raise InputError('initial fast HSE bundle changed')
        if _digest(self.executable) != self.sha256:
            raise InputError('fast HSE worker changed after verification')
        if self.work_root.is_relative_to(initial.directory):
            raise InputError('work root must be outside the initial bundle')
        tag = hashlib.sha256(str(point_id).encode()).hexdigest()[:12]
        work = Path(tempfile.mkdtemp(prefix='fast-hse-' + tag + '-', dir=self.work_root))
        bundle = work / 'bundle'
        if _copy_bundle(initial.directory, bundle) != initial.manifest:
            raise InputError('initial bundle changed during copy')
        groups, propagator = _input_text((bundle / initial.input_file).read_text())
        # Appending overrides only the two explicitly selected time-grid inputs;
        # Fortran namelists use the final assignment. Preserve all other settings.
        groups['tgrid'] = groups.get('tgrid', '') + '\n dt=' + format(dt, '.17g') + ', nt=' + str(initial.max_steps) + '\n'
        (bundle / 'inputfile').write_text('\n'.join('&' + k + v + '/' for k, v in groups.items()) + '\n')
        return _State(self, bundle, dt, initial.max_steps, propagator)


class _State:
    checkpoint_protocol = 'salmonpy.direct-state/1'
    checkpoint_format = 'salmonpy.fast-hse-worker-native/2'
    checkpoint_modes = ('save','restore')
    checkpoint_decomposition = (1,1,1,1,1)

    def __init__(self, backend, work, dt, capacity, propagator):
        self._execution_mode=backend.execution_mode
        self.work = work
        self.dt = dt
        self.capacity = capacity
        self.steps = 0
        self.A = np.zeros(3)
        self.current = np.zeros(3)
        self.closed = False
        self._failed = False
        self._timeout = backend.timeout
        self._lock = threading.Lock()
        self._messages = queue.Queue()
        self._process = None
        self._reader = None
        self._log = (work.parent / 'worker.log').open('w', encoding='utf-8')
        self.timing = _response_timing(dt, 'fast-PR HSE ' + propagator + '; native charge current j')
        try:
            self._process = subprocess.Popen([str(backend.executable)], cwd=work,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', errors='replace', bufsize=1,
                env={**os.environ, 'OMP_NUM_THREADS': str(backend.threads), 'OPENBLAS_NUM_THREADS': '1'},
                start_new_session=True)
            self._reader = threading.Thread(target=self._read_stdout, name='salmonpy-fast-hse-reader', daemon=True)
            self._reader.start()
            parts = self._response('READY')
            if len(parts) != 4 or self._step_number(parts[0]) != 0:
                raise RuntimeFailure('invalid fast HSE READY response')
            self.current = self._numbers(parts[1:4])
            self.initial_current = self.current.copy()
            self.energy = None
            if self.execution_mode=='performance':
                version=self._request('MODE_PERFORMANCE','PERFORMANCE')
                if version!=['1']:raise RuntimeFailure('unsupported fast HSE performance protocol')
        except BaseException:
            self._failed = True
            self._terminate()
            raise

    @property
    def execution_mode(self):return self._execution_mode

    def _read_stdout(self):
        try:
            for line in self._process.stdout:
                self._log.write(line)
                self._log.flush()
                stripped = line.strip()
                if stripped.startswith('SALMONPY_HSE_'):
                    self._messages.put(stripped.split())
        except BaseException as exc:
            self._messages.put(exc)
        finally:
            self._messages.put(None)

    def _response(self, expected):
        try:
            parts = self._messages.get(timeout=self._timeout)
        except queue.Empty as exc:
            raise RuntimeFailure('fast HSE worker timed out; see ' + str(self._log.name)) from exc
        if parts is None:
            raise RuntimeFailure('fast HSE worker ended without ' + expected + '; see ' + str(self._log.name))
        if isinstance(parts, BaseException):
            raise RuntimeFailure('fast HSE worker output failed') from parts
        if not parts or parts[0] != 'SALMONPY_HSE_' + expected:
            raise RuntimeFailure('unexpected fast HSE response: ' + ' '.join(parts))
        return parts[1:]

    @staticmethod
    def _step_number(text):
        try:
            result = int(text)
        except ValueError as exc:
            raise RuntimeFailure('invalid native step number') from exc
        if result < 0:
            raise RuntimeFailure('negative native step number')
        return result

    @staticmethod
    def _numbers(parts, *, checked=True):
        try:
            result = np.array([float(p.replace('D', 'E').replace('d', 'e')) for p in parts])
        except ValueError as exc:
            raise RuntimeFailure('invalid native numeric response') from exc
        if checked and not np.isfinite(result).all():
            raise RuntimeFailure('nonfinite native response')
        return result

    def _request(self, command, expected):
        if self.closed or self._failed:
            raise InputError('fast HSE state is closed or failed')
        try:
            self._process.stdin.write(command + '\n')
            self._process.stdin.flush()
            return self._response(expected)
        except BaseException:
            self._failed = True
            self._terminate()
            raise

    def _enter(self):
        if not self._lock.acquire(blocking=False):
            raise InputError('concurrent calls on one fast HSE state are unsupported')

    def advance(self, old, new):
        self._enter()
        try:
            if self.execution_mode=='checked':
                old=_vector(old,'old A');new=_vector(new,'new A')
            else:
                from .ld import _array_unchecked
                old=_array_unchecked(old,'old A',(3,));new=_array_unchecked(new,'new A',(3,))
            if self.closed or self._failed or self.steps >= self.capacity:
                raise InputError('fast HSE state is closed, failed or exhausted')
            if self.execution_mode=="checked" and not np.array_equal(old, self.A):
                raise InputError('old A differs from accepted fast HSE endpoint')
            parts = self._request('STEP ' + ' '.join(format(x, '.17g') for x in new), 'STEP')
            try:
                if len(parts) != 5 or self._step_number(parts[0]) != self.steps + 1:
                    raise RuntimeFailure('fast HSE worker step mismatch')
                values = self._numbers(parts[1:],checked=self.execution_mode=="checked")
            except BaseException:
                self._failed = True; self._terminate(); raise
            self.steps += 1
            self.A = new
            self.current = values[:3].copy()
            self.energy = float(values[3])
            return self.current.copy()
        finally:
            self._lock.release()

    @property
    def time(self):
        return self.steps * self.dt

    def step(self, A):
        return self.advance(self.A.copy(), A)

    def checkpoint(self, path, restore):
        self._enter()
        try:
            if not isinstance(restore, (bool, np.bool_)):
                raise InputTypeError('restore must be boolean')
            destination = Path(path).expanduser().resolve()
            spelling = str(destination)
            if any(c in spelling for c in ('\x00', '\n', '\r')) or len(os.fsencode(spelling)) > 4000:
                raise InputError('invalid native checkpoint path')
            if restore and self.steps != 0:
                raise InputError('restore requires a fresh fast HSE state')
            if restore and not destination.is_file():
                raise InputError('checkpoint file is missing')
            if not restore and (destination.exists() or not destination.parent.is_dir()):
                raise InputError('new checkpoint file in an existing directory required')
            parts = self._request(('LOAD ' if restore else 'SAVE ') + spelling, 'LOADED' if restore else 'SAVED')
            try:
                expected = 4 if restore else 1
                if len(parts) != expected:
                    raise RuntimeFailure('invalid fast HSE checkpoint response')
                step = self._step_number(parts[0])
                if step > self.capacity or (not restore and step != self.steps):
                    raise RuntimeFailure('fast HSE checkpoint step mismatch')
                if restore:
                    current = self._numbers(parts[1:])
            except BaseException:
                self._failed = True; self._terminate(); raise
            if restore:
                self.steps = step
                self.current = current
                self.energy = None
                # A is supplied by the caller's authenticated checkpoint metadata.
        finally:
            self._lock.release()

    def _terminate(self):
        process = self._process
        if process is not None:
            if process.poll() is None:
                try:
                    process.terminate()
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=2.)
                except subprocess.TimeoutExpired:
                    try:
                        process.kill()
                    except ProcessLookupError:
                        pass
                    process.wait(timeout=2.)
            if process.stdin is not None:
                try:
                    process.stdin.close()
                except OSError:
                    pass
            if self._reader is not None:
                self._reader.join(timeout=2.)
            if process.stdout is not None:
                process.stdout.close()
        if not self._log.closed:
            self._log.close()
        self.closed = True

    def close(self):
        self._enter()
        try:
            if self.closed:
                return
            try:
                if not self._failed:
                    parts = self._request('CLOSE', 'CLOSED')
                    if len(parts) != 1 or self._step_number(parts[0]) != self.steps:
                        raise RuntimeFailure('invalid fast HSE close response')
                    code = self._process.wait(timeout=min(self._timeout, 5.))
                    if code:
                        raise RuntimeFailure('fast HSE worker failed during close')
            finally:
                self._terminate()
        finally:
            self._lock.release()
