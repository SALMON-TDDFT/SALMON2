"""Experimental validation-only E-worker driver; not a supported public API.

Explicit prepared directory, one native child, bounded reply waits. Backend
failures close the process; this driver makes no atomic rollback guarantee.
"""
import os
import queue
import subprocess
import threading
import numpy as np


class IsolatedWorker:
    def __init__(self, executable, directory, input_text, *, dt, capacity, timeout=30.):
        self.dt = dt
        self.time = 0.
        self.steps = 0
        self.capacity = capacity
        self.timeout = timeout
        self._events = queue.Queue()
        env = dict(os.environ)
        env.update(SALMONPY_TDDFT_WORKER='E3', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        self._proc = subprocess.Popen([str(executable)], cwd=directory, env=env,
                                      stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, text=True, bufsize=1)
        def read():
            try:
                with (directory / 'worker.log').open('w') as log:
                    for line in self._proc.stdout:
                        log.write(line)
                        log.flush()
                        if line.startswith('SALMONPY_TDDFT_E_V3 '):
                            self._events.put(line)
            finally:
                self._proc.stdout.close()
                self._events.put(None)
        self._reader = threading.Thread(target=read, daemon=True)
        self._reader.start()
        try:
            self._send(input_text+'\n!SALMONPY_END_INPUT\n'+f'INIT_E {dt:.17e} {capacity} 0 0 0')
            self.dipoles = self._receive(0)
        except BaseException:
            self.close()
            raise

    def _send(self, text):
        done = queue.Queue()
        def write():
            try:
                self._proc.stdin.write(text+'\n')
                self._proc.stdin.flush()
                done.put(None)
            except BaseException as error:
                done.put(error)
        writer = threading.Thread(target=write, daemon=True)
        writer.start()
        try:
            error = done.get(timeout=self.timeout)
        except queue.Empty:
            # Kill before closing stdin so a blocked writer releases its lock.
            self._proc.kill()
            self._proc.wait(timeout=5)
            writer.join(timeout=5)
            raise TimeoutError('E-worker input write timed out') from None
        if error is not None:
            raise error

    def _receive(self, step):
        line = self._events.get(timeout=self.timeout)
        if line is None:
            raise RuntimeError('E-worker stopped; inspect worker.log')
        fields = line.split()
        if len(fields) != 11 or int(fields[1]) != step:
            raise RuntimeError('invalid E-worker reply')
        values = np.asarray([float(v) for v in fields[2:]])
        if not np.isfinite(values).all():
            raise RuntimeError('nonfinite dipole')
        values.setflags(write=False)
        return values

    def step_fields(self, fields):
        if self.steps >= self.capacity:
            raise ValueError('worker capacity exhausted')
        try:
            self._send('STEP_E '+' '.join(f'{v:.17e}' for v in fields.E.value))
            values = self._receive(self.steps+1)
        except BaseException:
            self.close()
            raise
        self.steps += 1
        self.time = self.steps*self.dt
        self.dipoles = values
        return values.copy()

    def close(self):
        if self._proc.poll() is None:
            self._proc.terminate()
            try:
                self._proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._proc.kill()
                self._proc.wait(timeout=5)
        self._proc.stdin.close()
        self._reader.join(timeout=5)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
