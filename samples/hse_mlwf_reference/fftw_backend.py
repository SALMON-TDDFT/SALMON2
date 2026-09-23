"""Cached single-thread FFTW complex 3-D periodic convolution.

Loads only libfftw3, never the threaded library. Each instance owns fixed,
FFTW-aligned storage and two in-place plans. Instances are not reentrant:
use a distinct instance for each worker. Returned arrays are independent
copies. Planning time is exposed separately from convolution execution.
"""
import ctypes as ct
import os
import time
import numpy as np

DEFAULT_FFTW = '/opt/homebrew/opt/fftw/lib/libfftw3.dylib'

class FFTWConvolution:
    def __init__(self, multiplier, library=None):
        multiplier = np.asarray(multiplier)
        if multiplier.ndim != 3 or min(multiplier.shape) < 1 or not np.isfinite(multiplier).all():
            raise ValueError('A finite nonempty three-dimensional multiplier is required')
        self.shape = tuple(multiplier.shape)
        self.multiplier = np.array(multiplier, dtype=np.complex128, order='C', copy=True)
        self._scaled_multiplier = self.multiplier / self.multiplier.size
        self._memory = self._forward = self._backward = None
        self._buffer = None
        self.lib = ct.CDLL(library or os.environ.get('FFTW_LIBRARY', DEFAULT_FFTW))
        lib = self.lib
        lib.fftw_malloc.argtypes = [ct.c_size_t]
        lib.fftw_malloc.restype = ct.c_void_p
        lib.fftw_free.argtypes = [ct.c_void_p]
        lib.fftw_free.restype = None
        lib.fftw_plan_dft_3d.argtypes = [ct.c_int]*3 + [ct.c_void_p]*2 + [ct.c_int, ct.c_uint]
        lib.fftw_plan_dft_3d.restype = ct.c_void_p
        lib.fftw_execute.argtypes = [ct.c_void_p]
        lib.fftw_execute.restype = None
        lib.fftw_destroy_plan.argtypes = [ct.c_void_p]
        lib.fftw_destroy_plan.restype = None
        try:
            self._memory = lib.fftw_malloc(self.multiplier.size * 16)
            if not self._memory:
                raise MemoryError('fftw_malloc failed')
            raw = (ct.c_double * (2*self.multiplier.size)).from_address(self._memory)
            self._buffer = np.ctypeslib.as_array(raw).view(np.complex128).reshape(self.shape)
            start = time.perf_counter()
            # FFTW_MEASURE=0; FFTW_FORWARD=-1, FFTW_BACKWARD=+1. Planning
            # may overwrite storage; every convolve fills it from fresh input.
            self._forward = lib.fftw_plan_dft_3d(*self.shape, self._memory, self._memory, -1, 0)
            self._backward = lib.fftw_plan_dft_3d(*self.shape, self._memory, self._memory, 1, 0)
            self.plan_seconds = time.perf_counter()-start
            if not self._forward or not self._backward:
                raise RuntimeError('FFTW could not construct convolution plans')
        except Exception:
            self.close()
            raise

    def convolve(self, rho):
        if self._memory is None:
            raise RuntimeError('FFTW convolution is closed')
        rho = np.asarray(rho)
        if rho.shape != self.shape:
            raise ValueError('Density shape differs from convolution grid')
        np.copyto(self._buffer, rho, casting='same_kind')
        self.lib.fftw_execute(self._forward)
        self._buffer *= self._scaled_multiplier
        self.lib.fftw_execute(self._backward)
        return self._buffer.copy()

    def close(self):
        for name in ('_forward', '_backward'):
            plan = getattr(self, name)
            if plan:
                self.lib.fftw_destroy_plan(plan)
                setattr(self, name, None)
        self._buffer = None
        if self._memory:
            self.lib.fftw_free(self._memory)
            self._memory = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
