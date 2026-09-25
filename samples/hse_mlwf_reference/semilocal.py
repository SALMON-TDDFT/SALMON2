"""Unpolarized PBE and HSE06 semilocal remainder through Libxc (atomic units).

HSE06 here is *only* the weighted semilocal remainder, PBE correlation plus
WPBEH exchange at omega=0 minus 25% short-range WPBEH exchange.
This is Libxc HSE06's exchange-hole model; it differs slightly from GGA_X_PBE
at finite gradient. Add 25% short-range Fock
exchange separately; do not multiply this remainder by 0.75 again.

Inputs are total density rho and sigma=|grad rho|**2. Returned eps is energy
per electron; vrho and vsigma differentiate energy per volume rho*eps.
The multiplicative GGA potential still requires vrho - 2 div(vsigma grad rho).
"""
import ctypes as ct
import os
import numpy as np

DEFAULT_LIBXC = '/opt/homebrew/opt/libxc/lib/libxc.dylib'

class Semilocal:
    def __init__(self, name='hse06', library=None):
        self.name = name.lower()
        ids = {'hse06': (428,), 'pbe': (101, 130)}
        if self.name not in ids:
            raise ValueError('Supported functionals: hse06, pbe')
        self.lib = ct.CDLL(library or os.environ.get('LIBXC_LIBRARY', DEFAULT_LIBXC))
        self._funcs = []
        lib = self.lib
        ptr = ct.POINTER(ct.c_double)
        lib.xc_func_alloc.argtypes = []
        lib.xc_func_alloc.restype = ct.c_void_p
        lib.xc_func_init.argtypes = [ct.c_void_p, ct.c_int, ct.c_int]
        lib.xc_func_init.restype = ct.c_int
        for name_ in ('xc_func_end', 'xc_func_free'):
            getattr(lib, name_).argtypes = [ct.c_void_p]
            getattr(lib, name_).restype = None
        lib.xc_gga_exc_vxc.argtypes = [ct.c_void_p, ct.c_size_t] + [ptr]*5
        lib.xc_gga_exc_vxc.restype = None
        lib.xc_hyb_cam_coef.argtypes = [ct.c_void_p, ptr, ptr, ptr]
        lib.xc_hyb_cam_coef.restype = None
        try:
            for functional in ids[self.name]:
                f = lib.xc_func_alloc()
                if not f:
                    raise MemoryError('xc_func_alloc failed')
                if lib.xc_func_init(f, functional, 1):
                    lib.xc_func_free(f)
                    raise RuntimeError(f'Libxc could not initialize functional {functional}')
                self._funcs.append(f)
            omega, alpha, beta = ct.c_double(), ct.c_double(), ct.c_double()
            if self.name == 'hse06':
                lib.xc_hyb_cam_coef(self._funcs[0], ct.byref(omega), ct.byref(alpha), ct.byref(beta))
            # Libxc CAM convention: alpha is LR exact fraction; alpha+beta is SR.
            self.coefficients = dict(omega=omega.value,
                                     long_range_exact=alpha.value,
                                     short_range_exact=alpha.value+beta.value)
        except Exception:
            self.close()
            raise

    def evaluate(self, rho, sigma):
        """Return (eps, vrho, vsigma), arrays of the same shape as rho."""
        if not self._funcs:
            raise RuntimeError('Semilocal functional is closed')
        rho, sigma = np.asarray(rho, dtype=np.float64), np.asarray(sigma, dtype=np.float64)
        if rho.shape != sigma.shape or not np.isfinite(rho).all() or not np.isfinite(sigma).all():
            raise ValueError('rho and sigma must have matching shapes and finite values')
        if (rho < 0).any() or (sigma < 0).any():
            raise ValueError('rho and sigma must be nonnegative')
        shape = rho.shape
        r, s = np.ascontiguousarray(rho).reshape(-1), np.ascontiguousarray(sigma).reshape(-1)
        result = np.zeros((3, r.size), dtype=np.float64)
        part = np.empty_like(result)
        ptr = ct.POINTER(ct.c_double)
        def address(a):
            return a.ctypes.data_as(ptr)
        if r.size:
            for f in self._funcs:
                self.lib.xc_gga_exc_vxc(f, r.size, address(r), address(s),
                                      address(part[0]), address(part[1]), address(part[2]))
                result += part
        # Exact vacuum is defined as zero; Libxc applies its native tiny-density
        # threshold internally (also to derivatives). No artificial density floor.
        result[:, r == 0] = 0.
        if not np.isfinite(result).all():
            raise FloatingPointError('Libxc returned a nonfinite semilocal result')
        return tuple(a.reshape(shape) for a in result)

    def close(self):
        for f in self._funcs:
            self.lib.xc_func_end(f)
            self.lib.xc_func_free(f)
        self._funcs.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
