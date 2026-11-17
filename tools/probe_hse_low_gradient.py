"""Explicit diagnostic of Libxc HSE gradient derivatives; does not change physics.

BUILD_RECORD OUTPUT_JSON. Loads only the hash-verified selected native library.
This is a sensitivity diagnostic, not a physical validation or a library patch.
"""
import argparse
import ctypes as c
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record = json.loads(args.build.read_text())
    library = Path(record['library']).resolve(strict=True)
    with library.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    if digest != record['sha256']:
        raise ValueError('native library checksum mismatch')
    lib = c.CDLL(str(library))
    lib.xc_version_string.restype = c.c_char_p
    lib.xc_func_alloc.restype = c.c_void_p
    lib.xc_func_init.argtypes = [c.c_void_p, c.c_int, c.c_int]
    lib.xc_func_set_ext_params.argtypes = [c.c_void_p, c.c_void_p]
    lib.xc_gga_exc_vxc.argtypes = [c.c_void_p, c.c_size_t, *([c.c_void_p]*5)]
    lib.xc_func_end.argtypes = [c.c_void_p]
    lib.xc_func_free.argtypes = [c.c_void_p]
    func = lib.xc_func_alloc()
    if not func:
        raise RuntimeError('Libxc allocation failed')
    initialized = False
    try:
        if lib.xc_func_init(func, 428, 1):
            raise RuntimeError('Libxc HSE initialization failed')
        initialized = True
        params = np.array([.25, .11, .11], dtype=np.float64)
        lib.xc_func_set_ext_params(func, params.ctypes.data)
        gradient = np.r_[0., np.logspace(-24, -4, 21)]
        rho = np.full(len(gradient), .008)
        sigma = gradient**2
        energy = np.empty_like(rho)
        vrho = np.empty_like(rho)
        vsigma = np.empty_like(rho)
        lib.xc_gga_exc_vxc(func, len(rho), *[a.ctypes.data for a in (rho, sigma, energy, vrho, vsigma)])
        derivative = 2*gradient*vsigma
        finite = all(np.isfinite(a).all() for a in (energy, vrho, vsigma, derivative))
        if not finite:
            raise RuntimeError('nonfinite HSE diagnostic')
        report = dict(schema='salmonpy.hse-gradient-diagnostic/1', status='diagnostic',
                      library_sha256=digest, libxc_version=lib.xc_version_string().decode(),
                      rho=.008, omega=.11, gradient=gradient.tolist(), energy_per_particle=energy.tolist(),
                      vrho=vrho.tolist(), vsigma=vsigma.tolist(), gradient_derivative=derivative.tolist(),
                      low_gradient_sensitivity_detected=bool(np.max(abs(derivative[gradient < 1e-12])) > 1e-6),
                      caveat='No smoothing/clamping, no physical validation, no source/library modification.')
        with args.output.open('x') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    finally:
        if initialized:
            lib.xc_func_end(func)
        lib.xc_func_free(func)
    print(args.output)


if __name__ == '__main__':
    main()
