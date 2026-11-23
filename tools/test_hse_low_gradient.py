"""Verify the native HSE small-gradient probe against 80-digit energy differences.

BUILD_RECORD PROBE NEW_OUTPUT_DIRECTORY. Explicit execution only. The probe is
compiled from tests/native/test_hse_low_gradient.f90 against the selected library.
The Libxc components are accessed through its documented C ABI. No MPI/build on
import. This verifies numerical evaluation, not physical convergence.
"""
import argparse
import ctypes as c
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

PI = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
GAMMA = D('0.577215664901532860606512090082402431042159335939923598805767234884867726778')


def exchange_energy(rho, sigma):
    """Independent high-precision energy; derivatives use numerical differences."""
    one = D(1)
    a, b, c0, d, e0 = map(D, ('1.0161144', '-.37170836', '-.077215461', '.57786348', '-.051955731'))
    k = (3*PI*PI*rho)**(one/3)
    u = sigma/(4*k*k*rho*rho)
    h = (D('.00979681')*u+D('.0410834')*u*u)/(one+D('.187440')*u*u+D('.00120824')*u**D('2.5')+D('.0347188')*u**3)
    y = u*h
    z = d+y
    factor = 9/(4*a)
    x = factor*y
    if x:
        series = D(0)
        term = one
        for j in range(1, 16):
            term *= -x/j
            series += term/j
        t1 = a/2*((x.exp()-one)*(-GAMMA-x.ln())-x.exp()*series-GAMMA-factor.ln()-z.ln())
    else:
        t1 = a/2*(-GAMMA-factor.ln()-z.ln())
    f = D('6.4753871')*h+D('.47965830')
    eg = D('-.02628417880')-D('.07117647788')*u+D('.08534541323')*u*u
    t2 = (z*z*b+z*c0+2*e0+z*u*c0*f+2*u*eg)/(2*z**3)
    enhancement = -D(8)/9*(t1+t2)
    return -D(3)/4*(3/PI)**(one/3)*rho**(D(4)/3)*enhancement


def exchange_reference(rho, sigma):
    with localcontext() as ctx:
        ctx.prec = 80
        r, s = D(str(rho)), D(str(sigma))
        dr = r*D('1e-20')
        energy = exchange_energy(r, s)
        vrho = (exchange_energy(r+dr, s)-exchange_energy(r-dr, s))/(2*dr)
        if s:
            ds = s*D('1e-10')
            vsigma = (exchange_energy(r, s+ds)-exchange_energy(r, s-ds))/(2*ds)
        else:
            ds = 4*(3*PI*PI*r)**(D(2)/3)*r*r*D('1e-25')
            vsigma = (exchange_energy(r, ds)-energy)/ds
        return np.array([float(energy/r), float(vrho), float(vsigma)])


def screened_limit(rho, omega):
    # Independent 80-digit ENERGY difference. At u=0 the screened term1
    # depends on u only through y=u*H(sqrt(u))=O(u²), so its first derivative
    # vanishes. Evaluate terms2--5 without using the production F2 formula.
    with localcontext() as ctx:
        ctx.prec = 80
        r, w0 = D(str(rho)), D(str(omega))
        k = (3*PI*PI*r)**(D(1)/3)
        w = w0/k
        def energy(u):
            d, c0, b, e0 = map(D, ('.57786348', '-.077215461', '-.37170836', '-.051955731'))
            h = (D('.00979681')*u+D('.0410834')*u*u)/(1+D('.187440')*u*u+D('.00120824')*u**D('2.5')+D('.0347188')*u**3)
            z = d+u*h
            q = z+w*w
            f = D('6.4753871')*h+D('.47965830')
            g = D('-.02628417880')-D('.07117647788')*u+D('.08534541323')*u*u
            t2 = (z*z*b+z*c0+2*e0+z*u*c0*f+2*u*g)/(2*z**3)
            t3 = -w*(4*q*q*b+6*q*c0+15*e0+6*q*u*c0*f+15*u*g)/(8*z*q**D('2.5'))
            t4 = -w**3*(q*c0+5*e0+q*u*c0*f+5*u*g)/(2*z*z*q**D('2.5'))
            t5 = -w**5*(e0+u*g)/(z**3*q**D('2.5'))
            return -D(3)/4*(3/PI)**(D(1)/3)*r**(D(4)/3)*(-D(8)/9)*(t2+t3+t4+t5)
        du = D('1e-25')
        return float((energy(du)-energy(D(0)))/(du*4*k*k*r*r))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('probe', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record = json.loads(args.build.read_text())
    library = Path(record['library']).resolve(strict=True)
    digest = hashlib.sha256(library.read_bytes()).hexdigest()
    if digest != record['sha256']:
        raise ValueError('native library checksum mismatch')
    args.output.mkdir(parents=True, exist_ok=False)
    data = args.output/'sweep.txt'
    with (args.output/'probe.log').open('w') as log:
        subprocess.run([str(args.probe.resolve()), str(data.resolve())], stdout=log, stderr=subprocess.STDOUT, check=True)
    lib = c.CDLL(str(library))
    lib.xc_func_alloc.restype = c.c_void_p
    lib.xc_func_init.argtypes = [c.c_void_p, c.c_int, c.c_int]
    lib.xc_func_set_ext_params.argtypes = [c.c_void_p, c.c_void_p]
    lib.xc_gga_exc_vxc.argtypes = [c.c_void_p, c.c_size_t, *([c.c_void_p]*5)]
    lib.xc_func_end.argtypes = lib.xc_func_free.argtypes = [c.c_void_p]

    def raw(identifier, rho, sigma, params=None):
        func = lib.xc_func_alloc()
        if not func:
            raise RuntimeError('Libxc allocation failed')
        ready = False
        try:
            if lib.xc_func_init(func, identifier, 1):
                raise RuntimeError('Libxc initialization failed')
            ready = True
            if params is not None:
                parameters = np.array(params, dtype=np.float64)
                lib.xc_func_set_ext_params(func, parameters.ctypes.data)
            arrays = [np.array([rho]), np.array([sigma]), *[np.empty(1) for _ in range(3)]]
            lib.xc_gga_exc_vxc(func, 1, *[a.ctypes.data for a in arrays])
            return np.array([a[0] for a in arrays[2:]])
        finally:
            if ready:
                lib.xc_func_end(func)
            lib.xc_func_free(func)

    rows = np.loadtxt(data)
    errors = []
    relative_errors = []
    unchanged = 0
    screened_density = 0
    for rho, omega, reduced, sigma, energy, vrho, vsigma in rows:
        actual = np.array([energy, vrho, vsigma])
        if reduced > 1e-4 or rho <= 2e-15:
            expected = raw(428, rho, sigma, [.25, omega, omega])
            if not np.array_equal(actual.view(np.uint64), expected.view(np.uint64)):
                raise AssertionError(('ordinary-gradient bytes changed', rho, omega, reduced, actual, expected))
            unchanged += 1
            screened_density += int(reduced <= 1e-4)
            continue
        screened = raw(524, rho, sigma, [omega])
        if reduced**2 <= np.finfo(float).eps:
            screened[2] = screened_limit(rho, omega)
        expected = exchange_reference(rho, sigma)-.25*screened+raw(130, rho, sigma)
        error = abs(actual-expected)
        if not np.allclose(actual, expected, rtol=2e-12, atol=2e-13):
            raise AssertionError(('80-digit energy derivative mismatch', rho, omega, reduced, actual, expected, error))
        errors.append(error)
        relative_errors.append(error/np.maximum(abs(expected), np.finfo(float).tiny))
    report = dict(schema='salmonpy.hse-low-gradient-verification/1', status='pass',
                  library_sha256=digest, cases=len(rows), high_precision_cases=len(errors),
                  preserved_byte_equal_cases=unchanged, ordinary_gradient_byte_equal_cases=unchanged-screened_density,
                  density_screening_byte_equal_cases=screened_density, decimal_precision=80,
                  max_absolute_errors=np.max(errors, axis=0).tolist(),
                  max_relative_errors=np.max(relative_errors, axis=0).tolist(),
                  caveat='Numerical evaluation check; not a physical convergence study.')
    (args.output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
