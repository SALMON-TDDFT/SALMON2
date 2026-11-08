"""Independent compensated reference for original-H iterative refinement.

Split real operands into fixed-point slices with exact dot products, then add
all slice products with error-free summation. This is validation-only Python;
production SALMON uses the Fortran compensated kernel.
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def _slices(a, axis, bits):
    rest = np.asarray(a, dtype=float).copy()
    result = []
    for _ in range(6):
        exponent = np.frexp(np.max(abs(rest), axis=axis, keepdims=True))[1]
        part = np.ldexp(np.rint(np.ldexp(rest, bits-exponent)), exponent-bits)
        result.append(part)
        rest -= part
    return result


def precise_product(a, b):
    """About double-double precision for finite, normally scaled test matrices."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    bits = (51-int(np.ceil(np.log2(max(1,a.shape[1])))))//2
    if bits < 8:
        raise ValueError('reference matrix too large for fixed-point splitting')
    ar, ai = _slices(a.real, 1, bits), _slices(a.imag, 1, bits)
    br, bi = _slices(b.real, 0, bits), _slices(b.imag, 0, bits)
    total = np.zeros((a.shape[0],b.shape[1]),complex)
    correction = np.zeros_like(total)
    for order in range(11):
        for i in range(6):
            j = order-i
            if not 0 <= j < 6:
                continue
            for term in (ar[i]@br[j], -(ai[i]@bi[j]), 1j*(ar[i]@bi[j]), 1j*(ai[i]@br[j])):
                updated = total+term
                z = updated-total
                correction += (total-(updated-z))+(term-z)
                total = updated
    return total, correction


def _precise_diagonal_product(diagonal, x):
    """Independent fixed-point elementwise product, retaining its low part."""
    a = np.asarray(diagonal, complex)[:, None]
    # No inner summation: two 26-bit fixed-point factors multiply exactly.
    ar, ai = _slices(a.real, 1, 26), _slices(a.imag, 1, 26)
    br, bi = _slices(x.real, 1, 26), _slices(x.imag, 1, 26)
    total = np.zeros_like(x)
    correction = np.zeros_like(x)
    for order in range(11):
        for i in range(6):
            j = order-i
            if not 0 <= j < 6:
                continue
            for term in (ar[i]*br[j], -(ai[i]*bi[j]), 1j*(ar[i]*bi[j]), 1j*(ai[i]*br[j])):
                updated = total+term
                z = updated-total
                correction += (total-(updated-z))+(term-z)
                total = updated
    return total, correction


def original_residual(h, energy, eta, x, b, xlow=None):
    # Form E-Hii before multiplication, but retain subtraction roundoff so
    # the operation remains the original H even far from resonance.
    diagonal = h.diagonal().real
    off_diagonal = h.copy()
    off_diagonal[np.diag_indices(len(h))] -= diagonal
    shift = energy-diagonal
    z = shift-energy
    shift_low = (energy-(shift-z))+(-diagonal-z)
    total = b.copy()
    error = np.zeros_like(total)
    for solution in (x, xlow):
        if solution is None or not np.any(solution):
            continue
        high, low = precise_product(off_diagonal, solution)
        sh, sl = _precise_diagonal_product(shift+1j*eta, solution)
        lh, ll = _precise_diagonal_product(shift_low, solution)
        for term in (-high, -low, sh, sl, lh, ll):
            updated = total+term
            z = updated-total
            error += (total-(updated-z))+(term-z)
            total = updated
    return total+error


def original_gram(h, b, energy, eta):
    """Dense-LU reference; no shared native solver or compensated implementation."""
    op = h-(energy+1j*eta)*np.eye(len(h))
    factors = lu_factor(op)
    x = lu_solve(factors,b)
    xlow = np.zeros_like(x)
    previous = None
    for iteration in range(10):
        dx = lu_solve(factors, original_residual(h,energy,eta,x,b,xlow))
        updated = x+dx
        z = updated-x
        xlow += (x-(updated-z))+(dx-z)
        x = updated
        updated = x+xlow
        z = updated-x
        xlow = (x-(updated-z))+(xlow-z)
        x = updated
        high, low = precise_product(b.conj().T,x)
        lh, ll = precise_product(b.conj().T,xlow)
        gram = .5*(high+high.conj().T)+.5*(low+low.conj().T)
        gram += .5*(lh+lh.conj().T)+.5*(ll+ll.conj().T)
        if previous is not None and np.all(abs(gram-previous) <= 1e-12*(1+abs(gram))):
            break
        previous = gram
    else:
        raise ValueError('original-H response refinement did not converge')
    residual = original_residual(h,energy,eta,x,b,xlow)
    err = np.max(np.linalg.norm(residual,axis=0)/np.maximum(np.linalg.norm(b,axis=0),1e-14))
    if not np.isfinite(err) or err > 1e-8 or not np.isfinite(gram).all():
        raise ValueError(f'original-H refinement residual {err}')
    return gram, err
