"""Finite-mesh, one-spin density-matrix distance to a free-electron ensemble.

Equal spin multiplicities cancel in the normalized distance. All momenta are
canonical atomic units; comparisons use a common vector-potential gauge.
"""
import itertools
import numpy as np


def momentum_cube(n, spacing, origin=None):
    if origin is None:
        origin = -(n-1)/2
    axis = (np.arange(n)+origin)*spacing
    return np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1)


def global_indices(k, ngrid, length, mesh):
    spacing = 2*np.pi/(mesh*length)
    size = ngrid*mesh
    origin = -ngrid*mesh/2-(mesh-1)/2
    q = momentum_cube(size,spacing,origin)
    gaxis = np.fft.fftfreq(ngrid)*ngrid*2*np.pi/length
    g = np.stack(np.meshgrid(gaxis,gaxis,gaxis,indexing='ij'),axis=-1).reshape(-1,3,order='F')
    coordinates = (k[:,None,:]+g[None,:,:])/spacing-origin
    ints = np.rint(coordinates).astype(int)
    if np.max(abs(coordinates-ints))>1e-9:
        raise ValueError('Unexpected shifted k mesh')
    indices = np.ravel_multi_index(tuple(ints.transpose(2,0,1)),(size,)*3)
    return indices,q


def fourier_coefficients(u, ngrid, dv):
    return np.array([np.fft.fftn(v.reshape((ngrid,)*3+(u.shape[-1],),order='F'),
                    axes=(0,1,2),norm='ortho').reshape(-1,u.shape[-1],order='F')*np.sqrt(dv) for v in u])


def fermi_reference(q, number):
    energy = np.sum(q*q,axis=-1)/2
    if not 0<number<energy.size:
        raise ValueError('Reference number outside finite basis')
    threshold = np.sort(energy.ravel())[int(np.ceil(number))-1]
    shell = np.isclose(energy,threshold,rtol=0,atol=1e-12)
    lower = (energy<threshold)&~shell
    f = lower.astype(float)
    f[shell] = (number-lower.sum())/shell.sum()
    return f


def mean_momentum(f,q):
    return np.sum(f[...,None]*q,axis=(0,1,2))/f.sum()


def flow_reference(base, target, spacing):
    """Convex combination of the eight adjacent integer boosts of a resting sea.

Matches mean canonical momentum exactly, but fractional boosts mix ensembles.
No cyclic momentum wrapping is permitted.
    """
    displacement = np.asarray(target)/spacing
    lo = np.floor(displacement).astype(int)
    fraction = displacement-lo
    result = np.zeros_like(base)
    support = np.argwhere(base>0)
    for bits in itertools.product((0,1),repeat=3):
        shift = lo+bits
        weight = np.prod(np.where(bits,fraction,1-fraction))
        if weight<1e-16:
            continue
        if np.any(support.min(axis=0)+shift<0) or np.any(support.max(axis=0)+shift>=base.shape):
            raise ValueError('Boost reaches finite momentum boundary')
        result += weight*np.roll(base,shift,axis=(0,1,2))
    return result


def hs_distance(c,f):
    """Mean over k of Tr[(c c† - diag f)^2], including actual orbital norm."""
    gram = c.conj().transpose(0,2,1)@c
    purity = np.sum(abs(gram)**2)/len(c)
    diagonal = np.sum(abs(c)**2,axis=2)
    value = purity+np.sum(f*f)/len(c)-2*np.sum(f*diagonal)/len(c)
    if value < -1e-9:
        raise ValueError('Negative squared distance')
    return max(0.,float(value))


def rank_lower_bound(f,rank):
    """Best possible distance for rank-rank projectors separately at every k."""
    return float(rank+np.sum(f*f)/len(f)-2*np.sum(np.sort(f,axis=1)[:,-rank:])/len(f))


def alpha_from_distance(distance, initial):
    if initial<=0:
        raise ValueError('Initial reference distance must be positive')
    return .2*np.asarray(distance)/initial
