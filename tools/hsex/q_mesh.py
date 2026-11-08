"""Commensurate screening meshes in x-fastest order."""
from itertools import product

import numpy as np


def signed_indices(mesh):
    indices = np.indices(tuple(mesh[::-1])).reshape(3, -1)[::-1].T
    return np.where(2 * indices >= mesh, indices - mesh, indices)


def interpolate_boxes(coarse, qmesh, mesh, long_range):
    """Interpolate the short-range box kernel periodically, restoring the Coulomb tail."""
    sites = signed_indices(qmesh)
    for j, site in enumerate(sites):
        coarse[j] -= long_range(site / qmesh)
    result = np.empty((int(np.prod(mesh)), *coarse.shape[1:]), complex)
    for j, site in enumerate(signed_indices(mesh)):
        point = site * qmesh / mesh
        lower = np.floor(point).astype(int)
        fraction = point - lower
        value = long_range(site / mesh)
        for delta in product(range(2), repeat=3):
            weight = np.prod(np.where(delta, fraction, 1-fraction))
            if weight:
                index = np.ravel_multi_index(tuple((lower + delta) % qmesh), qmesh, order='F')
                value += weight * coarse[index]
        result[j] = value
    return result
