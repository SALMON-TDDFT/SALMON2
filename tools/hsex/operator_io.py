"""SALMON column-major dense (v1) and packed Hermitian lower triangle (v2)."""
from pathlib import Path

import numpy as np


def read_operator(path, offset, ng, packed):
    count = ng * (ng + 1) // 2 if packed else ng * ng
    if Path(path).stat().st_size != offset + 16 * count:
        raise ValueError(f'incomplete grid Hamiltonian: {path}')
    data = np.memmap(path, '=c16', 'r', offset=offset, shape=(count,))
    if not packed:
        return data.reshape((ng, ng), order='F')
    matrix = np.empty((ng, ng), dtype=complex, order='F')
    start = 0
    for col in range(ng):
        column = data[start:start + ng - col]
        matrix[col:, col] = column
        matrix[col, col+1:] = column[1:].conj()
        start += len(column)
    return matrix
