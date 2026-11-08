"""Compare a spin-degenerate metallic response with its spinor representation."""
from pathlib import Path
import sys

import numpy as np

from mpi_test_support import run_mpi, shared_fixture


def main():
    executable = Path(sys.argv[1])
    ranks = int(sys.argv[2])
    threads = int(sys.argv[3])
    mesh = np.array([6] * 3)
    n = np.array([3] * 3)
    nk, nreal = int(np.prod(mesh)), int(np.prod(n))
    sites = np.indices(tuple(mesh)).reshape(3, -1)[::-1].T
    k = ((sites - 2.5) / mesh) * (2 * np.pi / 9)
    variation = 1e-7 * np.random.default_rng(544).uniform(-1, 1, nk)
    signs = np.r_[-np.ones(nreal // 2), 0, np.ones(nreal // 2)]
    levels = np.linspace(-1e-4, 1e-4, nreal) + variation[:, None] * signs
    occupations = 2 / nreal - 100 * levels

    with shared_fixture() as path:
        results = []
        for components in (1, 2):
            ng = nreal * components
            energies = np.repeat(levels, components, axis=1)
            populations = np.repeat(occupations / components, components, axis=1)
            states = np.eye(ng, dtype='c16')
            if components == 2:
                # A unitary spin rotation introduces both components in each state.
                rotation = (0.5 + 0.5j) * np.array([[1, 1j], [1j, 1]])
                states = np.kron(rotation, np.eye(nreal))[:, np.arange(ng).reshape(2, -1).T.ravel()]
            with (path / 'pipeline.bin').open('wb') as stream:
                np.r_[n, mesh, mesh, 1, components, ng].astype('i4').tofile(stream)
                np.r_[.2, 1e-13, (9 * np.eye(3)).ravel()].astype('f8').tofile(stream)
                for array in (k, energies, populations):
                    array.tofile(stream)
                for _ in range(nk):
                    states.T.tofile(stream)
                # The Hamiltonian commutes with the spin rotation. Write it
                # directly so roundoff cannot split its exact degeneracies.
                for energy in levels:
                    np.diag(np.tile(energy, components)).astype('c16').tofile(stream)
            run_mpi(executable, ranks, threads, path, f'components-{components}',
                    timeout=600, extra_args=(str(components),))
            results.append(np.fromfile(path / 'pipeline-output.bin', dtype='f8'))
        np.testing.assert_allclose(results[1], results[0], atol=2e-12, rtol=2e-11)
        print('SPINOR_METAL_W_MAX_ERROR', np.max(np.abs(results[1] - results[0])))


if __name__ == '__main__':
    main()
