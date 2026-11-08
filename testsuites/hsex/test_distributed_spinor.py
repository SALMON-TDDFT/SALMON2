"""Compare SOC screening samples across replicated and uneven k ownership."""
from pathlib import Path
import sys

import numpy as np

from mpi_test_support import run_mpi, shared_fixture


def main():
    executable = Path(sys.argv[1]).resolve()
    mesh = np.array([6] * 3)
    grid = np.array([3] * 3)
    nk, ng, nb = int(np.prod(mesh)), 2 * int(np.prod(grid)), 3
    frac = (np.indices(tuple(mesh)).reshape(3, -1)[::-1].T - 2.5) / 6
    variation = np.cos(2 * np.pi * frac[:, 0])
    filled = .55 + .01 * variation
    occupations = np.column_stack((filled, 1 - filled, np.zeros(nk)))
    energies = np.column_stack((.002 * np.log((1 - filled) / filled),
                                .002 * np.log(filled / (1 - filled)), np.full(nk, .5)))
    # The source orbitals vary with k, so an incorrect owner or packing order
    # changes the response. The empty complement is gapped and cheap to solve.
    orbitals = []
    for angle in .005 * variation:
        basis = np.eye(ng, dtype='c16')
        basis[np.ix_([0, 3], [0, 3])] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        orbitals.append(basis)

    with shared_fixture() as path:
        with (path / 'pipeline.bin').open('wb') as stream:
            np.r_[grid, mesh, mesh, 1, 1, nb].astype('i4').tofile(stream)
            np.r_[.2, .001, (9 * np.eye(3)).ravel()].astype('f8').tofile(stream)
            for array in (frac * (2 * np.pi / 9), energies, occupations):
                array.tofile(stream)
            for basis in orbitals:
                basis[:, :nb].T.tofile(stream)
            for basis, levels in zip(orbitals, energies):
                ((basis * np.r_[levels, np.ones(ng - nb)]) @ basis.conj().T).T.tofile(stream)
        outputs = []
        for mode in ('replicated', 'distributed'):
            # This fixture deliberately fails the smoothness check. The
            # finite-q samples must still agree before the head is rejected.
            run_mpi(executable, 5, 1, path, mode, expect_failure=True, extra_args=('2', mode))
            transcript = (path / mode / 'stdout.rank-0').read_text()
            assert 'HSEX_W_REJECT stage=70 error=4' in transcript, transcript
            sample_file = path / 'hsex-spinor-head-samples.bin'
            outputs.append(sample_file.read_bytes())
            sample_file.rename(path / f'{mode}-samples.bin')
        assert outputs[0] == outputs[1], 'finite-q samples differ under uneven k ownership'
        print('DISTRIBUTED_SPINOR_SAMPLES_BITWISE_PASS')


if __name__ == '__main__':
    main()
