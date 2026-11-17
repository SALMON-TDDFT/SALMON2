"""Two independent k-only groups compare equivalent full/reduced fast-PR HSE.

Run with mpiexec -n 2 (group1) or -n 4 --group-size 2 (group2).
Arguments: BUILD FULL_BUNDLE REDUCED_BUNDLE NEW_OUTPUT.
Use --drive NPZ --point-index 1 to replay the original failing local FDTD A/c.
Without --drive, use a weak z field starting at zero. No fixture conversion,
build, download, tolerance relaxation or field projection is performed.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def main():
    from mpi4py import MPI
    from salmonpy.fast_hse import FastHSEMaterialProvider
    from salmonpy.fast_hse_mpi import NativeFastHSEBackend

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('full', type=Path)
    parser.add_argument('reduced', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--drive', type=Path)
    parser.add_argument('--group-size', type=int, choices=(1,2), default=1)
    parser.add_argument('--point-index', type=int, default=1)
    parser.add_argument('--dt', type=float, default=.003)
    parser.add_argument('--steps', type=int, default=32)
    args = parser.parse_args()
    world = MPI.COMM_WORLD
    if world.size != 2*args.group_size:
        raise ValueError('require exactly two solver groups of --group-size ranks')
    group = world.rank//args.group_size
    solver = world.Split(group,world.rank)
    try:
        if not np.isfinite(args.dt) or args.dt <= 0 or args.steps <= 0:
            raise ValueError('positive finite dt and positive step count required')
        if args.drive:
            with np.load(args.drive, allow_pickle=False) as data:
                drive = np.array(data['A_over_c'][:, args.point_index, :], copy=True)
        else:
            drive = np.zeros((args.steps, 3))
            drive[1:, 2] = np.geomspace(1e-23, 1e-8, args.steps-1)
        if drive.shape != (args.steps, 3) or not np.isfinite(drive).all():
            raise ValueError('drive must contain one finite A/c vector per step')
        folder = args.output/str(world.rank)
        folder.mkdir(parents=True, exist_ok=False)
        context = SimpleNamespace(solver_comm=solver, point_ids=(0,))
        bundle = args.full if group == 0 else args.reduced
        provider = FastHSEMaterialProvider(bundle, max_steps=args.steps)
        state = NativeFastHSEBackend(args.build, folder, allow_symmetry=True).create_state(
            provider.prepare(context), 0, args.dt, context)
        currents, energies = [], []
        try:
            for A in drive:
                currents.append(state.step(A))
                energies.append(state.energy)
        finally:
            state.close()
        np.savez(folder/'result.npz', current=currents, energy=energies, A_over_c=drive)
        world.Barrier()
        if world.rank == 0:
            with np.load(args.output/str(args.group_size)/'result.npz', allow_pickle=False) as reduced:
                full = np.asarray(currents)
                energy = np.asarray(energies)
                np.testing.assert_allclose(full, reduced['current'], rtol=2e-10, atol=2e-12)
                np.testing.assert_allclose(energy, reduced['energy'], rtol=2e-12, atol=2e-11)
                report = dict(schema='salmonpy.hse-low-gradient-mpi/1', status='pass',
                              dt=args.dt, steps=args.steps, solver_group_sizes=[args.group_size, args.group_size],
                              max_current_difference=float(np.max(abs(full-reduced['current']))),
                              max_energy_difference=float(np.max(abs(energy-reduced['energy']))),
                              caveat='Synthetic numerical regression, not physical convergence.')
                (args.output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
                print(json.dumps(report))
        world.Barrier()
    except BaseException:
        # A rank-local failure must not leave another rank hung at a barrier.
        import traceback
        traceback.print_exc()
        world.Abort(1)
        raise


if __name__ == '__main__':
    main()
