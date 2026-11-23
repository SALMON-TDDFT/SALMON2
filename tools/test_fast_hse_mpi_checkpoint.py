"""Run with mpiexec -n 2; invoke save, then fresh jobs for each other mode.

Tests direct native rank shards, deliberately bypassing manifest authentication.
Use an external launcher timeout to detect collective hangs. Input files remain
unchanged: negative modes copy shards into their own new work directory.
"""
import argparse
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np


def main():
    from mpi4py import MPI
    from salmonpy.fast_hse import FastHSEMaterialProvider
    from salmonpy.fast_hse_mpi import NativeFastHSEBackend
    from salmonpy.errors import RuntimeFailure

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'bundle', 'workroot'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--mode', choices=('save', 'resume', 'truncate', 'swap'), required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--input-file', default='rt.inp')
    parser.add_argument('--dt', type=float, default=.02)
    parser.add_argument('--split', type=int, default=7)
    parser.add_argument('--steps', type=int, default=12)
    args = parser.parse_args()
    comm = MPI.COMM_WORLD
    if comm.size != 2:
        raise ValueError('exactly two MPI ranks required')
    if not 0 < args.split < args.steps:
        raise ValueError('require 0 < split < steps')
    work = args.workroot.resolve()
    checkpoint = args.checkpoint.resolve()
    error = None
    if comm.rank == 0:
        try:
            work.mkdir(parents=True, exist_ok=False)
            if args.mode == 'save':
                checkpoint.mkdir(parents=True, exist_ok=False)
            elif not checkpoint.is_dir():
                raise ValueError('checkpoint directory does not exist')
        except Exception as exc:
            error = str(exc)
    error = comm.bcast(error, root=0)
    if error:
        raise ValueError(error)
    context = SimpleNamespace(solver_comm=comm, point_ids=(0,))
    provider = FastHSEMaterialProvider(args.bundle, args.input_file, max_steps=args.steps)
    state = NativeFastHSEBackend(args.build, work).create_state(provider.prepare(context), 0, args.dt, context)
    shard = checkpoint / f'rank-{comm.rank}.bin'
    reference = checkpoint / f'rank-{comm.rank}.npz'

    def endpoint(step):
        return np.array([.001*np.sin(.17*step), .0003*np.cos(.11*step), 0.])

    try:
        if args.mode == 'save':
            rows = []
            for step in range(1, args.steps+1):
                state.step(endpoint(step))
                if step == args.split:
                    state.checkpoint(shard, False)
                    saved_current = state.current.copy()
                if step > args.split:
                    rows.append(np.r_[state.current, state.energy])
            np.savez(reference, trace=np.array(rows), current=saved_current,
                     A=endpoint(args.split), split=args.split, steps=args.steps, dt=args.dt)
        elif args.mode == 'resume':
            with np.load(reference) as expected:
                assert int(expected['split']) == args.split and int(expected['steps']) == args.steps
                assert float(expected['dt']) == args.dt
                state.checkpoint(shard, True)
                assert state.steps == args.split
                np.testing.assert_array_equal(state.A, expected['A'])
                np.testing.assert_array_equal(state.current, expected['current'])
                rows = []
                for step in range(args.split+1, args.steps+1):
                    state.step(endpoint(step))
                    rows.append(np.r_[state.current, state.energy])
                np.testing.assert_array_equal(np.array(rows), expected['trace'])
        else:
            damaged = work / f'negative-{comm.rank}.bin'
            source = checkpoint / f'rank-{1-comm.rank if args.mode == "swap" else comm.rank}.bin'
            shutil.copyfile(source, damaged)
            if args.mode == 'truncate' and comm.rank == 1:
                with damaged.open('r+b') as stream:
                    stream.truncate(damaged.stat().st_size//2)
            comm.Barrier()
            initial_current = state.current.copy()
            rejected = False
            try:
                state.checkpoint(damaged, True)
            except RuntimeFailure:
                rejected = True
            results = comm.allgather((rejected, state.steps, state.A.tolist(), state.current.tolist()))
            assert all(row[0] and row[1] == 0 and row[2] == [0., 0., 0.] for row in results), results
            np.testing.assert_array_equal(state.current, initial_current)
        comm.Barrier()
        if comm.rank == 0:
            print('SALMONPY_HSE_MPI_CHECKPOINT_PASS ' + json.dumps({'mode': args.mode, 'ranks': comm.size}), flush=True)
    finally:
        state.close()


if __name__ == '__main__':
    main()
