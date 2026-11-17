"""Explicit public symmetry opt-in and direct restart on independent MPI groups.

mpiexec -n 1|2 python TOOL BUILD BUNDLE NEW_OUTPUT --mode save|resume --checkpoint DIR
Use --group-size2 for two internal k ranks per independent group.
The deny-old-k2 mode targets an old library lacking the k-MPI capability. Input states are synthetic
when using prepare_fast_hse_symmetry_fixture.py; no physical convergence claim.
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
    from salmonpy.errors import InputError
    from salmonpy import fast_hse_mpi, embedded_tddft

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'bundle', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--mode', choices=('save', 'resume', 'deny-default', 'deny-old-k2', 'deny-old'), required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=32)
    parser.add_argument('--group-size', type=int, choices=(1,2), default=1)
    parser.add_argument('--admission-fault', choices=('truncated','nonfinite','nonclosed'))
    parser.add_argument('--split', type=int, default=7)
    args = parser.parse_args()
    if not 3 < args.split < args.steps:
        raise ValueError('require 3 < split < steps')
    world = MPI.COMM_WORLD
    if world.size % args.group_size:
        raise ValueError('world size must be divisible by group size')
    comm = world if args.mode == 'deny-old-k2' else world.Split(world.rank//args.group_size,world.rank)
    point = 0 if args.mode == 'deny-old-k2' else world.rank//args.group_size
    root = args.output.resolve() / f'rank-{world.rank}'
    root.mkdir(parents=True, exist_ok=False)
    bundle = root / 'initial'
    shutil.copytree(args.bundle, bundle)
    text = (bundle / 'rt.inp').read_text().replace("xc='hse06'", f"xc='hse06',hse_omega={.11+.01*point:.17g}")
    (bundle / 'rt.inp').write_text(text)
    provider = FastHSEMaterialProvider(bundle, max_steps=args.steps)
    context = SimpleNamespace(solver_comm=comm, point_ids=(point,))
    backend = NativeFastHSEBackend(args.build, root, allow_symmetry=args.mode != 'deny-default')
    if args.mode.startswith('deny-'):
        try:
            state = backend.create_state(provider.prepare(context), point, .001, context)
        except InputError:
            assert not fast_hse_mpi._used and not embedded_tddft._used
        else:
            state.close()
            raise AssertionError('unsupported/default symmetry admitted')
        world.Barrier()
        if world.rank == 0:
            print('Public symmetry rejection before lifetime PASS', args.mode, flush=True)
        return
    if args.admission_fault:
        path = bundle/'sym.dat'
        original = path.read_text()
        if comm.rank == comm.size-1:
            if args.admission_fault == 'truncated':
                damaged = original+'1 0 0 0\n'
            elif args.admission_fault == 'nonfinite':
                damaged = original.replace('1.0','NaN',1)
            else:
                damaged = '1 0 0 0\n0 1 0 0\n0 0 1 0\n0 -1 0 0\n1 0 0 0\n0 0 1 0\n'
            path.write_text(damaged)
        try:
            invalid_provider = FastHSEMaterialProvider(bundle,max_steps=args.steps)
            rejected = backend.create_state(invalid_provider.prepare(context), point, .001, context)
        except InputError as error:
            assert 'sym.dat' in str(error)
            assert not fast_hse_mpi._used and not embedded_tddft._used
        else:
            rejected.close()
            raise AssertionError('malformed rank-local symmetry admitted')
        path.write_text(original)
        comm.Barrier()
    state = backend.create_state(provider.prepare(context), point, .001, context)
    target = args.checkpoint.resolve() / f'rank-{world.rank}'
    if args.mode == 'save':
        target.mkdir(parents=True, exist_ok=False)
    shard = target / 'checkpoint.bin'
    trace = target / 'trace.npz'
    endpoint = lambda step: np.array([0., 0., .0001*np.sin(.3*step)])
    try:
        first = 0
        if args.mode == 'resume':
            state.checkpoint(shard, True)
            first = state.steps
            assert first == args.split
            with np.load(trace) as expected:
                np.testing.assert_array_equal(state.current, expected['current'])
                np.testing.assert_array_equal(state.A, expected['A'])
        rows = []
        for step in range(first+1, args.steps+1):
            if step == 4:
                before = root / 'before-field.bin'
                after = root / 'after-field.bin'
                state.checkpoint(before, False)
                old = (state.steps, state.A.tobytes(), state.current.tobytes(), state.energy)
                try:
                    state.step([1e-250, 0, 0])
                except InputError as error:
                    assert 'little group' in str(error)
                else:
                    raise AssertionError('forbidden field accepted')
                assert not state.failed
                assert old == (state.steps, state.A.tobytes(), state.current.tobytes(), state.energy)
                state.checkpoint(after, False)
                assert before.read_bytes() == after.read_bytes()
            state.step(endpoint(step))
            if step == args.split and args.mode == 'save':
                state.checkpoint(shard, False)
                saved_current = state.current.copy()
            if step > args.split:
                rows.append(np.r_[state.current, state.energy])
        if args.mode == 'save':
            np.savez(trace, trace=np.array(rows), current=saved_current, A=endpoint(args.split))
            state.checkpoint(target / 'final.bin', False)
        else:
            with np.load(trace) as expected:
                np.testing.assert_array_equal(np.array(rows), expected['trace'])
            state.checkpoint(root / 'final.bin', False)
            assert (root / 'final.bin').read_bytes() == (target / 'final.bin').read_bytes()
        (root / 'result.json').write_text(json.dumps(dict(status='passed', mode=args.mode,
            point_id=point, internal_ranks=comm.size, steps=args.steps, resumed_calls=args.steps-first,
            elapsed_steps_replayed=0, field_rejection_preserved=True,
            admission_fault=args.admission_fault, admission_retry_preserved=True if args.admission_fault else None))+'\n')
        world.Barrier()
        if world.rank == 0:
            print('Public symmetry direct checkpoint/field retry PASS', args.mode, flush=True)
    finally:
        state.close()


if __name__ == '__main__':
    main()
