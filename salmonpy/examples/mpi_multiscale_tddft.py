"""One Python Multiscale1D rank coupled to two independent TDDFT groups.

mpiexec -n 5 python -m salmonpy.examples.mpi_multiscale_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
Use --material-ranks 1 with three MPI ranks. Set OpenMP/BLAS threads explicitly.
The host-MPI library and GS must already exist. Nothing is built on import.
This is a connection example with artificial pulse/settings, not converged optics.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI

from salmonpy.coupled import CoupledStepper
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.embedded_tddft import NativeTDDFTBackend, NativeTDDFTMaterialProvider
from salmonpy.fields import AcosPulse
from salmonpy.material_response import MaterialBinding
from salmonpy.maxwell import Multiscale1D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build_record')
    parser.add_argument('bundle')
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--material-ranks', type=int, choices=(1, 2), default=2)
    parser.add_argument('--steps', type=int, default=400)
    args = parser.parse_args()
    world = MPI.COMM_WORLD
    size = args.material_ranks
    if world.size != 1 + 2 * size:
        raise ValueError('launch 1 + 2*material-ranks MPI ranks')
    if args.steps < 1:
        raise ValueError('positive steps required')
    output = Path(args.output).resolve()
    if world.rank == 0:
        output.mkdir(parents=True, exist_ok=False)
        print(f'One field rank; two material groups of {size}; two electronic states.')
        print('Memory unmeasured. Start with a small supported GS; no automatic memory limit.')
    world.Barrier()
    plan = ParallelPlan((
        SolverGroup('field', (0,), role='field'),
        SolverGroup('left', tuple(range(1, size + 1)), point_ids=(1,)),
        SolverGroup('right', tuple(range(size + 1, 2 * size + 1)), point_ids=(2,)),
    ))
    # Same immutable GS, independent mutable electronic state in each group.
    # Replace provider/backend bindings to use other materials or current solvers.
    provider = NativeTDDFTMaterialProvider(args.bundle, args.input, max_steps=args.steps)
    backend = NativeTDDFTBackend(
        args.build_record, output, decomposition=(1, size, 1, 1, 1))
    binding = MaterialBinding('same-GS', provider, backend)
    dt = .02
    with ParallelContext.create(world, plan) as context:
        field = None
        if context.group.role == 'field':
            field = Multiscale1D(
                2, spacing=5., vacuum=(300, 50),
                incident=AcosPulse(.02, .5, 6., polarization=(0., 0., 1.)),
            ).new_state(dt=dt, steps=args.steps)
        coupled = CoupledStepper.from_multiscale(context, {1: binding, 2: binding}, field)
        try:
            history = []
            grids = [field.grid_A.copy()] if field is not None else []
            for _ in range(args.steps):
                frame = coupled.step()
                if frame is not None:
                    history.append({
                        key: value.tolist() if isinstance(value, np.ndarray) else value
                        for key, value in frame.items()
                    })
                    grids.append(field.grid_A.copy())
            communication = coupled.response.communication_report()
        finally:
            coupled.response.close()
            if coupled.field is not None:
                coupled.field.close()
    if world.rank == 0:
        np.savez(output / 'field.npz', x=field.x, grid_A=np.asarray(grids),
                 time=np.arange(args.steps + 1) * dt)
        (output / 'response.json').write_text(json.dumps({
            'units': 'atomic', 'point_ids': plan.point_ids, 'timing': coupled.timing,
            'history': history, 'initial_current': coupled.response.initial_current,
            'note': 'History starts at dt; initial feedback is zero, not a measured GS current.',
        }, indent=2) + '\n')
        (output / 'communication.json').write_text(json.dumps(communication, indent=2) + '\n')
        print(output / 'response.json')


if __name__ == '__main__':
    try:
        main()
    except BaseException:
        # No cross-group rollback. A launcher timeout also bounds native hangs.
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
