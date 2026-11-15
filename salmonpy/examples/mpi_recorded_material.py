"""Three-rank recorded custom/SBE model example, exact serial/MPI comparison.

mpiexec -n 3 python -m salmonpy.examples.mpi_recorded_material NEW_OUTPUT
Optional --catalog materials.json --material ID selects a configured SBE model.
"""
import argparse
from dataclasses import dataclass
import hashlib
from pathlib import Path
import numpy as np
from salmonpy.examples.custom_material import ExampleMaterial, waveform


@dataclass(frozen=True)
class RecordedExample(ExampleMaterial):
    def record(self):
        return {'material_id': 'custom-film',
                'definition_digest': hashlib.sha256(self.identity.encode()).hexdigest()}


def main():
    from salmonpy.distributed import ParallelPlan, ParallelContext, SolverGroup
    from salmonpy.materials import bind_recorded_material, response_identity
    from salmonpy.material_response import ResponseSession
    from salmonpy.crystal import respond
    from salmonpy.point_runs import open_point_arrays
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--catalog', type=Path)
    parser.add_argument('--material')
    args = parser.parse_args()
    from mpi4py import MPI
    if bool(args.catalog) != bool(args.material):
        parser.error('provide both catalog and material')
    world = MPI.COMM_WORLD
    if world.size != 3:
        raise ValueError('exactly three ranks required')
    model = None
    error = None
    try:
        if args.catalog:
            from salmonpy.config import load_materials
            from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
            catalog = load_materials(args.catalog, registry=builtin_material_registry())
            model = material_from_catalog(catalog, args.material, base_directory=args.catalog.resolve().parent)
        else:
            model = RecordedExample()
        binding = bind_recorded_material(model, key='prepared-material')
        bindings = {10: binding, 90: binding}
        mapping = {10: 'right', 90: 'left'}
        record = response_identity(bindings, mapping, run_id='recorded-example')
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    errors = world.allgather(error)
    if any(errors):
        raise ValueError(next(e for e in errors if e))
    plan = ParallelPlan((SolverGroup('field', (0,), role='field'),
                         SolverGroup('right-group', (1,), point_ids=(10,)),
                         SolverGroup('left-group', (2,), point_ids=(90,))))
    drives = [waveform(point) for point in plan.point_ids]
    with ParallelContext.create(world, plan) as context:
        session = ResponseSession(context, bindings, dt=.02, point_record=record, point_map=mapping)
        history = []
        try:
            for n in range(1, len(drives[0].time)):
                A = np.stack([drive.values[n] for drive in drives]) if world.rank == 0 else None
                value = session.step(A)
                if world.rank == 0:
                    history.append(value.copy())
            # Finish collective operations before root-only serial reference calculation.
            archive = session.save_arrays(args.directory if world.rank == 0 else None,
                {'current': history, 'time': drives[0].time[1:]} if world.rank == 0 else None,
                point_axes={'current': 1} if world.rank == 0 else None,
                metadata={'units':'atomic', 'drive':'A/c', 'current':'charge j',
                          'scope':'synthetic/interface regression'} if world.rank == 0 else None)
        finally:
            session.close()
    if world.rank == 0:
        reference = np.stack([respond(model, drive).current[1:] for drive in drives], axis=1)
        np.testing.assert_array_equal(history, reference)
        loaded = open_point_arrays(archive.directory)
        assert loaded.metadata['point_run'] == record
        print('PASS recorded model, serial/MPI equality and saved IDs:', args.directory)


if __name__ == '__main__':
    main()
