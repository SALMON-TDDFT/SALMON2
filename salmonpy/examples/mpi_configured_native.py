"""Three-rank configured LD endpoint or managed TDDFT response regression.

Supply local catalog, material/initial-state ID, waveform NPZ (A and scalar dt),
new output directory; --mode is explicit. TDDFT also requires --runtime.
"""
import argparse
import hashlib
from pathlib import Path
import numpy as np


def main():
    from mpi4py import MPI
    from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
    from salmonpy.config import load_materials
    from salmonpy.initial_states import tddft_from_catalog
    from salmonpy.runtime import load_runtime
    from salmonpy.materials import bind_configured_ld_endpoints, bind_configured_tddft, response_identity
    from salmonpy.material_response import ResponseSession
    from salmonpy.distributed import ParallelPlan, ParallelContext, SolverGroup
    from salmonpy.point_runs import open_point_arrays
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog', type=Path)
    parser.add_argument('identity')
    parser.add_argument('waveform', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--mode', choices=('ld-endpoints', 'tddft'), required=True)
    parser.add_argument('--runtime', type=Path)
    args = parser.parse_args()
    if (args.mode == 'tddft') != (args.runtime is not None):
        parser.error('--runtime is required only for tddft')
    world = MPI.COMM_WORLD
    if world.size != 3:
        raise ValueError('exactly three MPI ranks required')
    model = None
    error = None
    try:
        with np.load(args.waveform, allow_pickle=False) as data:
            A = np.asarray(data['A'])
            raw_dt = np.asarray(data['dt'])
            if A.dtype.kind not in 'iuf' or raw_dt.dtype.kind not in 'iuf' or raw_dt.shape != ():
                raise ValueError('real numeric A and scalar dt required')
            A = A.astype(float)
            dt = float(raw_dt)
        if A.ndim != 2 or A.shape[1] != 3 or len(A) < 2 or not np.isfinite(A).all() or np.any(A[0]) or not np.isfinite(dt) or dt <= 0:
            raise ValueError('finite (steps+1,3) A with A[0]=0 and positive dt required')
        catalog = load_materials(args.catalog, registry=builtin_material_registry())
        base = args.catalog.resolve().parent
        if args.mode == 'tddft':
            engine = load_runtime(args.runtime, required=['tddft', 'single-rank'])
            model = tddft_from_catalog(catalog, args.identity, base_directory=base, executable=engine,
                                       max_steps=len(A)-1, threads=1)
            binding = bind_configured_tddft(model, key='configured', runtime=args.runtime)
        else:
            model = material_from_catalog(catalog, args.identity, base_directory=base)
            binding = bind_configured_ld_endpoints(model, key='configured')
        bindings = {10:binding, 90:binding}
        mapping = {10:'positive', 90:'negative'}
        record = response_identity(bindings, mapping, run_id='configured-native-example')
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    errors = world.allgather(error)
    if any(errors):
        if model is not None: model.close()
        raise ValueError(next(e for e in errors if e))
    signature = (dt, A.shape, hashlib.sha256(A.tobytes()).hexdigest())
    signatures = world.allgather(signature)
    if any(item != signature for item in signatures):
        model.close()
        raise ValueError('waveform differs across MPI ranks')
    plan = ParallelPlan((SolverGroup('field',(0,),role='field'),
                         SolverGroup('positive',(1,),point_ids=(10,)),
                         SolverGroup('negative',(2,),point_ids=(90,))))
    try:
        with ParallelContext.create(world, plan) as context:
            session = ResponseSession(context, bindings, dt=dt, point_record=record, point_map=mapping)
            history = []
            try:
                for field in A[1:]:
                    value = session.step(np.stack([field,-field]) if world.rank == 0 else None)
                    if world.rank == 0: history.append(value.copy())
                archive = session.save_arrays(args.directory if world.rank == 0 else None,
                    {'current':history,'A_over_c':np.stack([A[1:],-A[1:]],axis=1),'time':np.arange(1,len(A))*dt} if world.rank == 0 else None,
                    point_axes={'current':1,'A_over_c':1} if world.rank == 0 else None,
                    metadata={'units':'atomic','mode':args.mode,'scope':'interface regression'} if world.rank == 0 else None)
            finally:
                session.close()
        if world.rank == 0:
            reference = []
            for sign in (1.,-1.):
                with model.new_state(dt=dt) as state:
                    values = []
                    for old,new in zip(A[:-1],A[1:]):
                        if args.mode == 'tddft': current = state.step(sign*new)
                        else: current = state.step(-sign*(new-old)/dt)[1]
                        values.append(current.copy())
                    reference.append(values)
            np.testing.assert_array_equal(history,np.stack(reference,axis=1))
            loaded = open_point_arrays(archive.directory)
            assert loaded.metadata['point_run'] == record
            print('PASS direct/MPI currents and material/initial/point/run IDs:', args.mode, args.directory)
    finally:
        model.close()


if __name__ == '__main__':
    main()
