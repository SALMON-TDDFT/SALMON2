"""Native E/H with field-local LD, SBE and internally parallel TDDFT.

mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD RUNTIME GS_BUNDLE NEW_OUTPUT
Use --tddft-ranks 1 with 3 ranks for a decomposition comparison.
Use --field-grid 1 2 1 --tddft-ranks 1 with 4 ranks for distributed fields.
Tiny 3D connection fixture, not recommended production TDDFT geometry.
"""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from salmonpy.coupled import NativeEHField, CoupledStepper
from salmonpy.distributed import ParallelPlan, ParallelContext, SolverGroup
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider, NativeTDDFTBackend
from salmonpy.material_response import MaterialBinding, ResponseSession
from salmonpy.materials import bind_material
from salmonpy.examples.mpi_mixed_multiscale import models
from salmonpy.fdtd import FDTD, FDTDMPI, PlaneWave, SalmonExecutable
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.ld import LorentzDrude
from salmonpy.fields import AcosPulse
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'runtime', 'bundle', 'output'): parser.add_argument(name, type=Path)
    parser.add_argument('--fast-hse-build', type=Path, help='explicit fast-PR serial worker build record for hse-point')
    parser.add_argument('--fast-hse-mpi-build', type=Path, help='explicit native fast-PR host-MPI build record')
    parser.add_argument('--hse-ranks', type=int, choices=(1, 2), default=1)
    parser.add_argument('--dt', type=float, default=.003)
    parser.add_argument('--tddft-ranks', type=int, choices=(1, 2), default=2)
    parser.add_argument('--layers', type=Path, help='native-layers/1 JSON referencing sbe-point, si-tddft, hse-point, drude, lorentz')
    parser.add_argument('--memory-budget-gib', type=float, help='aggregate advisory memory budget across MPI ranks; warning only')
    parser.add_argument('--rank-budget', type=int, help='planning cap; defaults to the launched world size')
    parser.add_argument('--field-grid', type=int, nargs=3, default=(1, 1, 1), metavar=('PX', 'PY', 'PZ'))
    parser.add_argument('--compare', type=Path, help='compare against an existing same-physics run with another MPI decomposition')
    parser.add_argument('--steps', type=int, default=128)
    parser.add_argument('--without-ld', action='store_true', help='vacuum control for local LD regions')
    parser.add_argument('--save-at', type=int, help='stop after this step and save a replay checkpoint')
    parser.add_argument('--resume', type=Path, help='reconstruct and continue a replay checkpoint')
    parser.add_argument('--direct-save-at', type=int, help='save internal states and stop at this step')
    parser.add_argument('--direct-resume', type=Path, help='protected internal-state restore without replay; candidate-capable runtime required')
    parser.add_argument('--unprotected-direct-restore', action='store_true', help='explicit legacy restore with no failure-preservation guarantee')
    args = parser.parse_args()
    if args.unprotected_direct_restore and args.direct_resume is None:
        parser.error('unprotected-direct-restore requires direct-resume')
    if args.direct_save_at is not None and not 0 < args.direct_save_at < args.steps:
        parser.error('direct-save-at must be between zero and steps')
    if (args.direct_save_at is not None or args.direct_resume is not None) and (args.save_at is not None or args.resume is not None):
        parser.error('direct and replay checkpoint options cannot be mixed')
    if args.save_at is not None and not 0 < args.save_at < args.steps:
        parser.error('save-at must be between zero and steps')
    if args.layers is not None and args.without_ld: parser.error('express vacuum in the layout; do not combine --layers with --without-ld')
    if not np.isfinite(args.dt) or args.dt <= 0: parser.error('positive finite dt required')
    if (args.fast_hse_build is not None or args.fast_hse_mpi_build is not None) and args.layers is None:
        parser.error('fast HSE requires an explicit layer layout')
    if args.hse_ranks != 1 and args.fast_hse_mpi_build is None:
        parser.error('hse-ranks > 1 requires --fast-hse-mpi-build')
    if args.memory_budget_gib is not None and (not np.isfinite(args.memory_budget_gib) or args.memory_budget_gib <= 0):
        parser.error('positive finite memory-budget-gib required')
    if args.steps < 1: parser.error('positive steps required')
    if any(v < 1 for v in args.field_grid): parser.error('positive field-grid entries required')
    field_ranks = int(np.prod(args.field_grid))
    from mpi4py import MPI
    world = MPI.COMM_WORLD
    try:
        if args.layers is None and world.size != field_ranks + 1 + args.tddft_ranks:
            raise ValueError('launch product(field-grid) + 1 + tddft-ranks MPI ranks')
        output = args.output.resolve()
        if world.rank == 0: output.mkdir(parents=True, exist_ok=False)
        world.Barrier()
        dt = args.dt
        grid = Grid((40, 40, 40), (2., 2., 2.))
        points = {10: (19, 20, 20), 20: (20, 20, 20)}
        shapes = (Shape('box', (2., 4., 4.), (7., 1., 1.), material=1),
                  Shape('box', (2., 4., 4.), (9., 1., 1.), material=2))
        materials = (LorentzDrude(plasma_frequency=20., strength=[1.], resonance=[0.], damping=[2.]),
                     LorentzDrude(plasma_frequency=15., strength=[.7, .3], resonance=[25., 35.], damping=[2., 3.]))
        geometry = Geometry(grid, () if args.without_ld else shapes)
        observations = tuple(tuple(grid.coordinates[d][p[d]] for d in range(3)) for p in points.values())
        calculation = FDTD(geometry, () if args.without_ld else materials, dt, args.steps,
                           source=PlaneWave(AcosPulse(.001, 20., .15), 0, -16.), observations=observations)
        _, sbe = models()
        layout = json.loads(args.layers.read_text()) if args.layers is not None else None
        selected = {layer['material_id'] for layer in layout.get('layers', [])} if layout is not None else {'sbe-point','si-tddft'}
        available = {'drude': materials[0], 'lorentz': materials[1]}
        if 'sbe-point' in selected:
            available['sbe-point'] = bind_material(sbe, key='sbe-point', identity='synthetic-two-band-v1')
        if 'si-tddft' in selected:
            available['si-tddft'] = MaterialBinding('si-tddft',
                NativeTDDFTMaterialProvider(args.bundle, 'rt.inp', max_steps=args.steps),
                NativeTDDFTBackend(args.build, output, decomposition=(1, args.tddft_ranks, 1, 1, 1)))
        if 'hse-point' in selected:
            if args.fast_hse_build is None: raise ValueError('hse-point requires --fast-hse-build')
            from salmonpy.fast_hse import FastHSEMaterialProvider, FastHSEResponseBackend
            from salmonpy.fast_hse_mpi import NativeFastHSEBackend
            available['hse-point'] = MaterialBinding('hse-point',
                FastHSEMaterialProvider(args.bundle, 'rt.inp', max_steps=args.steps),
                (FastHSEResponseBackend(args.fast_hse_build, output) if args.fast_hse_mpi_build is None
                 else NativeFastHSEBackend(args.fast_hse_mpi_build, output)))
        assembly = None
        if layout is None:
            bindings = {10: available['sbe-point'], 20: available['si-tddft']}
            plan = ParallelPlan((SolverGroup('field', tuple(range(field_ranks)), role='field'),
                                 SolverGroup('sbe', (field_ranks,), point_ids=(10,)),
                                 SolverGroup('tddft', tuple(range(field_ranks+1, world.size)), point_ids=(20,))))
        else:
            from salmonpy.layers import assemble_native_layers
            rank_settings = {key: count for key, count in {'sbe-point': 1, 'si-tddft': args.tddft_ranks, 'hse-point': args.hse_ranks}.items() if key in selected}
            assembly = assemble_native_layers(layout, available,
                grid=grid, dt=dt, steps=args.steps, source=calculation.source, field_grid=tuple(args.field_grid),
                ranks_per_point=rank_settings, rank_budget=args.rank_budget if args.rank_budget is not None else world.size)
            if assembly.record()['required_ranks'] != world.size:
                raise ValueError(f"launch exactly {assembly.record()['required_ranks']} ranks for this layout")
            calculation = assembly.calculation; geometry = calculation.geometry
            points = assembly.points; bindings = assembly.bindings; plan = assembly.parallel_plan
        # Replay checkpoints deliberately reconstruct native hidden histories from zero.
        # Bind to input bytes and the native library, excluding output-directory identities.
        def digest_file(path):
            h = hashlib.sha256()
            with Path(path).open('rb') as stream:
                for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
            return h.hexdigest()
        build = json.loads(args.build.read_text())
        assets = {}
        if 'si-tddft' in selected or 'hse-point' in selected:
            assets = {str(p.relative_to(args.bundle)): digest_file(p)
                      for p in sorted(args.bundle.rglob('*')) if p.is_file()}
        signature = {'input': FDTDMPI(tuple(args.field_grid)).input_text(calculation),
                     'points': {str(k): list(v) for k,v in points.items()},
                     'materials': {str(k): v.key for k,v in bindings.items()},
                     'library': digest_file(build['library']), 'assets': assets,
                     'field_grid': list(args.field_grid), 'tddft_ranks': args.tddft_ranks,
                     'layout': layout, 'steps': args.steps,
                     'sbe_model': hashlib.sha256(b''.join(np.asarray(getattr(sbe, k)).tobytes()
                         for k in ('eigen','momentum','nonlocal_velocity','weights','lattice'))
                         + repr((sbe.electrons,sbe.correction_order,sbe.nonlocal_correction)).encode()).hexdigest()}
        if 'hse-point' in selected:
            signature['fast_hse_backend'] = available['hse-point'].backend.identity
            signature['hse_ranks'] = args.hse_ranks
        restored = None
        if args.resume is not None:
            metadata = json.loads((args.resume/'checkpoint.json').read_text())
            if metadata.get('schema') != 'salmonpy.native-replay/1' or metadata['signature'] != signature:
                raise ValueError('replay checkpoint input/library/topology differs')
            if digest_file(args.resume/'history.npz') != metadata['history_sha256']:
                raise ValueError('replay checkpoint checksum mismatch')
            with np.load(args.resume/'history.npz', allow_pickle=False) as archive:
                restored = {k: archive[k].copy() for k in ('A_over_c','E_at_points','J_charge')}
            if not 0 < metadata['step'] < args.steps or any(v.shape != (metadata['step'],len(points),3) for v in restored.values()):
                raise ValueError('invalid replay checkpoint step/shape')
        history = []
        context = ParallelContext.create(world, plan)
        budget=None if args.memory_budget_gib is None else int(args.memory_budget_gib*2**30)
        memory_enabled=any(value is not None for value in world.allgather(budget))
        memory_reports=[]
        memory_path=output/'memory-advice.json'
        if memory_enabled:
            from salmonpy._memory_advice import initial_advice
            initial_advice(world,context,calculation,args.bundle,tuple(args.field_grid),sbe,budget,memory_reports,memory_path,
                           kind="field" if context.group.role=="field" else bindings[context.group.point_ids[0]].key)
        field = None
        if context.group.role == 'field':
            field = NativeEHField(context, calculation, args.build, output/'field', points=points, grid=tuple(args.field_grid))
        session = ResponseSession(context, bindings, dt=dt)
        for point,binding in bindings.items():
            if binding.key=='hse-point' and np.max(np.abs(session.initial_current[point]))>1e-12:
                raise ValueError('zero-initial FDTD requires negligible initial HSE current (<=1e-12 au); prepare a zero-current initial state')
        coupled = CoupledStepper(session, field)
        coupled.memory_budget_bytes=budget
        coupled.memory_report_path=memory_path
        coupled.memory_reports=memory_reports if memory_enabled else None
        if memory_enabled:
            from salmonpy._memory_advice import live_advice
            live_advice(world,budget,memory_reports,memory_path,
                        children=any(hasattr(state,"_process") for state in session.states.values()))
        start = 0
        if args.direct_resume is not None:
            from salmonpy.examples._direct_checkpoint import checkpoint
            start = checkpoint(world, coupled, args.direct_resume, signature, restore=True, protected=not args.unprotected_direct_restore)
            if world.rank == 0:
                history = (json.loads((args.direct_resume/'frames.json').read_text()) if args.unprotected_direct_restore
                           else coupled.checkpoint_frames)
                for frame in history:
                    for key, value in frame.items():
                        if isinstance(value, list): frame[key] = np.asarray(value)
                if len(history) != start: raise ValueError('checkpoint frame count mismatch')
        stop = args.direct_save_at or args.save_at or args.steps
        if args.direct_resume is not None and stop <= start:
            raise ValueError('direct save/stop must advance beyond restored step')
        if restored is not None and stop <= metadata['step']:
            raise ValueError('save-at must advance beyond the resumed checkpoint')
        for index in range(start, stop):
            frame = coupled.step()
            if frame is not None:
                if restored is not None and index < metadata['step']:
                    for key, values in restored.items():
                        np.testing.assert_array_equal(frame[key], values[index])
                history.append(frame)
        if args.direct_save_at is not None:
            from salmonpy.examples._direct_checkpoint import checkpoint, digest
            destination = output/'direct-checkpoint'
            checkpoint(world, coupled, destination, signature, frames=history)
        communication = session.communication_report()
        if memory_enabled:
            live_advice(world,budget,memory_reports,memory_path,
                        children=any(hasattr(state,"_process") for state in session.states.values()),
                        label='mixed_segment_end')
        session.close()
        if field is not None: field.close()
        context.close()
        if args.direct_save_at is not None:
            if world.rank == 0: print(output/'direct-checkpoint')
            return
        if args.save_at is not None:
            if world.rank == 0:
                checkpoint = output/'checkpoint'
                checkpoint.mkdir()
                np.savez(checkpoint/'history.npz', **{key: np.asarray([h[key] for h in history])
                         for key in ('A_over_c','E_at_points','J_charge')})
                (checkpoint/'checkpoint.json').write_text(json.dumps({
                    'schema': 'salmonpy.native-replay/1', 'step': stop,
                    'signature': signature, 'history_sha256': digest_file(checkpoint/'history.npz')}, indent=2)+'\n')
                print(checkpoint)
            world.Barrier()
            return
        # Native lifetimes are finished; replay material drives in independent serial states.
        if world.rank == 0:
            A = np.asarray([h['A_over_c'] for h in history])
            j = np.asarray([h['J_charge'] for h in history])
            E = np.asarray([h['E_at_points'] for h in history])
            peaks = np.max(abs(j), axis=(0, 2))
            if not np.all(peaks > 0): raise AssertionError('nonzero response required for both electronic materials')
            replay = np.empty_like(j)
            initial = None
            try:
                if any(b.key == 'si-tddft' for b in bindings.values()):
                    initial = TDDFTInitialState.from_salmon(args.bundle, 'rt.inp', executable=(SalmonExecutable(args.runtime, digest_file(args.runtime))
                                                                      if args.runtime.is_file() else load_runtime(args.runtime)),
                                                          max_steps=args.steps, threads=1)
                for column, point in enumerate(points):
                    is_sbe = bindings[point].key == 'sbe-point'
                    is_hse = bindings[point].key == 'hse-point'
                    if is_hse:
                        from types import SimpleNamespace
                        from salmonpy.material_response import SolverContext
                        local = SolverContext(SimpleNamespace(rank=0,size=1),'serial-reference',(point,))
                        binding = bindings[point]
                        from salmonpy.fast_hse import FastHSEResponseBackend
                        reference = FastHSEResponseBackend(args.fast_hse_build, output)
                        state = reference.create_state(binding.provider.prepare(local),point,dt,local)
                    else: state = sbe.new_state(dt=dt) if is_sbe else initial.new_state(dt=dt)
                    try:
                        old = np.zeros(3)
                        for n in range(args.steps):
                            replay[n, column] = state.advance(old,A[n,column]) if is_hse else state.step(A[n,column])
                            old = A[n,column].copy()
                    finally: state.close()
                    if is_sbe: np.testing.assert_array_equal(j[:, column], replay[:, column])
                    else: np.testing.assert_allclose(j[:, column], replay[:, column], rtol=2e-10, atol=2e-12)
            finally:
                if initial is not None: initial.close()
            previous = np.zeros((len(points), 3)); old_A = np.zeros((len(points), 3))
            for n, frame in enumerate(history, 1):
                np.testing.assert_array_equal(frame['J_previous'], previous)
                np.testing.assert_allclose(frame['A_over_c'], old_A-dt*frame['E_at_points'], rtol=0, atol=1e-18)
                if frame['electronic_field_clock'] != (n+.5)*dt: raise AssertionError('electronic clock shifted')
                previous = frame['J_charge']; old_A = frame['A_over_c']
            np.savez(output/'response.npz', point_ids=list(points), A_over_c=A, E=E, j=j, serial_replay=replay)
            ids, counts = np.unique(geometry.material_ids(), return_counts=True)
            report = {'status': 'passed', 'steps': args.steps, 'dt': dt, 'tddft_ranks': args.tddft_ranks if 'si-tddft' in selected else 0, 'field_grid': list(args.field_grid),
                      'point_materials': {str(p): bindings[p].key for p in points},
                      'initial_charge_current': session.initial_current,
                      'outer_mpi_ranks':world.size,
                      'serial_hse_workers':sum(b.key=='hse-point' for b in bindings.values()) if args.fast_hse_mpi_build is None else 0,
                      'hse_internal_mpi':args.fast_hse_mpi_build is not None,
                      'hse_k_ranks':args.hse_ranks,
                      'constitutive_cell_counts': dict(zip(map(str, ids.tolist()), counts.tolist())),
                      'ld_materials': {} if args.without_ld else {'1': 'drude', '2': 'two-lorentz-poles'},
                      'maximum_abs_current': peaks.tolist(), 'maximum_replay_error': np.max(abs(j-replay), axis=(0, 2)).tolist(),
                      'timing': coupled.timing, 'communication': communication,
                      'native_steps_this_job': stop-start, 'restored_step': start,
                      'checkpoint_validation': getattr(coupled, 'checkpoint_validation', None),
                      'scope': 'tiny native E/H coupling fixture; no physical convergence or whole-session restart claim'}
            if assembly is not None:
                report['assembly'] = assembly.record()
                report['ld_materials'] = {str(v): k for k, v in assembly.record()['constitutive_material_ids'].items()}
            if args.compare is not None:
                expected = json.loads((args.compare/'result.json').read_text())
                if expected.get('status') != 'passed': raise ValueError('successful comparison run required')
                for key in ('steps', 'dt', 'point_materials', 'constitutive_cell_counts', 'ld_materials'):
                    if expected.get(key) != report[key]: raise ValueError('comparison configuration differs: '+key)
                differences = {}
                with np.load(args.compare/'response.npz', allow_pickle=False) as reference:
                    np.testing.assert_array_equal(reference['point_ids'], list(points))
                    for key, actual in (('A_over_c', A), ('E', E), ('j', j)):
                        if actual.shape != reference[key].shape or not np.isfinite(reference[key]).all():
                            raise ValueError('invalid comparison history: '+key)
                        np.testing.assert_allclose(actual, reference[key], rtol=2e-10, atol=2e-12)
                        differences[key] = float(np.max(abs(actual-reference[key])))
                report['decomposition_comparison'] = {'maximum_absolute_difference': differences,
                                                       'rtol': 2e-10, 'atol': 2e-12}
            (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
            print(output)
    except BaseException:
        # Native errors are not recoverable: avoid stranding peers in cleanup collectives.
        import traceback
        traceback.print_exc()
        world.Abort(1)


if __name__ == '__main__': main()
