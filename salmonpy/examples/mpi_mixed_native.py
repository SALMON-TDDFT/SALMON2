"""Native E/H with field-local LD, SBE and internally parallel TDDFT.

mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native BUILD RUNTIME GS_BUNDLE NEW_OUTPUT
Use --tddft-ranks 1 with 3 ranks for a decomposition comparison.
Use --field-grid 1 2 1 --tddft-ranks 1 with 4 ranks for distributed fields.
3D connection fixture; --field-shape changes Maxwell workload, not TDDFT precision.
"""
import argparse
import json
import hashlib
import time
import math
from dataclasses import asdict
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


def _field_layout(shape):
    shape=tuple(shape)
    if len(shape)!=3 or any(type(n) is not int or n<40 or n%2 for n in shape):
        raise ValueError('sample field-shape requires three even integers >=40')
    grid=Grid(shape,(2.,2.,2.))
    center=tuple(n//2 for n in shape)
    return grid,{10:(center[0]-1,*center[1:]),20:center}


def _clock():
    try:
        value=time.perf_counter()
        return value if math.isfinite(value) else None
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'runtime', 'bundle', 'output'): parser.add_argument(name, type=Path)
    parser.add_argument('--fast-hse-build', type=Path, help='explicit fast-PR serial worker build record for hse-point')
    parser.add_argument('--fast-hse-mpi-build', type=Path, help='explicit native fast-PR host-MPI build record')
    parser.add_argument('--hse-bundle', type=Path, help='optional separate fast-PR HSE GS bundle; default is BUNDLE')
    parser.add_argument('--hse-symmetry', action='store_true', help='explicit fast-PR symmetry opt-in; native HSE group1/2 with matching mode capability')
    parser.add_argument('--hse-reference-bundle', type=Path, help='explicit equivalent full-k HSE bundle for serial comparison')
    parser.add_argument('--tddft-points-per-group',type=int,default=1,help='explicit shared standard TDDFT points per native group')
    parser.add_argument('--hse-points-per-group',type=int,default=1,help='explicit shared same-model HSE points per native group')
    parser.add_argument('--hse-ranks', type=int, choices=(1, 2), default=1)
    parser.add_argument('--dt', type=float, default=.003)
    parser.add_argument('--tddft-ranks', type=int, choices=(1, 2), default=2)
    parser.add_argument('--tddft-materials',type=Path,help='offline native-tddft-materials/1 JSON: per-ID bundle, input_file, decomposition, points_per_group')
    parser.add_argument('--layers', type=Path, help='native-layers/1 JSON referencing sbe-point, si-tddft, hse-point, drude, lorentz')
    parser.add_argument('--run-id', default='mixed-example', help='explicit experiment/run identity; fixture default mixed-example')
    parser.add_argument('--legacy-checkpoint-identity', action='store_true', help='explicit unidentified session for old checkpoint readers; cannot bypass identified archives')
    parser.add_argument('--memory-budget-gib', type=float, help='aggregate advisory memory budget across MPI ranks; warning only')
    parser.add_argument('--rank-budget', type=int, help='advisory planning budget; defaults to the launched world size')
    parser.add_argument('--field-shape', type=int, nargs=3, default=(40,40,40), help='even Maxwell grid counts >=40; electronic GS grid is unchanged')
    parser.add_argument('--field-grid', type=int, nargs=3, default=(1, 1, 1), metavar=('PX', 'PY', 'PZ'))
    parser.add_argument('--compare', type=Path, help='compare against an existing same-physics run with another MPI decomposition')
    parser.add_argument('--steps', type=int, default=128)
    parser.add_argument('--protect-steps', action='store_true', help='explicit whole-step numerical/ordinary-output rollback; requires current native transaction ABI')
    parser.add_argument('--without-ld', action='store_true', help='vacuum control for local LD regions')
    parser.add_argument('--save-at', type=int, help='stop after this step and save a replay checkpoint')
    parser.add_argument('--resume', type=Path, help='reconstruct and continue a replay checkpoint')
    parser.add_argument('--direct-save-at', type=int, help='save internal states and stop at this step')
    parser.add_argument('--direct-resume', type=Path, help='protected internal-state restore without replay; candidate-capable runtime required')
    parser.add_argument('--unprotected-direct-restore', action='store_true', help='explicit legacy restore with no failure-preservation guarantee')
    args = parser.parse_args()
    try:
        grid,points=_field_layout(args.field_shape)
    except ValueError as exc:
        parser.error(str(exc))
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
    if args.hse_bundle is not None and args.fast_hse_build is None:
        parser.error('--hse-bundle requires --fast-hse-build')
    if args.hse_symmetry and (args.fast_hse_mpi_build is None or args.hse_reference_bundle is None):
        parser.error('--hse-symmetry requires native HSE and --hse-reference-bundle')
    if args.hse_reference_bundle is not None and not args.hse_symmetry:
        parser.error('--hse-reference-bundle requires --hse-symmetry')
    if args.tddft_materials is not None and args.layers is None:parser.error('tddft-materials requires --layers')
    if args.tddft_points_per_group<1:parser.error('positive tddft-points-per-group required')
    if args.tddft_points_per_group>1 and args.layers is None:parser.error('tddft grouping requires --layers')
    if args.hse_points_per_group<1 or (args.hse_points_per_group>1 and args.fast_hse_mpi_build is None):
        parser.error('positive hse-points-per-group; grouping requires fast-hse-mpi-build')
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
        if args.hse_symmetry and 'hse-point' not in selected:
            raise ValueError('--hse-symmetry requires an hse-point layer')
        if args.hse_bundle is not None and 'hse-point' not in selected:
            raise ValueError('--hse-bundle requires an hse-point layer')
        available = {'drude': materials[0], 'lorentz': materials[1]}
        if 'sbe-point' in selected:
            available['sbe-point'] = bind_material(sbe, key='sbe-point', identity='synthetic-two-band-v1')
        from salmonpy.examples._tddft_materials import MaterialInput,load_materials
        standard_inputs={'si-tddft':MaterialInput(args.bundle.resolve(),'rt.inp',
                         (1,args.tddft_ranks,1,1,1),args.tddft_points_per_group)}
        if args.tddft_materials is not None:standard_inputs.update(load_materials(args.tddft_materials))
        standard_inputs={key:value for key,value in standard_inputs.items() if key in selected}
        for key,spec in standard_inputs.items():
            available[key]=MaterialBinding(key,
                NativeTDDFTMaterialProvider(spec.bundle,spec.input_file,max_steps=args.steps),
                NativeTDDFTBackend(args.build,output,decomposition=spec.decomposition,
                                   shared_model=spec.points_per_group>1))
        if 'hse-point' in selected:
            if args.fast_hse_build is None: raise ValueError('hse-point requires --fast-hse-build')
            from salmonpy.fast_hse import FastHSEMaterialProvider, FastHSEResponseBackend
            from salmonpy.fast_hse_mpi import NativeFastHSEBackend
            available['hse-point'] = MaterialBinding('hse-point',
                FastHSEMaterialProvider(args.hse_bundle or args.bundle, 'rt.inp', max_steps=args.steps),
                (FastHSEResponseBackend(args.fast_hse_build, output) if args.fast_hse_mpi_build is None
                 else NativeFastHSEBackend(args.fast_hse_mpi_build, output, allow_symmetry=args.hse_symmetry, shared_model=args.hse_points_per_group>1)))
        assembly = None
        if layout is None:
            bindings = {10: available['sbe-point'], 20: available['si-tddft']}
            plan = ParallelPlan((SolverGroup('field', tuple(range(field_ranks)), role='field'),
                                 SolverGroup('sbe', (field_ranks,), point_ids=(10,)),
                                 SolverGroup('tddft', tuple(range(field_ranks+1, world.size)), point_ids=(20,))))
        else:
            from salmonpy.layers import assemble_native_layers
            rank_settings = {key: count for key, count in {'sbe-point': 1, 'hse-point': args.hse_ranks,
                             **{key:spec.ranks for key,spec in standard_inputs.items()}}.items() if key in selected}
            assembly = assemble_native_layers(layout, available,
                grid=grid, dt=dt, steps=args.steps, source=calculation.source, field_grid=tuple(args.field_grid),
                ranks_per_point=rank_settings, rank_budget=args.rank_budget if args.rank_budget is not None else world.size,
                points_per_group=({**{key:spec.points_per_group for key,spec in standard_inputs.items() if spec.points_per_group>1},
                                   **({'hse-point':args.hse_points_per_group} if args.hse_points_per_group>1 else {})} or None))
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
        if (args.tddft_materials is None and 'si-tddft' in selected) or ('hse-point' in selected and args.hse_bundle is None):
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
        if args.tddft_points_per_group>1:signature['tddft_points_per_group']=args.tddft_points_per_group
        if args.tddft_materials is not None:
            signature.pop('tddft_ranks',None)
            signature.pop('tddft_points_per_group',None)
            signature['standard_material_inputs']={key:{'input_file':spec.input_file,
                'decomposition':list(spec.decomposition),'points_per_group':spec.points_per_group,
                'assets':{str(p.relative_to(spec.bundle)):digest_file(p)
                    for p in sorted(spec.bundle.rglob('*')) if p.is_file()}}
                for key,spec in standard_inputs.items()}
        if 'hse-point' in selected:
            signature['fast_hse_backend'] = available['hse-point'].backend.identity
            signature['hse_ranks'] = args.hse_ranks
            if args.hse_points_per_group>1:signature['hse_points_per_group']=args.hse_points_per_group
            if args.hse_symmetry:
                signature['hse_symmetry'] = True
                signature['hse_reference_assets'] = {str(p.relative_to(args.hse_reference_bundle)): digest_file(p)
                    for p in sorted(args.hse_reference_bundle.rglob('*')) if p.is_file()}
            if args.hse_bundle is not None:
                signature['hse_assets'] = {str(p.relative_to(args.hse_bundle)): digest_file(p)
                    for p in sorted(args.hse_bundle.rglob('*')) if p.is_file()}
        restored = None
        if args.resume is not None:
            metadata = json.loads((args.resume/'checkpoint.json').read_text())
            if metadata.get('schema') not in ('salmonpy.native-replay/1','salmonpy.native-replay/2') or metadata['signature'] != signature:
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
            material_key=bindings[context.group.point_ids[0]].key if context.group.point_ids else None
            advice_bundle=(standard_inputs[material_key].bundle if material_key in standard_inputs else
                           (args.hse_bundle or args.bundle) if material_key=='hse-point' else args.bundle)
            initial_advice(world,context,calculation,advice_bundle,
                    tuple(args.field_grid),sbe,budget,memory_reports,memory_path,
                           kind="field" if context.group.role=="field" else bindings[context.group.point_ids[0]].key,
                           protected_steps=args.protect_steps)
        field = None
        if context.group.role == 'field':
            field = NativeEHField(context, calculation, args.build, output/'field', points=points, grid=tuple(args.field_grid))
        point_record=point_map=constitutive_identity=None
        if not args.legacy_checkpoint_identity:
            from salmonpy.examples._mixed_identity import sample_records
            ld_records={}
            if assembly is None:
                ld_definitions={} if args.without_ld else {mid:asdict(model.config) for mid,model in zip(('drude','lorentz'),materials)}
                cell_ids={} if args.without_ld else {'1':'drude','2':'lorentz'}
            else:
                ar=assembly.record()
                ld_definitions={mid:entry['parameters'] for mid,entry in ar['materials'].items() if entry['kind']=='ld'}
                ld_records={mid:entry.get('material_record') for mid,entry in ar['materials'].items() if entry['kind']=='ld'}
                cell_ids={str(number):mid for mid,number in ar['constitutive_material_ids'].items()}
            point_record,point_map,constitutive_identity=sample_records(args.run_id,signature,points,bindings,
                constitutive=ld_definitions,cell_ids=cell_ids,assembly=assembly,constitutive_records=ld_records)
        session = ResponseSession(context, bindings, dt=dt,point_record=point_record,point_map=point_map)
        for point,binding in bindings.items():
            if binding.key=='hse-point' and np.max(np.abs(session.initial_current[point]))>1e-12:
                raise ValueError('zero-initial FDTD requires negligible initial HSE current (<=1e-12 au); prepare a zero-current initial state')
        coupled = CoupledStepper(session, field)
        if constitutive_identity is not None:
            coupled.constitutive_identity=constitutive_identity
        from salmonpy.examples._mixed_identity import checkpoint_identity
        identity=checkpoint_identity(coupled)
        if args.resume is not None:
            from salmonpy.examples._mixed_identity import check_identity
            check_identity(identity,metadata.get('identity'),metadata['schema'])
        error=None
        if world.rank==0:
            try:(output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
            except Exception as exc:error=str(exc)
        session.agree(error)
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
        step_targets=None
        if args.protect_steps:
            from salmonpy.native_steps import native_step_targets
            step_targets=native_step_targets(coupled)
        loop_start=_clock()
        for index in range(start, stop):
            frame = coupled.step_protected(step_targets) if step_targets is not None else coupled.step()
            if frame is not None:
                if restored is not None and index < metadata['step']:
                    for key, values in restored.items():
                        np.testing.assert_array_equal(frame[key], values[index])
                history.append(frame)
        loop_end=_clock()
        elapsed=None if loop_start is None or loop_end is None else loop_end-loop_start
        loop_seconds=world.allgather(elapsed)
        if args.direct_save_at is not None:
            from salmonpy.examples._direct_checkpoint import checkpoint, digest
            destination = output/'direct-checkpoint'
            checkpoint(world, coupled, destination, signature, frames=history)
        if identity is not None and args.direct_save_at is None and args.save_at is None:
            arrays=destination=axes=info=None
            if world.rank==0:
                destination=output/'point-response'
                arrays={key:np.asarray([frame[source] for frame in history])
                        for key,source in (('A_over_c','A_over_c'),('E','E_at_points'),('j','J_charge'))}
                axes={key:1 for key in arrays}
                info={'mixed_identity':identity,'units':{'A_over_c':'A/c:au','E':'E:au','j':'charge-current-density:au'}}
            session.save_arrays(destination,arrays,point_axes=axes,metadata=info)
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
                    'schema': ('salmonpy.native-replay/2' if identity is not None else 'salmonpy.native-replay/1'), 'identity': identity, 'step': stop,
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
            initials = {}
            try:
                for key,spec in standard_inputs.items():
                    initials[key]=TDDFTInitialState.from_salmon(spec.bundle,spec.input_file,
                        executable=(SalmonExecutable(args.runtime,digest_file(args.runtime))
                            if args.runtime.is_file() else load_runtime(args.runtime)),
                        max_steps=args.steps,threads=1)
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
                        if args.hse_reference_bundle is not None:
                            from salmonpy.fast_hse import FastHSEMaterialProvider
                            reference_provider = FastHSEMaterialProvider(args.hse_reference_bundle, max_steps=args.steps)
                        else:
                            reference_provider = binding.provider
                        state = reference.create_state(reference_provider.prepare(local),point,dt,local)
                    else: state = sbe.new_state(dt=dt) if is_sbe else initials[bindings[point].key].new_state(dt=dt)
                    try:
                        old = np.zeros(3)
                        for n in range(args.steps):
                            replay[n, column] = state.advance(old,A[n,column]) if is_hse else state.step(A[n,column])
                            old = A[n,column].copy()
                    finally: state.close()
                    if is_sbe: np.testing.assert_array_equal(j[:, column], replay[:, column])
                    else: np.testing.assert_allclose(j[:, column], replay[:, column], rtol=2e-10, atol=2e-12)
            finally:
                for initial in initials.values():initial.close()
            previous = np.zeros((len(points), 3)); old_A = np.zeros((len(points), 3))
            for n, frame in enumerate(history, 1):
                np.testing.assert_array_equal(frame['J_previous'], previous)
                np.testing.assert_allclose(frame['A_over_c'], old_A-dt*frame['E_at_points'], rtol=0, atol=1e-18)
                if frame['electronic_field_clock'] != (n+.5)*dt: raise AssertionError('electronic clock shifted')
                previous = frame['J_charge']; old_A = frame['A_over_c']
            np.savez(output/'response.npz', point_ids=list(points), A_over_c=A, E=E, j=j, serial_replay=replay)
            ids, counts = np.unique(geometry.material_ids(), return_counts=True)
            report = {'status': 'passed', 'protected_steps': args.protect_steps, 'steps': args.steps, 'dt': dt, 'tddft_ranks': standard_inputs['si-tddft'].ranks if 'si-tddft' in standard_inputs else 0, 'field_grid': list(args.field_grid),
                      'point_materials': {str(p): bindings[p].key for p in points},
                      'initial_charge_current': session.initial_current,
                      'outer_mpi_ranks':world.size,
                      'field_shape':list(grid.counts),
                      'coupled_loop_seconds_by_rank':loop_seconds,
                      'coupled_loop_seconds_max':None if any(v is None for v in loop_seconds) else max(loop_seconds),
                      'serial_hse_workers':sum(b.key=='hse-point' for b in bindings.values()) if args.fast_hse_mpi_build is None else 0,
                      'hse_internal_mpi':args.fast_hse_mpi_build is not None,
                      'hse_k_ranks':args.hse_ranks,
                      'hse_symmetry':args.hse_symmetry,
                      'constitutive_cell_counts': dict(zip(map(str, ids.tolist()), counts.tolist())),
                      'ld_materials': {} if args.without_ld else {'1': 'drude', '2': 'two-lorentz-poles'},
                      'maximum_abs_current': peaks.tolist(), 'maximum_replay_error': np.max(abs(j-replay), axis=(0, 2)).tolist(),
                      'timing': coupled.timing, 'communication': communication,
                      'native_steps_this_job': stop-start, 'restored_step': start,
                      'checkpoint_validation': getattr(coupled, 'checkpoint_validation', None),
                      'identity':identity,
                      'scope': 'tiny native E/H fixture; fixed-binary/input/topology restart; no physical convergence claim'}
            if args.tddft_materials is not None:
                report['standard_material_inputs'] = signature['standard_material_inputs']
            if assembly is not None:
                report['assembly'] = assembly.record()
                report['ld_materials'] = {str(v): k for k, v in assembly.record()['constitutive_material_ids'].items()}
            if args.compare is not None:
                expected = json.loads((args.compare/'result.json').read_text())
                if expected.get('status') != 'passed': raise ValueError('successful comparison run required')
                for key in ('steps', 'dt', 'point_materials', 'constitutive_cell_counts', 'ld_materials'):
                    if expected.get(key) != report[key]: raise ValueError('comparison configuration differs: '+key)
                if expected.get('field_shape',[40,40,40]) != report['field_shape']: raise ValueError('comparison field shape differs')
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
        # Unhandled errors remain fatal; --protect-steps handles only declared healthy-group rollback.
        import traceback
        traceback.print_exc()
        world.Abort(1)


if __name__ == '__main__': main()
