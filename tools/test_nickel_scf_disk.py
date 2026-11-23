"""Explicit Ni SOC fresh-process SCF restart test; saves electronic checkpoint data.

BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF [--size 2]
Requires explicit authorization to retain orbital/density/checkpoint data.
No download, build or RT continuation is performed.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--potential', type=Path, required=True)
    parser.add_argument('--size', type=int, choices=(1, 2), default=1)
    parser.add_argument('--restore-build',type=Path)
    parser.add_argument('--extend-limit',type=int)
    parser.add_argument('--check-mixing',action='store_true')
    parser.add_argument('--adaptive-mixing',action='store_true')
    parser.add_argument('--spin-field',action='store_true',help='validate auxiliary spin-field stages, independence and direct restart')
    parser.add_argument('--mixing-method',choices=('simple_dm','broyden'),default='simple_dm')
    parser.add_argument('--compatible-runtime',choices=('declared',),default=False)
    parser.add_argument('--worker', choices=('continuous', 'save', 'resume'))
    args = parser.parse_args()
    if args.worker is None:
        args.output.mkdir(parents=True, exist_ok=False)
        for phase in ('continuous', 'save', 'resume'):
            command = ['mpiexec', '--bind-to', 'none', '-n', str(args.size),
                       sys.executable, str(Path(__file__).resolve()),
                       str((args.restore_build if phase=='resume' and args.restore_build else args.build).resolve()), str(args.output.resolve()),
                       '--potential', str(args.potential.resolve()),
                       '--size', str(args.size), '--worker', phase]
            command+=['--mixing-method',args.mixing_method]
            if args.check_mixing:command+=['--check-mixing']
            if args.adaptive_mixing:command+=['--adaptive-mixing']
            if args.spin_field:command+=['--spin-field']
            if args.extend_limit is not None:command+=['--extend-limit',str(args.extend_limit)]
            if phase=='resume' and args.restore_build:command+=['--compatible-runtime','declared']
            with (args.output / (phase + '.log')).open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT,
                               env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
        result = {'status': 'passed', 'solver_ranks': args.size, 'fresh_processes': 3,
                  'electrons': 18, 'states': 36, 'electronic_temperature_k': 300,
                  'kmesh': [2, 2, 2], 'restored_iterations': 4 if args.spin_field else 2,
                  'spin_field_stages':[[.001,-.0005,.002],[.0005,-.00025,.001],[0.,0.,0.]] if args.spin_field else None,
                  'restore_executed_iterations': 0, 'arrays_exact': True,
                  'adaptive_mixing':args.adaptive_mixing,'extended_limit':args.extend_limit,'source_transition':bool(args.restore_build),
                  'scope': 'short unconverged SOC SCF; not physical convergence or TDDFT'}
        (args.output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        print('Ni SOC fresh-process SCF restart PASS', args.output / 'result.json')
        return

    import numpy as np
    from mpi4py import MPI
    from types import SimpleNamespace
    from salmonpy.dft import Pseudopotential, SCFSettings
    from salmonpy.native_gs import NativeGSBackend
    from salmonpy.pseudopotentials import catalog_reference
    from test_nickel_soc import nickel_model
    from test_nickel_scf_mpi import InMemoryDFT, frame, compare
    from dataclasses import fields

    comm = MPI.COMM_WORLD
    if comm.size != args.size:
        raise ValueError('MPI size mismatch')
    work = args.output / args.worker
    if comm.rank == 0:
        work.mkdir(exist_ok=False)
    comm.Barrier()
    pseudo = Pseudopotential.from_file(args.potential, atomic_number=28, local_channel=2,
                                      reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni'))
    base = nickel_model(pseudo)
    model = InMemoryDFT(**{f.name: getattr(base, f.name) for f in fields(type(base))})
    scf = SCFSettings(max_iterations=12, tolerance=1e-12, mixing_method=args.mixing_method,
                      mixing_rate=.2, adaptive_mixing=args.adaptive_mixing, initial_wavefunction='random', random_seed=0)
    backend = NativeGSBackend(args.build, work, decomposition=(args.size, 1, 1, 1, 1))
    context = SimpleNamespace(solver_comm=comm, closed=False, _dependents=0)
    checkpoint = args.output / 'checkpoint'
    rankfile = f'rank-{comm.rank:08d}.npz'
    def final_info(state):
        from dataclasses import asdict
        info=asdict(state.info)
        if args.worker=='resume':
            if info!=json.loads((args.output/'continuous-info.json').read_text()):
                raise AssertionError('resumed scalar energy/residual/status mismatch')
        elif args.worker=='continuous' and comm.rank==0:
            (args.output/'continuous-info.json').write_text(json.dumps(info))
    with backend.open(context) as session:
        with session.run(model, scf=scf, material_id='ni-soc',
                         steps=None if args.worker == 'resume' else 2 if args.worker == 'save' or args.spin_field else 4,
                         restart=checkpoint if args.worker == 'resume' else None,
                         compatible_runtime=args.compatible_runtime) as state:
            if args.spin_field:
                stages=([.001,-.0005,.002],[.0005,-.00025,.001],[0.,0.,0.])
                if args.worker=='resume':
                    if state.iterations!=4 or state.scf_status['executed_iterations']!=0:
                        raise AssertionError('spin-field restore replayed SCF')
                    with np.load(args.output/'save'/rankfile) as expected:compare(dict(expected),frame(state))
                    expected=json.loads((args.output/'spin-field-policy.json').read_text())
                    if state.scf_spin_field!=expected:raise AssertionError('field or spin observables not restored')
                    for field in stages[1:]:
                        state.set_scf_spin_field(field);state.continue_scf(steps=2)
                    with np.load(args.output/'continuous'/rankfile) as expected:compare(dict(expected),frame(state))
                else:
                    boundary=frame(state);counts=state.scf_status
                    state.set_scf_spin_field([0,0,0]);compare(boundary,frame(state))
                    if state.scf_status!=counts:raise AssertionError('zero-field no-op changed status')
                    zero_energy=state.info.energy_hartree
                    state.set_scf_spin_field(stages[0]);compare(boundary,frame(state))
                    spin=state.scf_spin_field
                    if abs(spin['zeeman_energy_hartree']+.5*np.dot(spin['field_au'],spin['spin_magnetic_moment_muB']))>2e-12:
                        raise AssertionError('Zeeman energy/moment sign or units disagree')
                    from salmonpy._gs_tddft import open_model
                    from salmonpy.errors import InputError
                    try:open_model(state,.001,1)
                    except InputError as exc:
                        if 'auxiliary spin field' not in str(exc):raise
                    else:raise AssertionError('TDDFT silently inherited auxiliary field')
                    if abs(state.scf_spin_field['energy_without_auxiliary_hartree']-zero_energy)>2e-12:
                        raise AssertionError('field energy double counted')
                    if comm.size>1:
                        bad=np.asarray(stages[0],dtype=np.float64);bad[0]+=.001*comm.rank
                        fn=session.lib.sp_gs_scf_spin_field
                        if state.session._gs_call(state.directory,fn,state._token,bad.ctypes.data)!=11:
                            raise AssertionError('different rank fields accepted')
                        if state.scf_spin_field['field_au']!=stages[0]:raise AssertionError('rejected field mutated state')
                    saved=state.snapshot()
                    state.set_scf_spin_field(stages[1]);state.restore(saved)
                    from salmonpy.errors import NumericalError
                    stable=frame(state);stable_control=state.scf_spin_field
                    try:
                        extreme=-np.full(3,np.finfo(np.float64).max)
                        state.set_scf_spin_field(extreme)
                    except NumericalError:
                        compare(stable,frame(state))
                        if state.scf_spin_field!=stable_control:raise AssertionError('rejected large field changed state')
                    else:
                        before_extreme=state.scf_spin_field;frozen=frame(state)
                        try:state.set_scf_spin_field(-extreme)
                        except NumericalError:pass
                        else:raise AssertionError('overflow in field delta was not rejected')
                        compare(frozen,frame(state))
                        if state.scf_spin_field!=before_extreme:raise AssertionError('failed field energy changed state')
                    state.restore(saved);saved.close()
                    if state.scf_spin_field['field_au']!=stages[0]:raise AssertionError('memory restore lost field')
                    if args.worker=='continuous':
                        with session.run(model,scf=scf,steps=2,material_id='ni-zero-field-control') as sibling:
                            sibling_start=frame(sibling)
                            state.continue_scf(steps=2)
                            compare(sibling_start,frame(sibling))
                            if any(sibling.scf_spin_field['field_au']):raise AssertionError('field leaked to sibling')
                            sibling.continue_scf(steps=2)
                            if np.array_equal(frame(state)['orbitals'],frame(sibling)['orbitals']):
                                raise AssertionError('field did not affect orbital propagation')
                        for field in stages[1:]:
                            state.set_scf_spin_field(field);state.continue_scf(steps=2)
                        np.savez(work/rankfile,**frame(state))
                    else:
                        state.continue_scf(steps=2);state.save_scf(checkpoint)
                        policy=state.scf_spin_field
                        if comm.rank==0:(args.output/'spin-field-policy.json').write_text(json.dumps(policy))
                        np.savez(work/rankfile,**frame(state))
                if args.worker!='save' and any(state.scf_spin_field['field_au']):
                    raise AssertionError('final field not zero')
                if args.worker!='save':final_info(state)
                return
            if args.worker == 'save':
                state.save_scf(checkpoint)
                if args.adaptive_mixing:
                    policy=state.scf_mixing
                    if comm.rank==0:(args.output/'mixing-policy.json').write_text(json.dumps(policy))
                np.savez(work / rankfile, **frame(state))
            elif args.worker == 'continuous':
                np.savez(work / rankfile, **frame(state))
                final_info(state)
            else:
                if state.iterations != 2 or state.scf_status['executed_iterations'] != 0:
                    raise AssertionError('restart recalculated SCF or loaded wrong iteration')
                with np.load(args.output / 'save' / rankfile) as expected:
                    compare(dict(expected), frame(state))
                if args.adaptive_mixing:
                    if state.scf_mixing!=json.loads((args.output/'mixing-policy.json').read_text()):
                        raise AssertionError('adaptive controller not restored from disk')
                    before=frame(state);policy=state.scf_mixing
                    with state.snapshot() as saved:
                        state.set_scf_adaptive_mixing(min_rate=.02,growth=1.1,patience=2,cooldown=4)
                        compare(before,frame(state))
                        if comm.size>1:
                            fn=session.lib.sp_gs_scf_adaptive_mixing
                            if state.session._gs_call(state.directory,fn,state._token,1,.02,1.1,2+comm.rank,4)!=11:
                                raise AssertionError('different per-rank adaptive policies accepted')
                        state.restore(saved)
                        if state.scf_mixing!=policy:raise AssertionError('snapshot did not restore adaptive policy')
                if args.check_mixing:
                    before=frame(state);counts=state.scf_status
                    with state.snapshot() as saved:
                        counts=state.scf_status
                        state.set_scf_mixing(.1)
                        compare(before,frame(state))
                        if state.scf_status!=counts:raise AssertionError('mixing change ran SCF')
                        if comm.size>1:
                            fn=session.lib.sp_gs_scf_set_mixing
                            if state.session._gs_call(state.directory,fn,state._token,.1+.01*comm.rank)!=11:
                                raise AssertionError('different per-rank mixing rates accepted')
                        state.restore(saved)
                if args.extend_limit is not None:
                    before=frame(state);counts=state.scf_status
                    state.extend_scf_limit(args.extend_limit)
                    compare(before,frame(state))
                    if state.scf_status['executed_iterations']!=counts['executed_iterations']:
                        raise AssertionError('extension recalculated SCF')
                    if comm.size>1:
                        fn=session.lib.sp_gs_scf_extend_limit
                        if state.session._gs_call(state.directory,fn,state._token,args.extend_limit+comm.rank+1)!=11:
                            raise AssertionError('mixed per-rank limits accepted')
                        if state.scf_status['max_iterations']!=args.extend_limit:
                            raise AssertionError('rejected extension changed limit')
                state.continue_scf(steps=2)
                final_info(state)
                with np.load(args.output / 'continuous' / rankfile) as expected:
                    compare(dict(expected), frame(state))
                if args.extend_limit is not None:
                    from dataclasses import replace
                    from salmonpy.errors import SCFReviewRequired
                    with session.run(model,scf=replace(scf,max_iterations=2),steps=1,
                                     material_id='terminal-cap') as capped:
                        try:capped.continue_scf(steps=1)
                        except SCFReviewRequired:pass
                        if capped.iterations!=2:raise AssertionError('terminal cap not reached')
                        boundary=frame(capped);count=capped.scf_status['executed_iterations']
                        capped.extend_scf_limit(args.extend_limit)
                        compare(boundary,frame(capped))
                        if capped.scf_status['executed_iterations']!=count:raise AssertionError('terminal extension ran SCF')
                        capped.continue_scf(steps=2)
                        with np.load(args.output/'continuous'/rankfile) as expected:compare(dict(expected),frame(capped))
                        if state.scf_status['max_iterations']!=args.extend_limit:raise AssertionError('sibling cap was affected')


if __name__ == '__main__':
    main()
