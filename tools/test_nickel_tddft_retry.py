"""Ni SOC handled-exception rollback: current/energy/spin streams and sibling isolation.

Launch SIZE ranks with an audited runtime. GS stays in memory; only normal
logs/streams and scalar evidence are retained. Unconverged GS is software-only.
"""
import argparse
from dataclasses import fields
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend
from salmonpy.pseudopotentials import catalog_reference
from salmonpy.native_steps import _NativeBackup
from salmonpy.stepping import StepTarget,guarded_step,StepRecoveryError
from test_nickel_soc import nickel_model
from test_nickel_scf_mpi import InMemoryDFT


class SpinProbeDFT(InMemoryDFT):
    def input_text(self,scf=SCFSettings()):
        # Explicit observation frequency for this test, not a global default.
        return super().input_text(scf)+"\n&analysis\n out_rt_spin_step=1\n/\n"



def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--potential',type=Path,required=True)
    p.add_argument('--size',type=int,choices=(1,2),default=1)
    p.add_argument('--execution-mode',choices=('checked','performance'),default='checked')
    p.add_argument('--deep-faults',action='store_true',help='requires explicit SALMONPY_RECOVERY_TESTING build')
    args=p.parse_args();comm=MPI.COMM_WORLD
    if args.deep_faults and args.execution_mode!='checked':p.error('deep faults require checked mode')
    if comm.size!=args.size:raise ValueError('MPI size mismatch')
    if comm.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    pp=Pseudopotential.from_file(args.potential,atomic_number=28,local_channel=2,
                               reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni'))
    base=nickel_model(pp);calc=SpinProbeDFT(**{f.name:getattr(base,f.name) for f in fields(DFT)})
    scf=SCFSettings(max_iterations=12,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,
                   initial_wavefunction='random',random_seed=0)
    backend=NativeGSBackend(args.build,args.output,decomposition=(args.size,1,1,1,1))
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    def drive(i):return np.array([1e-5*np.sin(.3*i),7e-6*(np.cos(.2*i)-1),1e-5*np.sin(.17*i)])
    def target(state):
        b=_NativeBackup(state,comm,'sp_tddft_step_transaction',state.work,
                        ('A','steps','failed','closed','current'))
        return StepTarget('ni-soc',state,b.capture,b.rollback,b.release)
    def data(state):
        return {str(p):p.read_bytes() for p in sorted(state.work.glob('*.data'))
                if str(p.relative_to(state.work)) in state.output_files}
    def rows(path):return [line for line in path.read_text().splitlines() if line.strip() and not line.startswith('#')]
    with backend.open(context) as session:
        with session.run(calc,scf=scf,steps=4,material_id='ni-soc-retry-probe') as gs:
            with gs.tddft_model(dt=.001,max_steps=8) as model:
                control=model.new_state('control');retry=model.new_state('retry',execution_mode=args.execution_mode)
                sibling=model.new_state('sibling');sibling_control=model.new_state('sibling-control')
                try:
                    targets=[target(retry)]
                    guarded_step(lambda:None,targets)
                    before=data(retry) if comm.rank==0 else None
                    old_A=retry.A.copy();old_j=retry.current.copy();old_n=retry.electron_number
                    expected_first=control.step(drive(1))
                    cases=[]
                    if args.deep_faults:
                        import ctypes as c
                        from salmonpy.errors import NumericalError
                        hook=session.lib.sp_fast_hse_test_fault
                        hook.argtypes=[c.c_int,c.c_int];hook.restype=c.c_int
                        fast_probe=model.new_state('unchecked-probe',execution_mode='performance')
                        try:
                            if hook(20,comm.size-1)!=0:raise AssertionError('performance probe setup failed')
                            np.testing.assert_array_equal(fast_probe.step(drive(1)),expected_first)
                            # The unconsumed hook must still fire in the next checked step.
                            try:retry.step(drive(1))
                            except NumericalError:pass
                            else:raise AssertionError('performance unexpectedly consumed a stage check')
                            if retry.steps!=0 or retry.failed:raise AssertionError('performance probe changed checked state')
                            if comm.rank==0 and data(retry)!=before:raise AssertionError('checked probe wrote failed output')
                        finally:fast_probe.close()
                        for stage in (None,20,21,22,23,24):
                            token=targets[0].capture()
                            try:
                                if stage is not None and hook(stage,comm.size-1)!=0:
                                    raise AssertionError('fault setup failed')
                                try:retry.step(np.array([1e100,0.,0.]) if stage is None else drive(1))
                                except NumericalError:pass
                                else:raise AssertionError('deep failure was not returned')
                                if retry.steps!=0 or retry.failed:raise AssertionError('deep failure changed state flags')
                                np.testing.assert_array_equal(retry.A,old_A)
                                np.testing.assert_array_equal(retry.current,old_j)
                                if retry.electron_number!=old_n:raise AssertionError('deep failure changed electron count')
                                if comm.rank==0 and data(retry)!=before:raise AssertionError('deep failure wrote output')
                                np.testing.assert_array_equal(retry.step(drive(1)),expected_first)
                                cases.append('overflow' if stage is None else 'stage-'+str(stage))
                            finally:
                                targets[0].rollback(token)
                                targets[0].release(token)
                    class InjectedFailure(RuntimeError):pass
                    def fail():
                        retry.step(drive(1));sibling.step(-drive(1))
                        raise InjectedFailure('handled host failure after accepted SOC step')
                    try:guarded_step(fail,targets)
                    except StepRecoveryError as exc:
                        if not exc.original_preserved:raise AssertionError('rollback not preserved')
                    else:raise AssertionError('injected failure was ignored')
                    if retry.steps!=0 or retry.failed:raise AssertionError('bad restored state flags')
                    np.testing.assert_array_equal(retry.A,old_A)
                    np.testing.assert_array_equal(retry.current,old_j)
                    if retry.electron_number!=old_n:raise AssertionError('electron number changed')
                    if comm.rank==0 and data(retry)!=before:raise AssertionError('SOC output bytes changed after rollback')
                    np.testing.assert_array_equal(sibling.current,sibling_control.step(-drive(1)))
                    for i in range(1,5):
                        expected=expected_first if i==1 else control.step(drive(i))
                        np.testing.assert_array_equal(guarded_step(lambda:retry.step(drive(i)),targets),expected)
                        if i>1:np.testing.assert_array_equal(sibling.step(-drive(i)),sibling_control.step(-drive(i)))
                    for state in (control,retry,sibling,sibling_control):guarded_step(lambda:None,[target(state)])
                    if comm.rank==0:
                        if not any(n.endswith('_rt_spin.data') for n in retry.output_files):
                            raise AssertionError('spin stream is not owned by rollback')
                        from salmonpy.outputs import read_spin_dynamics
                        for a,b in ((retry,control),(sibling,sibling_control)):
                            for x,y in zip(a.output_files,b.output_files):
                                if rows(a.work/x)!=rows(b.work/y):raise AssertionError('stream differs: '+x)
                                if x.endswith('_rt_spin.data'):
                                    sa=read_spin_dynamics(a.work/x);sb=read_spin_dynamics(b.work/y)
                                    np.testing.assert_array_equal(sa.magnetization,sb.magnetization)
                                    np.testing.assert_array_equal(sa.spin_current,sb.spin_current)
                                    if sa.magnetization.shape!=(5,3):raise AssertionError('missing per-step spin observations')
                finally:
                    for state in (sibling_control,sibling,retry,control):state.close()
    if comm.rank==0:
        result={'status':'passed','MPI_ranks':args.size,'state_preserved':True,
                'all_three_streams_byte_preserved':True,'retry_currents_and_streams_exact':True,
                'sibling_unaffected':True,'GS_exported':False,'execution_mode':args.execution_mode,
                'deep_cases':cases,'performance_skips_stage_checks':args.deep_faults,
                'scope':'checked automatic numerical recovery or explicit performance rollback; not crash/MPI failure recovery'}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Ni SOC rollback/retry PASS',args.output/'result.json')


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
