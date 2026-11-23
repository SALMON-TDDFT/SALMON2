"""Ni SOC handled-exception rollback: current/energy/spin streams and sibling isolation.

Launch SIZE ranks with an audited runtime. GS stays in memory; only normal
logs/streams and scalar evidence are retained. Unconverged GS is software-only.
"""
import argparse
from dataclasses import fields,replace
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



def check_pauli_zero(state):
    """Pauli identity-channel current must equal native matter current Jm."""
    from salmonpy.outputs import read_spin_dynamics
    from salmonpy.analysis import read_table
    spin=read_spin_dynamics(state.work/next(n for n in state.output_files if n.endswith('_rt_spin.data')))
    current=read_table(state.work/next(n for n in state.output_files if n.endswith('_rt.data')))
    times,si,ci=np.intersect1d(spin.time,current.column('Time'),return_indices=True)
    if len(times)!=len(current.values):raise AssertionError('missing matching spin/current times')
    values=np.column_stack([current.column('Jm_'+axis) for axis in ('x','y','z')])[ci]
    np.testing.assert_allclose(spin.pauli_current[si,0,:],values,rtol=2e-12,atol=2e-14,equal_nan=False)
    return float(np.max(np.abs(spin.pauli_current[si,0,:]-values)))


def close_output_probe(gs,comm):
    """Read completed spin streams while the calling program remains attached."""
    from salmonpy.outputs import read_spin_dynamics
    with gs.tddft_model(dt=.001,max_steps=4) as model:
        first=model.new_state('close-first');sibling=model.new_state('close-sibling')
        traces=[]
        try:
            for i in range(1,5):traces.append(first.step(np.array([1e-5*i,0.,0.])).tolist())
            first.close();error=None;before={}
            if comm.rank==0:
                try:
                    for name in first.output_files:before[name]=(first.work/name).read_bytes()
                    name=next(n for n in first.output_files if n.endswith('_rt_spin.data'))
                    spin=read_spin_dynamics(first.work/name)
                    if len(spin.time)!=5:raise AssertionError('missing final spin observations after close')
                    check_pauli_zero(first)
                except Exception as exc:error=type(exc).__name__+': '+str(exc)
            errors=comm.allgather(error)
            if any(errors):raise AssertionError('SOC close output probe: '+str(errors))
            for i in range(1,5):sibling.step(np.array([-1e-5*i,0.,0.]))
            sibling.close();error=None
            if comm.rank==0:
                try:
                    if any((first.work/n).read_bytes()!=data for n,data in before.items()):
                        raise AssertionError('closed sibling output changed')
                    name=next(n for n in sibling.output_files if n.endswith('_rt_spin.data'))
                    if len(read_spin_dynamics(sibling.work/name).time)!=5:raise AssertionError('sibling final spin observations missing')
                    check_pauli_zero(sibling)
                except Exception as exc:error=type(exc).__name__+': '+str(exc)
            errors=comm.allgather(error)
            if any(errors):raise AssertionError('SOC sibling close probe: '+str(errors))
            return {'status':'passed','spin_rows_per_state':5,'closed_sibling_files_unchanged':True,
                    'streams_read_while_host_attached':True,'explicit_flush_or_transaction_used':False,
                    'GS_exported':False,'current_trace':traces}
        finally:
            for state in (first,sibling):
                if not state.closed:state.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--potential',type=Path,required=True)
    p.add_argument('--size',type=int,default=1)
    p.add_argument('--execution-mode',choices=('checked','performance'),default='checked')
    p.add_argument('--decomposition',type=int,nargs=5,metavar=('K','ORBITAL','X','Y','Z'),help='explicit internal MPI split; product must equal size')
    p.add_argument('--expect-spatial-rejection',action='store_true',help='verify unsupported SOC spatial TDDFT is rejected while MPI and retained GS survive')
    p.add_argument('--conventional-cell',action='store_true',help='cubic fcc Ni: 4 atoms, 72 electrons, 144 spinors for spatial MPI tests')
    p.add_argument('--numerical-overflow',action='store_true',help='checked production build: recover from an overflowing finite A')
    p.add_argument('--close-outputs',action='store_true',help='normal close must finalize SOC output without transaction or explicit flush')
    p.add_argument('--deep-faults',action='store_true',help='requires explicit SALMONPY_RECOVERY_TESTING build')
    args=p.parse_args();comm=MPI.COMM_WORLD
    if (args.deep_faults or args.numerical_overflow) and args.execution_mode!='checked':p.error('numerical failure probes require checked mode')
    from math import prod
    decomposition=tuple(args.decomposition) if args.decomposition is not None else (args.size,1,1,1,1)
    if min(decomposition)<1 or prod(decomposition)!=args.size:p.error('positive MPI decomposition with product equal to size required')
    if comm.size!=args.size:raise ValueError('MPI size mismatch')
    if prod(decomposition[2:])>1 and not args.conventional_cell:p.error('spatial test requires --conventional-cell; SALMON does not support spatial splitting of the primitive nonorthogonal cell')
    if comm.rank==0:args.output.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    pp=Pseudopotential.from_file(args.potential,atomic_number=28,local_channel=2,
                               reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni'))
    base=nickel_model(pp)
    if args.conventional_cell:
        base=replace(base,numbers=[28]*4,cell=np.eye(3)*(2*base.cell[0,1]),
                     fractional_positions=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]],electrons=72,states=144)
    calc=SpinProbeDFT(**{f.name:getattr(base,f.name) for f in fields(DFT)})
    scf=SCFSettings(max_iterations=12,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,
                   initial_wavefunction='random',random_seed=0)
    backend=NativeGSBackend(args.build,args.output,decomposition=decomposition)
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
            if args.expect_spatial_rejection:
                from salmonpy.errors import InputError
                before=gs.scf_status.copy();density=gs.array('density').values.copy();message=None
                try:
                    with gs.tddft_model(dt=.001,max_steps=8):pass
                except InputError as exc:
                    if 'SOC' not in str(exc) or 'spatial MPI' not in str(exc):raise
                    message=str(exc)
                else:raise AssertionError('unsupported spatial SOC TDDFT was not rejected')
                if comm.allreduce(comm.rank+1)!=comm.size*(comm.size+1)//2:raise AssertionError('caller MPI unavailable')
                if before!=gs.scf_status:raise AssertionError('rejection changed GS counters')
                np.testing.assert_array_equal(density,gs.array('density').values)
                if comm.rank==0:
                    result={'status':'passed','expected_rejection':message,'SOC_spatial_TDDFT_supported':False,
                            'decomposition':list(decomposition),'MPI_caller_survives':True,'GS_preserved':True,'GS_exported':False}
                    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
                return
            if args.close_outputs:
                result=close_output_probe(gs,comm)
                if comm.rank==0:(args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
                return
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
                    if args.numerical_overflow:
                        from salmonpy.errors import NumericalError
                        token=targets[0].capture()
                        try:
                            try:retry.step(np.array([1e100,0.,0.]))
                            except NumericalError:pass
                            else:raise AssertionError('overflow did not return a numerical error')
                            if retry.steps!=0 or retry.failed:raise AssertionError('overflow changed state flags')
                            np.testing.assert_array_equal(retry.A,old_A)
                            np.testing.assert_array_equal(retry.current,old_j)
                            if retry.electron_number!=old_n:raise AssertionError('overflow changed electron count')
                            if comm.rank==0 and data(retry)!=before:raise AssertionError('overflow changed output')
                            np.testing.assert_array_equal(retry.step(drive(1)),expected_first)
                            cases.append('production-overflow')
                        finally:
                            targets[0].rollback(token);targets[0].release(token)
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
                            check_pauli_zero(a);check_pauli_zero(b)
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
        result={'status':'passed','MPI_ranks':args.size,'decomposition':list(decomposition),'atoms':len(calc.numbers),'electrons':calc.electrons,'states':calc.states,
                'kmesh':list(map(int,calc.kmesh)),'temperature_k':300,'conventional_cell':args.conventional_cell,'state_preserved':True,
                'all_three_streams_byte_preserved':True,'retry_currents_and_streams_exact':True,
                'sibling_unaffected':True,'Pauli0_matches_matter_current':True,'GS_exported':False,'execution_mode':args.execution_mode,
                'deep_cases':cases,'performance_skips_stage_checks':args.deep_faults,
                'scope':'checked automatic numerical recovery or explicit performance rollback; not crash/MPI failure recovery'}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Ni SOC rollback/retry PASS',args.output/'result.json')


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
