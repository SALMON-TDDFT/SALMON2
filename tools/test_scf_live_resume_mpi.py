"""Three MPI groups: paused/restored SCF vs continuous new/previous runtimes.

BUILD BASELINE NEW_OUTPUT --size2. No GS checkpoint/density/orbital files collected;
all numerical comparison arrays remain in RAM. Only scalar evidence is written.
"""
import argparse,json,warnings
import ctypes as c
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
from mpi4py import MPI
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend
from salmonpy.errors import InputError


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','baseline','output'):p.add_argument(name,type=Path)
    p.add_argument('--size',type=int,default=2)
    p.add_argument('--split',choices=('orbital','spatial','k'),default='orbital')
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','hse06'),default='pz')
    p.add_argument('--initial-no-mix',type=int,default=3)
    p.add_argument('--current-steps',type=int,default=0)
    p.add_argument('--pre-scf-threshold',type=float,default=0.)
    p.add_argument('--pre-scf-steps',type=int,default=4)
    p.add_argument('--gamma',action='store_true')
    p.add_argument('--tolerance',type=float,default=1e-30)
    args=p.parse_args();world=MPI.COMM_WORLD
    if world.size!=3*args.size:raise ValueError('launch exactly3*size ranks')
    g=world.rank//args.size;comm=world.Split(g,world.rank)
    work=args.output/str(g)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    dec=(1,args.size,1,1,1) if args.split=='orbital' else (1,1,args.size,1,1) if args.split=='spatial' else (args.size,1,1,1,1)
    pp=Path(__file__).resolve().parents[1]/'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pp.exists():pp=Path(__file__).resolve().parents[1].parent/'examples/si_intro_data/Si_rps.dat'
    pseudo=Pseudopotential.from_file(pp,atomic_number=14,local_channel=2)
    class InMemoryDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            text=super().input_text(scf).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
            if self.xc=='hse06':
                text=text.replace('exx_pre_scf_threshold=0',f'exx_pre_scf_threshold={args.pre_scf_threshold},exx_pre_scf_steps={args.pre_scf_steps}')
            return text.replace('&scf',f'&scf\n nscf_init_mix_zero={args.initial_no_mix}',1)
    a=InMemoryDFT((14,),np.diag([8.]*3),[[.3,.4,.45]],(pseudo,),4,2,(8,8,8),((1,1,1) if args.gamma else (2,1,1)),xc=args.xc)
    b=InMemoryDFT((14,14),np.diag([9.]*3),[[.3,.4,.45],[.7,.6,.55]],(pseudo,),8,4,(12,12,12),((1,1,1) if args.gamma else (1,2,1)),xc='hse06' if args.xc=='hse06' else 'libxc_pbe')
    assert "write_gs_restart_data='no'" in a.input_text()
    scf=SCFSettings(max_iterations=8,tolerance=args.tolerance)
    backend=NativeGSBackend(args.baseline if g==2 else args.build,work,decomposition=dec)
    def frame(state,final=False):
        return {q:state.array(q).values for q in (('density','occupations','eigenvalues','orbitals','forces') if final else ('density','occupations','eigenvalues','orbitals'))}
    with warnings.catch_warnings():
      warnings.simplefilter('ignore',RuntimeWarning)
      with backend.open(ctx) as session:
        if g==0:
            state=session.run(a,scf=scf,material_id='A',steps=3)
            assert state.scf_status['paused'] and state.scf_status['initializations']==1
            assert state.scf_status['executed_iterations']==3
            cut3=frame(state);snap3=state.snapshot()
            sibling=session.run(b,scf=scf,material_id='B',steps=2)
            sibling_cut=frame(sibling);sibling_snapshot=sibling.snapshot()
            # Native preflight must reject bad/mixed tokens without any mutation.
            before=state.scf_status
            bad_token=snap3._token+comm.rank if comm.size>1 else sibling_snapshot._token
            code=session.lib.sp_gs_scf_restore(state._token,bad_token)
            assert comm.allgather(code)==[11]*comm.size
            assert state.scf_status==before
            if comm.size>1:
                rejected=c.c_int64()
                code=session._gs_call(state.directory,session.lib.sp_gs_start,b'gs.inp',3+comm.rank,c.byref(rejected))
                assert comm.allgather(code)==[11]*comm.size and rejected.value==0
                assert state.scf_status==before
            try:state.restore(sibling_snapshot)
            except InputError:pass
            else:raise AssertionError('cross-owner snapshot accepted')
            state.continue_scf(steps=2);assert state.iterations==5
            snap5=state.snapshot()
            state.continue_scf(steps=3)
            assert state.iterations==8 and not state.scf_status['paused']
            linear=state.scf_status
            assert linear['executed_iterations']==8 and linear['initializations']==1
            expected=frame(state,True);info=state.info
            for saved,n,remaining in ((snap3,3,5),(snap5,5,3)):
                before=state.scf_status['executed_iterations']
                state.restore(saved)
                assert state.iterations==n and state.scf_status['executed_iterations']==before
                assert state.scf_status['initializations']==1
                if n==3:
                    for q,value in frame(state).items():np.testing.assert_array_equal(value,cut3[q])
                state.continue_scf(steps=remaining)
                assert state.info==info
                for q,value in frame(state,True).items():np.testing.assert_array_equal(value,expected[q])
            for q,value in frame(sibling).items():np.testing.assert_array_equal(value,sibling_cut[q])
            sibling.continue_scf(steps=1)
            sibling_future=frame(sibling);sibling_work=sibling.scf_status['executed_iterations']
            sibling.restore(sibling_snapshot)
            assert sibling.iterations==2 and sibling.scf_status['executed_iterations']==sibling_work
            assert sibling.scf_status['initializations']==1
            sibling.continue_scf(steps=1)
            for q,value in frame(sibling).items():np.testing.assert_array_equal(value,sibling_future[q])
            for q,value in frame(state,True).items():np.testing.assert_array_equal(value,expected[q])
            snap3.close();snap5.close();sibling_snapshot.close();sibling.close()
            assert state.scf_status['snapshots']==0
            for q in ('wfn.bin','occupation.bin'):assert not list(work.rglob(q))
        else:
            state=session.run(a,scf=scf,material_id='A')
            expected=frame(state,True);info=state.info
            linear=None
        values=world.allgather((expected,info))
        for rank in range(args.size):
            for other in (rank+args.size,rank+2*args.size):
                assert values[rank][1]==values[other][1]
                for q,array in values[rank][0].items():np.testing.assert_array_equal(array,values[other][0][q])
        if args.current_steps:
            counters=state.scf_status if g==0 else None
            model=state.tddft_model(dt=.001,max_steps=args.current_steps+1)
            response=model.new_state('continued-gs')
            try:
                current=np.asarray([response.step(1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)]))
                                    for i in range(1,args.current_steps+1)])
                traces=world.allgather(current)
                for rank in range(args.size):
                    for other in (rank+args.size,rank+2*args.size):
                        np.testing.assert_array_equal(traces[rank],traces[other])
                if g==0:assert state.scf_status==counters
            finally:response.close();model.close()
        state.close()
    if world.rank==0:
        result={'status':'pass','size':args.size,'split':args.split,'xc':args.xc,'initial_no_mix':args.initial_no_mix,
                'continuous_arrays_exact':True,
                'continuous_new_and_previous_exact':args.build.resolve()!=args.baseline.resolve(),'snapshot_restore_numerical_arrays_exact':True,
                'mixed_history_continuation_exact':True,'sibling_preserved':True,'restore_executed_iterations':0,
                'linear_initializations':linear['initializations'],'linear_executed_iterations':linear['executed_iterations'],
                'final_info':asdict(info),'GS_datasets_collected':False,
                'pre_scf_threshold':args.pre_scf_threshold,'pre_scf_steps':args.pre_scf_steps,'gamma':args.gamma,
                'continued_GS_to_RT_steps':args.current_steps,
                'continued_GS_to_RT_exact':bool(args.current_steps)}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback;traceback.print_exc();MPI.COMM_WORLD.Abort(1)
