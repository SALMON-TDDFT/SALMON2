"""Three MPI groups: paused/restored SCF vs continuous new/previous runtimes.

BUILD BASELINE NEW_OUTPUT --size2. No GS checkpoint/density/orbital files collected;
all numerical comparison arrays remain in RAM. Only scalar evidence is written.
"""
import argparse,json,warnings,copy
import ctypes as c
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict,replace
import numpy as np
from mpi4py import MPI
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,SymmetryOperations
from salmonpy.native_gs import NativeGSBackend
from salmonpy.errors import InputError


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','baseline','output'):p.add_argument(name,type=Path)
    p.add_argument('--size',type=int,default=2)
    p.add_argument('--split',choices=('orbital','spatial','k','mixed'),default='orbital')
    p.add_argument('--decomposition',type=int,nargs=5)
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','hse06'),default='pz')
    p.add_argument('--initial-no-mix',type=int,default=3)
    p.add_argument('--current-steps',type=int,default=0)
    p.add_argument('--pre-scf-threshold',type=float,default=0.)
    p.add_argument('--pre-scf-steps',type=int,default=4)
    p.add_argument('--gamma',action='store_true')
    p.add_argument('--symmetry',choices=('none','identity','mirror-x','c2z'),default='none')
    p.add_argument('--tolerance',type=float,default=1e-30)
    p.add_argument('--common-api',action='store_true')
    p.add_argument('--protected-api',action='store_true')
    p.add_argument('--full-k-control',action='store_true',help='also compare SCF-converged diamond-Si full/reduced k states')
    args=p.parse_args();args.common_api=args.common_api or args.protected_api;world=MPI.COMM_WORLD
    if world.size!=3*args.size:raise ValueError('launch exactly3*size ranks')
    g=world.rank//args.size;comm=world.Split(g,world.rank)
    work=args.output/str(g)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    dec=(1,args.size,1,1,1) if args.split=='orbital' else (1,1,args.size,1,1) if args.split=='spatial' else (args.size,1,1,1,1)
    if args.decomposition:dec=tuple(args.decomposition)
    pp=Path(__file__).resolve().parents[1]/'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pp.exists():pp=Path(__file__).resolve().parents[1].parent/'examples/si_intro_data/Si_rps.dat'
    pseudo=Pseudopotential.from_file(pp,atomic_number=14,local_channel=2)
    class InMemoryDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            text=super().input_text(scf).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
            if self.xc=='hse06':
                text=text.replace('exx_pre_scf_threshold=0',f'exx_pre_scf_threshold={args.pre_scf_threshold},exx_pre_scf_steps={args.pre_scf_steps}')
            return text.replace('&scf',f'&scf\n nscf_init_mix_zero={args.initial_no_mix}',1)
    symmetry=None
    if args.symmetry!='none':
        if args.gamma:raise ValueError('symmetry control uses a non-Gamma mesh')
        operations=np.zeros((1 if args.symmetry=='identity' else 2,3,4));operations[:,:,:3]=np.eye(3)
        if len(operations)>1:operations[1,0,0]=-1
        if args.symmetry=='c2z':operations[1,1,1]=-1
        symmetry=SymmetryOperations(operations,directions='yyn')
    mesh=(2,2,2) if symmetry is not None else (1,1,1) if args.gamma else (2,1,1)
    a=InMemoryDFT((14,),np.diag([8.]*3),[[0,0,0]] if symmetry is not None else [[.3,.4,.45]],
        (pseudo,),4,2,(8,8,8),mesh,xc=args.xc,symmetry=symmetry)
    b=InMemoryDFT((14,14),np.diag([9.]*3),[[0,0,0],[.5,.5,.5]] if symmetry is not None else [[.3,.4,.45],[.7,.6,.55]],
        (pseudo,),8,4,(12,12,12),mesh if symmetry is not None else (1,1,1) if args.gamma else (1,2,1),
        xc='hse06' if args.xc=='hse06' else 'libxc_pbe',symmetry=symmetry)
    assert "write_gs_restart_data='no'" in a.input_text()
    scf=SCFSettings(max_iterations=8,tolerance=args.tolerance)
    backend=NativeGSBackend(args.baseline if g==2 else args.build,work,decomposition=dec,hse_spatial=args.xc=='hse06' and args.split in ('spatial','mixed'))
    def frame(state,final=False):
        values={}
        for q in (('density','occupations','eigenvalues','orbitals','forces') if final else ('density','occupations','eigenvalues','orbitals')):
            array=state.array(q)
            if q=='forces' and g!=2:assert array.unit=='hartree/bohr'
            values[q]=array.values
        return values
    def full_k_control(session):
        if not args.full_k_control or g!=0:return None
        if args.xc!='hse06':raise ValueError('full-k control requires fast HSE')
        operations=np.zeros((2,3,4));operations[:,:,:3]=np.eye(3)
        operations[1,:,:3]=[[0,1,0],[1,0,0],[0,0,1]]
        model=InMemoryDFT((14,14),np.array([[0,1,1],[1,0,1],[1,1,0]])*5.13,
            [[0,0,0],[.25,.25,.25]],(pseudo,),8,4,(8,8,8),(2,2,2),xc='hse06',
            symmetry=SymmetryOperations(operations,directions='yyn'))
        settings=SCFSettings(max_iterations=128,tolerance=1e-8)
        with session.run(model,scf=settings,material_id='diamond-reduced') as reduced:
          with session.run(replace(model,symmetry=None),scf=settings,material_id='diamond-full') as full:
            first,second=reduced.info,full.info
            assert first.converged and second.converged,'full-k comparison requires SCF convergence (solver results still retained)'
            np.testing.assert_allclose(first.energy_hartree,second.energy_hartree,atol=1e-6,rtol=0)
            np.testing.assert_allclose(reduced.array('density').values,full.array('density').values,atol=1e-7,rtol=1e-4)
            reduced_array=reduced.array('forces');full_array=full.array('forces')
            assert reduced_array.unit==full_array.unit=='hartree/bohr'
            reduced_force=reduced_array.values;full_force=full_array.values
            np.testing.assert_allclose(reduced_force,full_force,atol=1e-6,rtol=0)
            force_agrees=True
            return {'SCF_converged':True,'reduced_kpoints':reduced.array('occupations').values.shape[1],
                'full_kpoints':full.array('occupations').values.shape[1],
                'reduced_info':asdict(first),'full_info':asdict(second),
                'energy_atol_hartree':1e-6,'density_rtol':1e-4,'density_atol_bohr3':1e-7,
                'force_atol_hartree_per_bohr':1e-6,'density_within_tolerance':True,
                'force_within_tolerance':force_agrees,
                'force_max_absolute_difference':float(np.max(abs(reduced_force-full_force))),
                'status':'pass' if force_agrees else 'energy-density-pass-force-mismatch',
                'scope':'same fixed cell/grid/k sampling, not grid/k convergence certification'}
    with warnings.catch_warnings():
      warnings.simplefilter('ignore',RuntimeWarning)
      with backend.open(ctx) as session:
        if g==0:
            state=session.run(a,scf=scf,material_id='A',steps=3)
            assert state.scf_status['paused'] and state.scf_status['initializations']==1
            assert state.scf_status['executed_iterations']==3
            cut3=frame(state);snap3=state.snapshot()
            sibling=session.run(b,scf=scf,material_id='B',steps=2)
            sibling_cut=frame(sibling)
            if args.common_api:
                from salmonpy.checkpointing import execute_checkpoint,restore_checkpoint,CheckpointState
                groups=[{'name':'scf','ranks':list(range(comm.size)),'point_ids':['A','B']}]
                checkpoint=execute_checkpoint(comm,groups,states={
                    'A':state.checkpoint_target(),'B':sibling.checkpoint_target()})
                snap3.close();snap3=checkpoint['states']['A'];sibling_snapshot=checkpoint['states']['B']
                targets={'A':state.checkpoint_target(snap3),'B':sibling.checkpoint_target(sibling_snapshot)}
                counters=state.scf_status,sibling.scf_status
                if targets['A'].stage_restore is None:
                    if args.protected_api:raise AssertionError('SCF isolated candidates missing')
                    try:
                        restore_checkpoint(comm,groups,states=targets,saved=checkpoint['compatibility'],
                            live=CheckpointState(states={'A':state,'B':sibling}))
                    except InputError:pass
                    else:raise AssertionError('direct-only SCF target accepted as protected restore')
                incompatible=copy.deepcopy(checkpoint['compatibility'])
                incompatible['participants'][0]['states']['A']['identity']='wrong owner or runtime'
                try:execute_checkpoint(comm,groups,states=targets,operation='restore',saved=incompatible)
                except InputError:pass
                else:raise AssertionError('incompatible SCF target admitted')
                assert (state.scf_status,sibling.scf_status)==counters
                for q,value in frame(state).items():np.testing.assert_array_equal(value,cut3[q])
                for q,value in frame(sibling).items():np.testing.assert_array_equal(value,sibling_cut[q])
            else:sibling_snapshot=sibling.snapshot()
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
                if args.protected_api:
                    from salmonpy.checkpointing import CheckpointTarget,AtomicRestoreError
                    targets={'A':state.checkpoint_target(saved),'B':sibling.checkpoint_target(sibling_snapshot)}
                    if n==3:
                        sibling.continue_scf(steps=1)
                        state.restore(snap5)
                        originals=frame(state),frame(sibling)
                        counts=state.scf_status,sibling.scf_status
                        def fail_late():
                            candidate=targets['B'].stage_restore()
                            if comm.rank==comm.size-1:raise RuntimeError('injected late B candidate failure')
                            return candidate
                        broken={'A':targets['A'],'B':CheckpointTarget(targets['B'].capability,
                            stage_restore=fail_late,discard_restore=targets['B'].discard_restore)}
                        try:
                            restore_checkpoint(comm,groups,states=broken,saved=checkpoint['compatibility'],
                                live=CheckpointState(states={'A':state,'B':sibling}))
                        except AtomicRestoreError as exc:assert exc.original_preserved and not exc.cleanup_errors
                        else:raise AssertionError('late candidate failure accepted')
                        assert (state.scf_status,sibling.scf_status)==counts
                        for owner,frame0,final in ((state,originals[0],False),(sibling,originals[1],False)):
                            for q,value in frame(owner,final).items():np.testing.assert_array_equal(value,frame0[q])
                        # Continue the preserved mixing/HSE state before retrying the restore.
                        state.continue_scf(steps=3)
                        assert state.info==info
                        for q,value in frame(state,True).items():np.testing.assert_array_equal(value,expected[q])
                        before=state.scf_status['executed_iterations']
                        # New candidates are required after the failed attempt was discarded.
                        targets={'A':state.checkpoint_target(saved),'B':sibling.checkpoint_target(sibling_snapshot)}
                    holder=CheckpointState(states={'A':state,'B':sibling})
                    staged=restore_checkpoint(comm,groups,states=targets,saved=checkpoint['compatibility'],live=holder)
                    assert state.iterations==8  # only detached payloads have changed
                    reserved_counts=state.scf_status
                    assert session.lib.sp_gs_continue(state._token,1)==11
                    if comm.size>1:
                        mode=3 if comm.rank==0 else 4
                        assert session.lib.sp_gs_scf_candidate(state._token,0,mode)==11
                    assert state.scf_status==reserved_counts
                    for candidate in staged['states'].values():candidate.publish()
                elif args.common_api:
                    execute_checkpoint(comm,groups,states={'A':state.checkpoint_target(saved),
                        'B':sibling.checkpoint_target(sibling_snapshot)},operation='restore',saved=checkpoint['compatibility'])
                else:state.restore(saved)
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
        query=getattr(session.lib,'sp_gs_cartesian_force_protocol',None)
        if query is not None:query.argtypes=[];query.restype=c.c_int
        force_cartesian=query is not None and query()==1
        values=world.allgather((expected,info,force_cartesian))
        for rank in range(args.size):
            for other in (rank+args.size,rank+2*args.size):
                assert values[rank][1]==values[other][1]
                for q,array in values[rank][0].items():
                    previous=values[other][0][q]
                    if q=='forces' and args.symmetry!='none' and other>=2*args.size and not values[other][2]:
                        # Old sym_vector_force_xyz returned direct-lattice components.
                        # Only this intentional units correction differs from the old runtime.
                        previous=a.cell.T@previous
                    np.testing.assert_array_equal(array,previous)
        if args.current_steps:
            counters=state.scf_status if g==0 else None
            model=state.tddft_model(dt=.001,max_steps=args.current_steps+1)
            response=model.new_state('continued-gs')
            try:
                current=np.asarray([response.step(1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)]) if args.symmetry=='none' else np.array([0.,0.,1e-4*np.sin(.17*i)]))
                                    for i in range(1,args.current_steps+1)])
                traces=world.allgather(current)
                for rank in range(args.size):
                    for other in (rank+args.size,rank+2*args.size):
                        np.testing.assert_array_equal(traces[rank],traces[other])
                if g==0:assert state.scf_status==counters
            finally:response.close();model.close()
        state.close()
        full_check=full_k_control(session)
    if world.rank==0:
        result={'status':'pass','size':args.size,'split':args.split,'xc':args.xc,'initial_no_mix':args.initial_no_mix,
                'continuous_arrays_exact':True,
                'continuous_new_and_previous_exact':args.build.resolve()!=args.baseline.resolve() and (args.symmetry=='none' or values[2*args.size][2]),
                'previous_equivalent_with_cartesian_force_units':args.build.resolve()!=args.baseline.resolve(),
                'old_force_lattice_components_converted':args.symmetry!='none' and not values[2*args.size][2],'snapshot_restore_numerical_arrays_exact':True,
                'mixed_history_continuation_exact':True,'sibling_preserved':True,'restore_executed_iterations':0,
                'linear_initializations':linear['initializations'],'linear_executed_iterations':linear['executed_iterations'],
                'final_info':asdict(info),'GS_datasets_collected':False,
                'pre_scf_threshold':args.pre_scf_threshold,'pre_scf_steps':args.pre_scf_steps,'gamma':args.gamma,'symmetry':args.symmetry,'retained_kpoints':expected['occupations'].shape[1],
                'cartesian_force_basis_protocol':force_cartesian,
                'common_direct_checkpoint_api':args.common_api,'protected_SCF_candidates':args.protected_api,
                'late_candidate_failure_preserves_originals':args.protected_api,
                'original_continues_before_retry_exact':args.protected_api,
                'common_incompatible_rejected_before_mutation':args.common_api,
                'continued_GS_to_RT_steps':args.current_steps,
                'continued_GS_to_RT_exact':bool(args.current_steps),'full_k_control':full_check}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback;traceback.print_exc();MPI.COMM_WORLD.Abort(1)
