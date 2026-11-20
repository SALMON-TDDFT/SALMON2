"""Explicit MPI/fresh-process SCF disk regression; writes requested GS data.

Launch parent with BUILD NEW_OUTPUT. Fresh MPI hosts run continuous/save/resume
and compare all local numerical arrays. Never compile or download implicitly.
"""
import argparse,json,os,subprocess,sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--size',type=int,default=2)
    p.add_argument('--restore-size',type=int,help='pure k/orbital SCF checkpoint migration to a different rank count')
    p.add_argument('--repartition',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','hse06'),default='libxc_pbe')
    p.add_argument('--split',choices=('k','orbital','spatial','mixed'),default='k')
    p.add_argument('--symmetry',action='store_true')
    p.add_argument('--decomposition',type=int,nargs=5)
    p.add_argument('--restore-decomposition',type=int,nargs=5)
    p.add_argument('--spatial-axis',choices=('x','y','z'),default='x')
    p.add_argument('--restore-axis',choices=('x','y','z'))
    p.add_argument('--rgrid',type=int,nargs=3)
    p.add_argument('--restore-rgrid',type=int,nargs=3)
    p.add_argument('--grid',type=int,default=8)
    p.add_argument('--two-atoms',action='store_true',help='four occupied bands for uneven orbital splits')
    p.add_argument('--gamma',action='store_true',help='Gamma-only ordinary-XC orbital restart control')
    p.add_argument('--pre-scf-threshold',type=float,default=0.)
    p.add_argument('--pre-scf-steps',type=int,default=4)
    p.add_argument('--without-plot',action='store_true',help='verify unavailable optional plotting does not fail saving')
    p.add_argument('--worker',choices=('continuous','save','resume','reject','protected','cleanup'))
    args=p.parse_args()
    if args.worker is None:
        args.output.mkdir(parents=True,exist_ok=False)
        for phase in ('continuous','save','resume','reject'):
            ranks=(args.restore_size or args.size) if phase in ('resume','reject') else args.size
            command=['mpiexec','-n',str(ranks),sys.executable,str(Path(__file__).resolve()),str(args.build.resolve()),str(args.output.resolve()),
                     '--size',str(ranks),'--xc',args.xc,'--split',args.split,'--worker',phase,
                     '--pre-scf-threshold',str(args.pre_scf_threshold),'--pre-scf-steps',str(args.pre_scf_steps)]
            command+=['--grid',str(args.grid),'--spatial-axis',args.restore_axis if phase in ('resume','reject') and args.restore_axis else args.spatial_axis]
            layout=args.restore_rgrid if phase in ('resume','reject') and args.restore_rgrid else args.rgrid
            if layout:command+=['--rgrid',*map(str,layout)]
            full=args.restore_decomposition if phase in ('resume','reject') and args.restore_decomposition else args.decomposition
            if full:command+=['--decomposition',*map(str,full)]
            if (args.restore_size or args.restore_axis or args.restore_rgrid or args.restore_decomposition) and phase in ('resume','reject'):command+=['--repartition']
            if args.symmetry:command+=['--symmetry']
            if args.gamma:command+=['--gamma']
            if args.two_atoms:command+=['--two-atoms']
            if args.without_plot:command+=['--without-plot']
            with (args.output/(phase+'.log')).open('w') as log:
                subprocess.run(command,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT,check=True)
        result={'status':'pass','fresh_processes':4,'MPI_size':args.size,'xc':args.xc,'split':args.split,
                'symmetry':args.symmetry,'continuous_resume_arrays_exact':not bool(args.restore_size or args.restore_axis or args.restore_rgrid or args.restore_decomposition),'restore_MPI_size':args.restore_size or args.size,
                'comparison_rtol':1e-9 if args.restore_size or args.restore_axis or args.restore_rgrid or args.restore_decomposition else 0.,'comparison_atol':1e-10 if args.restore_size or args.restore_axis or args.restore_rgrid or args.restore_decomposition else 0.,
                'restored_accepted_iterations':3,'bootstrap_executed_iterations':0,
                'resume_executed_iterations':5,'restored_k_slices_exact':args.split=='k','restored_orbital_slices_exact':args.split=='orbital','restored_spatial_interiors_exact':args.split in ('spatial','mixed'),'restored_mixed_orbital_k_interiors_exact':args.split=='mixed','corrupt_rank_rejected_before_creation':True,
                'native_read_failure_preserves_live_state':True,'pre_scf_threshold':args.pre_scf_threshold,
                'GS_datasets_saved':True}
        report=json.loads((args.output/'checkpoint/diagnostics/report.json').read_text())
        assert report['plot_status']==('unavailable' if args.without_plot else 'saved')
        result['residual_plot_status']=report['plot_status']
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result);return
    import numpy as np
    from mpi4py import MPI
    from types import SimpleNamespace
    from salmonpy.dft import DFT,Pseudopotential,SCFSettings,SymmetryOperations
    from salmonpy.native_gs import NativeGSBackend
    from salmonpy.errors import InputError
    comm=MPI.COMM_WORLD
    if comm.size!=args.size:raise ValueError('MPI size mismatch')
    work=args.output/args.worker
    if comm.rank==0:work.mkdir()
    comm.Barrier()
    ctx=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    from salmonpy import dft
    pp=Pseudopotential.from_file(Path(dft.__file__).parent/'examples/si_intro_data/Si_rps.dat',atomic_number=14,local_channel=2)
    symmetry=None
    if args.symmetry:
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:2,:2]=-np.eye(2)
        symmetry=SymmetryOperations(ops,directions='yyn')
    mesh=(1,1,1) if args.gamma else (2,2,2) if args.symmetry else (2,1,1)
    class DiskDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            text=super().input_text(scf)
            if self.xc=='hse06':text=text.replace('exx_pre_scf_threshold=0',f'exx_pre_scf_threshold={args.pre_scf_threshold},exx_pre_scf_steps={args.pre_scf_steps}')
            return text
    calc=DiskDFT((14,),np.diag([8.]*3),[[0,0,0]] if args.symmetry else [[.3,.4,.45]],(pp,),4,2,(args.grid,)*3,mesh,xc=args.xc,symmetry=symmetry)
    if args.two_atoms:
        from dataclasses import replace
        calc=replace(calc,numbers=(14,14),fractional_positions=[[.3,.4,.45],[.7,.6,.55]],electrons=8,states=4)
    scf=SCFSettings(max_iterations=8,tolerance=1e-30)
    dec=(args.size,1,1,1,1) if args.split=='k' else (1,args.size,1,1,1) if args.split=='orbital' else (1,1,args.size,1,1)
    if args.split=='spatial':
        axes=[1,1,1];axes['xyz'.index(args.spatial_axis)]=args.size;dec=(1,1,*(args.rgrid or axes))
    if args.decomposition:dec=tuple(args.decomposition)
    backend=NativeGSBackend(args.build,work,decomposition=dec,hse_orbital=args.xc=='hse06' and args.split=='orbital',hse_spatial=args.xc=='hse06' and args.split in ('spatial','mixed'))
    checkpoint=args.output/'checkpoint'
    def spatial_expected(directory,array,lower):
        wave=array.ndim==7
        result=np.empty_like(array);covered=np.zeros(array.shape[:3]+(array.shape[4:6] if wave else ()),dtype=bool)
        for file in sorted(directory.glob('rank-*.npz')):
            with np.load(file) as old:
                box=old['grid_bounds'];orig=old['orbital_lower'] if wave else old['density_lower']
                if not wave and (old['orbital_lower'][4]!=1 or old['orbital_lower'][5]!=1):continue
                if wave:
                    ol=max(lower[4],orig[4]);oh=min(lower[4]+array.shape[4],orig[4]+old['orbitals'].shape[4])
                    kl=max(lower[5],orig[5]);kh=min(lower[5]+array.shape[5],orig[5]+old['orbitals'].shape[5])
                    if ol>=oh or kl>=kh:continue
                    to=slice(ol-lower[4],oh-lower[4]);tk=slice(kl-lower[5],kh-lower[5])
                    so=slice(ol-orig[4],oh-orig[4]);sk=slice(kl-orig[5],kh-orig[5])
                for z in range(array.shape[2]):
                    for y in range(array.shape[1]):
                        for x in range(array.shape[0]):
                            g=(np.array([x,y,z])+np.asarray(lower[:3])-1)%args.grid+1
                            if np.all(g>=box[:,0]) and np.all(g<=box[:,1]):
                                q=g-orig[:3]
                                if wave:
                                    result[x,y,z,:,to,tk,:]=old['orbitals'][q[0],q[1],q[2],:,so,sk,:]
                                    covered[x,y,z,to,tk]=True
                                else:
                                    result[x,y,z]=old['density'][tuple(q)];covered[x,y,z]=True
        assert covered.all()
        return result
    def interior(value,state):
        rho=state.array('density');orb=state.array('orbitals')
        start=np.asarray(rho.lower_indices)-np.asarray(orb.lower_indices[:3])
        return value[tuple(slice(int(a),int(a+n)) for a,n in zip(start,rho.values.shape))]
    def expected_orbitals(directory,orbital):
        if args.split in ('spatial','mixed'):return spatial_expected(directory,orbital.values,orbital.lower_indices)
        # Test-only global-k assembly, never used by the production loader.
        dim=4 if args.split=='orbital' else 5
        combined=np.empty_like(orbital.values);covered=np.zeros(combined.shape[dim],dtype=bool)
        klo=int(orbital.lower_indices[dim]);khi=klo+combined.shape[dim]-1
        for file in sorted(directory.glob('rank-*.npz')):
            with np.load(file) as old:
                lo=int(old['orbital_lower'][dim]);hi=lo+old['orbitals'].shape[dim]-1
                a=max(lo,klo);b=min(hi,khi)
                if a<=b:
                    dst=[slice(None)]*7;src=list(dst)
                    dst[dim]=slice(a-klo,b-klo+1);src[dim]=slice(a-lo,b-lo+1)
                    combined[tuple(dst)]=old['orbitals'][tuple(src)]
                    covered[a-klo:b-klo+1]=True
        assert covered.all()
        return combined
    with backend.open(ctx) as session:
        if args.worker=='cleanup':
            from salmonpy.native_gs import _scf_candidate_function
            from salmonpy.errors import RuntimeFailure
            class HostFailure(RuntimeError):pass
            def values(state):return {q:state.array(q).values for q in ('density','orbitals','occupations','eigenvalues','scf_history')}
            state=session.run(calc,scf=scf,material_id='memory',steps=3)
            saved=state.snapshot();counts=state.scf_status;before=values(state)
            candidate=state.checkpoint_target(saved).stage_restore()
            try:
                with candidate:raise HostFailure('host exception after staging')
            except HostFailure:pass
            assert state.scf_status==counts
            for q,value in before.items():np.testing.assert_array_equal(state.array(q).values,value)
            saved.close()
            try:
                with state.snapshot() as scoped:
                    scoped_candidate=state.checkpoint_target(scoped).stage_restore()
                    raise HostFailure('snapshot scope must preserve the original host exception')
            except HostFailure:pass
            assert scoped_candidate.finished and scoped.closed
            for q,value in before.items():np.testing.assert_array_equal(state.array(q).values,value)
            saved=state.snapshot()
            candidate=state.checkpoint_target(saved).stage_restore()
            state.close()  # the live owner deterministically owns pending numerical data
            assert candidate.finished and state.closed
            candidate.close();saved.close()
            assert session.resources()['states']==session.resources()['live_child_communicators']==0
            try:
                with session.run(calc,scf=scf,material_id='disk',steps=3) as disk_state:
                    disk_state.save_scf(work/'checkpoint')
                    disk_candidate=disk_state.checkpoint_target(directory=work/'checkpoint').stage_restore()
                    raise HostFailure('exception leaves an unpublished disk candidate')
            except HostFailure:pass
            assert disk_state.closed and disk_candidate.finished
            assert session.resources()['states']==session.resources()['live_child_communicators']==0
            # A reservation created by an external ABI host is not Python-owned.
            with session.run(calc,scf=scf,material_id='external',steps=3) as external:
                snap=external.snapshot();fn=_scf_candidate_function(session.lib)
                assert fn(external._token,snap._token,2)==0
                try:
                    try:external.close()
                    except RuntimeFailure:pass
                    else:raise AssertionError('unowned native candidate was silently discarded')
                    assert not external.closed
                finally:
                    assert fn(external._token,0,4)==0
                    snap.close()
            assert session.resources()['states']==session.resources()['live_child_communicators']==0
            # The same session can still create and evolve a new independent owner.
            with session.run(calc,scf=scf,material_id='after-cleanup',steps=3) as fresh:
                for q,value in before.items():np.testing.assert_array_equal(fresh.array(q).values,value)
                fresh.continue_scf(steps=5)
                expected=values(fresh)
            with session.run(calc,scf=scf,material_id='continuous-reference') as reference:
                for q,value in expected.items():np.testing.assert_array_equal(reference.array(q).values,value)
            if comm.rank==0:(work/'result.json').write_text(json.dumps({'status':'pass','memory_and_disk_candidates_released':True,'original_exception_preserved':True,'external_reservation_not_discarded':True,'GS_and_child_communicators_after_cleanup':0,'normal_results_exact':True})+'\n')
            return
        if args.worker=='protected':
            from dataclasses import replace
            from salmonpy.checkpointing import execute_checkpoint,restore_checkpoint,CheckpointState,AtomicRestoreError
            from unittest.mock import patch
            other=replace(calc,cell=np.diag([9.]*3),grid=(12,12,12))
            with session.run(calc,scf=scf,material_id='A',steps=3) as a, session.run(other,scf=scf,material_id='B',steps=3) as b:
                groups=[{'name':'scf','ranks':list(range(comm.size)),'point_ids':['A','B']}]
                targets={'A':a.checkpoint_target(directory=work/'A-checkpoint'),
                         'B':b.checkpoint_target(directory=work/'B-checkpoint')}
                saved=execute_checkpoint(comm,groups,states=targets)
                a.continue_scf(steps=2);b.continue_scf(steps=1)
                def values(state):return {q:state.array(q).values for q in ('density','orbitals','occupations','eigenvalues','scf_history')}
                before=values(a),values(b);counts=a.scf_status,b.scf_status
                real=targets['B'].stage_restore
                def fail_after_stage():
                    value=real()
                    if comm.rank==comm.size-1:raise RuntimeError('injected late disk candidate failure')
                    return value
                from dataclasses import replace as replace_target
                bad=dict(targets);bad['B']=replace_target(targets['B'],stage_restore=fail_after_stage)
                live=CheckpointState(states={'A':a,'B':b})
                try:restore_checkpoint(comm,groups,states=bad,saved=saved['compatibility'],live=live)
                except AtomicRestoreError:pass
                else:raise AssertionError('late disk candidate failure not reported')
                assert (a.scf_status,b.scf_status)==counts
                for state,expected in zip((a,b),before):
                    for q,value in expected.items():np.testing.assert_array_equal(state.array(q).values,value)
                # Original states still run before retry; compare to independent full runs.
                a.continue_scf(steps=3);b.continue_scf(steps=4)
                expected=values(a),values(b)
                for model,old in zip((calc,other),expected):
                    with session.run(model,scf=scf,material_id='reference') as reference:
                        for q,value in old.items():np.testing.assert_array_equal(reference.array(q).values,value)
                # Create fresh wrappers after the failed staged wrappers were retired.
                targets={'A':a.checkpoint_target(directory=work/'A-checkpoint'),
                         'B':b.checkpoint_target(directory=work/'B-checkpoint')}
                counts=a.scf_status,b.scf_status
                prepared=restore_checkpoint(comm,groups,states=targets,saved=saved['compatibility'],live=live)
                for candidate in prepared['states'].values():candidate.publish()
                assert a.iterations==b.iterations==3
                assert a.scf_status['executed_iterations']==counts[0]['executed_iterations']
                assert b.scf_status['executed_iterations']==counts[1]['executed_iterations']
                a.continue_scf(steps=5);b.continue_scf(steps=5)
                for state,old in zip((a,b),expected):
                    for q,value in old.items():np.testing.assert_array_equal(state.array(q).values,value)
                # Direct one-owner disk restore uses the same staged implementation.
                b_before=values(b)
                a.restore_scf(work/'A-checkpoint')
                assert a.iterations==3 and b.iterations==8
                for q,value in b_before.items():np.testing.assert_array_equal(b.array(q).values,value)
                a.continue_scf(steps=5)
                for q,value in expected[0].items():np.testing.assert_array_equal(a.array(q).values,value)
            if comm.rank==0:(work/'result.json').write_text(json.dumps({'status':'pass','late_failure_preserves_both':True,'originals_continue_before_retry':True,'retry_matches_continuous':True,'restore_SCF_iterations':0})+'\n')
            return
        if args.worker=='reject':
            path=checkpoint/json.loads((checkpoint/'manifest.json').read_text())['files'][-1]['name']
            if comm.rank==comm.size-1:
                with path.open('r+b') as f:
                    f.seek(-1,2);v=f.read(1);f.seek(-1,2);f.write(bytes([v[0]^1]))
            comm.Barrier()
            try:
                try:session.run(calc,scf=scf,material_id='Si',restart=checkpoint,**({'repartition':True} if args.repartition else {}))
                except InputError:pass
                else:raise AssertionError('corrupt disk payload accepted')
            finally:
                comm.Barrier()
                if comm.rank==comm.size-1:
                    with path.open('r+b') as f:f.seek(-1,2);f.write(v)
                comm.Barrier()
            assert not session._states
            return
        with session.run(calc,scf=scf,material_id='Si',steps=None if args.worker=='resume' else 3 if args.worker=='save' else 8,
                         restart=checkpoint if args.worker=='resume' else None,**({'repartition':True} if args.repartition else {})) as state:
            if args.worker=='save':
                assert state.scf_status['executed_iterations']==3
                if args.without_plot:
                    from unittest.mock import patch
                    with patch.dict(sys.modules,{'matplotlib.figure':None}):state.save_scf(checkpoint)
                else:state.save_scf(checkpoint)
                cut={q:state.array(q).values for q in ('density','orbitals','occupations','eigenvalues','scf_history')}
                cut['orbital_lower']=np.asarray(state.array('orbitals').lower_indices,dtype=np.int64)
                rho=state.array('density')
                cut['density_lower']=np.asarray(rho.lower_indices,dtype=np.int64)
                cut['grid_bounds']=np.asarray([[lo,lo+n-1] for lo,n in zip(rho.lower_indices,rho.values.shape)],dtype=np.int64)
                np.savez(work/f'rank-{comm.rank:08d}.npz',**cut)
                if comm.rank==0:
                    import csv
                    rows=list(csv.reader((checkpoint/'diagnostics/residuals.csv').open()))
                    assert len(rows)==4  # header and three accepted iterations
                return
            if args.worker=='resume':
                assert state.iterations==3
                if args.repartition:
                    with np.load(args.output/'save/rank-00000000.npz') as cut:
                        for q in ('density','occupations','eigenvalues','scf_history'):
                            if args.split in ('spatial','mixed') and q=='density':
                                density=state.array(q)
                                np.testing.assert_array_equal(density.values,spatial_expected(args.output/'save',density.values,density.lower_indices))
                            else:np.testing.assert_array_equal(state.array(q).values,cut[q],err_msg='restored cut '+q)
                    orbital=state.array('orbitals')
                    np.testing.assert_array_equal(orbital.values,expected_orbitals(args.output/'save',orbital))
                assert state.scf_status['executed_iterations']==0,'bootstrap replayed SCF'
                if args.xc=='hse06' and args.split=='orbital' and not args.repartition:
                    # Round-trip every serialized numerical field, including
                    # spatial exchange workspace and ACE, before any continuation.
                    roundtrip=work/'restored-checkpoint'
                    state.save_scf(roundtrip)
                    original=(checkpoint/f'rank-{comm.rank:08d}.bin').read_bytes()
                    restored=(roundtrip/f'rank-{comm.rank:08d}.bin').read_bytes()
                    # Five header integers and two execution counters are host data.
                    import ctypes
                    offset=7*ctypes.sizeof(ctypes.c_int)
                    assert restored[offset:]==original[offset:],'HSE numerical payload changed on restore'

                before={q:state.array(q).values for q in ('density','orbitals','occupations','eigenvalues')}
                bad_path=work/'missing.bin' if comm.rank==comm.size-1 else checkpoint/f'rank-{comm.rank:08d}.bin'
                code=session.lib.sp_gs_scf_disk(state._token,os.fsencode(bad_path),2)
                assert comm.allgather(code)==[4]*comm.size
                assert state.scf_status['executed_iterations']==0
                for q,value in before.items():np.testing.assert_array_equal(state.array(q).values,value)
                if args.repartition:
                    from salmonpy.native_gs import _scf_o_repartition_api,_scf_k_repartition_api,_scf_r_repartition_api,_scf_m_repartition_api,_scf_saved_boxes,_scf_repartition_kind
                    saved=json.loads((checkpoint/'manifest.json').read_text())
                    kind=_scf_repartition_kind(saved['decomposition'],session.backend.decomposition)
                    if kind=='m':
                        bounds=np.asarray([_scf_saved_boxes(entry,True) for entry in saved['files']],dtype=np.int32)
                        migrate=_scf_m_repartition_api(session.lib)
                    else:
                        key={'r':'grid_bounds','o':'orbital_bounds','k':'k_bounds'}[kind]
                        bounds=np.asarray([entry[key] for entry in saved['files']],dtype=np.int32)
                        migrate=(_scf_r_repartition_api(session.lib,hse=args.xc=='hse06') if kind=='r' else
                                 _scf_o_repartition_api(session.lib,hse=args.xc=='hse06') if kind=='o' else
                                 _scf_k_repartition_api(session.lib,hse=args.xc=='hse06'))
                    last=checkpoint/saved['files'][-1]['name'];withheld=last.with_suffix('.withheld')
                    if comm.rank==0:last.rename(withheld)
                    comm.Barrier()
                    try:
                        code=migrate(state._token,os.fsencode(checkpoint),saved['ranks'],bounds.ctypes.data)
                        assert comm.allgather(code)==[4]*comm.size
                    finally:
                        comm.Barrier()
                        if comm.rank==0:withheld.rename(last)
                        comm.Barrier()
                    for q,value in before.items():np.testing.assert_array_equal(state.array(q).values,value)
                    assert state.scf_status['executed_iterations']==0

                state.continue_scf(steps=5)
                assert state.scf_status['executed_iterations']==5
                state.write_scf_report(work/'diagnostics')
            arrays={q:state.array(q).values for q in ('density','orbitals','occupations','eigenvalues','forces','scf_history')}
            orbital=state.array('orbitals')
            arrays['info']=np.array([state.info.energy_hartree,state.info.residual,state.info.iterations])
            arrays['orbital_lower']=np.asarray(orbital.lower_indices,dtype=np.int64)
            rho=state.array('density')
            arrays['density_lower']=np.asarray(rho.lower_indices,dtype=np.int64)
            arrays['grid_bounds']=np.asarray([[lo,lo+n-1] for lo,n in zip(rho.lower_indices,rho.values.shape)],dtype=np.int64)
            name=f'rank-{comm.rank:08d}.npz'
            if args.worker=='continuous':np.savez(work/name,**arrays)
            elif not args.repartition:
                with np.load(args.output/'continuous'/name) as expected:
                    for key,value in arrays.items():np.testing.assert_array_equal(value,expected[key],err_msg=key)
            else:
                # Reference retains the source layout/initial random stream, not
                # a separately initialized target layout with different guesses.
                with np.load(args.output/'continuous/rank-00000000.npz') as expected:
                    for key,value in arrays.items():
                        if key not in ('orbitals','orbital_lower','density_lower','grid_bounds') and not (args.split in ('spatial','mixed') and key=='density'):
                            np.testing.assert_allclose(value,expected[key],rtol=1e-9,atol=1e-10,err_msg=key)
                combined=expected_orbitals(args.output/'continuous',orbital)
                if args.split in ('spatial','mixed'):
                    np.testing.assert_allclose(arrays['density'],spatial_expected(args.output/'continuous',rho.values,rho.lower_indices),rtol=1e-9,atol=1e-10)
                    np.testing.assert_allclose(interior(arrays['orbitals'],state),interior(combined,state),rtol=1e-9,atol=1e-10)
                else:np.testing.assert_allclose(arrays['orbitals'],combined,rtol=1e-9,atol=1e-10,err_msg='migrated orbitals')
if __name__=='__main__':main()
