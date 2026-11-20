"""Existing GS→memory RT vs existing-file RT. No SCF; use 2*size ranks."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.native_gs import NativeGSBackend
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider
def fast_controls(a):
    import shutil
    from salmonpy.fast_hse import FastHSEMaterialProvider
    from salmonpy.fast_hse_mpi import NativeFastHSEBackend
    w=MPI.COMM_WORLD;assert w.size==2*a.size
    g=w.rank//a.size;comm=w.Split(g,w.rank);root=Path(a.output)/str(g)
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier();ctx=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    fixture=Path(a.fixture);trace=[]
    anchor=Path.cwd()
    from unittest.mock import patch
    guard=patch('os.chdir',side_effect=AssertionError('GS to HSE host chdir forbidden')) if a.no_chdir and g==0 else None
    if guard:guard.start()
    from salmonpy.bands import _input_groups,_scalar
    groups=_input_groups((fixture/'rt.inp').read_text())
    symmetric='y' in _scalar(groups.get('system',''),'yn_symmetry','n')
    def drive(i):
        if symmetric:return 1e-4*np.array([0.,0.,np.sin(.23*i)])
        return 1e-4*np.array([np.sin(.23*i),np.cos(.17*i)-1,np.sin(.13*i)])
    session=None
    if g==0:
        staged=root/'initial'
        if comm.rank==0:shutil.copytree(fixture,staged)
        comm.Barrier()
        session=NativeGSBackend(a.build,root,decomposition=(a.size,1,1,1,1)).open(ctx)
        gs=session.load(staged,material_id='silicon')
        model=gs.tddft_model(dt=.001,max_steps=32)
        unused=gs.tddft_model(dt=.001,max_steps=32);unused.close()
        variable=gs.tddft_model(dt=.001,max_steps=32)
        different_dt=variable.new_state('different-dt',dt=.0005,capacity=16)
        different_dt.close();variable.close()
        state=model.new_state('first');sibling=model.new_state('second')
        # A later GS import must not overwrite the RT model's active input view.
        other=session.load(fixture,material_id='other')
        other.close();gs.close()
        if comm.rank==0:shutil.rmtree(staged)
        comm.Barrier()
    else:
        initial=FastHSEMaterialProvider(fixture,max_steps=32).prepare(ctx)
        model=NativeFastHSEBackend(a.build,root,allow_symmetry=True).open_model(initial,.001,ctx)
        state=model.new_state('first');sibling=model.new_state('second')
    try:
        for i in range(1,17):
            trace.append(np.r_[state.step(drive(i)),state.energy,sibling.step(-drive(i)),sibling.energy])
            if i==8:
                path=root/f'state-{comm.rank}.bin';state.checkpoint(path,False)
                restored=model.new_state('restored');restored.checkpoint(path,True)
            if i>8:
                np.testing.assert_array_equal(restored.step(drive(i)),state.current)
                assert restored.energy==state.energy
        if g==0:
            assert not (model.work/'in-memory-no-gs').exists()
            assert not list(model.work.rglob('wfn.bin'))
        fresh=model.new_state('fresh')
        np.testing.assert_array_equal(fresh.current,state.initial_current)
        fresh.close();restored.close()
        model.close()
        state.step(drive(17));sibling.step(-drive(17))
    finally:
        state.close();sibling.close();model.close()
        if session is not None:session.close()
    assert Path.cwd()==anchor
    if guard:guard.stop()
    rows=comm.gather(np.array(trace),root=0)
    gathered=w.gather(rows if comm.rank==0 else None,root=0)
    if w.rank==0:
        for first,second in zip(gathered[0],gathered[a.size]):np.testing.assert_allclose(first,second,rtol=3e-13,atol=3e-15)
        error=np.abs(gathered[0][0]-gathered[a.size][0])
        result={'status':'pass','ranks':w.size,'memory_vs_file_current_energy':True,'new_runtime_no_chdir':a.no_chdir,
                'max_current_difference_au':float(error[:,[0,1,2,4,5,6]].max()),
                'max_energy_difference_hartree':float(error[:,[3,7]].max()),
                'independent_states':True,'source_GS_closed_and_removed':True,
                'interleaved_GS_import':True,'full_k_three_component_drive':not symmetric,
                'direct_resume_exact':True,'SCF_calls':0,'replayed_steps':0}
        (Path(a.output)/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('build');p.add_argument('fixture');p.add_argument('output');p.add_argument('--size',type=int,default=1);p.add_argument('--split',default='orbital');p.add_argument('--external-k',action='store_true');p.add_argument('--describe-input',action='store_true');p.add_argument('--fast-hse',action='store_true');p.add_argument('--control-build');p.add_argument('--interleave-fixture');p.add_argument('--measure-norm',action='store_true');p.add_argument('--no-chdir',action='store_true');a=p.parse_args()
    if a.fast_hse:
        try:fast_controls(a)
        except BaseException:
            import traceback;traceback.print_exc();MPI.COMM_WORLD.Abort(1)
        return
    w=MPI.COMM_WORLD;assert w.size==2*a.size
    g=w.rank//a.size;comm=w.Split(g,w.rank);root=Path(a.output)/str(g)
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier();ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,),closed=False,_dependents=0)
    dec=((1,a.size,1,1,1) if a.split=='orbital' else (a.size,1,1,1,1) if a.split=='k' else (1,1,a.size,1,1))
    def A(i):return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])
    anchor=Path.cwd()
    from unittest.mock import patch
    guard=patch('os.chdir',side_effect=AssertionError('GS to TDDFT host chdir forbidden')) if a.no_chdir and g==0 else None
    if guard:guard.start()
    trace=[]
    try:
     if g==0:
      fixture=Path(a.fixture)
      if a.external_k:
       import re,shutil
       from salmonpy.bands import _input_groups,_scalar
       staged=root/'external-k-input';error=None
       if comm.rank==0:
        try:
         groups=_input_groups((fixture/'rt.inp').read_text())
         kbody=groups['kgrid']
         match=re.search(r'num_kgrid\s*=\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)',kbody,re.I)
         if match is None or _scalar(kbody,'file_kw','none')!='none' or _scalar(kbody,'yn_gamma_centered','n')!='n':
          raise ValueError('external-k comparison requires an unshifted regular full source mesh')
         if re.search(r'\bdk_shift\b',kbody,re.I):
          raise ValueError('external-k comparison uses source meshes without explicit dk_shift')
         mx,my,mz=map(int,match.groups())
         points=[[(ix+.5)/mx-.5,(iy+.5)/my-.5,(iz+.5)/mz-.5]
                 for iz in range(mz) for iy in range(my) for ix in range(mx)]
         shutil.copytree(fixture,staged)
         (staged/'tables').mkdir(exist_ok=True)
         lines=[str(len(points))]+[str(i+1)+' '+' '.join(format(float(v),'.17e') for v in (*point,1/len(points)))
                                   for i,point in enumerate(points)]
         (staged/'tables/kpoints.dat').write_text('\n'.join(lines)+'\n')
         groups['kgrid']="\n file_kw='tables/kpoints.dat'\n"
         text=''.join('&'+key+body+'/\n' for key,body in groups.items())
         (staged/'rt.inp').write_text(text)
         (staged/'rt-missing.inp').write_text(text.replace('tables/kpoints.dat','tables/missing.dat'))
        except Exception as exc:error=str(exc)
       error=comm.bcast(error,root=0)
       if error:raise ValueError(error)
       fixture=staged
      with NativeGSBackend(a.build,root,decomposition=dec).open(ctx) as session:
       gs=session.load(fixture,material_id='silicon')
       if a.external_k:
        assert (gs.directory/'tables/kpoints.dat').read_bytes()==(fixture/'tables/kpoints.dat').read_bytes()
        assert any(name=='tables/kpoints.dat' for name,_ in gs._assets)
       if a.describe_input:
        metadata_error=None
        try:
         description=gs.input_metadata()
         record=description.as_dict()
         assert record['material_id']=='silicon'
         assert record['evidence']['origin']=='imported' and record['evidence']['scf_converged'] is None
         assert record['evidence']['checkpoint_provenance']=='unknown'
         assert description.settings['cell_bohr'] is not None
         assert description.settings['fractional_positions'] is not None
         if a.external_k:assert description.settings['k_sampling']['mode']=='file'
        except Exception as exc:metadata_error=str(exc) or type(exc).__name__
        errors=comm.allgather(metadata_error)
        if any(errors):raise RuntimeError('input description failed: '+str(errors))
        assert comm.allgather(description.source_id)==[description.source_id]*comm.size
       assert gs.info.converged is None
       from salmonpy.errors import InputError
       try:gs.array('forces')
       except InputError:pass
       else:raise AssertionError('uncomputed import forces exposed')
       model=gs.tddft_model(dt=.003,max_steps=32)
       state=model.new_state('a');sibling=model.new_state('b')
       norm_trace=[]
       if a.measure_norm:
        initial_norm=state.electron_number;sibling_norm=sibling.electron_number
        assert initial_norm>0 and initial_norm==sibling_norm
       for i in range(1,9):
        trace.append(state.step(A(i)))
        if a.measure_norm:
         before=state.current.copy();norm_trace.append(state.electron_number)
         np.testing.assert_array_equal(before,state.current)
         if i==1:assert sibling.electron_number==sibling_norm
        np.testing.assert_array_equal(sibling.step(-A(i)),sibling.current)
       if getattr(session.backend,'_ordinary_rt',False):
        import shutil
        from salmonpy.bands import _input_groups
        bad=root/'unsupported-xc'
        if comm.rank==0:
         shutil.copytree(fixture,bad)
         groups=_input_groups((bad/'rt.inp').read_text())
         groups['functional']="\n xc='hse06'\n"
         (bad/'rt.inp').write_text(''.join('&'+key+body+'/\n' for key,body in groups.items()))
        comm.Barrier();before=session.resources();current=state.current.copy()
        try:session.load(bad,material_id='unsupported-hybrid')
        except InputError as exc:assert 'HSE-disabled' in str(exc)
        else:raise AssertionError('HSE reached an ordinary-only native runtime')
        assert session.resources()==before
        np.testing.assert_array_equal(state.current,current)
       state.checkpoint(root/f'state-{comm.rank}.bin')
       if a.external_k:
        # Admission failure must precede native initialization and preserve RT/GS.
        before=session.resources();current=state.current.copy();caught=False
        try:session.load(fixture,'rt-missing.inp',material_id='missing-table')
        except InputError:caught=True
        assert all(comm.allgather(caught))
        assert session.resources()==before
        np.testing.assert_array_equal(state.current,current)
       if a.measure_norm:saved_norm=state.electron_number
       # An existing GS load while RT is active must preserve its inputs/caches.
       loaded=session.load(a.interleave_fixture or fixture,material_id='second');loaded.close();gs.close()
       if a.describe_input:assert gs.input_metadata().as_dict()==description.as_dict()
       for i in range(9,17):trace.append(state.step(A(i)))
       fresh=model.new_state('restore-after-GS-close')
       fresh.checkpoint(root/f'state-{comm.rank}.bin',restore=True)
       if a.measure_norm:assert fresh.electron_number==saved_norm
       for i in range(9,17):np.testing.assert_array_equal(fresh.step(A(i)),trace[i-1])
       independent=model.new_state('independent')
       for i in range(1,17):np.testing.assert_array_equal(independent.step(A(i)),trace[i-1])
       if a.measure_norm:
        assert state.electron_number==fresh.electron_number==independent.electron_number
       for s in (state,sibling,fresh,independent):s.close()
       model.close()
       assert not list(root.glob('gs-*/**/wfn.bin'))
     else:
      initial=NativeTDDFTMaterialProvider(a.fixture,'rt.inp',max_steps=32).prepare(ctx)
      state=NativeTDDFTBackend(a.control_build or a.build,root,decomposition=dec).create_state(initial,1,.003,ctx)
      try:
       for i in range(1,17):trace.append(state.advance(state.A.copy(),A(i)))
      finally:state.close()
     assert Path.cwd()==anchor
     if guard:guard.stop()
     values=w.allgather(np.asarray(trace))
     difference=float(np.max(abs(values[0]-values[a.size])))
     exact=bool(np.array_equal(values[0],values[a.size]))
     if a.control_build:
      np.testing.assert_allclose(values[0],values[a.size],rtol=3e-13,atol=3e-15)
     else:np.testing.assert_array_equal(values[0],values[a.size])
     if w.rank==0:
      result={'status':'pass','size':a.size,'split':a.split,'file_route_exact':exact,'new_runtime_no_chdir':a.no_chdir,'file_route_within_tolerance':True,'max_current_difference_au':difference,'independent_states_exact':True,'resume_exact':True,'SCF_calls':0,'replayed_steps':0,'input_metadata':a.describe_input,'external_k_table':a.external_k,'missing_asset_rejected_before_native':a.external_k,'measured_norm_checked':a.measure_norm}
      (Path(a.output)/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
    except BaseException:
     import traceback;traceback.print_exc();w.Abort(1)


if __name__=="__main__":
    main()
