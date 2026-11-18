"""Existing GS→memory RT vs existing-file RT. No SCF; use 2*size ranks."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.native_gs import NativeGSBackend
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider
def main():
    p=argparse.ArgumentParser();p.add_argument('build');p.add_argument('fixture');p.add_argument('output');p.add_argument('--size',type=int,default=1);p.add_argument('--split',default='orbital');a=p.parse_args()
    w=MPI.COMM_WORLD;assert w.size==2*a.size
    g=w.rank//a.size;comm=w.Split(g,w.rank);root=Path(a.output)/str(g)
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier();ctx=SimpleNamespace(solver_comm=comm,point_ids=(1,),closed=False,_dependents=0)
    dec=(1,a.size,1,1,1) if a.split=='orbital' else (1,1,a.size,1,1)
    def A(i):return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])
    trace=[]
    try:
     if g==0:
      with NativeGSBackend(a.build,root,decomposition=dec).open(ctx) as session:
       gs=session.load(a.fixture,material_id='silicon')
       assert gs.info.converged is None
       from salmonpy.errors import InputError
       try:gs.array('forces')
       except InputError:pass
       else:raise AssertionError('uncomputed import forces exposed')
       model=gs.tddft_model(dt=.003,max_steps=32)
       state=model.new_state('a');sibling=model.new_state('b')
       for i in range(1,9):
        trace.append(state.step(A(i)))
        np.testing.assert_array_equal(sibling.step(-A(i)),sibling.current)
       state.checkpoint(root/f'state-{comm.rank}.bin')
       # An existing GS load while RT is active must preserve its inputs/caches.
       loaded=session.load(a.fixture,material_id='second');loaded.close();gs.close()
       for i in range(9,17):trace.append(state.step(A(i)))
       fresh=model.new_state('restore-after-GS-close')
       fresh.checkpoint(root/f'state-{comm.rank}.bin',restore=True)
       for i in range(9,17):np.testing.assert_array_equal(fresh.step(A(i)),trace[i-1])
       independent=model.new_state('independent')
       for i in range(1,17):np.testing.assert_array_equal(independent.step(A(i)),trace[i-1])
       for s in (state,sibling,fresh,independent):s.close()
       model.close()
       assert not list(root.glob('gs-*/**/wfn.bin'))
     else:
      initial=NativeTDDFTMaterialProvider(a.fixture,'rt.inp',max_steps=32).prepare(ctx)
      state=NativeTDDFTBackend(a.build,root,decomposition=dec).create_state(initial,1,.003,ctx)
      try:
       for i in range(1,17):trace.append(state.advance(state.A.copy(),A(i)))
      finally:state.close()
     values=w.allgather(np.asarray(trace));np.testing.assert_array_equal(values[0],values[a.size])
     if w.rank==0:
      result={'status':'pass','size':a.size,'split':a.split,'file_route_exact':True,'independent_states_exact':True,'resume_exact':True,'SCF_calls':0,'replayed_steps':0}
      (Path(a.output)/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
    except BaseException:
     import traceback;traceback.print_exc();w.Abort(1)


if __name__=="__main__":
    main()
