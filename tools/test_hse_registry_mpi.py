"""Different full/reduced HSE models vs independent old-runtime controls, no SCF."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import ctypes as c
import numpy as np
from mpi4py import MPI
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend,_history_info
from salmonpy.errors import InputError,RuntimeFailure
def main():
    p=argparse.ArgumentParser();p.add_argument('build');p.add_argument('baseline');p.add_argument('fixture_a');p.add_argument('fixture_b');p.add_argument('output');p.add_argument('--size',type=int,default=1);a=p.parse_args()
    w=MPI.COMM_WORLD;assert w.size==3*a.size
    g=w.rank//a.size;comm=w.Split(g,w.rank);out=Path(a.output)/str(g)
    if comm.rank==0:out.mkdir(parents=True,exist_ok=False)
    comm.Barrier();context=SimpleNamespace(solver_comm=comm,point_ids=(1,2) if not g else (g,),closed=False,_dependents=0)
    fixtures={1:a.fixture_a,2:a.fixture_b};dts={1:.001,2:.002}
    traces={};hist={};states={};models={};registry=None
    initial={k:FastHSEMaterialProvider(v,max_steps=32).prepare(context) for k,v in fixtures.items()}
    def drive(k,i):return np.array([0.,0.,(1 if k==1 else -.7)*1e-4*np.sin(.3*i)])
    def step(s,k,i):return np.r_[s.step(drive(k,i)),s.energy]
    try:
     if g:
      s=NativeFastHSEBackend(a.baseline,out,allow_symmetry=True).create_state(initial[g],g,dts[g],context)
      try:
       traces[g]=[step(s,g,i) for i in range(1,26)];hist[g]=s._call(_history_info,s.lib).tolist()
      finally:s.close()
     else:
      registry=NativeFastHSEBackend(a.build,out,allow_symmetry=True).open_registry(context)
      for k in (1,2):
       models[k]=registry.open_model(initial[k],dts[k]);states[k]=models[k].new_state(str(k));traces[k]=[]
      assert context._dependents==1
      # Different model tokens cannot select different workspaces on different ranks.
      if comm.size>1:
       try:models[1 if comm.rank==0 else 2].new_state('wrong',dt=.001,capacity=32)
       except InputError:pass
       else:raise AssertionError('mixed public model selection accepted')
       token=models[1 if comm.rank==0 else 2].lib.token
       handle=c.c_int64();current=np.zeros(3)
       assert registry.lib.sp_fast_hse_model_spawn(token,.001,32,c.byref(handle),current)==11
      for i in range(1,25):
       for k in (1,2):
        traces[k].append(step(states[k],k,i))
        if i==7:states[k].checkpoint(out/f'{k}-{comm.rank}.bin',False)
      for k in (1,2):
       fresh=models[k].new_state('restored-'+str(k))
       try:
        fresh.checkpoint(out/f'{k}-{comm.rank}.bin',True)
        for i in range(8,25):np.testing.assert_array_equal(step(fresh,k,i),traces[k][i-1])
       finally:fresh.close()
       independent=models[k].new_state('fresh-'+str(k))
       try:
        for i in range(1,25):np.testing.assert_array_equal(step(independent,k,i),traces[k][i-1])
       finally:independent.close()
      # Closing models/registry stops creation but cannot retire live numerical children.
      registry.close();assert context._dependents==1
      for k in (1,2):
       try:models[k].new_state('must-reject')
       except InputError:pass
       else:raise AssertionError('closed model allowed a new state')
       traces[k].append(step(states[k],k,25));hist[k]=states[k]._call(_history_info,states[k].lib).tolist()
      for k in (1,2):states[k].close()
      states.clear();assert context._dependents==0
     traces={k:np.array(v) for k,v in traces.items()}
     all_results=w.allgather((traces,hist))
     for k in (1,2):
      np.testing.assert_array_equal(all_results[0][0][k],all_results[k*a.size][0][k])
      assert all_results[0][1][k]==all_results[k*a.size][1][k]
     if w.rank==0:
      result={'status':'pass','size':a.size,'old_independent_current_energy_history_exact':True,'fresh_state_exact':True,'direct_resume_exact':True,'registry_close_preserves_states':True,'mixed_model_selection_rejected':a.size>1,'replayed_steps':0,'SCF_calls':0}
      (Path(a.output)/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
    except BaseException:
     import traceback;traceback.print_exc();w.Abort(1)


if __name__=="__main__":main()
