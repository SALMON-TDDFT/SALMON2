"""Compare two Python shared TDDFT states to independent controls, then direct resume.
Three groups; --size is ranks per group; --rgrid chooses spatial splitting.
"""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
p=argparse.ArgumentParser();p.add_argument('build');p.add_argument('fixture');p.add_argument('out');p.add_argument('--size',type=int,default=2);p.add_argument('--rgrid',action='store_true');p.add_argument('--baseline');p.add_argument('--no-chdir',action='store_true');a=p.parse_args()
world=MPI.COMM_WORLD
try:
 assert world.size==3*a.size
 g=world.rank//a.size;c=world.Split(g,world.rank);root=Path(a.out)/str(g)
 if c.rank==0:root.mkdir(parents=True,exist_ok=False)
 c.Barrier();ctx=SimpleNamespace(solver_comm=c,point_ids=(1,2) if g==0 else (g,),closed=False,_dependents=0)
 anchor=Path.cwd()
 from unittest.mock import patch
 guard=patch('os.chdir',side_effect=AssertionError('shared TDDFT host chdir forbidden')) if a.no_chdir and g==0 else None
 if guard:guard.start()
 provider=NativeTDDFTMaterialProvider(a.fixture,'rt.inp',max_steps=32);initial=provider.prepare(ctx)
 backend=NativeTDDFTBackend(a.baseline if g and a.baseline else a.build,root,decomposition=(1,1,a.size,1,1) if a.rgrid else (1,a.size,1,1,1))
 states=[];model=None;traces={}
 def drive(i,k):return np.array([1e-4*np.sin(.3*i),7e-5*(np.cos(.2*i)-1),1e-4*np.sin(.17*i)])*(1 if k==1 else -.7)
 try:
  if g==0:
   model=backend.open_model(initial,.003,ctx);states=[model.new_state('a'),model.new_state('b')]
   assert len(list(root.glob('*/bundle/data_for_restart')))==1
  else:states=[backend.create_state(initial,g,.003,ctx)]
  for offset,s in enumerate(states):traces[offset+1 if g==0 else g]=[]
  for i in range(1,33):
   for offset,s in enumerate(states):
    k=offset+1 if g==0 else g
    if g==0 and i==1:
     from salmonpy.errors import InputError
     for bad in ([float('nan'),0,0],[float('inf'),0,0],[1j,0,0]):
      try:s.step(bad)
      except InputError:pass
      else:raise AssertionError('invalid A accepted')
     assert s.steps==0
    j=s.step(drive(i,k)) if g==0 else s.advance(s.A.copy(),drive(i,k))
    traces[k].append(j)
    if g==0 and i in (8,9) and k==i-7:s.checkpoint(root/f'{k}-{c.rank}.bin')
  if g==0:
   from salmonpy.errors import RuntimeFailure
   collision=model.work/model.base/'native_tddft_h3_g1_rt.data'
   if c.rank==0:collision.write_text('keep this output')
   c.Barrier()
   try:model.new_state('collision')
   except RuntimeFailure:pass
   else:raise AssertionError('existing output overwritten')
   assert not model.failed
   if c.rank==0:
    assert collision.read_text()=='keep this output';collision.unlink()
   c.Barrier()
   for k,split in ((1,8),(2,9)):
    fresh=model.new_state(k);states.append(fresh);fresh.checkpoint(root/f'{k}-{c.rank}.bin',True)
    for i in range(split+1,33):np.testing.assert_array_equal(fresh.step(drive(i,k)),traces[k][i-1])
   if a.size>1:
    mixed=model.new_state('bad-shards');states.append(mixed)
    path=root/(f'1-{c.rank}.bin' if c.rank==0 else f'2-{c.rank}.bin')
    try:mixed.checkpoint(path,True)
    except RuntimeFailure:pass
    else:raise AssertionError('mixed checkpoint epochs accepted')
    assert mixed.failed
    mixed.close();states.remove(mixed)
   for s in reversed(states):s.close()
   states=[]
   fresh=model.new_state('idle');states.append(fresh);model.close()
   np.testing.assert_array_equal(fresh.step(drive(1,1)),traces[1][0])
  np.savez(root/f'trace-{c.rank}.npz',**{f'p{k}':v for k,v in traces.items()})
 finally:
  for s in reversed(states):s.close()
  if model:model.close()
 assert ctx._dependents==0
 assert Path.cwd()==anchor
 if guard:guard.stop()
 gathered=world.allgather(traces)
 for k in (1,2):np.testing.assert_array_equal(gathered[0][k],gathered[k*a.size][k])
 if world.rank==0:(Path(a.out)/'result.json').write_text(json.dumps({'status':'pass','group_size':a.size,'rgrid':a.rgrid,'charge_current_matches_independent':True,'new_runtime_no_chdir':a.no_chdir,'previous_runtime_compared':bool(a.baseline),'resume_exact':True,'no_replay':True})+'\n')
 c.Free()
except BaseException:
 import traceback;traceback.print_exc();world.Abort(1)
