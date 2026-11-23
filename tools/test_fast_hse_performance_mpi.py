"""Single-state fast-PR HSE performance regression against checked/previous runtime.

mpiexec -n 6 python tools/test_fast_hse_performance_mpi.py BUILD BASELINE FIXTURE NEW_OUT
Use --symmetry for a prepared yyn reduced-k fixture (z drive), otherwise xyz drive.
FIXTURE contains inputfile, pseudo1_rps.dat, sym.dat and existing data_for_restart.
No GS calculation/copy or electronic disk save. Use OMP/BLAS=1 and MPI timeout.
Synthetic fixture results establish update parity, not material convergence.
"""
import argparse,ctypes as c,hashlib,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext
from salmonpy.fast_hse_mpi import _State,_configure,_history_info
from salmonpy._native_group import _bind_performance_directory
from salmonpy.errors import InputError
from salmonpy.fast_hse import _input_text
p=argparse.ArgumentParser(description=__doc__)
for name in ('build','baseline','fixture','output'):p.add_argument(name)
p.add_argument('--symmetry',action='store_true');a=p.parse_args()
world=MPI.COMM_WORLD
if world.size<3 or world.size%3:raise ValueError('three equally sized solver groups required')
n=world.size//3;g=world.rank//n;comm=world.Split(g,world.rank)
root=Path(a.output)/str(g);fixture=Path(a.fixture).resolve()
if comm.rank==0:
 root.mkdir(parents=True,exist_ok=False)
 groups,_=_input_text((fixture/'inputfile').read_text())
 groups['parallel']=groups.get('parallel','')+'\n nproc_k='+str(n)+', nproc_ob=1, nproc_rgrid=1,1,1\n'
 (root/'inputfile').write_text('\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
 for name in ('pseudo1_rps.dat','sym.dat'):(root/name).write_bytes((fixture/name).read_bytes())
 (root/'data_for_restart').symlink_to(fixture/'data_for_restart',target_is_directory=True)
comm.Barrier()
record=json.loads(Path(a.baseline if g==0 else a.build).read_text())
assert hashlib.sha256(Path(record['library']).read_bytes()).hexdigest()==record['sha256']
lib=c.CDLL(record['library']);_configure(lib)
native=NativeParallelContext(SimpleNamespace(solver_comm=comm,closed=False,_dependents=0),record['library'],record['sha256'])
mode='performance' if g==2 else 'checked'
s=_State(native,lib,root,.001,16,'probe',execution_mode=mode)
if g==2:
 lib.sp_fast_hse_step_performance.argtypes=lib.sp_fast_hse_step.argtypes;lib.sp_fast_hse_step_performance.restype=c.c_int
 s._performance_enter,s._performance_leave=_bind_performance_directory(lib)
initializer=lib.sp_fast_hse_init_symmetry if a.symmetry else lib.sp_fast_hse_init
initializer.argtypes=[c.c_char_p,c.c_double,c.c_int,lib.sp_fast_hse_step.argtypes[0]];initializer.restype=c.c_int
j=np.zeros(3);assert s._call(initializer,b'inputfile',.001,16,j)==0
s.initialized=True;s.current=j.copy();s.initial_current=j.copy()
if a.symmetry:
 before=s.current.copy()
 try:s.advance(s.A.copy(),[.001,0.,0.])
 except InputError:pass
 else:raise AssertionError('symmetry-breaking drive accepted')
 assert not s.failed and s.steps==0;np.testing.assert_array_equal(s.current,before)
lib.sp_fast_hse_step_transaction.argtypes=[c.c_int];lib.sp_fast_hse_step_transaction.restype=c.c_int
trace=[];energies=[];histories=[]
def drive(i):
 return np.array([0.,0.,.001*np.sin(.17*i)]) if a.symmetry else np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)])
for i in range(1,17):
 trace.append(s.step(drive(i)));energies.append(s.energy)
 histories.append(s._call(lambda:_history_info(lib)))
 if i==8:
  assert s._call(lib.sp_fast_hse_step_transaction,0)==0
  saved=(s.steps,s.A.copy(),s.current.copy(),s.energy)
assert s._call(lib.sp_fast_hse_step_transaction,1)==0
s.steps,s.A,s.current,s.energy=saved
assert s._call(lib.sp_fast_hse_step_transaction,3)==0
for i in range(9,17):
 np.testing.assert_array_equal(s.step(drive(i)),trace[i-1]);assert s.energy==energies[i-1]
 np.testing.assert_array_equal(s._call(lambda:_history_info(lib)),histories[i-1])
assert s._call(lib.sp_fast_hse_step_transaction,2)==0
s.close()
rows=world.gather((g,comm.rank,np.array(trace),np.array(energies),np.array(histories)),root=0)
if world.rank==0:
 for row in rows:
  for index in (2,3,4):np.testing.assert_array_equal(row[index],rows[0][index])
 report={'ranks_per_group':n,'symmetry':a.symmetry,'steps':16,'current_energy_history_all_modes_array_equal':True,
 'explicit_memory_restore_array_equal':True,'symmetry_breaking_A_rejected_before_update':a.symmetry,
 'GS':'existing read symlink; no new GS/copies/electronic disk saves'}
 np.savez(Path(a.output)/'traces.npz',current=rows[0][2],energy=rows[0][3],history=rows[0][4])
 (Path(a.output)/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
