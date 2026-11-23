"""Synthetic single-state performance comparison; no new GS or GS copies.

mpiexec -n 6 python tools/test_tddft_performance_mpi.py BUILD BASELINE FIXTURE NEW_OUT XC
Three independent groups, ranks/group=world.size/3. FIXTURE must contain rt.inp,
pseudo1_rps.dat and existing data_for_restart for the normalized small synthetic
periodic GS. XC is pz/tbmbj/libxc_pz/libxc_pbe. Changing XC here compares update
paths on identical synthetic initial orbitals, not a converged material response.
Only input/pseudopotential bytes are copied; GS is linked read-only by convention.
Checkpoint restoration is explicit and in memory; no electron-state files saved.
Run with an MPI timeout and one OMP/BLAS thread. Reuse binaries only after SHA.
"""
import ctypes as c,hashlib,json,os,re,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext
from salmonpy.embedded_tddft import _NativeState
from salmonpy.api import _bind_state_operations
from salmonpy._native_group import _bind_performance_directory
build,baseline,fixture,out,xc=sys.argv[1:6];world=MPI.COMM_WORLD;n=world.size//3;g=world.rank//n;comm=world.Split(g,world.rank)
assert world.size>=3 and world.size%3==0
assert xc in ('pz','tbmbj','libxc_pz','libxc_pbe')
mode=('baseline','checked','performance')[g];root=Path(out)/mode
fixture=Path(fixture).resolve()
if comm.rank==0:
 root.mkdir(parents=True,exist_ok=False)
 text=(fixture/'rt.inp').read_text();text=re.sub("xc='pz'", "xc='"+xc+"'",text)
 (root/'rt.inp').write_text(text)
 (root/'pseudo1_rps.dat').write_bytes((fixture/'pseudo1_rps.dat').read_bytes())
 (root/'data_for_restart').symlink_to(fixture/'data_for_restart',target_is_directory=True)
comm.Barrier()
record=json.loads(Path(baseline if mode=='baseline' else build).read_text());assert hashlib.sha256(Path(record['library']).read_bytes()).hexdigest()==record['sha256']
lib=c.CDLL(record['library']);_bind_state_operations(lib,'tddft',route='single')
native=NativeParallelContext(SimpleNamespace(solver_comm=comm,closed=False,_dependents=0),record['library'],record['sha256'])
s=_NativeState(native,lib,root,16,execution_mode='performance' if g==2 else 'checked')
if g==2:
 lib.sp_tddft_step_performance.argtypes=lib.sp_tddft_step.argtypes;lib.sp_tddft_step_performance.restype=c.c_int
 s._performance_enter,s._performance_leave=_bind_performance_directory(lib)
J=np.zeros(3);layout=np.array([n,1,1,1,1],np.int32)
assert s._call(lib.sp_tddft_init,b'rt.inp',.001,16,np.zeros(3),layout,J)==0
s.initialized=True
lib.sp_tddft_step_transaction.argtypes=[c.c_int];lib.sp_tddft_step_transaction.restype=c.c_int
trace=[];norm=[]
def drive(i):return np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)])
for i in range(1,17):
 trace.append(s.advance(s.A.copy(),drive(i)));norm.append(s.electron_number)
 if i==8:
  assert s._call(lib.sp_tddft_step_transaction,0)==0;saved_A=s.A.copy()
assert s._call(lib.sp_tddft_step_transaction,1)==0
s.steps=8;s.A=saved_A
assert s._call(lib.sp_tddft_step_transaction,3)==0
replayed=[]
for i in range(9,17):replayed.append(s.advance(s.A.copy(),drive(i)))
np.testing.assert_array_equal(replayed,trace[8:])
assert s._call(lib.sp_tddft_step_transaction,2)==0
s.close()
rows=world.gather((g,comm.rank,np.array(trace),np.array(norm)),root=0)
if world.rank==0:
 for group,rank,j,ne in rows:
  np.testing.assert_array_equal(j,rows[0][2]);np.testing.assert_array_equal(ne,rows[0][3])
 report={'xc':xc,'group_ranks':n,'steps':16,'current_all_modes_array_equal':True,'electron_number_all_modes_array_equal':True,'explicit_memory_restore_array_equal':True,'new_GS_collected':False,'GS_access':'read existing symlink only','scope':'synthetic propagation, not physical convergence'}
 np.savez(Path(out)/'traces.npz',current=rows[0][2],electron_number=rows[0][3])
 (Path(out)/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
