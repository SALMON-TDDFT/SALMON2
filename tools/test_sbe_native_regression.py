"""Regenerate the SBE oracle fixture from the pinned pre-extraction SALMON source.
Run from a checkout: PYTHONPATH=. FC=gfortran python tools/test_sbe_native_regression.py
Build the current native package first. Requires LP64 BLAS (Accelerate on macOS).
The old routines are compiled separately; no current propagation code is reused.
On the same compiler/BLAS, all outputs must be bitwise equal before saving.
"""
from pathlib import Path
import subprocess
import os
import sys
import tempfile
root=Path(__file__).resolve().parents[1]
os.chdir(root)
work=tempfile.TemporaryDirectory(prefix='sbe-oracle-'); d=Path(work.name)
procs=(root/'tests/native/sbe_legacy_procedures.f90').read_text()
(d/'oracle.f90').write_text('''module salmon_global
 integer :: norder_correction=0
 character(32) :: gauge_sbe='velocity_gauge'
end module
module communication
 interface comm_summation
 module procedure sum_r,sum_c
 end interface
contains
 subroutine sum_r(a,b,i)
 real(8),intent(in)::a
 real(8),intent(out)::b
 integer::i
 b=a
 end subroutine
 subroutine sum_c(a,b,n,i)
 integer::n,i
 complex(8),intent(in)::a(n)
 complex(8),intent(out)::b(n)
 b=a
 end subroutine
end module
module oracle
 use sbe_types
 use communication
 implicit none
contains
'''+procs+'\nend module\n')
c=(root/'salmonpy/sbe_c_api.f90').read_text().replace('use sbe_velocity_kernel','use oracle\n use salmon_global').replace('module sbe_c_api','module oracle_c_api').replace('salmonpy_sbe','oracle_sbe').replace('call velocity_evolve','call dt_evolve_bloch')
a=c.index('  call velocity_current_local');z=c.index('  rho=sbe%rho',a)
c=c[:a]+'''  norder_correction=order
  call calc_current_bloch(sbe,gs,anew,j,0)
'''+c[z:]
(d/'api.f90').write_text(c)
subprocess.run([os.environ.get('FC','gfortran'),'-shared','-fPIC','-O3','-J'+str(d),'-I'+str(d),str(root/'native/sbe_types.f90'),str(d/'oracle.f90'),str(d/'api.f90'),*(['-framework','Accelerate'] if sys.platform=='darwin' else ['-lopenblas']),'-o',str(d/'oracle.dylib')],check=True)

import ctypes as ct
from pathlib import Path
import numpy as np
from salmonpy.sbe import SBEInitialState,_native
rng=np.random.default_rng(38219); nb=3;nk=2
p=rng.normal(size=(nb,nb,3,nk))+1j*rng.normal(size=(nb,nb,3,nk));p=(p+p.swapaxes(0,1).conj())*.05
v=p*.13
e=np.array([[-.5,-.49],[.2,.23],[.7,.74]]);w=np.array([1.,2.]); lattice=np.diag([4.,5.,6.]);A=rng.normal(size=(41,3))*.08;dt=.013
lib=ct.CDLL(str(d/'oracle.dylib'));fn=lib.oracle_sbe;fn.argtypes=_native().argtypes;fn.restype=None
results=[];densities=[]
for order in range(4):
 for vnl in (False,True):
  initial=SBEInitialState(e,p,v,w,lattice,2,order,vnl); state=initial.new_state(dt=dt,initial_A=A[0]);rho=state.density.copy(order='F');js=[];rs=[]
  for it in range(len(A)):
   old=np.ascontiguousarray(A[max(it-1,0)]);new=np.ascontiguousarray(A[it]);j=np.empty(3)
   arrays=(initial.eigen,initial.weights,initial.momentum,initial.nonlocal_velocity,rho,old,new)
   fn(nb,nk,2,order,int(vnl),initial.volume,dt,*(x.ctypes.data for x in arrays),int(it>0),j.ctypes.data)
   actual=state.current if it==0 else state.step(new)
   np.testing.assert_array_equal(actual,j)
   np.testing.assert_array_equal(state.density,rho)
   js.append(j.copy());rs.append(rho.copy())
  results.append(js);densities.append(rs)
np.savez('salmonpy/tests/sbe_reference.npz',eigen=e,momentum=p,nonlocal_velocity=v,weights=w,lattice=lattice,A=A,dt=dt,current=results,density=densities)
print('Legacy/native bitwise equality: 8 configurations x 41 endpoints, J and density')
