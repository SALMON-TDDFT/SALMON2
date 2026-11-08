"""The refinement must use the same lower-triangle Hermitian H as ZHETRD('L')."""
from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from pump_screen import solve_pump_target
rng=np.random.default_rng(182)
n=np.array([2]*3);ng=8;cell=np.eye(3)*2;k=np.zeros((1,3));mesh=np.ones(3,dtype=int)
u=np.linalg.qr(rng.normal(size=(ng,ng))+1j*rng.normal(size=(ng,ng)))[0]
levels=np.arange(ng)*.25;occ=np.r_[1.8,.2,np.zeros(ng-2)]
h=(u*levels)@u.conj().T
with TemporaryDirectory() as tmp:
    path=Path(tmp)/'hsex-pump-000001.bin'
    def write(matrix):
        with path.open('wb') as f:
            f.write(b'HSP1');np.r_[n,1,1,ng,80].astype('i4').tofile(f)
            for a in (cell,k[0],np.zeros(3),levels,occ): np.asarray(a,'f8').tofile(f)
            u.T.astype('c16').tofile(f);matrix.T.astype('c16').tofile(f)
    job=(tmp,n,k,cell,k,np.zeros((1,3),int),mesh,np.zeros((1,3),int),[np.zeros((1,3),int)],0,[.001,.0005])
    write(h)
    original,_=solve_pump_target(job)
    perturbed=h.copy();perturbed[np.triu_indices(ng,1)]+=1e-11*(1+2j)
    perturbed[np.diag_indices(ng)]+=1e-12j
    write(perturbed)
    changed,_=solve_pump_target(job)
    np.testing.assert_array_equal(original,changed)
    perturbed[0,-1]+=1e-7j
    write(perturbed)
    try:
        solve_pump_target(job)
    except AssertionError:
        pass
    else:
        raise AssertionError('non-Hermitian H outside existing tolerance was accepted')
print('PUMP_LOWER_TRIANGLE_PASS')
