"""Full reduced-W replicated/distributed uneven and empty-k MPI layouts."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tools/hsex'),str(Path(__file__).parent)]
from test_prepare_distributed import fixture
from sternheimer import equilibrium_operator
from mpi_test_support import shared_fixture,run_mpi
exe=Path(sys.argv[1])
with shared_fixture() as path:
    fixture(path);n=np.array([4]*3);ng=64;nb=4;cell=np.eye(3)*9
    frac=(np.indices((4,)*3).reshape(3,-1)[::-1].T-1.5)/4;k=frac*2*np.pi/9
    psi=np.fromfile(path/'data_for_restart/wfn.bin','c16').reshape(64,ng,ng)*np.sqrt(729/ng)
    h=np.array([equilibrium_operator(path,j,n,k,cell,np.array([0.]),np.zeros(ng),b'KLI1') for j in range(64)])
    for label,ranks,selected,mesh in [('uneven',5,np.arange(64),[4,4,4]),('empty',4,np.array([0,2]),[2,1,1])]:
        nk=len(selected);eig=np.tile([-.8,-.5,.2,.2459016393442623],(nk,1))
        # Match the fixture's actual empty levels; empty target values do not
        # enter occupied-only screening but its positive gap must be present.
        eig[:,2:]=[.2,.2459016393442623]
        pop=np.tile([2.,2.,0.,0.],(nk,1))
        with (path/'pipeline.bin').open('wb') as f:
            np.r_[n,mesh,mesh,2,2,nb].astype('i4').tofile(f)
            np.r_[.2,0.,cell.T.ravel()].astype('f8').tofile(f)
            for a in [k[selected],eig,pop]:np.asarray(a,dtype='f8').tofile(f)
            psi[selected,:nb].astype('c16').tofile(f)
            h[selected].transpose(0,2,1).astype('c16').tofile(f)
        run_mpi(exe,ranks,12,path,label+'-reference',extra_args=('1','insulating-replicated'))
        reference=np.fromfile(path/'pipeline-output.bin','f8')
        for mode in ['reduced-replicated','reduced-distributed']:
            run_mpi(exe,ranks,12,path,label+'-'+mode,timeout=120,extra_args=('1',mode))
            actual=np.fromfile(path/'pipeline-output.bin','f8')
            np.testing.assert_allclose(actual,reference,atol=3e-12,rtol=3e-11)
            print('REDUCED_W_LAYOUT_PASS',label,ranks,mode,np.max(abs(actual-reference)),flush=True)
print('REDUCED_W_UNEVEN_EMPTY_MPI_PASS')
