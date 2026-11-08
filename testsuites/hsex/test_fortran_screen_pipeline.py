"""End-to-end Fortran screening vs independent Python preparation."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tools/hsex'),str(Path(__file__).parent)]
from test_prepare_distributed import fixture
from prepare_w import prepare
from sternheimer import equilibrium_operator
from mpi_test_support import shared_fixture,run_mpi

def main():
    exe=Path(sys.argv[1]);nr=int(sys.argv[2]) if len(sys.argv)>2 else 1
    threads=int(sys.argv[3]) if len(sys.argv)>3 else 1
    boxes=int(sys.argv[4]) if len(sys.argv)>4 else 2
    with shared_fixture() as path:
        fixture(path)
        n=np.array([4]*3);mesh=np.array([4]*3);nk=64;ng=64;cell=np.eye(3)*9
        frac=(np.indices((4,)*3).reshape(3,-1)[::-1].T-1.5)/4;k=frac*2*np.pi/9
        eig=np.tile(np.r_[-.8,-.5,np.linspace(.2,3,ng-2)],(nk,1))
        psi=np.fromfile(path/'data_for_restart/wfn.bin','c16').reshape(nk,ng,ng)*np.sqrt(729/ng)
        h=np.array([equilibrium_operator(path,j,n,k,cell,np.array([0.]),np.zeros(ng),b'KLI1') for j in range(nk)])
        for pump in (False,True):
            nb=3 if pump else 2;eta=.001 if pump else 0.
            occ=np.tile(np.r_[2.,2.,np.zeros(nb-2)],(nk,1))
            for qm in (2,4):
                ref=path/f'ref-{pump}-{qm}'
                reference_input = {'pump': path, 'eta': eta} if pump else {'hamiltonians': path}
                prepare(path, ref, mesh, boxes, .2, qmesh=[qm]*3, workers=nr, **reference_input)
                with (path/'pipeline.bin').open('wb') as f:
                    np.r_[n,mesh,[qm]*3,boxes,2,nb].astype('i4').tofile(f)
                    np.r_[.2,eta,cell.T.ravel()].astype('f8').tofile(f)
                    for a in (k,eig[:,:nb],occ):np.asarray(a,dtype='f8').tofile(f)
                    psi[:,:nb].astype('c16').tofile(f);h.transpose(0,2,1).astype('c16').tofile(f)
                run_mpi(exe,nr,threads,path,f'case-{pump}-{qm}')
                actual=np.fromfile(path/'pipeline-output.bin','f8')
                run_mpi(exe, nr, threads, path, f'distributed-{pump}-{qm}',
                        extra_args=('1', 'dense-distributed'))
                distributed = np.fromfile(path/'pipeline-output.bin', 'f8')
                np.testing.assert_array_equal(distributed, actual)
                offset=40+8*(9+3*nk+3+ng+1+3*nk)
                expected=np.frombuffer(ref.read_bytes(),dtype='f8',count=nk*boxes**6,offset=offset)
                np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=2e-11)
                print('PIPELINE',pump,qm,nr,threads,'MAX_ERROR',np.max(abs(actual-expected)),flush=True)
        original=(path/'pipeline.bin').read_bytes()
        # One-rank invalid callback and metadata errors must terminate all MPI ranks.
        for invalid in ('nan','infinity','nan-last-state','fractional-small-q'):
            data=bytearray(original)
            headbytes=12*4+11*8
            if invalid=='fractional-small-q':
                offset=headbytes+3*nk*8+nb*nk*8
                np.frombuffer(data,dtype='f8',count=nb*nk,offset=offset)[:3]=[2.,1.9,.1]
            else:
                offset=headbytes+3*nk*8+2*nb*nk*8
                index=2*ng*nb*nk-1 if invalid=='nan-last-state' else 0
                np.frombuffer(data,dtype='f8',offset=offset)[index]=np.inf if invalid=='infinity' else np.nan
            (path/'pipeline.bin').write_bytes(data)
            run_mpi(exe,nr,threads,path,'reject-'+invalid,timeout=30,expect_failure=True)
            print('COLLECTIVE_REJECTION',invalid,nr,threads,flush=True)
    print('FORTRAN_SCREEN_PIPELINE_PASS')
if __name__=='__main__':main()
