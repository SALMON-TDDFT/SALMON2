"""Compare scalar matrix-free W with the dense resolvent, including the optical cell."""
from pathlib import Path
import sys,re
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tools/hsex'),str(Path(__file__).parent)]
from test_prepare_distributed import fixture
from sternheimer import equilibrium_operator
from mpi_test_support import shared_fixture,run_mpi


def main():
    exe=Path(sys.argv[1]);ranks=int(sys.argv[2]) if len(sys.argv)>2 else 1
    threads=int(sys.argv[3]) if len(sys.argv)>3 else 1
    with shared_fixture() as path:
        fixture(path)
        n=np.array([4]*3);mesh=np.array([4]*3);nk=64;ng=64;nb=4;cell=np.eye(3)*9
        frac=(np.indices((4,)*3).reshape(3,-1)[::-1].T-1.5)/4;k=frac*2*np.pi/9
        eig=np.tile(np.r_[-.8,-.5,np.linspace(.2,3,ng-2)],(nk,1))
        psi=np.fromfile(path/'data_for_restart/wfn.bin','c16').reshape(nk,ng,ng)*np.sqrt(729/ng)
        h=np.array([equilibrium_operator(path,j,n,k,cell,np.array([0.]),np.zeros(ng),b'KLI1') for j in range(nk)])
        occ=np.tile([2.,2.,0.,0.],(nk,1))
        # Deliberately rotate unused empty orbitals away from eigenvectors.
        # Integer occupied-projector screening must not depend on their
        # arbitrary empty-space gauge or convergence residuals.
        a=psi[:,2].copy();b=psi[:,3].copy()
        psi[:,2]=.6*a+.8*b;psi[:,3]=-.8*a+.6*b
        for qm in [2,4]:
            with (path/'pipeline.bin').open('wb') as f:
                np.r_[n,mesh,[qm]*3,2,2,nb].astype('i4').tofile(f)
                np.r_[.2,0.,cell.T.ravel()].astype('f8').tofile(f)
                for a in (k,eig[:,:nb],occ):np.asarray(a,dtype='f8').tofile(f)
                psi[:,:nb].astype('c16').tofile(f);h.transpose(0,2,1).astype('c16').tofile(f)
            run_mpi(exe,ranks,threads,path,f'dense-{qm}',extra_args=('1','dense-optical-distributed'))
            expected=np.fromfile(path/'pipeline-output.bin','f8')
            run_mpi(exe,ranks,threads,path,f'matrix-free-{qm}',extra_args=('1','insulating-distributed'))
            actual=np.fromfile(path/'pipeline-output.bin','f8')
            np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=2e-11)
            assert np.isfinite(actual).all()
            for reduced_mode in ['reduced-distributed','reduced-sym-distributed','reduced-replicated']:
                run_mpi(exe,ranks,threads,path,f'{reduced_mode}-{qm}',extra_args=('1',reduced_mode))
                compressed=np.fromfile(path/'pipeline-output.bin','f8')
                np.testing.assert_allclose(compressed,expected,atol=3e-12,rtol=3e-11)
                assert np.isfinite(compressed).all()
                print('REDUCED_FULL_W_EQUIVALENCE_PASS',qm,ranks,reduced_mode,
                      np.max(abs(compressed-expected)),flush=True)
            run_mpi(exe,ranks,threads,path,f'reduced-truncated-{qm}',
                    extra_args=('1','reduced-truncated-distributed'))
            truncated=np.fromfile(path/'pipeline-output.bin','f8')
            assert np.isfinite(truncated).all() and np.linalg.norm(truncated-expected)>1e-8
            print('REDUCED_FULL_W_TRUNCATION_PASS',qm,ranks,np.linalg.norm(truncated-expected),flush=True)
            for label in [f'dense-{qm}',f'matrix-free-{qm}']:
                lines=(path/label/'stdout.rank-0').read_text().splitlines()
                timing=[line for line in lines if line.startswith('HSEX_W_PHASE_SECONDS ')]
                assert len(timing)==1, 'Missing per-phase screening timers'
                seconds=np.array(timing[0].split(':')[1].split(),float)
                assert seconds.shape==(11,) and np.isfinite(seconds).all() and (seconds>=0).all()
                if label.startswith('matrix-free'):
                    assert np.array_equal(seconds[:2],[0.,0.]), 'Matrix-free route constructed/factored a dense H'
            print('SCALAR_MATRIX_FREE_W_MAX_ERROR',qm,ranks,np.max(abs(actual-expected)),flush=True)
            label=f'hybrid-{qm}'
            run_mpi(exe,ranks,threads,path,label,extra_args=('1','hybrid-distributed'))
            hybrid=np.fromfile(path/'pipeline-output.bin','f8')
            np.testing.assert_allclose(hybrid,expected,atol=2e-12,rtol=2e-11)
            transcript=(path/label/'stdout.rank-0').read_text()
            parallel=re.search(r'HSEX_W_Q_PARALLEL regions=(\d+) batches=(\d+)',transcript)
            assert parallel and int(parallel[1])==(nk+ranks-1)//ranks
            assert int(parallel[2])>=int(parallel[1])
            if qm==4 and threads<32:
                assert int(parallel[2])>int(parallel[1]), 'q batch fork/join was not amortized'
            calls=re.search(r'SCREEN_REALSPACE_ACTION_CALLS (\d+)',transcript)
            assert calls and int(calls[1])>0, 'Hybrid optical response did not use the real-space action'
            print('SCALAR_HYBRID_W_MAX_ERROR',qm,ranks,np.max(abs(hybrid-expected)),flush=True)
        original=(path/'pipeline.bin').read_bytes()
        # Explicitly reject unsupported fractional scalar screening and invalid
        # orbital data collectively, without using the dense callback.
        for invalid in ['fractional','nan-last-state','infinity']:
            data=bytearray(original);head=12*4+11*8
            if invalid=='fractional':
                occupations=np.frombuffer(data,'f8',count=nb*nk,offset=head+3*nk*8+nb*nk*8).reshape(nk,nb)
                occupations[:,1:3]=[1.9,.1]
            else:
                orbitals=np.frombuffer(data,'f8',count=2*ng*nb*nk,offset=head+3*nk*8+2*nb*nk*8)
                orbitals[-1 if invalid=='nan-last-state' else 0]=np.nan if invalid=='nan-last-state' else np.inf
            (path/'pipeline.bin').write_bytes(data)
            run_mpi(exe,ranks,threads,path,'reject-scalar-'+invalid,
                    extra_args=('1','insulating-distributed'),expect_failure=True)
            print('SCALAR_COLLECTIVE_REJECTION_PASS',invalid,ranks,flush=True)
        # Reject a vertex outside the finite-grid Nyquist range from inside
        # the persistent team, then terminate every rank without a barrier hang.
        data=bytearray(original)
        np.frombuffer(data,'f8',count=1,offset=12*4)[0]=3.0
        (path/'pipeline.bin').write_bytes(data)
        run_mpi(exe,ranks,threads,path,'reject-outer-vertex',timeout=30,
                extra_args=('1','hybrid-distributed'),expect_failure=True)
        print('OUTER_VERTEX_COLLECTIVE_REJECTION_PASS',ranks,threads,flush=True)
        (path/'pipeline.bin').write_bytes(original)
        if threads>1:
            run_mpi(exe,ranks,threads,path,'reject-mpi-thread-single',timeout=30,
                    extra_args=('1','hybrid-distributed','mpi-single'),expect_failure=True)
            print('MPI_THREAD_SINGLE_REJECTION_PASS',ranks,threads,flush=True)
        run_mpi(exe,ranks,threads,path,'reject-short-team',timeout=30,
                extra_args=('1','hybrid-distributed','short-team'),expect_failure=True)
        print('OMP_SHORT_TEAM_COLLECTIVE_REJECTION_PASS',ranks,threads,flush=True)
        print('SCALAR_MATRIX_FREE_SCREEN_PASS')


if __name__=='__main__':main()
