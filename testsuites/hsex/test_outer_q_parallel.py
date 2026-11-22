"""Offline matched-backend comparison; retain no synthetic/physical GS payloads."""
from pathlib import Path
import sys,re,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tools/hsex'),str(Path(__file__).parent)]
from test_prepare_distributed import fixture
from sternheimer import equilibrium_operator
from mpi_test_support import shared_fixture,run_mpi

old,new=map(lambda x:Path(x).resolve(),sys.argv[1:3])
ranks=int(sys.argv[3]);threads=int(sys.argv[4]);output=Path(sys.argv[5]).resolve()
hashes={str(exe):hashlib.sha256(exe.read_bytes()).hexdigest() for exe in [old,new]}
records={}
with shared_fixture() as path:
    fixture(path)
    side=2 if '--small-mesh' in sys.argv[6:] else 4
    n=np.array([4]*3);mesh=np.array([side]*3);nk=side**3;ng=64;nb=3;cell=np.eye(3)*9
    frac=(np.indices((side,)*3).reshape(3,-1)[::-1].T-(side-1)/2)/side;k=frac*2*np.pi/9
    eig=np.tile(np.r_[-.8,-.5,np.linspace(.2,3,ng-2)],(nk,1))
    psi=np.fromfile(path/'data_for_restart/wfn.bin','c16').reshape(64,ng,ng)[:nk]*np.sqrt(729/ng)
    reference_k=(np.indices((4,)*3).reshape(3,-1)[::-1].T-1.5)/4*2*np.pi/9
    h=np.array([equilibrium_operator(path,j,n,reference_k,cell,np.array([0.]),np.zeros(ng),b'KLI1') for j in range(nk)])
    for pump in [False,True]:
        eta=.001 if pump else 0.;occ=np.tile([2.,2.,0.],(nk,1))
        with (path/'pipeline.bin').open('wb') as f:
            np.r_[n,mesh,mesh,2,2,nb].astype('i4').tofile(f)
            np.r_[.2,eta,cell.T.ravel()].astype('f8').tofile(f)
            for a in [k,eig[:,:nb],occ]:np.asarray(a,dtype='f8').tofile(f)
            psi[:,:nb].astype('c16').tofile(f);h.transpose(0,2,1).astype('c16').tofile(f)
        reference=None;case={}
        for trial in range(3):
            for mode,exe in [('before',old),('outer',new)]:
                assert hashlib.sha256(exe.read_bytes()).hexdigest()==hashes[str(exe)]
                label=f'{pump}-{trial}-{mode}'
                run_mpi(exe,ranks,threads,path,label,timeout=180,extra_args=('1','dense-distributed'))
                w=np.fromfile(path/'pipeline-output.bin','f8');assert np.isfinite(w).all()
                if reference is None:reference=w
                else:np.testing.assert_allclose(w,reference,atol=2e-12,rtol=2e-11)
                log=(path/label/'stdout.rank-0').read_text()
                seconds=list(map(float,re.search(r'HSEX_W_SECONDS[^:]+:\s+([^\n]+)',log)[1].split()))
                regions=re.search(r'HSEX_W_Q_PARALLEL regions=(\d+) batches=(\d+)',log)
                if mode=='outer':
                    assert regions and int(regions[1])==(nk+ranks-1)//ranks
                    assert int(regions[2])>int(regions[1])
                case[label]={'seconds':seconds,'whole_W_seconds':sum(seconds),
                             'max_W_error':float(np.max(abs(w-reference))),
                             'regions_batches':list(map(int,regions.groups())) if regions else None}
        before=np.median([v['whole_W_seconds'] for name,v in case.items() if name.endswith('before')])
        outer=np.median([v['whole_W_seconds'] for name,v in case.items() if name.endswith('outer')])
        case['median_whole_W_speedup']=float(before/outer)
        records[str(pump)]=case
records['executable_sha256']=hashes;records['MPI']=ranks;records['OMP']=threads
output.write_text(json.dumps(records,indent=2)+'\n')
print('OUTER_Q_MATCHED_BENCHMARK_PASS',output,flush=True)
