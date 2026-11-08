"""Whole-file screening comparison across q ownership and pair-distributed I/O."""
from pathlib import Path
import argparse,json,os,subprocess,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/hsex'))
from prepare_w import prepare


def fixture(path,n=4,kmesh=4):
    restart=path/'data_for_restart';restart.mkdir()
    nk=kmesh**3;ng=n**3;nv=2;cell=np.eye(3)*9
    frac=(np.indices((kmesh,)*3).reshape(3,-1)[::-1].T-(kmesh-1)/2)/kmesh
    k=frac*2*np.pi/9
    eig=np.tile(np.r_[-.8,-.5,np.linspace(.2,3,ng-2)],(nk,1))
    occ=np.tile(np.r_[2.,2.,np.zeros(ng-2)],(nk,1))
    with (restart/'kli-v1.bin').open('wb') as f:
        f.write(b'KLI1');np.array([n,n,n,1,nk,ng,1],dtype='i4').tofile(f)
        for a in (cell,k,occ,np.zeros(3),np.zeros(ng),[0.]):np.asarray(a,dtype='f8').tofile(f)
    rng=np.random.default_rng(5202)
    with (restart/'wfn.bin').open('wb') as wf:
        for j in range(nk):
            u=np.linalg.qr(rng.normal(size=(ng,ng))+1j*rng.normal(size=(ng,ng)))[0]
            states=u.T*np.sqrt(ng/729)
            states.tofile(wf);h=(u*eig[j])@u.conj().T
            with (path/f'kli-hamiltonian-{j+1:06d}.bin').open('wb') as f:
                f.write(b'KSH1');np.array([n,n,n,j+1],dtype='i4').tofile(f)
                for a in (k[j],cell,[0.],np.zeros(ng)):np.asarray(a,dtype='f8').tofile(f)
                h.T.tofile(f)
            with (path/f'hsex-pump-{j+1:06d}.bin').open('wb') as f:
                f.write(b'HSP1');np.array([n,n,n,nk,j+1,3,80],dtype='i4').tofile(f)
                for a in (cell,k[j],np.zeros(3),eig[j,:3],occ[j,:3]):np.asarray(a,dtype='f8').tofile(f)
                states[:3].tofile(f);h.T.tofile(f)
    np.savetxt(path/'Si_eigen.data',np.c_[np.tile(np.arange(ng)+1,nk),eig.ravel(),occ.ravel()],fmt=['%d','%.17e','%.1f'])


def run(comm=None,baseline=None):
    rank=0 if comm is None else comm.rank
    temporary=tempfile.TemporaryDirectory(prefix='w-distribution-check-') if rank==0 else None
    path=Path(temporary.name) if rank==0 else None
    if rank==0:
        fixture(path)
        for kind in ('equilibrium','pump'):
            for q in (2,4):
                output=path/f'ref-{kind}-{q}'
                if baseline is None:
                    args=dict(hamiltonians=path) if kind=='equilibrium' else dict(pump=path,eta=.001)
                    prepare(path,output,[4]*3,2,.2,qmesh=[q]*3,**args)
                else:
                    args=['--hamiltonians',str(path)] if kind=='equilibrium' else ['--pump',str(path),'--eta','.001']
                    subprocess.run([sys.executable,str(baseline/'prepare_w.py'),str(path),str(output),
                        '--mesh','4','4','4','--q-mesh',*[str(q)]*3,'--boxes','2','--cutoff','.2',*args],check=True)
    if comm is not None:path=Path(comm.bcast(str(path) if rank==0 else None,root=0))
    for kind in ('equilibrium','pump'):
        for q in (2,4):
            output=path/f'actual-{kind}-{q}'
            args=dict(hamiltonians=path) if kind=='equilibrium' else dict(pump=path,eta=.001)
            prepare(path,output,[4]*3,2,.2,qmesh=[q]*3,comm=comm,**args,**(dict(q_kernel=os.environ["HSEX_Q_KERNEL"],q_threads=int(os.environ.get("HSEX_Q_THREADS","2"))) if os.environ.get("HSEX_Q_KERNEL") else {}))
            if rank==0:
                a=(path/f'ref-{kind}-{q}').read_bytes();b=output.read_bytes()
                assert a[:40]==b[:40] and len(a)==len(b)
                x,y=[np.frombuffer(v,dtype='f8',offset=40) for v in (a,b)]
                np.testing.assert_allclose(x,y,rtol=2e-11,atol=2e-12)
                print(kind,q,'whole-file maxdiff',np.max(abs(x-y)),flush=True)
            if comm is not None:comm.Barrier()
    # Root-only population validation must fail collectively, not strand peers.
    if rank==0:
        with (path/'hsex-pump-000001.bin').open('r+b') as f:
            f.seek(176);np.asarray([2.,1.9,.1],dtype='f8').tofile(f)
    if comm is not None:comm.Barrier()
    rejected=False
    try:prepare(path,path/'must-not-exist',[4]*3,2,.2,qmesh=[2]*3,comm=comm,pump=path,eta=.001)
    except ValueError as exc:
        rejected='demonstrated insulating gap' in str(exc)
    assert rejected and not (path/'must-not-exist').exists()
    if comm is not None:comm.Barrier()
    if rank==0:temporary.cleanup()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mpi',action='store_true');p.add_argument('--baseline',type=Path)
    a=p.parse_args();comm=None
    if a.mpi:
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
    try:run(comm,a.baseline)
    except BaseException:
        import traceback
        traceback.print_exc()
        if comm is not None:comm.Abort(1)
        raise
