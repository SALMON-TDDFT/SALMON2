"""RPA regression on an inhomogeneous, insulating synthetic grid checkpoint."""
from pathlib import Path
import sys,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from prepare_w import prepare


def test():
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp); restart=root/'data_for_restart';restart.mkdir()
        n=6; nk=8; nb=3; cell=8*np.eye(3);ng=n**3
        frac=(np.indices((2,2,2)).reshape(3,-1).T-.5)/2
        k=2*np.pi*frac/8
        with (restart/'kli-v1.bin').open('wb') as f:
            f.write(b'KLI1');np.array([n,n,n,1,nk,nb,1],dtype='i4').tofile(f)
            for a in [cell,k,np.tile([2,0,0],(nk,1)),np.zeros(3),np.zeros(ng),[0]]:
                np.asarray(a,dtype='f8').tofile(f)
        with (root/'Si_eigen.data').open('w') as f:
            for ik in range(nk):
                for ib,energy in enumerate([-0.3,0.1,0.2]):
                    print(ib+1,energy,2*(ib==0),file=f)
        rng=np.random.default_rng(322)
        grid=np.indices((n,n,n)).reshape(3,-1)
        envelope=np.exp(-.3*((grid-1)**2).sum(0))
        psi=[]
        for ik in range(nk):
            # Time-reversed symmetry is not needed to make static W real.
            v=(rng.normal(size=(ng,nb))+1j*rng.normal(size=(ng,nb)))*envelope[:,None]
            psi.append(np.linalg.qr(v)[0].T*np.sqrt(ng/512))
        psi=np.array(psi)
        p=restart/'wfn.bin'
        def run(name,phase):
            (psi*phase).tofile(p);file=root/name
            prepare(root,file,[2]*3,2,.50,2)
            return np.fromfile(file,dtype='f8',offset=file.stat().st_size-8*8**2*nk)
        first=run('W',1.); second=run('W-gauge',np.exp(2j*np.pi*rng.random((nk,nb,1))))
        np.testing.assert_allclose(first,second,atol=1e-13,rtol=0)
        assert np.ptp(first)>0.01

if __name__=='__main__':test()
