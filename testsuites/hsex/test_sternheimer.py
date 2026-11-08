"""Compare both microscopic screening paths on a complete synthetic band sum."""
from pathlib import Path
import sys,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from prepare_w import prepare


def run():
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp);restart=p/'data_for_restart';restart.mkdir()
        n=4;nk=8;ng=n**3;nv=2;cell=np.eye(3)*7
        k=(np.indices((2,2,2)).reshape(3,-1).T-.5)*np.pi/7
        rng=np.random.default_rng(4291)
        eig=np.tile(np.linspace(-.9,2,ng),(nk,1))+.005*rng.random((nk,ng))
        occ=np.tile(np.r_[2,2,np.zeros(ng-2)],(nk,1))
        with (restart/'kli-v1.bin').open('wb') as f:
            f.write(b'KLI1');np.array([n,n,n,1,nk,ng,1],'i4').tofile(f)
            for a in [cell,k,occ,np.zeros(3),np.zeros(ng),[0]]: np.asarray(a,'f8').tofile(f)
        psi=[]
        for j in range(nk):
            u=np.linalg.qr(rng.normal(size=(ng,ng))+1j*rng.normal(size=(ng,ng)))[0]
            psi.append(u.T*np.sqrt(ng/343))
            with (p/f'kli-hamiltonian-{j+1:06d}.bin').open('wb') as f:
                f.write(b'KSH1');np.array([n,n,n,j+1],'i4').tofile(f)
                for a in [k[j],cell,[0],np.zeros(ng)]:np.asarray(a,'f8').tofile(f)
                ((u*eig[j])@u.conj().T).T.tofile(f)
        psi=np.array(psi);psi.tofile(restart/'wfn.bin')
        np.savetxt(p/'Si_eigen.data',np.c_[np.tile(np.arange(ng)+1,nk),eig.ravel(),occ.ravel()],fmt=['%d','%.17e','%.1f'])
        prepare(p,p/'sum',[2]*3,2,.48,ng-nv)
        # Non-occupied grid states may not be accessed by the Sternheimer backend.
        psi[:,nv:]=np.nan;psi.tofile(restart/'wfn.bin')
        prepare(p,p/'stern',[2]*3,2,.48,hamiltonians=p)
        count=2**6*nk
        a,b=[np.frombuffer((p/f).read_bytes()[-count*8:],'f8') for f in ['sum','stern']]
        np.testing.assert_allclose(a,b,rtol=3e-12,atol=3e-14)
        print('max |W(stern)-W(all bands)|:',abs(a-b).max())
        # Identical synthetic grid operator with a TB-mBJ reference must give
        # the same screening values, but a distinct RT compatibility tag.
        old=restart/'kli-v1.bin'
        old.rename(restart/'tbmbj-wref.bin')
        for path,tag in [(restart/'tbmbj-wref.bin',b'TBM1')]+[(f,b'TSH1') for f in p.glob('kli-hamiltonian-*.bin')]:
            with path.open('r+b') as stream:stream.write(tag)
        prepare(p,p/'tb-stern',[2]*3,2,.48,hamiltonians=p)
        tb=(p/'tb-stern').read_bytes()
        assert tb[:4]==b'HXT1'
        np.testing.assert_allclose(np.frombuffer(tb[-count*8:],'f8'),a,rtol=3e-12,atol=3e-14)
        # Reuse must reject a genuinely TR-breaking reference.
        try:
            prepare(p,p/'wrong',[2]*3,2,.48,hamiltonians=p,time_reversal=True)
        except ValueError as error:
            assert 'time reversal' in str(error)
        else:
            raise AssertionError('accepted TR-breaking occupied manifold')
        # Construct exact conjugate pairs with a different occupied gauge.
        # Phases should not affect reduction; the q-grid includes Nyquist faces.
        for j in range(nk//2):
            jm=nk-1-j
            psi[jm,:nv]=psi[j,:nv].conj()*np.exp(1j*rng.random((nv,1)))
            eig[jm]=eig[j]
            h=p/f'kli-hamiltonian-{j+1:06d}.bin'
            hm=p/f'kli-hamiltonian-{jm+1:06d}.bin'
            offset=h.stat().st_size-16*ng*ng
            with hm.open('r+b') as f:
                f.seek(offset);np.fromfile(h,'c16',offset=offset).conj().tofile(f)
        np.savetxt(p/'Si_eigen.data',np.c_[np.tile(np.arange(ng)+1,nk),eig.ravel(),occ.ravel()],fmt=['%d','%.17e','%.1f'])
        psi.tofile(restart/'wfn.bin')
        for name,tr in [('full',False),('paired',True)]:
            prepare(p,p/name,[2]*3,2,.48,hamiltonians=p,time_reversal=tr)
        full,paired=[np.frombuffer((p/name).read_bytes()[-count*8:],'f8') for name in ('full','paired')]
        np.testing.assert_allclose(full,paired,rtol=3e-12,atol=3e-14)




if __name__=='__main__':run()
