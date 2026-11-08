"""Independent Python comparisons for the native Fortran screening pipeline."""
from pathlib import Path
import sys,subprocess,tempfile,os
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/hsex'))
from box_factors import box_forms
from metal_head import metallic_cell,check_eta
from q_mesh import signed_indices

def main():
    exe=Path(sys.argv[1])
    rng=np.random.default_rng(345)
    with tempfile.TemporaryDirectory() as td:
        d=Path(td);n=np.array([6,6,6]);boxes=3;ng=5
        q=rng.normal(size=(ng,3));q[0]=0
        rec=np.eye(3)*.7;vol=(2*np.pi/.7)**3;nk=8
        a=rng.normal(size=(ng,ng))+1j*rng.normal(size=(ng,ng));chi=a@a.conj().T
        with (d/'input.bin').open('wb') as f:
            np.r_[n,boxes,ng,nk].astype('i4').tofile(f)
            np.r_[vol,rec.T.ravel(),q.ravel()].astype('f8').tofile(f)
            chi.T.astype('c16').tofile(f)
        subprocess.run([str(exe),'assembly'],cwd=d,check=True)
        actual=np.fromfile(d/'output.bin','c16').reshape((boxes**3,boxes**3),order='F')
        q2=np.sum((q@rec)**2,axis=1);bare=np.divide(4*np.pi,q2,out=np.zeros(ng),where=q2>1e-14)
        b=np.sqrt(bare);eps=np.eye(ng)+2/(nk*vol)*b[:,None]*chi*b[None,:]
        forms=box_forms(n,boxes,q);expected=forms.conj()@(b[:,None]*np.linalg.inv(eps)*b[None,:])@forms.T/vol
        np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=2e-12)
        print('ASSEMBLY_MAX_ERROR',np.max(abs(actual-expected)))
        invalid=bytearray((d/'input.bin').read_bytes())
        np.frombuffer(invalid,dtype='i4',count=6)[3]=0
        (d/'input.bin').write_bytes(invalid)
        rejected=subprocess.run([str(exe),'assembly'],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        assert rejected.returncode==0, (rejected.returncode,rejected.stdout,rejected.stderr)
        assert b'SCREEN_REJECTED status=1' in rejected.stdout.splitlines(), (rejected.stdout,rejected.stderr)

    print('FORTRAN_SCREEN_PASS')
if __name__=='__main__':main()
