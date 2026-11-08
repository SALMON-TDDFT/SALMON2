"""Compile the actual projected HSEX velocity and compare an energy derivative."""
from pathlib import Path
import os,subprocess,tempfile,shlex
import numpy as np
root=Path(__file__).resolve().parents[2]

def run():
    rng=np.random.default_rng(4051); na,nb=5,3
    rand=lambda shape:rng.normal(size=shape)+1j*rng.normal(size=shape)
    p=rand((na,nb));rp=rand((na,nb,3));s=rand((na,na));s=s+s.conj().T
    ds=rand((na,na,3));ds=ds+ds.transpose(1,0,2).conj()
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        for name,v in [('p',p),('rp',rp),('s',s),('ds',ds)]:v.ravel(order='F').astype('c16').tofile(tmp/name)
        driver=tmp/'test.f90'
        driver.write_text('''program test
use hsex_velocity, only: hsex_velocity_diagonal
implicit none
complex(8)::p(5,3),rp(5,3,3),s(5,5),ds(5,5,3)
real(8)::v(3,3)
open(10,file='p',access='stream',form='unformatted');read(10)p;close(10)
open(10,file='rp',access='stream',form='unformatted');read(10)rp;close(10)
open(10,file='s',access='stream',form='unformatted');read(10)s;close(10)
open(10,file='ds',access='stream',form='unformatted');read(10)ds;close(10)
call hsex_velocity_diagonal(p,rp,s,ds,v)
open(10,file='v',access='stream',form='unformatted');write(10)v;close(10)
end program
''')
        checker=root/'tools/check_forbidden_eigensolvers.py'
        subprocess.run(['python3',str(checker),'--source-root',str(root/'src')],check=True)
        source=root/'src/xc/hsex_velocity.f90'
        assert source.exists(), 'missing projected HSEX velocity implementation'
        fc=os.environ.get('FC','gfortran')
        subprocess.run([fc,'-c','-fcheck=all',str(source),str(driver)],cwd=tmp,check=True)
        subprocess.run(['python3',str(checker),'--object-tree',str(tmp)],check=True)
        subprocess.run([fc,'hsex_velocity.o','test.o','-o','test',*shlex.split(os.environ.get('HSEX_BLAS_FLAGS','-lblas'))],cwd=tmp,check=True)
        subprocess.run([str(tmp/'test')],cwd=tmp,check=True)
        actual=np.fromfile(tmp/'v').reshape((3,nb),order='F')
        expected=np.empty_like(actual)
        for axis in range(3):
            def energy(h):
                ph=p+1j*h*rp[:,:,axis]
                return np.einsum('ib,ij,jb->b',ph.conj(),s+h*ds[:,:,axis],ph).real
            expected[axis]=(energy(1e-6)-energy(-1e-6))/2e-6
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=2e-8)
        # NaN/Infinity must reject, not silently contaminate projected current.
        for invalid in (np.nan,np.inf):
            bad=p.copy();bad[-1,-1]=invalid;bad.ravel(order='F').tofile(tmp/'p')
            result=subprocess.run([str(tmp/'test')],cwd=tmp,capture_output=True)
            assert result.returncode != 0 and b'nonfinite' in result.stderr
        print('projected HSEX velocity energy derivative / nonfinite rejection PASS')

if __name__=='__main__':run()
