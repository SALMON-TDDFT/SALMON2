"""Metallic full-matrix head and rejection checks against Python reference."""
from pathlib import Path
import sys,tempfile,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/hsex'))
from metal_head import metallic_cell
from q_mesh import signed_indices

def main():
 exe=Path(sys.argv[1]);rng=np.random.default_rng(55)
 qm=np.array([8]*3);mesh=np.array([16]*3);na=4
 rec=np.array([[.7,.04,.02],[0,.8,.03],[0,0,.9]])
 q=signed_indices(qm)/qm@rec.T
 a=rng.normal(size=(na,na));base=a@a.T+np.eye(na)*5
 v=np.array([base+(z@z)*np.eye(na)*.03+z[0]*z[1]*np.ones((na,na))*.01 for z in q],complex)
 with tempfile.TemporaryDirectory() as td:
  p=Path(td)
  for kind in ('valid','nan','infinity','nonsmooth','singular-reciprocal'):
   values=v.copy();test_rec=rec.copy()
   if kind=='singular-reciprocal':
    test_rec=np.eye(3)*1e-7
    values[:]=base
   if kind=='nan':values[1,0,0]=np.nan
   if kind=='infinity':values[1,0,0]=np.inf
   if kind=='nonsmooth':values[1]*=3
   with (p/'head.bin').open('wb') as f:
    np.r_[qm,mesh,na].astype('i4').tofile(f);test_rec.T.astype('f8').tofile(f);np.array([.05]).tofile(f)
    values.transpose(0,2,1).astype('c16').tofile(f)
   subprocess.run([str(exe)],cwd=p,check=True)
   data=(p/'head-output.bin').read_bytes();status=np.frombuffer(data,'i4',count=1)[0]
   if kind=='valid':
    expected=metallic_cell(values,qm,rec.T,mesh)[:2]
    assert status==0
    actual=np.frombuffer(data,'c16',offset=4).reshape(2,na,na).transpose(0,2,1)
    np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=2e-12)
    print('HEAD_MAX_ERROR',np.max(abs(actual-np.asarray(expected))))
   else:assert status!=0,kind
 print('FORTRAN_SCREEN_HEAD_PASS')
if __name__=='__main__':main()
