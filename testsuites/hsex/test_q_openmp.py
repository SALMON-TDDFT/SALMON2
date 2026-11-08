"""Optional native q backend: compare independent SciPy reference kernels."""
import os
import sys
import subprocess
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from q_openmp import QKernel
from sternheimer import resolvent_gram
from pump_screen import broadened_gram

@unittest.skipUnless(os.environ.get('HSEX_Q_KERNEL'),'set HSEX_Q_KERNEL to an audited native library')
class QOpenMP(unittest.TestCase):
    @unittest.skipIf(int(os.environ.get('HSEX_Q_MIN_THREADS','1'))>=12,'Fugaku minimum-thread policy')
    def test_reject_reduced_team(self):
        code="""
import os
import numpy as np
from q_openmp import QKernel
k=QKernel(os.environ['HSEX_Q_KERNEL'])
try:
    k.solve(np.eye(2),np.ones(2),np.zeros(1),np.empty((2,0)),
            np.ones((2,1,1)),np.ones((2,1,1)),[1],np.zeros((1,1)),
            np.ones((1,1)),[0.],12,False)
except ValueError as exc:
    assert 'differs from requested 12' in str(exc), str(exc)
else:
    raise AssertionError('reduced team must be rejected')
"""
        env=os.environ.copy();env['OMP_THREAD_LIMIT']='2';env['OMP_DYNAMIC']='FALSE'
        env['PYTHONPATH']=str(Path(__file__).resolve().parents[2]/'tools/hsex')+os.pathsep+env.get('PYTHONPATH','')
        subprocess.run([sys.executable,'-c',code],env=env,check=True)

    def test_responses(self):
        rng=np.random.default_rng(982)
        n,nb,nq,mg=24,3,25,5
        adj=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
        u=np.linalg.qr(rng.normal(size=(n,2))+1j*rng.normal(size=(n,2)))[0]
        psi=rng.normal(size=(n,nb,nq))+1j*rng.normal(size=(n,nb,nq))
        phase=np.exp(1j*rng.normal(size=(n,mg,nq)))
        counts=np.resize(np.array([1,5,3,2,4,5,1]),nq)
        levels=-rng.random((nb,nq));occ=rng.random((nb,nq))*2;occ[1,2]=0;occ[2,1]=-1e-12
        d=np.linspace(1,3,n); e=np.full(n-1,.15)
        kernel=QKernel(os.environ['HSEX_Q_KERNEL'])
        for pump in (False,True):
            if pump:levels=np.tile(d[[5,11,18],None],(1,nq))
            etas=np.array([.0002,.0001]) if pump else np.array([0.])
            projector=np.empty((n,0),complex) if pump else u
            expected=[]
            for q in range(nq):
                g=counts[q];total=np.zeros((len(etas),g,g),complex)
                for band in range(nb):
                    b=psi[:,band,q,None]*phase[:,:g,q]
                    if not pump:b-=u@(u.conj().T@b)
                    b=adj@b
                    for c,eta in enumerate(etas):
                        term,_=broadened_gram(b,d,e,levels[band,q],eta) if pump else resolvent_gram(b,d,e,levels[band,q])
                        total[c]+=occ[band,q]*.5*term
                expected.append(total)
            for threads in sorted(set((max(1,int(os.environ.get('HSEX_Q_MIN_THREADS','1'))),12))):
                got,res,team=kernel.solve(adj,d,e,projector,psi,phase,counts,levels,occ,etas,threads,pump)
                self.assertEqual(team,threads)
                self.assertLess(res,1e-8)
                for a,b in zip(got,expected):np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-12)
        # Indefinite equilibrium operator must fail, like solveh_banded.
        with self.assertRaises(ValueError):
            kernel.solve(adj,-d,e,u,psi,phase,counts,levels,occ,np.array([0.]),max(2,int(os.environ.get("HSEX_Q_MIN_THREADS","1"))),False)
        from q_openmp import response_batches
        loads=[]
        def load(q):
            loads.append(q)
            return psi[:,:,q],phase[:,:counts[q],q],levels[:,q],occ[:,q]
        batched,res,team=response_batches(kernel,adj,d,e,np.empty((n,0),complex),nq,load,max(4,int(os.environ.get("HSEX_Q_MIN_THREADS","1"))),True,[.0002,.0001],budget=1)
        self.assertEqual(loads,list(range(nq)))
        for a,b in zip(batched,expected):np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-12)
        for invalid in (np.nan,np.inf):
            bad=psi.copy();bad[0,0,0]=invalid
            with self.assertRaises(ValueError):kernel.solve(adj,d,e,u,bad,phase,counts,levels,occ,np.array([0.]),max(2,int(os.environ.get("HSEX_Q_MIN_THREADS","1"))),False)

if __name__=='__main__':unittest.main()
