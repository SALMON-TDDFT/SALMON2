"""Fixed-target folding/compression and owner-only collective failure test."""
from pathlib import Path
import re,sys
import itertools
import numpy as np
from mpi_test_support import shared_fixture,run_mpi
exe=Path(sys.argv[1]);ranks=int(sys.argv[2]);threads=int(sys.argv[3])
with shared_fixture() as work:
    for mode in ['valid','fail','tiny-cutoff']:
        run_mpi(exe,ranks,threads,work,mode,timeout=120,extra_args=(mode,))
        texts=[(work/mode/f'stdout.rank-{rank}').read_text() for rank in range(ranks)]
        for rank,text in enumerate(texts):
            assert f'REDUCED_Q_MPI_PASS rank={rank} case={mode}' in text
        if mode=='valid':
            joined='\n'.join(texts)
            assert 'converged=1' in joined
            errors=re.findall(r'HSEX_REDUCED_TARGET_CHI_ERROR q/relative/norm\s+\d+\s+(\S+)',joined)
            assert len(errors)==4 and max(map(float,errors))<1e-9
            # Independent Fourier-diagonal model checks actual source/folding,
            # rather than comparing two solvers with shared vertex construction.
            rows=re.findall(r'HSEX_REDUCED_Q target/q/source/NG/qfrac\s+(.*)',joined)
            norms=re.findall(r'HSEX_REDUCED_TARGET_CHI_ERROR q/relative/norm\s+\d+\s+\S+\s+(\S+)',joined)
            wrapped=False
            for row,actual_norm in zip(rows,norms):
                fields=row.split();target,qid,source,ng=map(int,fields[:4]);q=np.array(fields[4:],float)
                order=np.array([(target-1)%4,((target-1)//4)%4,(target-1)//16])
                source_order=(order-np.rint(4*q).astype(int))%4
                expected_source=1+source_order[0]+4*source_order[1]+16*source_order[2]
                assert source==expected_source
                centered=lambda x:np.where(2*x>=4,x-4,x)/4
                wrap=np.rint(centered(source_order)+q-centered(order)).astype(int)
                wrapped |= bool(np.any(wrap))
                vectors=[np.array(g) for g in itertools.product(range(-2,3),repeat=3)
                         if np.dot(np.array(g)+q,np.array(g)+q)*(2*np.pi/10)**2/2<=.12]
                assert len(vectors)==ng
                diagonal=[]
                for g in vectors:
                    total=0.
                    for band,level in enumerate([-2.,-1.]):
                        mode=(g+wrap+np.array([band,0,0]))%6
                        index=1+mode[0]+6*mode[1]+36*mode[2]
                        if index in [1,2]:continue
                        total+=1/(1+.01*index+.001*target-level-.001*source)
                    diagonal.append(total)
                np.testing.assert_allclose(float(actual_norm),np.linalg.norm(diagonal),rtol=2e-10,atol=2e-12)
            assert wrapped, 'fixture did not exercise k+q folding'
            print('REDUCED_Q_INDEPENDENT_FOLD_RESPONSE_PASS')
            counts=re.search(r'HSEX_REDUCED_COUNTS CG_RHS/CG_H/shared_rank/shared_H/validation_H\s+(.*)',joined)
            rhs,cg,rank,shared,validation=map(int,counts[1].split())
            assert shared==rank+validation and validation==rhs
            print('REDUCED_Q_MPI_COUNTS',rhs,cg,rank,shared,validation)
print('REDUCED_Q_MPI_COLLECTIVE_PASS',ranks,threads)
