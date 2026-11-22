"""Little-group symmetry plus common response space against all-q CG."""
from pathlib import Path
import re,sys
from mpi_test_support import shared_fixture,run_mpi
exe=Path(sys.argv[1]);ranks=int(sys.argv[2]);threads=int(sys.argv[3])
with shared_fixture() as work:
    for mode in ['valid','broken-h','broken-atoms','fail']:
        run_mpi(exe,ranks,threads,work,mode,timeout=120,extra_args=(mode,))
        texts=[(work/mode/f'stdout.rank-{rank}').read_text() for rank in range(ranks)]
        for rank,text in enumerate(texts):assert f'RESPONSE_STAR_MPI_PASS rank={rank} case={mode}' in text
        if mode=='fail':continue
        joined='\n'.join(texts)
        sizes=re.search(r'HSEX_REDUCED_SYMMETRY_COUNTS operations/q/full_RHS/rep_q/rep_RHS\s+(.*)',joined)
        operations,q,rhs,rep_q,rep_rhs=map(int,sizes[1].split())
        assert (q,rhs)==(6,12)
        expected={'valid':(6,2,4),'broken-h':(2,4,8),'broken-atoms':(1,6,12)}[mode]
        assert (operations,rep_q,rep_rhs)==expected,(mode,operations,rep_q,rep_rhs)
        errors=re.findall(r'HSEX_REDUCED_TARGET_CHI_ERROR q/relative/norm\s+\d+\s+(\S+)',joined)
        exact=re.findall(r'HSEX_REDUCED_SYMMETRY_CHI_ERROR q/relative\s+\d+\s+(\S+)',joined)
        assert len(errors)==6 and max(map(float,errors))<1e-9
        assert len(exact)==6 and max(map(float,exact))<1e-9
        assert 'converged=1' in joined
        print('RESPONSE_STAR_CASE_PASS',mode,operations,rep_q,rep_rhs,max(map(float,errors)))
print('RESPONSE_STAR_MPI_COLLECTIVE_PASS',ranks,threads)
