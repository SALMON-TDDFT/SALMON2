"""Fortran finite-q PDEP collective failures and nv=1 bounds regression."""
from pathlib import Path
import hashlib
import sys
from mpi_test_support import shared_fixture,run_mpi

exe=Path(sys.argv[1])
ranks=int(sys.argv[2]) if len(sys.argv)>2 else 4
threads=int(sys.argv[3]) if len(sys.argv)>3 else 12
with shared_fixture() as path:
    reference=path/'reference.bin'
    cases=[('reference',-1,'valid',reference),('rank1',1,'valid',reference)]
    cases += [(case,-1,case,reference) for case in ('bad-gap','bad-gram','nan-state','inf-state','fail-cg','io-fail','gamma')]
    cases += [('fail-cg-pdep',1,'fail-cg',reference)]
    for label,rank,mode,ref in cases:
        before=hashlib.sha256(reference.read_bytes()).hexdigest() if reference.exists() else None
        run_mpi(exe,ranks,threads,path,label,timeout=30,extra_args=(str(ref),str(rank),mode))
        for process in range(ranks):
            log=(path/label/f'stdout.rank-{process}').read_text()
            assert f'PDEP_MPI_RETURN_PASS rank={process} case={mode}' in log
        if before is not None:
            assert hashlib.sha256(reference.read_bytes()).hexdigest()==before
    bad=path/'bad-reference.bin'
    data=bytearray(reference.read_bytes());data[:4]=b'FAIL';bad.write_bytes(data)
    run_mpi(exe,ranks,threads,path,'bad-reference',timeout=30,
            extra_args=(str(bad),'1','bad-reference'))
    for process in range(ranks):
        assert 'PDEP_MPI_RETURN_PASS' in (path/'bad-reference'/f'stdout.rank-{process}').read_text()
print('PDEP_Q_COLLECTIVE_IEEE_NV1_PASS',ranks,threads,flush=True)
