"""Native Gamma FFT/optical equivalence, cache and collective-failure regression."""
from pathlib import Path
import hashlib,sys
from mpi_test_support import shared_fixture,run_mpi
exe=Path(sys.argv[1]);ranks=int(sys.argv[2]) if len(sys.argv)>2 else 4
threads=int(sys.argv[3]) if len(sys.argv)>3 else 12
with shared_fixture() as work:
    cache=work/'basis.bin'
    for mode in ['reference','cold','warm','changed','write-exists','bad-cache-rank','bad-cache-cell','nan-state','inf-state','bad-gram','fail-cg','maxiter']:
        before=hashlib.sha256(cache.read_bytes()).hexdigest() if cache.exists() else None
        run_mpi(exe,ranks,threads,work,mode,timeout=120,extra_args=(str(cache),mode))
        for rank in range(ranks):
            assert f'GAMMA_PDEP_MPI_PASS rank={rank} case={mode}' in (work/mode/f'stdout.rank-{rank}').read_text()
        if before is not None:assert hashlib.sha256(cache.read_bytes()).hexdigest()==before
    corrupt=work/'corrupt.bin';data=bytearray(cache.read_bytes());data[:4]=b'FAIL';corrupt.write_bytes(data)
    run_mpi(exe,ranks,threads,work,'bad-cache',timeout=120,extra_args=(str(corrupt),'bad-cache'))
    for rank in range(ranks):assert 'GAMMA_PDEP_MPI_PASS' in (work/'bad-cache'/f'stdout.rank-{rank}').read_text()
print('GAMMA_PDEP_FFT_HEAD_WINGS_CACHE_COLLECTIVE_PASS',ranks,threads)
