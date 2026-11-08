"""Portable MPI test launches and failure-preserving shared scratch fixtures."""
from contextlib import contextmanager
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile


@contextmanager
def shared_fixture():
    # Fugaku jobs must set TMPDIR to a directory on their shared filesystem.
    path=Path(tempfile.mkdtemp(prefix='hsex-native-',dir=os.environ.get('TMPDIR'))).resolve()
    print('MPI_TEST_FIXTURE',path,flush=True)
    try:
        yield path
    except BaseException:
        print('MPI_TEST_FAILURE_FIXTURE_PRESERVED',path,flush=True)
        raise
    else:
        if os.environ.get('HSEX_KEEP_TEST_FIXTURES')=='1':
            print('MPI_TEST_FIXTURE_PRESERVED',path,flush=True)
        else:
            shutil.rmtree(path)


def run_mpi(executable,ranks,threads,path,label,timeout=180,expect_failure=False,extra_args=()):
    path=Path(path).resolve();executable=Path(executable).resolve()
    logs=path/label;logs.mkdir()
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),PARALLEL=str(threads),
             OMP_DYNAMIC='FALSE',OMP_THREAD_LIMIT=str(threads),FLIB_FASTOMP='FALSE',
             OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    script='''set -eu
cd "$1"
rank=${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-}}}
if test -z "$rank"; then echo 'MPI rank environment unavailable' >&2; exit 2; fi
exe=$2
threads=$3
logs=$4
shift 4
exec "$exe" "$threads" "$@" > "$logs/stdout.rank-$rank" 2> "$logs/stderr.rank-$rank"
'''
    command=['mpiexec','-n',str(ranks),'bash','-c',script,'hsex-test',str(path),str(executable),str(threads),str(logs),*extra_args]
    with (logs/'launcher.log').open('w') as out,(logs/'launcher.stderr').open('w') as err:
        result=subprocess.run(command,cwd=path,env=env,stdout=out,stderr=err,timeout=timeout)
    stdout='\n'.join((logs/f'stdout.rank-{rank}').read_text(errors='replace')
                     for rank in range(ranks) if (logs/f'stdout.rank-{rank}').exists())
    if expect_failure:
        if result.returncode==0 or 'SCREEN_FAILURE' not in stdout:
            raise AssertionError(f'Expected collective numerical rejection, exit={result.returncode}; inspect {logs}')
    else:
        if result.returncode!=0:
            raise RuntimeError(f'MPI numerical test failed, exit={result.returncode}; inspect {logs}')
        missing=[rank for rank in range(ranks) if not (logs/f'stdout.rank-{rank}').exists()]
        if missing:raise AssertionError(f'Missing rank logs {missing}: {logs}')
        if not re.search(r'HSEX_W_OMP team=\s*'+str(threads)+r'(?:\s|$)',stdout):
            raise AssertionError(f'Actual OpenMP team {threads} not observed: {logs}')
        print(stdout.strip(),flush=True)
    return result
