"""Launch one MPI rank and record the SALMON child's lifetime high-water RSS."""
import json,os,sys,time
from pathlib import Path
rank=os.environ.get('OMPI_COMM_WORLD_RANK',os.environ.get('PMI_RANK'))
if rank is None:raise RuntimeError('MPI rank identity is unavailable')
start=time.perf_counter()
pid=os.posix_spawn(sys.argv[1],[sys.argv[1]],os.environ)
_,status,usage=os.wait4(pid,0)
code=os.waitstatus_to_exitcode(status)
Path(f'rank-{int(rank)}.json').write_text(json.dumps(dict(rank=int(rank),returncode=code,
    wall_seconds=time.perf_counter()-start,peak_rss_bytes=int(usage.ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))+'\n')
sys.exit(code if code>=0 else 128-code)
