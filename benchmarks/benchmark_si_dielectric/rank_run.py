"""Record per-SALMON-rank child peak RSS; do not use launcher RSS as rank memory."""
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

rank = os.environ.get('OMPI_COMM_WORLD_RANK', os.environ.get('PMI_RANK'))
if rank is None:
    raise RuntimeError('MPI rank environment unavailable')
started = time.monotonic()
run = subprocess.run(sys.argv[1:])
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
rss_bytes = int(usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024))
Path(f'rank-{int(rank):04d}.json').write_text(json.dumps(dict(
    rank=int(rank), returncode=run.returncode, wall_seconds=time.monotonic()-started,
    peak_rss_bytes=rss_bytes), indent=2)+'\n')
sys.exit(run.returncode)
