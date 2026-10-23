from pathlib import Path
import os,sys,subprocess,resource,json,socket,time
rank=int(os.environ.get('OMPI_COMM_WORLD_RANK',os.environ.get('PMIX_RANK','-1')))
assert rank>=0
p=Path.cwd();t=time.monotonic()
with (p/('stdout.rank-%d'%rank)).open('x') as out,(p/('stderr.rank-%d'%rank)).open('x') as err:
 q=subprocess.run([sys.argv[1]],stdin=sys.stdin,stdout=out,stderr=err)
r=resource.getrusage(resource.RUSAGE_CHILDREN)
rec=dict(rank=rank,exit=q.returncode,host=socket.gethostname(),OMP=os.environ['OMP_NUM_THREADS'],wall=time.monotonic()-t,maxrss=r.ru_maxrss,maxrss_unit='bytes' if sys.platform=='darwin' else 'KiB',hugetlb='not available on Darwin' if sys.platform=='darwin' else 'not sampled')
with (p/('rank-%d.json'%rank)).open('x') as f:json.dump(rec,f)
sys.exit(q.returncode if q.returncode>=0 else 128-q.returncode)
