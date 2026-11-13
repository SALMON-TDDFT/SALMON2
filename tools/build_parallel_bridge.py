"""Explicitly build the subgroup ABI against a prepared runtime and host MPI.
Usage: python tools/build_parallel_bridge.py PREPARED_SOURCE NEW_DIRECTORY [--fc mpifort]
No MPI library is bundled; use the same MPI installation as mpi4py.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source');p.add_argument('directory');p.add_argument('--fc',default='mpifort')
a=p.parse_args();source=Path(a.source).resolve();out=Path(a.directory).resolve()
provenance=json.loads((source/'salmonpy-runtime-source.json').read_text())
paths=[source/'src/misc/nvtx_wrapper.f90',source/'src/parallel/communication.f90',source/'src/salmonpy/parallel_bridge.f90']
if not all(p.is_file() for p in paths): raise ValueError('prepared runtime lacks MPI attachment sources')
if out.is_relative_to(source): raise ValueError('build outside prepared source')
out.mkdir(parents=True,exist_ok=False)
lib=out/('libsalmonpy_parallel.dylib' if sys.platform=='darwin' else 'libsalmonpy_parallel.so')
checker=Path(__file__).resolve().with_name('check_forbidden_eigensolvers.py')
subprocess.run([sys.executable,str(checker),'--source-root',str(source)],check=True)
cmd=[a.fc,'-cpp','-fPIC','-c','-O2','-fallow-argument-mismatch',*map(str,paths)]
with (out/'build.log').open('w') as log: subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
subprocess.run([a.fc,'-shared',*map(str,sorted(out.glob('*.o'))),'-o',str(lib)],cwd=out,check=True)
record={'schema':'salmonpy.parallel-bridge/1','sha256':hashlib.sha256(lib.read_bytes()).hexdigest(),
        'library':str(lib),'command':cmd,'source':provenance,
        'source_hashes':{str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
(out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
