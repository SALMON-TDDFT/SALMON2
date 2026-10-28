#!/usr/bin/env python3
"""Run numerical shared PT-CN regressions against an audited ScaLAPACK build.

Usage: python3 run.py BUILD --libs '-L/path/to/lib -lscalapack -lfftw3 -lopenblas'
Requires mpifort/mpiexec on PATH. Does not run any production trajectory.
"""
import argparse, os, shlex, subprocess, tempfile
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build',type=Path)
p.add_argument('--libs',required=True)
a=p.parse_args(); b=a.build.resolve(); root=Path(__file__).resolve().parents[2]
checker=root/'tools/check_forbidden_eigensolvers.py'
def run(cmd,**kw): subprocess.run(list(map(str,cmd)),check=True,**kw)
run(['python3',checker,'--source-root',root/'src'])
obj=b/'src/CMakeFiles/salmon.dir'
run(['python3',checker,'--object-tree',obj])
libs=shlex.split(a.libs)
with tempfile.TemporaryDirectory(prefix='ptcn-shared-regression-') as tmp:
 d=Path(tmp)
 cases={'occupation_probe':['ptcn_core'],
        'lowdin_owner_probe':['ptcn_distributed_matrix','ptcn_native_matrix'],
        'lowdin_probe':['ptcn_core'],
        'learned_exchange_probe':['ptcn_core'],
        'fresh_exchange_probe':['ptcn_core'],
        'kinetic_probe':['ptcn_rectangular_fft','ptcn_kinetic_preconditioner'],
        'spinor_probe':['ptcn_core','exx_ptcn_core','ptcn_distributed_matrix','ptcn_native_matrix']}
 for name,modules in cases.items():
  o=d/(name+'.o'); exe=d/name
  run(['mpifort','-fcheck=all','-ffree-line-length-none','-I'+str(b),'-c',Path(__file__).with_name(name+'.f90'),'-o',o],cwd=d)
  run(['python3',checker,'--object-tree',d])
  run(['mpifort',o,*[obj/'rt'/(m+'.f90.o') for m in modules],*libs,'-o',exe],cwd=d)
  env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
  if name=='lowdin_owner_probe':
   for ranks in [1,2,4]:run(['mpiexec','-n',ranks,exe,'16'],cwd=d,env=env)
  elif name=='spinor_probe':
   for ranks in [1,2]:
    for mode in [0,1,2]:run(['mpiexec','-n',ranks,exe,mode],cwd=d,env=env)
  else:run([exe],cwd=d,env=env)
print('SHARED_PT_CN_REGRESSIONS_PASS')
