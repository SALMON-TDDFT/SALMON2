#!/usr/bin/env python3
"""Verify frozen auxiliary-Pauli NLCC response and optionally full SALMON bridge."""
from pathlib import Path
import argparse
import os
import shlex
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--checker',type=Path,required=True)
p.add_argument('--build-dir',type=Path)
p.add_argument('--fc',default='gfortran')
p.add_argument('--libs',default='-lopenblas')
p.add_argument('--mpi',action='store_true')
p.add_argument('--launcher',default='mpiexec')
a=p.parse_args()
root=Path(__file__).resolve().parents[2]

def run(args,**kwargs):
    return subprocess.run([str(x) for x in args],check=True,**kwargs)

run(['python3',a.checker,'--source-root',root/'src'])
with tempfile.TemporaryDirectory(prefix='nc-mgga23-core-') as scratch:
    scratch=Path(scratch)
    if not a.mpi:
        flags=['-std=f2008','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-J',scratch,'-I',scratch]
        source_names=[root/'src/xc/nc_mgga23_kernel.f90',root/'src/xc/nc_mgga23_core.f90',
                      Path(__file__).with_name('test_core.f90')]
        objects=[]
        for i,source in enumerate(source_names):
            obj=scratch/f'pure_{i}.o'
            run([a.fc,*flags,'-c',source,'-o',obj]);objects.append(obj)
        run(['python3',a.checker,'--object-tree',scratch])
        exe=scratch/'pure'
        run([a.fc,*objects,'-o',exe]);run([exe])
    if a.build_dir:
        build=a.build_dir.resolve();objroot=build/'src/CMakeFiles/salmon.dir'
        run(['python3',a.checker,'--object-tree',objroot])
        objects=sorted(x for x in objroot.rglob('*.o') if x.name!='main.f90.o')
        if not objects: raise SystemExit('no SALMON application objects found')
        obj=scratch/'bridge.o';exe=scratch/'bridge'
        run([a.fc,'-cpp',*(['-DTEST_MPI'] if a.mpi else []),'-fopenmp','-ffree-line-length-none','-fcheck=all',
             '-I'+str(build),'-c',Path(__file__).with_name('test_core_bridge.f90'),'-o',obj])
        run(['python3',a.checker,'--object-tree',scratch])
        run([a.fc,'-fopenmp',obj,*objects,*shlex.split(a.libs),'-o',exe])
        launch=[a.launcher,'-n','2'] if a.mpi else []
        run([*launch,exe],cwd=scratch,env=dict(os.environ,OMP_NUM_THREADS='1'))
        for mode in ['--invalid-metadata','--bad-core','--missing-tau','--bad-layout']:
            rejected=subprocess.run([*map(str,launch),str(exe),mode],cwd=scratch,text=True,capture_output=True,
                                    env=dict(os.environ,OMP_NUM_THREADS='1'))
            expected='incomplete V2' if mode in ['--missing-tau','--bad-layout'] else 'consistent core tau/current model'
            if rejected.returncode==0 or expected not in rejected.stderr:
                raise SystemExit('invalid core accepted or wrong rejection: '+mode+'\n'+rejected.stderr)
            print('PASS invalid-core rejection:',mode)
