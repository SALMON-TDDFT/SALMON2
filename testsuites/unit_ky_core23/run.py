#!/usr/bin/env python3
from pathlib import Path
import argparse
import importlib.util
import math
import os
import shlex
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--build-dir',type=Path,required=True)
p.add_argument('--checker',type=Path,required=True)
p.add_argument('--fc',default='gfortran')
p.add_argument('--libs',default='-lopenblas')
p.add_argument('--mpi',action='store_true')
p.add_argument('--launcher',default='mpiexec')
p.add_argument('--pseudo',type=Path,help='optional actual V2 file, tested in addition to synthetic core')
a=p.parse_args()
root=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('core23_test',Path(__file__).with_name('test_reader.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
objroot=a.build_dir.resolve()/'src/CMakeFiles/salmon.dir'
subprocess.run(['python3',str(a.checker),'--source-root',str(root/'src')],check=True)
subprocess.run(['python3',str(a.checker),'--object-tree',str(objroot)],check=True)
objects=sorted(str(x) for x in objroot.rglob('*.o') if x.name!='main.f90.o')
with tempfile.TemporaryDirectory(prefix='ky-core23-reader-') as tmp:
    d=Path(tmp);exe=d/'probe';obj=d/'probe.o'
    subprocess.run([a.fc,'-cpp',*(['-DTEST_MPI'] if a.mpi else []),'-fopenmp','-ffree-line-length-none','-fcheck=all','-I'+str(a.build_dir.resolve()),
                    '-c',str(Path(__file__).with_name('driver.f90')),'-o',str(obj)],check=True)
    subprocess.run(['python3',str(a.checker),'--object-tree',str(d)],check=True)
    subprocess.run([a.fc,'-fopenmp',str(obj),*objects,*shlex.split(a.libs),'-o',str(exe)],check=True)
    v2=d/'synthetic_v2_rps.dat';v2.write_text(m.synthetic_v2())
    lines=v2.read_text().splitlines();v1=d/'synthetic_v1_rps.dat'
    v1.write_text('\n'.join(lines[:lines.index('CORE_AUX_PAULI_V1')]).replace('KY_SOC_NLCC_V2','KY_SOC_NLCC_V1')+'\nEND_KY_SOC_NLCC_V1\n')
    launch=[a.launcher,'-n','2'] if a.mpi else []
    def run(file,mode=[]):
        subprocess.run([*launch,str(exe),str(file.resolve()),*mode],cwd=d,check=True,
                       env=dict(os.environ,OMP_NUM_THREADS='1'),stdout=subprocess.DEVNULL)
        out=(d/'reader.rank0.data').read_text()
        if a.mpi and out!=(d/'reader.rank1.data').read_text():
            raise SystemExit('MPI core table/grid broadcasts differ')
        h,*rows=out.splitlines();rows=[list(map(float,x.split())) for x in rows]
        if not all(math.isfinite(v) for row in rows for v in row):
            raise SystemExit('nonfinite reader/grid output')
        return h.split(),rows
    h2,r2=run(v2);h1,r1=run(v1)
    if h2[:2]!=['T','T'] or h1[:2]!=['F','F']:
        raise SystemExit('wrong core completeness flags')
    n=int(h2[2])+1
    # Grid/projectors remain identical. V2 retains density tails that V1's
    # historical 1e-7 relative cutoff omitted; all retained V1 samples agree.
    for x,y in zip(r1[:n],r2[:n]):
        if x[0]!=y[0] or x[4:]!=y[4:] or (x[1]!=0 and x[1]!=y[1]):
            raise SystemExit('V2 changed retained V1 density/projector samples')
        if x[1]==0 and y[1]>1e-7*r2[0][1]:
            raise SystemExit('V2 differs outside documented legacy density tail')
    original=lines[lines.index('CORE_AUX_PAULI_V1')+5:-1]
    for actual,raw in zip(r2[:n],original):
        wanted=float(raw.split()[1])
        if abs(actual[1]-wanted)>2e-8*max(1.0,wanted):
            raise SystemExit('V2 independently truncated radial core density')
    # Check tensor orientation on +x: J_y^z=+C, J_z^y=-C; other entries zero.
    sample=r2[n+1];r=0.4;table=r2[:n]
    j=next(i for i in range(1,n) if table[i][0]>r)
    w=(r-table[j-1][0])/(table[j][0]-table[j-1][0])
    tau=(1-w)*table[j-1][2]+w*table[j][2]
    cur=(1-w)*table[j-1][3]+w*table[j][3]
    if abs(sample[1]-tau)>1e-12:
        raise SystemExit('core tau interpolation mismatch')
    tensor=sample[2:]
    expected=[0.0]*12;expected[7]=cur;expected[10]=-cur
    if max(abs(x-y) for x,y in zip(tensor,expected))>1e-12:
        raise SystemExit('core spin-current Cartesian orientation mismatch')
    _,periodic=run(v2,['--periodic'])
    def interp(r,col):
        h=table[1][0]-table[0][0]
        i=int(r/h)
        if i>=len(table)-1: return 0.0
        w=(r-table[i][0])/h
        return (1-w)*table[i][col]+w*table[i+1][col]
    rho_sum=tau_sum=old_rho_sum=0.0
    for x in range(-8,9):
      for y in range(-8,9):
       for z in range(-8,9):
        r=math.sqrt(x*x+y*y+z*z)
        rho_sum+=interp(r,1);tau_sum+=interp(r,2)
        if max(abs(x),abs(y),abs(z))<=2: old_rho_sum+=interp(r,1)
    if abs(periodic[n][0]-rho_sum)>2e-12 or abs(periodic[n][1]-tau_sum)>2e-12:
        raise SystemExit('adaptive periodic images omit core support')
    if rho_sum-old_rho_sum<1e-9:
        raise SystemExit('periodic regression does not exercise beyond two images')
    if r2[-1][1]<=0.0:
        raise SystemExit('loaded synthetic core never reached gKS action')
    print('PASS V1/V2 projector and retained density identity, joint core interpolation, tensor orientation, MPI broadcasts')
    print('PASS adaptive periodic images; omitted by legacy +/-2:',rho_sum-old_rho_sum)
    if a.pseudo:
        h,rows=run(a.pseudo)
        if h[:2]!=['T','T']:
            raise SystemExit('actual V2 core metadata incomplete')
        if rows[-1][1]<=0.0:
            raise SystemExit('actual V2 core never reached gKS action')
        print('PASS actual V2 production reader/grid/gKS:',a.pseudo,rows[-1])
