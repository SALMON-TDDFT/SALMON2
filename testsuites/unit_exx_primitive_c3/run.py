#!/usr/bin/env python3
"""Bounded geometry/SU2 test, no SALMON job or material calculation."""
from pathlib import Path
import argparse
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--checker',type=Path,required=True)
p.add_argument('--fc',default='gfortran')
a=p.parse_args()
root=Path(__file__).resolve().parents[2]
subprocess.run(['python3',str(a.checker),'--source-root',str(root/'src')],check=True)
with tempfile.TemporaryDirectory(prefix='primitive-c3-') as tmp:
    d=Path(tmp)
    for source,name in [(root/'src/xc/exx_symmetry.f90','sym'),(Path(__file__).with_name('probe.f90'),'probe')]:
        subprocess.run([a.fc,'-ffree-line-length-none','-fcheck=all','-J'+str(d),'-I'+str(d),
                        '-c',str(source),'-o',str(d/(name+'.o'))],check=True)
    subprocess.run(['python3',str(a.checker),'--object-tree',str(d)],check=True)
    subprocess.run([a.fc,str(d/'sym.o'),str(d/'probe.o'),'-o',str(d/'probe')],check=True)
    subprocess.run([str(d/'probe')],check=True)
