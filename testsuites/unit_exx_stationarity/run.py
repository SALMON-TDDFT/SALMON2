"""Tiny scalar-metric tests; no material job or eigensolver is executed."""
from pathlib import Path
import argparse, subprocess, tempfile
p=argparse.ArgumentParser();p.add_argument('--checker',type=Path,required=True);p.add_argument('--fc',default='gfortran');a=p.parse_args()
r=Path(__file__).resolve().parents[2]
subprocess.run(['python3',str(a.checker),'--source-root',str(r/'src')],check=True)
with tempfile.TemporaryDirectory(prefix='exx-stationarity-') as name:
 d=Path(name)
 for source in [r/'src/xc/exx_stationarity.f90',Path(__file__).with_name('probe.f90')]:
  subprocess.run([a.fc,'-fcheck=all','-ffree-line-length-none','-J'+str(d),'-I'+str(d),'-c',str(source),'-o',str(d/(source.stem+'.o'))],check=True)
 subprocess.run(['python3',str(a.checker),'--object-tree',str(d)],check=True)
 subprocess.run([a.fc,str(d/'exx_stationarity.o'),str(d/'probe.o'),'-o',str(d/'probe')],check=True)
 subprocess.run([str(d/'probe')],check=True)
