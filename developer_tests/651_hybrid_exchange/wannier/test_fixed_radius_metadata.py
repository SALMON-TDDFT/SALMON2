"""Fresh DC to RT support changes must not weaken physical/restart metadata checks."""
import argparse,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args()
b=a.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory(prefix='fixed-radius-') as tmp:
 exe=Path(tmp)/'probe'
 subprocess.run(['mpifort','-fopenmp','-I'+str(b),str(Path(__file__).with_name('fixed_radius_metadata_probe.f90')),
  str(obj/'io/salmon_global.f90.o'),str(obj/'xc/exx_functional.f90.o'),'-o',str(exe)],check=True)
 subprocess.run([str(exe)],cwd=tmp,check=True)
