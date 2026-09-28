"""Independent sparse candidate completeness and scaling checks under MPI."""
import argparse, os, subprocess, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(); p.add_argument('--build', type=Path, required=True); a = p.parse_args()
b = a.build.resolve(); obj = b/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory() as tmp:
    exe = Path(tmp)/'probe'
    subprocess.run(['mpifort', '-fopenmp', '-fcheck=all', '-ffree-line-length-none', '-I'+str(b),
        str(ROOT/'src/xc/exx_pair_candidates.f90'), str(Path(__file__).with_name('pair_candidates_probe.f90')),
        str(obj/'parallel/communication.f90.o'), str(obj/'misc/nvtx_wrapper.f90.o'), '-o', str(exe)], cwd=tmp, check=True)
    for ranks in (1, 2, 4):
        subprocess.run(['mpiexec', '-n', str(ranks), str(exe)], cwd=tmp, check=True,
            env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=30)
