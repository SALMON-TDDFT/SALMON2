from pathlib import Path
import os,subprocess,tempfile
here=Path(__file__).resolve().parent;root=here.parents[1]
with tempfile.TemporaryDirectory() as folder:
 p=Path(folder);(p/'config.h').write_text('')
 subprocess.run(['mpifort','-cpp','-ffree-line-length-none','-fallow-argument-mismatch','-DUSE_MPI','-I',folder,'-O2','-fcheck=all',str(root/'src/misc/nvtx_wrapper.f90'),str(root/'src/parallel/communication.f90'),str(root/'src/xc/lcfo_dist_rows.f90'),str(here/'gather_memory_probe.f90'),'-o',str(p/'probe')],cwd=p,check=True)
 for ranks in (2,4):subprocess.run(['mpirun','-np',str(ranks),str(p/'probe')],cwd=p,check=True)
