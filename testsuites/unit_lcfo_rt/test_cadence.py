from pathlib import Path
import os,subprocess,tempfile
root=Path(__file__).resolve().parents[2];here=Path(__file__).parent
with tempfile.TemporaryDirectory() as d:
    exe=Path(d)/'probe'
    (Path(d)/'config.h').write_text('')
    cmd=[os.environ.get('FC','gfortran'),'-I',d,'-cpp','-fopenmp','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',
         str(here/'transport_stubs.f90'),str(root/'src/xc/hse_wannier_gauge.f90'),
         *[str(root/'src/xc'/name) for name in ['hse_ace.f90','lcfo_dist_rows.f90','lcfo_dist_dense.f90','lcfo_seed.f90','lcfo_mlwf_links.f90','lcfo_wf_support.f90','lcfo_rt_wannier.f90']],str(here/'cadence_probe.f90'),
         '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)]
    subprocess.run(cmd,cwd=d,check=True)
    for threads in ['1','2','4']:
      for radius in ['0','1']:
        subprocess.run([str(exe),radius],cwd=d,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',
                       OMP_NUM_THREADS=threads))

