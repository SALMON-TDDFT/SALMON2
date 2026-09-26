from pathlib import Path
import os,subprocess,tempfile,struct
root=Path(__file__).resolve().parents[2];here=Path(__file__).parent
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 (Path(folder)/'config.h').write_text('')
 subprocess.run(['gfortran','-I',folder,'-cpp','-fopenmp','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all',
  str(here/'transport_stubs.f90'),str(root/'src/xc/hse_wannier_gauge.f90'),
  *[str(root/'src/xc'/name) for name in ['hse_ace.f90','lcfo_dist_rows.f90','lcfo_dist_dense.f90','lcfo_rt_wannier.f90']],str(here/'sphere_probe.f90'),
  '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],cwd=folder,check=True)
 for threads in ['1','2','4']:
  for radius,mode in [('1.51','normal'),('4.1','normal'),('0','normal'),('9','normal'),('1.51','weak')]:
   subprocess.run([str(exe),radius,mode],cwd=folder,check=True,env=dict(os.environ,
    SALMON_LCFO_RT_MLWF='1',SALMON_LCFO_RT_RADIUS=radius,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS=threads))

   data=(Path(folder)/'lcfo_mlwf_initial.bin').read_bytes()
   assert struct.unpack_from('=7i',data)==(16909060,2,1,1,8,8,8)
   assert len(data)==144
   centers=struct.unpack_from('=3d',data,84)
   for center,expected,length in zip(centers,[7,7.5,1.5],[8,10,12]):
    assert abs((center-expected+length/2)%length-length/2)<1e-12
   assert struct.unpack_from('=i',data,116)[0]==(mode=='weak')
