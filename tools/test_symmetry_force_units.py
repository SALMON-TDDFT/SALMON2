"""Explicit synthetic units/covariance test; no GS computation or datasets."""
from pathlib import Path
import subprocess,hashlib,json,shutil,sys
import argparse

def main():
 parser=argparse.ArgumentParser(description='Audit/compile actual force-symmetrization module; Cartesian identity/C2/skew controls.')
 parser.add_argument('output',type=Path)
 parser.add_argument('--source',type=Path)
 parser.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
 parser.add_argument('--compiler',default=shutil.which('gfortran-15') or shutil.which('gfortran'))
 args=parser.parse_args()
 root=Path(__file__).resolve().parents[1]
 source=args.source or root/'runtime/sources/solver/src/symmetry/symmetry_vector.f90'
 data=source.read_bytes()
 if not args.compiler:raise ValueError('explicit Fortran compiler required')
 out=args.output.resolve()
 if out.is_relative_to(root) or out.is_relative_to(source.resolve().parent):raise ValueError('new external output directory required')
 out.mkdir(exist_ok=False);src=out/'src';src.mkdir()
 (src/'symmetry_vector.f90').write_bytes(data)
 (src/'shim.f90').write_text('module sym_sub\nimplicit none\nreal(8)::Amat(3,3),Ainv(3,3)\nreal(8),allocatable::SymMatA(:,:,:),SymMatB(:,:,:)\nlogical::use_symmetry=.true.\nend module\n')
 (src/'probe.f90').write_text('''program probe
 use sym_sub
 use sym_vector_sub
 implicit none
 real(8)::f(3,2),original(3,2),r(3,2),fractional(3,2),wanted(3,2)
 integer::case_id,i
 fractional(:,1)=0d0;fractional(:,2)=.5d0
 original(:,1)=[1d0,2d0,3d0];original(:,2)=[-4d0,1d0,7d0]
 do case_id=1,3
  Amat=0d0;Ainv=0d0
  if(case_id==1.or.case_id==2)then
   do i=1,3
    Amat(i,i)=real(i+1,8);Ainv(i,i)=1d0/real(i+1,8)
   enddo
  else
   Amat=2d0;Ainv=.25d0
   do i=1,3
    Amat(i,i)=0d0;Ainv(i,i)=-.25d0
   enddo
  endif
  allocate(SymMatA(3,4,merge(1,2,case_id==1)));SymMatA=0d0
  do i=1,3
   SymMatA(i,i,1)=1d0
  enddo
  wanted=original
  if(case_id==2)then
   SymMatA(:,:,2)=SymMatA(:,:,1)
   SymMatA(1,1,2)=-1d0;SymMatA(2,2,2)=-1d0
   wanted(1:2,:)=0d0
  else if(case_id==3)then
   SymMatA(:,1,2)=[0d0,1d0,0d0]
   SymMatA(:,2,2)=[1d0,0d0,0d0]
   SymMatA(:,3,2)=[0d0,0d0,1d0]
   wanted(:,1)=[1.5d0,1.5d0,3d0];wanted(:,2)=[-1.5d0,-1.5d0,7d0]
  endif
  r=matmul(Amat,fractional);f=original
  call sym_vector_force_xyz(f,r)
  if(maxval(abs(f-wanted))>1d-13)then
   print *, 'force Cartesian units FAIL',case_id,f,wanted
   stop 1
  endif
  deallocate(SymMatA)
 enddo
 print *, 'force Cartesian units PASS'
 end program
 ''')
 checker=str(args.checker)
 subprocess.run([sys.executable,checker,'--source-root',str(src)],check=True)
 for file in ['shim','symmetry_vector','probe']:
  subprocess.run([args.compiler,'-I'+str(out),'-J'+str(out),'-c',str(src/(file+'.f90')),'-o',str(out/(file+'.o'))],check=True)
 subprocess.run([sys.executable,checker,'--object-tree',str(out)],check=True)
 subprocess.run([args.compiler,*[str(out/(n+'.o')) for n in ['shim','symmetry_vector','probe']],'-o',str(out/'probe')],check=True)
 subprocess.run([str(out/'probe')],check=True)
 (out/'result.json').write_text(json.dumps({'status':'pass','actual_module_sha256':hashlib.sha256((src/'symmetry_vector.f90').read_bytes()).hexdigest(),'executable_sha256':hashlib.sha256((out/'probe').read_bytes()).hexdigest(),'cases':['orthogonal identity','orthogonal C2','skew swap'],'source_objects_audited':True},indent=2)+'\n')

if __name__=="__main__":main()
