"""Compile the actual band path code against minimal types; test skew basis/padding."""
from pathlib import Path
import os,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[1]
compiler=shutil.which(os.environ.get('FC','gfortran'))
if compiler is None: raise SystemExit('gfortran required for native path test')
stubs='''
module parallelization
 integer :: nproc_id_global=0
end module
module communication
contains
 logical function comm_is_root(i)
 integer :: i
 comm_is_root=i==0
 end function
end module
module structures
 type s_dft_system
 integer :: nk=4,no=8,nspin=1
 real(8) :: primitive_b(3,3),primitive_a(3,3),wtk(4),vec_k(3,4)
 end type
 type s_band_dft
 integer :: nref_band,num_band_kpt
 real(8),allocatable :: band_kpt(:,:)
 logical,allocatable :: check_conv_esp(:,:,:)
 end type
 type s_parallel_info
 integer :: ik_s=1,ik_e=4
 end type
 type s_dft_energy
 real(8),allocatable :: esp(:,:,:)
 end type
end module
module salmon_global
 character(256) :: file_band_kpoints='none'
 character(3) :: lattice='non'
 integer :: nref_band=3,num_of_segments=2,ndiv_segment(10)
 real(8) :: kpt(3,11)
 character :: kpt_label(11)
end module
'''
program='''
program test
 use band_dft_sub
 use structures
 use salmon_global
 implicit none
 type(s_dft_system) :: sys
 type(s_band_dft) :: band
 type(s_parallel_info) :: info
 real(8) :: endpoint(3),midpoint(3)
 integer :: nref,i
 real(8),allocatable :: points(:,:)
 sys%primitive_b=reshape([2d0,0d0,0d0,1d0,3d0,0d0,0d0,.5d0,4d0],[3,3])
 kpt=0d0
 kpt(:,2)=[.5d0,.25d0,0d0]
 kpt(:,3)=[.25d0,0d0,.5d0]
 midpoint=[1.25d0,.75d0,0d0]
 endpoint=[.5d0,.25d0,2d0]
 ! 7 intervals + endpoint: exactly two batches.
 ndiv_segment(1:2)=[3,4]
 call get_band_kpt(points,nref,sys)
 if(size(points,2)/=8) error stop 'exact division size'
 if(maxval(abs(points(:,4)-midpoint))>1d-14) error stop 'skew midpoint'
 if(maxval(abs(points(:,8)-endpoint))>1d-14) error stop 'exact endpoint'
 deallocate(points)
 ! 5 intervals + endpoint: two copies pad the final batch.
 ndiv_segment(1:2)=[2,3]
 call get_band_kpt(points,nref,sys)
 if(size(points,2)/=8) error stop 'padding size'
 do i=6,8
 if(maxval(abs(points(:,i)-endpoint))>1d-14) error stop 'padded endpoint'
 end do
 call move_alloc(points,band%band_kpt)
 allocate(band%check_conv_esp(3,4,1))
 call calc_band_write(1,sys,band,info)
 if(maxval(abs(sys%vec_k-band%band_kpt(:,1:4)))>1d-14) error stop 'double transformation'
 open(91,file='points.txt',status='replace')
 write(91,*) 3
 write(91,*) 0d0,0d0,0d0
 write(91,*) .5d0,.25d0,0d0
 write(91,*) .25d0,0d0,.5d0
 close(91)
 file_band_kpoints='points.txt'
 call get_band_kpt(points,nref,sys)
 if(size(points,2)/=4) error stop 'explicit padding size'
 if(maxval(abs(points(:,2)-midpoint))>1d-14) error stop 'explicit exact point'
 do i=3,4
 if(maxval(abs(points(:,i)-endpoint))>1d-14) error stop 'explicit endpoint'
 end do
 print *, 'PASS native skew path, divisible endpoint, padding and single transformation'
end program
'''
with tempfile.TemporaryDirectory(prefix='salmonpy-band-path-') as tmp:
    tmp=Path(tmp);(tmp/'stubs.f90').write_text(stubs);(tmp/'test.f90').write_text(program)
    exe=tmp/('test.exe' if os.name=='nt' else 'test')
    subprocess.run([compiler,'-ffree-line-length-none','-fcheck=all','-g','stubs.f90',str(root/'src/gs/band_dft.f90'),'test.f90','-o',str(exe)],cwd=tmp,check=True)
    completed=subprocess.run([str(exe)],cwd=tmp,check=True,capture_output=True,text=True)
    print(completed.stdout.splitlines()[-1])
