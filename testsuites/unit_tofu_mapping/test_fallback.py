"""Exercise the production Fujitsu mapping routine with mocked topology APIs.
No Fujitsu runtime or MPI launch required; this verifies control flow only.
"""
from pathlib import Path
import subprocess,tempfile
repo=Path(__file__).resolve().parents[2]
s=(repo/'src/parallel/init_communicator.f90').read_text()
a=s.index('  subroutine tofu_network_oriented_mapping(iret)')
b=s.index('  end subroutine tofu_network_oriented_mapping',a)
routine=s[a:b]+ '  end subroutine tofu_network_oriented_mapping\n'
stubs='''module salmon_global
 integer :: nx_m=1,ny_m=1,nz_m=1,num_fragment(3)=[4,4,4]
 character(64) :: theory='dft'
 character :: yn_dc='y'
end module
module mpi_ext
 integer,parameter :: MPI_SUCCESS=0,FJMPI_LOGICAL=0
 integer :: dim=1,shape(3)=[16,1,1],ppn=1,rank_calls=0
contains
 subroutine FJMPI_Topology_get_dimension(d,ierr)
 integer :: d,ierr
 d=dim;ierr=0
 end subroutine
 subroutine FJMPI_Topology_get_shape(x,y,z,ierr)
 integer :: x,y,z,ierr
 x=shape(1);y=shape(2);z=shape(3);ierr=0
 end subroutine
 subroutine FJMPI_Topology_get_coords(comm,rank,mode,d,coords,ierr)
 integer :: comm,rank,mode,d,coords(*),ierr
 coords(1:d)=0;ierr=0
 end subroutine
 subroutine FJMPI_Topology_get_ranks(comm,mode,coords,maxppn,outppn,ranks,ierr)
 integer :: comm,mode,coords(*),maxppn,outppn,ranks(*),ierr
 rank_calls=rank_calls+1;outppn=ppn;ranks(1:ppn)=0;ierr=0
 end subroutine
end module
program probe
 use mpi_ext
 use salmon_global
 implicit none
 type info_type
 integer :: icomm_rko=0,id_rko=0,isize_rko=1,iaddress(5)=0
 integer :: imap(0:0,0:0,0:0,0:0,0:0)=-1
 end type
 type(info_type) :: info
 integer :: nproc_d_o(3)=[1,1,1],nproc_ob=1,nproc_k=1,myrank=0,nl,ix,iy,iz,iret,scenario
 character(32) :: process_allocation='grid_sequential',arg
 call get_command_argument(1,arg);read(arg,*)scenario
 select case(scenario)
 case(1) ! Reported DC fragment: global 64 ranks / 16 nodes, local comm size 1.
 case(2) ! Same bug can occur with a 3D allocation.
 dim=3;shape=[4,2,2]
 case(3) ! Total-system initialization also must skip unsupported topology early.
 yn_dc='t';info%isize_rko=64;ppn=4
 case(4) ! Retain supported ordinary 3D mapping.
 yn_dc='n';dim=3;shape=1
 end select
 call tofu_network_oriented_mapping(iret)
 if(scenario<4)then
 if(iret>=0.or.rank_calls/=0)error stop 'Expected fallback before node-rank query'
 else
 if(iret/=0.or.rank_calls/=1.or.any(info%imap/=0))error stop 'Supported mapping changed'
 endif
 print *,'PASS',scenario
contains
 logical function comm_is_root(rank)
 integer :: rank
 comm_is_root=rank==0
 end function
'''
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'probe.f90').write_text(stubs+routine+'end program\n')
 subprocess.run(['gfortran','-fcheck=all','-ffree-line-length-none',str(p/'probe.f90'),'-o',str(p/'probe')],cwd=d,check=True)
 for case in range(1,5):
  run=subprocess.run([str(p/'probe'),str(case)],capture_output=True,text=True,check=True)
  assert 'PASS' in run.stdout,(case,run.stdout,run.stderr)
  print('PASS scenario',case)
