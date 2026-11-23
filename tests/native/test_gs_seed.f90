! Analytic arrays only: no GS solve, checkpoint read or checkpoint write.
program test_gs_seed
use mpi
use structures
use salmonpy_gs_seed
use communication, only: comm_attach
implicit none
type(s_gs_seed) :: seed
type(s_dft_system) :: system
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_orbital) :: psi
integer :: rank,n,status,provided,o,mode
complex(8) :: expected
call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,status)
call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
call MPI_Comm_size(MPI_COMM_WORLD,n,status)
if(n/=2)stop 1
call comm_attach(MPI_COMM_WORLD,status)
if(status/=0)stop 2
! Source has 6 bands split 1:3/4:6, target has 4 split 1:2/3:4.
seed%orbital_start=rank*3+1;seed%orbital_end=(rank+1)*3
allocate(seed%occupations(6,1,1));seed%occupations=0d0;seed%occupations(:4,1,1)=2d0
system%no=4;system%nk=1;system%nspin=1
allocate(system%rocc(4,1,1))
mg%is=1;mg%ie=2
info%io_s=rank*2+1;info%io_e=(rank+1)*2
info%ik_s=1;info%ik_e=1;info%im_s=1;info%im_e=1;info%icomm_o=MPI_COMM_WORLD
allocate(psi%zwf(1:2,1:2,1:2,1,info%io_s:info%io_e,1,1))
do mode=1,2
 if(mode==1)then
  allocate(seed%psi%rwf(1:2,1:2,1:2,1,seed%orbital_start:seed%orbital_end,1,1))
  do o=seed%orbital_start,seed%orbital_end
   seed%psi%rwf(:,:,:,1,o,1,1)=real(o,8)
  enddo
 else
  deallocate(seed%psi%rwf)
  allocate(seed%psi%zwf(1:2,1:2,1:2,1,seed%orbital_start:seed%orbital_end,1,1))
  do o=seed%orbital_start,seed%orbital_end
   seed%psi%zwf(:,:,:,1,o,1,1)=cmplx(real(o,8),-real(o,8),8)
  enddo
 endif
 call load_gs_seed(seed,system,mg,info,psi)
 if(any(system%rocc/=2d0))stop 3
 do o=info%io_s,info%io_e
  expected=cmplx(real(o,8),0d0,8)
  if(mode==2)expected=cmplx(real(o,8),-real(o,8),8)
  if(any(psi%zwf(:,:,:,1,o,1,1)/=expected))stop 4
 enddo
enddo
if(rank==0)print *,'PASS analytic real/complex empty-band redistribution; SCF_calls=0'
call MPI_Finalize(status)
end program
