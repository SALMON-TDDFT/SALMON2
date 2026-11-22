program test_pdep_q_mpi
  use mpi
  use omp_lib, only: omp_get_num_threads
  use hsex_pdep_q, only: pdep_q_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer :: rank,size,ierr,ks,ke,j,status,pdep_rank,provided,offset,qp(3)
  real(8) :: cell(3,3),k(3,8),levels(2,8)
  complex(8),allocatable :: states(:,:,:)
  character(256) :: path,mode,arg
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
  offset=0
  if(command_argument_count()==4)offset=1
  call get_command_argument(1+offset,path);call get_command_argument(2+offset,arg);read(arg,*)pdep_rank
  call get_command_argument(3+offset,mode)
!$omp parallel
!$omp master
  if(rank==0)write(*,'(a,i0)')'HSEX_W_OMP team=',omp_get_num_threads()
!$omp end master
!$omp end parallel
  qp=[1,0,0]
  if(mode=='gamma')qp=0
  ks=1+8*rank/size;ke=8*(rank+1)/size
  allocate(states(64,2,max(0,ke-ks+1)));states=0d0
  if(ke>=ks)then
    states(1,1,:)=1d0;states(2,2,:)=1d0
  endif
  cell=0d0
  do j=1,3
    cell(j,j)=12d0
  enddo
  do j=0,7
    k(:,j+1)=dble([mod(j,2),mod(j/2,2),j/4])*acos(-1d0)/12d0
  enddo
  levels(1,:)=-.5d0;levels(2,:)=.2d0
  if(mode=='bad-gap')levels(2,:)=-.5d0
  if(mode=='bad-gram'.and.rank==1)states(1,1,:)=2d0
  if(mode=='nan-state'.and.rank==1)states(64,2,:)=ieee_value(0d0,ieee_quiet_nan)
  if(mode=='inf-state'.and.rank==1)states(64,2,:)=ieee_value(0d0,ieee_positive_inf)
  call pdep_q_test([4,4,4],[2,2,2],[2,2,2],qp,.05d0,cell,k,levels,states, &
    1,ks,ke,MPI_COMM_WORLD,rank,pdep_rank,trim(path),action,status)
  if(mode=='valid')then
    if(status/=0)error stop 'valid q diagnostic rejected'
  else
    if(status==0)error stop 'invalid q diagnostic accepted'
  endif
  write(*,'(a,i0,a,a)')'PDEP_MPI_RETURN_PASS rank=',rank,' case=',trim(mode)
  call MPI_Finalize(ierr)
contains
  subroutine action(index,x,y,code)
    implicit none
    integer,intent(in) :: index
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=0;y=.2d0*x;y(1)=-.5d0*x(1)
    if(index<ks.or.index>ke)code=1
    if(mode=='fail-cg'.and.rank==1.and.norm2(abs(x))==0d0)code=1
  end subroutine
end program
