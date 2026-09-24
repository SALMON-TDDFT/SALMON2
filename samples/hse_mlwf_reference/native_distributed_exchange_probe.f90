program probe
  use, intrinsic :: ieee_arithmetic
  use mpi
  use hse_exchange
  implicit none
  type(hse_kernel) :: op
  integer :: rank,np,status,n,m,no,nt,nk,ng,ik0,nloc,iu,ierr,p
  integer,allocatable :: starts(:),counts(:)
  real(8) :: h,omega
  real(8),allocatable :: k(:,:)
  complex(8),allocatable :: u(:,:,:),t(:,:,:),a(:,:,:),local_u(:,:,:),local_t(:,:,:),local_a(:,:,:)
  character(1024) :: input,output,invalid
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,np,status)
  call get_command_argument(1,input);call get_command_argument(2,output)
  open(newunit=iu,file=trim(input),access='stream',form='unformatted',status='old')
  read(iu)n,m,no,nt;read(iu)h,omega
  nk=m**3;ng=n**3
  allocate(k(3,nk),u(ng,no,nk),t(ng,nt,nk),a(ng,nt,nk),starts(np),counts(np))
  read(iu)k;read(iu)u;read(iu)t;close(iu)
  do p=0,np-1
    starts(p+1)=1+(p*nk)/np;counts(p+1)=((p+1)*nk)/np-(p*nk)/np
  enddo
  ik0=starts(rank+1);nloc=counts(rank+1)
  local_u=u(:,:,ik0:ik0+nloc-1);local_t=t(:,:,ik0:ik0+nloc-1)
  allocate(local_a(ng,nt,nloc))
  call hse_kernel_init(op,n,m,h,k,omega,2,ierr,ik0,nloc)
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,1,status)
  if(size(op%phase,2)/=nloc)call MPI_Abort(MPI_COMM_WORLD,2,status)
  call get_command_argument(3,invalid)
  if(trim(invalid)=='nan'.and.rank==0)local_t(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
  if(trim(invalid)=='shape'.and.rank==0)then
    deallocate(local_t);allocate(local_t(ng,nt,0))
  endif
  call hse_kernel_apply_distributed(op,local_u,local_t,local_a,starts,counts,rank,transpose_tiles,ierr)
  if(len_trim(invalid)>0)then
    if(ierr==0)call MPI_Abort(MPI_COMM_WORLD,5,status)
    call hse_kernel_destroy(op)
    call MPI_Finalize(status)
    stop
  endif
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,3,status)
  t=0;t(:,:,ik0:ik0+nloc-1)=local_a
  call MPI_Reduce(t,a,size(t),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,status)
  if(rank==0)then
    open(newunit=iu,file=trim(output),access='stream',form='unformatted',status='replace')
    write(iu)a;close(iu)
  endif
  call hse_kernel_destroy(op)
  call MPI_Finalize(status)
contains
  subroutine transpose_tiles(send,recv,count)
    complex(8),intent(in) :: send(:)
    complex(8),intent(out) :: recv(:)
    integer,intent(in) :: count
    integer :: e
    call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,e)
    if(e/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,4,e)
  end subroutine
end program
