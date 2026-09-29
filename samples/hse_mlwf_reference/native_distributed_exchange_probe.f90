program probe
  use, intrinsic :: ieee_arithmetic
  use omp_lib, only: omp_get_max_threads,omp_get_max_active_levels
  use mpi
  use exx_k_exchange
  implicit none
  type(exx_k_kernel) :: op
  integer :: rank,np,status,n,m,no,nt,nk,ng,ik0,nloc,iu,ierr,p,provided,expected_block,chosen_block,old_levels,callback_count
  integer,allocatable :: starts(:),counts(:)
  real(8) :: h,omega
  real(8),allocatable :: k(:,:)
  complex(8),allocatable :: u(:,:,:),t(:,:,:),a(:,:,:),local_u(:,:,:),local_t(:,:,:),local_a(:,:,:),callback_a(:,:,:)
  character(1024) :: input,output,invalid,setting
  character(64) :: layout
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,status)
  if(provided<MPI_THREAD_FUNNELED)call MPI_Abort(MPI_COMM_WORLD,7,status)
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
  call get_command_argument(3,invalid)
  chosen_block=2
  if(trim(invalid)=='block'.and.rank==0)chosen_block=3
  expected_block=0
  call get_command_argument(4,setting)
  if(len_trim(setting)>0)then
    read(setting,*,iostat=status)expected_block
    if(status/=0)call MPI_Abort(MPI_COMM_WORLD,20,status)
  endif
  layout='auto'
  call get_command_argument(5,setting)
  if(len_trim(setting)>0)layout=trim(setting)
  call exx_k_kernel_init(op,n,m,h,k,omega,chosen_block,ierr,ik0,nloc, &
    block_rows=expected_block,profile=.true.,fft_layout=layout)
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,1,status)
  if(expected_block==0)expected_block=chosen_block
  if(op%block/=expected_block)call MPI_Abort(MPI_COMM_WORLD,8,status)
  if(.not.op%profile)call MPI_Abort(MPI_COMM_WORLD,21,status)
  if(op%auto_fft.neqv.(trim(layout)=='auto'))call MPI_Abort(MPI_COMM_WORLD,14,status)
  if(trim(layout)/='auto')then
    if(op%contiguous_fft.neqv.(trim(layout)=='contiguous'))call MPI_Abort(MPI_COMM_WORLD,11,status)
  endif
  if(op%auto_fft)then
    if(any(op%fft_trial_seconds<0d0).or..not.all(ieee_is_finite(op%fft_trial_seconds)))then
      call MPI_Abort(MPI_COMM_WORLD,12,status)
    endif
    if(op%contiguous_fft.neqv.(minval(op%fft_trial_seconds)>1d-5.and. &
       op%fft_trial_seconds(2)<0.95d0*op%fft_trial_seconds(1))) &
      call MPI_Abort(MPI_COMM_WORLD,13,status)
  endif
  if(size(op%phase,2)/=nloc)call MPI_Abort(MPI_COMM_WORLD,2,status)
  call get_command_argument(3,invalid)
  if(trim(invalid)=='nan'.and.rank==0)local_t(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
  if(trim(invalid)=='shape'.and.rank==0)then
    deallocate(local_t);allocate(local_t(ng,nt,0))
  endif
  old_levels=omp_get_max_active_levels()
  call exx_k_kernel_apply_distributed(op,local_u,local_t,local_a,starts,counts,rank,transpose_tiles,ierr)
  if(len_trim(invalid)>0)then
    if(ierr==0)call MPI_Abort(MPI_COMM_WORLD,5,status)
    call exx_k_kernel_destroy(op)
    call MPI_Finalize(status)
    stop
  endif
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,3,status)
  if(omp_get_max_active_levels()/=old_levels)call MPI_Abort(MPI_COMM_WORLD,16,status)
  allocate(callback_a(ng,nt,nloc));callback_count=0
  call exx_k_kernel_apply_distributed(op,local_u,local_t,callback_a,starts,counts,rank,transpose_tiles, &
    ierr,fill_density)
  if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,17,status)
  if(maxval(abs(callback_a-local_a))>1d-10)call MPI_Abort(MPI_COMM_WORLD,18,status)
  if(callback_count/=nloc*((ng+np*op%block-1)/(np*op%block)))call MPI_Abort(MPI_COMM_WORLD,19,status)
  if(omp_get_max_active_levels()/=old_levels)call MPI_Abort(MPI_COMM_WORLD,16,status)
  if(op%profile)then
    if(any(op%seconds<0d0).or..not.all(ieee_is_finite(op%seconds)))call MPI_Abort(MPI_COMM_WORLD,9,status)
    if(sum(op%seconds)<=0d0)call MPI_Abort(MPI_COMM_WORLD,10,status)
  endif
  if(op%threads_used/=omp_get_max_threads())call MPI_Abort(MPI_COMM_WORLD,6,status)
  t=0;t(:,:,ik0:ik0+nloc-1)=local_a
  call MPI_Reduce(t,a,size(t),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,status)
  if(rank==0)then
    open(newunit=iu,file=trim(output),access='stream',form='unformatted',status='replace')
    write(iu)a;close(iu)
  endif
  call exx_k_kernel_destroy(op)
  call MPI_Finalize(status)
contains
  subroutine fill_density(j,lo,rows,density)
    integer,intent(in) :: j,lo,rows
    complex(8),intent(out) :: density(:,:)
    complex(8) :: phased(ng,no)
    integer :: c
    do c=1,no
      phased(:,c)=local_u(:,c,j)*op%phase(:,j)
    enddo
    density=0d0
    density(1:rows,:)=matmul(phased(lo:lo+rows-1,:),transpose(conjg(phased)))
    !$omp atomic update
    callback_count=callback_count+1
  end subroutine
  subroutine transpose_tiles(send,recv,count)
    complex(8),intent(in) :: send(:)
    complex(8),intent(out) :: recv(:)
    integer,intent(in) :: count
    integer :: e
    call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,e)
    if(e/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,4,e)
  end subroutine
end program
