program test_frozen_q_mpi
  use mpi
  use hsex_frozen_q, only: frozen_q_probe
  implicit none
  integer,parameter :: ng=216,nk=2
  integer :: rank,np,ierr,provided,ks,ke,status,i,unit,s,t,g,ref,rl,rr,actions,rows
  real(8) :: cell(3,3),k(3,nk),vals(7),seconds,density_error
  real(8) :: summed_error(2,3),trace_error
  complex(8),allocatable :: states(:,:,:)
  complex(8) :: probes(3)
  integer :: qp(3)
  character(256) :: line
  character(32) :: mode
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  ks=1+nk*rank/np;ke=nk*(rank+1)/np
  allocate(states(ng,1,max(0,ke-ks+1)));states=1d0/sqrt(dble(ng))
  mode='valid'
  if(command_argument_count()>0)call get_command_argument(1,mode)
  if(mode=='existing-output'.and.rank==0)then
    open(newunit=unit,file='hsex-frozen-q.dat',status='new')
    close(unit)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  cell=0d0
  do i=1,3
    cell(i,i)=12d0
  enddo
  k=0d0;k(1,:)=[-.25d0,.25d0]*2d0*acos(-1d0)/12d0
  qp=[1,0,0]
  if(mode=='frequency-q0')qp=0
  if(mode=='frequency'.or.mode=='frequency-q0')then
    probes=[cmplx(0d0,.05d0,8),cmplx(0d0,.025d0,8),cmplx(.12d0,.025d0,8)]
    call frozen_q_probe([6,6,6],[2,1,1],[2,1,1],qp,cell,k,states,1,ks,ke, &
      MPI_COMM_WORLD,rank,np,ng,1d-9,.05d0,.true.,action,status, &
      summed_error=summed_error,trace_error=trace_error,probe_shifts=probes)
  else
    call frozen_q_probe([6,6,6],[2,1,1],[2,1,1],[1,0,0],cell,k,states,1,ks,ke, &
      MPI_COMM_WORLD,rank,np,ng,1d-9,.05d0,.true.,action,status, &
      summed_error=summed_error,trace_error=trace_error)
  endif
  if(mode=='bad-action'.or.mode=='existing-output')then
    if(status==0)error stop 'remote H error accepted'
    if(rank==0)print *,'FROZEN_Q_MPI_REMOTE_FAILURE_REJECTED'
    call MPI_Finalize(ierr)
    stop
  endif
  if(status/=0)error stop 'distributed frozen probe rejected'
  if(maxval(summed_error)<0d0.or.maxval(summed_error)>1d-7)error stop 'k-summed density reference mismatch'
  if(trace_error<0d0.or.trace_error>1d-10)error stop 'k-fold or physical density normalization mismatch'
  if(rank==0)then
    if(mode=='frequency-q0')then
      open(newunit=unit,file='hsex-frozen-q0.dat',status='old')
    else
      open(newunit=unit,file='hsex-frozen-q.dat',status='old')
    endif
    read(unit,'(a)')line
    rows=0
    do
      read(unit,*,iostat=ierr)s,t,g,i,ref,vals,rl,rr,actions,seconds,density_error
      if(ierr/=0)exit
      rows=rows+1
      if(ref/=1.or.vals(5)>1d-7.or.vals(6)>1d-9)error stop 'reference error or false convergence'
      if(density_error>1d-7)error stop 'induced-density reference mismatch'
    enddo
    close(unit)
    if(rows/=merge(18,12,index(mode,'frequency')==1))error stop 'missing MPI probe rows'
    print *,'FROZEN_Q_MPI_OWNERS_AND_REFERENCE_PASS'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine action(index,x,y,code)
    implicit none
    integer,intent(in) :: index
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    integer :: j
    real(8) :: phase
    code=1
    if(index<ks.or.index>ke)return
    if(mode=='bad-action'.and.index==nk)return
    phase=.01d0*index
    do j=1,ng
      y(j)=.2d0*x(j)-.1d0*(exp(cmplx(0d0,phase,8))*x(1+mod(j,ng)) &
        +exp(cmplx(0d0,-phase,8))*x(1+modulo(j-2,ng)))
    enddo
    code=0
  end subroutine
end program
