program gauge_probe
  use hse_wannier_gauge
  implicit none
  complex(8),allocatable :: u(:,:,:),raw(:,:,:,:),d(:,:,:)
  real(8),allocatable :: b(:,:),w(:)
  integer,allocatable :: neighbors(:,:)
  real(8) :: f,variable
  integer :: iu,n,nk,nb,status
  character(1024) :: path
  call get_command_argument(1,path)
  open(newunit=iu,file=trim(path),form='unformatted',access='stream')
  read(iu)n,nk,nb
  allocate(u(n,n,nk),raw(n,n,nb,nk),d(n,n,nk),neighbors(nb,nk),b(3,nb),w(nb))
  read(iu)u,raw,neighbors,b,w
  close(iu)
  call gauge_functional(u,raw,neighbors,b,w,f,variable,d,status)
  if(status/=0)error stop 'functional failed'
  open(newunit=iu,file=trim(path)//'.out',form='unformatted',access='stream',status='replace')
  write(iu)f,variable,d
  close(iu)
  call check_phase_branch()
contains
  subroutine check_phase_branch()
    implicit none
    complex(8) :: v(1,1,2),links(1,1,2,2),derivative(1,1,2)
    real(8) :: vectors(3,2),weight(2),phases(1,2,2),trial_phases(1,2,2),base_value,trial_value,objective
    real(8) :: pi,epsilon
    integer :: neighbor(2,2),ierr
    pi=acos(-1d0);epsilon=1d-5;v=1d0
    vectors=0d0;vectors(1,:)=[1d0,-1d0];weight=0.5d0;neighbor=reshape([2,2,1,1],[2,2])
    links(:,:,1,:)=exp(cmplx(0d0,pi-epsilon,8))
    links(:,:,2,:)=conjg(links(:,:,1,:))
    call gauge_functional(v,links,neighbor,vectors,weight,base_value,objective,derivative,ierr,phase_out=phases)
    if(ierr/=0)error stop 'initial branch failed'
    v(:,:,2)=exp(cmplx(0d0,2*epsilon,8))
    call gauge_functional(v,links,neighbor,vectors,weight,trial_value,objective,derivative,ierr, &
                           phase_reference=phases,phase_out=trial_phases)
    if(ierr/=0.or.abs(trial_value-base_value)>1d-8)error stop 'phase branch discontinuity'
    if(abs(trial_phases(1,1,1)+trial_phases(1,2,2))>1d-12)error stop 'reverse link branch mismatch'
  end subroutine
end program
