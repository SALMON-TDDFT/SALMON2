! Fixed-H real-data propagation check only. NOT a TDHSE dielectric calculation.
program lcfo_frozen_probe
  use lcfo_rt_core
  use iso_fortran_env,only:int32
  implicit none
  complex(8),allocatable::h(:,:),c(:,:),next(:,:),initial(:,:),expected(:,:),gram(:,:)
  real(8),allocatable::eval(:)
  integer(int32)::header(3)
  integer::iu,n,no,step,status,j
  real(8)::dt,gram_err,phase_err,reversal_err
  character(1024)::path
  call get_command_argument(1,path)
  open(newunit=iu,file=trim(path),form='unformatted',access='stream',status='old')
  read(iu)header
  if(header(1)/=16909060.or.any(header(2:3)<1))stop 1
  n=header(2);no=header(3)
  allocate(h(n,n),c(n,no),next(n,no),initial(n,no),expected(n,no),eval(no),gram(no,no))
  read(iu)h,c,eval;close(iu)
  initial=c;dt=.02d0
  do step=1,4
    call lcfo_cayley_step(h,c,dt,next,status)
    if(status/=0)stop 2
    c=next
  end do
  expected=initial
  do j=1,no
    expected(:,j)=expected(:,j)*((1d0-cmplx(0d0,dt*eval(j)/2,8))/ &
                                       (1d0+cmplx(0d0,dt*eval(j)/2,8)))**4
  end do
  phase_err=maxval(abs(c-expected))
  gram=matmul(conjg(transpose(c)),c)
  do j=1,no
    gram(j,j)=gram(j,j)-1d0
  end do
  gram_err=maxval(abs(gram))
  do step=1,4
    call lcfo_cayley_step(h,c,-dt,next,status)
    if(status/=0)stop 3
    c=next
  end do
  reversal_err=maxval(abs(c-initial))
  write(*,'(a,i0)') 'dimension=',n
  write(*,'(a,i0)') 'states=',no
  write(*,'(a,es24.15)') 'analytic_cayley_phase_error=',phase_err
  write(*,'(a,es24.15)') 'orthogonality_error=',gram_err
  write(*,'(a,es24.15)') 'time_reversal_error=',reversal_err
  if(max(phase_err,gram_err,reversal_err)>1d-10)stop 4
  print *, 'PASS frozen Hamiltonian only; self-consistent RT not validated'
end program
