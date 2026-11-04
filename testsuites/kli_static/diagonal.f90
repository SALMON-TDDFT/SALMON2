program diagonal_probe
  use kli_diagonal
  use ieee_arithmetic
  implicit none
  real(8):: p(2,2),pq(2,2),d(2),v(2),r,c,e
  integer:: s
  ! Orthogonal eigenvectors (sqrt(.8),sqrt(.2)),(-sqrt(.2),sqrt(.8)).
  ! Orbital-dependent potentials (.5,0), (0,1); M12=.32, b1=.08.
  p=reshape([.8d0,.2d0,.2d0,.8d0],[2,2]);pq=reshape([.4d0,0d0,0d0,.8d0],[2,2])
  call run
  if(s/=0.or.maxval(abs(d-[.125d0,-.125d0]))>1d-12)error stop 1
  if(maxval(abs(v-[.475d0,.725d0]))>1d-12)error stop 2
  print *, 'D, row error, compatibility, residual:',d,r,c,e
  pq(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call run
  if(s==0)error stop 3
  pq(1,1)=ieee_value(0d0,ieee_positive_inf)
  call run
  if(s==0)error stop 4
  p=reshape([1d0,0d0,0d0,1d0],[2,2]);pq=0d0
  call run
  if(s==0)error stop 5
contains
  subroutine run
    call solve_kli(p,pq,1,1d0,1d-10,d,v,r,c,e,s)
  end subroutine
end program
