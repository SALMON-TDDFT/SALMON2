program diagonal_probe
  use kli_diagonal, only: solve_kli
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  real(8):: p(2,2),pq(2,2),d(2),v(2),r,c,e
  real(8):: pg(2,6),qg(2,6),dg(6),vg(2)
  real(8):: chain(2,3),qc(2,3),dc(3)
  integer:: s,i
  ! Orthogonal eigenvectors (sqrt(.8),sqrt(.2)),(-sqrt(.2),sqrt(.8)).
  ! Orbital-dependent potentials (.5,0), (0,1); M12=.32, b1=.08.
  p=reshape([.8d0,.2d0,.2d0,.8d0],[2,2]);pq=reshape([.4d0,0d0,0d0,.8d0],[2,2])
  call run
  if(s/=0.or.maxval(abs(d-[.125d0,-.125d0]))>1d-12)error stop 1
  if(maxval(abs(v-[.475d0,.725d0]))>1d-12)error stop 2
  print *, 'D, row error, compatibility, residual:',d,r,c,e
  ! Three k replicas exercise the smaller grid-space solve (6 > 2*2).
  do i=1,3
    pg(:,2*i-1:2*i)=p;qg(:,2*i-1:2*i)=pq
  enddo
  call grid
  if(s/=0.or.maxval(abs(vg-v))>1d-12)error stop 6
  do i=1,3
    if(maxval(abs(dg(2*i-1:2*i)-d))>1d-12)error stop 7
  enddo
  qg(2,6)=ieee_value(0d0,ieee_quiet_nan)
  call grid
  if(s==0)error stop 8
  qg(2,6)=ieee_value(0d0,ieee_positive_inf)
  call grid
  if(s==0)error stop 9
  pq(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call run
  if(s==0)error stop 3
  pq(1,1)=ieee_value(0d0,ieee_positive_inf)
  call run
  if(s==0)error stop 4
  p=reshape([1d0,0d0,0d0,1d0],[2,2]);pq=0d0
  call run
  if(s==0)error stop 5
  ! Zero pair overlap does not imply disconnection: the middle state
  ! connects the two endpoints. Exercise the spectral fallback.
  chain=reshape([1d0,0d0,.5d0,.5d0,0d0,1d0],[2,3])
  do i=1,3
    qc(:,i)=i*chain(:,i)
  enddo
  call solve_kli(chain,qc,1,1d0,1d-10,dc,v,r,c,e,s)
  if(s/=0.or.maxval(abs(dc-[1d0,0d0,-1d0]))>1d-12.or.maxval(abs(v-2d0))>1d-12)error stop 10
  ! Positive but numerically disconnected: retain the eigenvalue cutoff.
  p=reshape([1d0-1d-12,1d-12,1d-12,1d0-1d-12],[2,2])
  pq=0d0
  call run
  if(s==0)error stop 11
contains
  subroutine grid
    implicit none
    call solve_kli(pg,qg,3,1d0,1d-10,dg,vg,r,c,e,s)
  end subroutine
  subroutine run
    implicit none
    call solve_kli(p,pq,1,1d0,1d-10,d,v,r,c,e,s)
  end subroutine
end program
