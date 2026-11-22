program test_q2_head
  use hsex_screen_math, only: screen_sites,screen_q2_head
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter :: nq=8**3
  integer :: sites(3,nq),mesh(3),i,status,j
  real(8) :: rec(3,3),q(3),q2,qc,pi,eta,a,c,kappa,reference
  complex(8) :: values(2,2,nq),coefficient(2,2),average(2,2),expected(2,2)
  pi=acos(-1d0);mesh=8;rec=0d0
  do i=1,3
    rec(i,i)=.4d0
  enddo
  qc=(3d0/(4d0*pi)*abs(rec(1,1)*rec(2,2)*rec(3,3))/nq)**(1d0/3d0)
  call screen_sites(mesh,sites)
  ! Exact quadratic q^2 W including an imaginary odd wing. The wing must not
  ! pollute the spherical average; the body retains its regular contribution.
  values=0d0
  do i=1,nq
    q=matmul(rec,dble(sites(:,i))/mesh);q2=sum(q**2)
    if(q2==0d0)cycle
    values(1,1,i)=.2d0/q2+2d0
    values(2,2,i)=3d0
    values(1,2,i)=cmplx(0d0,.01d0*q(1)/q2,8)
    values(2,1,i)=conjg(values(1,2,i))
  enddo
  call screen_q2_head(values,mesh,rec,mesh,1d-8,qc,coefficient,average,status)
  if(status/=0)error stop 'q2 head fit failed'
  expected=0d0;expected(1,1)=.6d0/qc**2+2d0;expected(2,2)=3d0
  if(maxval(abs(average-expected))>1d-9)error stop 'q2 cell moment or wing error'
  if(abs(coefficient(1,1)-.2d0)>1d-12.or.abs(coefficient(2,2))>1d-12) &
    error stop 'q2 singular coefficient error'
  ! Small-q hydrodynamic metal: q^2 W=A(eta)+C(eta)q^2+O(q^4).
  ! As eta -> 0, A vanishes but C -> 4*pi/kappa^2. A head-only closure
  ! would incorrectly set the complete metallic cell contribution to zero.
  kappa=1.3d0
  do j=1,2
    eta=merge(.03d0,0d0,j==1)
    a=4d0*pi*eta**2/(eta**2+kappa**2)
    c=4d0*pi*kappa**2/(eta**2+kappa**2)**2
    do i=1,nq
      q=matmul(rec,dble(sites(:,i))/mesh);q2=sum(q**2)
      if(q2==0d0)cycle
      values(:,:,i)=0d0;values(1,1,i)=a/q2+c;values(2,2,i)=c
    enddo
    call screen_q2_head(values,mesh,rec,mesh,1d-8,qc,coefficient,average,status)
    if(status/=0)error stop 'metal q2 closure failed'
    reference=3d0*a/qc**2+c
    if(abs(average(1,1)-reference)>1d-9)error stop 'metal regular term lost'
    print '(a,3es20.10)','Q2_HEAD_METAL eta/cell/regular ',eta,real(average(1,1),8),c
  enddo
  call screen_q2_head(values,[2,2,2],rec,[2,2,2],.05d0,qc,coefficient,average,status)
  if(status==0)error stop 'insufficient small q sampling accepted'
  call screen_q2_head(values,mesh,rec,mesh,.05d0,ieee_value(0d0,ieee_quiet_nan),coefficient,average,status)
  if(status==0)error stop 'NaN cell radius accepted'
  do i=1,nq
    q=matmul(rec,dble(sites(:,i))/mesh);q2=sum(q**2)
    if(q2==0d0)cycle
    values(1,1,i)=values(1,1,i)-.01d0/q2
  enddo
  call screen_q2_head(values,mesh,rec,mesh,.05d0,qc,coefficient,average,status)
  if(status==0)error stop 'negative singular coefficient accepted'
  print *, 'Q2_HEAD_SPHERE_WINGS_METAL_CAP_IEEE_PASS'
end program
