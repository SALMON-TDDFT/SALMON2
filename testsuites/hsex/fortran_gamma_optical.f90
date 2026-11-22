program test_gamma_optical
  use hsex_longwave, only: longwave_cell,longwave_lowrank_cell
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  complex(8) :: vertices(4,5),raw(5,5),chi(5,5),basis(5,5),work(16)
  complex(8) :: point(8,8),average(8,8),p(8,8),a(8,8)
  real(8) :: rec(3,3),vol,weight(5),values(5),rwork(16),macro(3,3),m(3,3),signs(5)
  integer :: g(3,2),reverse(5),i,j,status,info
  vol=1000d0;rec=0d0
  do i=1,3
    rec(i,i)=2d0*acos(-1d0)/10d0
  enddo
  g(:,1)=[1,0,0];g(:,2)=-g(:,1);reverse=[1,2,3,5,4];signs=[-1d0,-1d0,-1d0,1d0,1d0]
  do j=1,5
    do i=1,4
      vertices(i,j)=cmplx(sin(dble(i*j)),cos(dble(i+2*j)),8)
    enddo
  enddo
  raw=matmul(conjg(transpose(vertices)),vertices)
  do j=1,5
    do i=1,5
      chi(i,j)=raw(i,j)+signs(i)*signs(j)*conjg(raw(reverse(i),reverse(j)))
    enddo
  enddo
  weight(:3)=sqrt(4d0*acos(-1d0));weight(4:)=weight(1)/rec(1,1)
  do j=1,5
    basis(:,j)=2d0/vol*weight*chi(:,j)*weight(j)
  enddo
  call zheev('V','L',5,basis,5,values,work,size(work),rwork,info)
  if(info/=0)error stop 'optical fixture eigensolver failure'
  call longwave_cell([4,4,4],2,g,rec,vol,1,chi,macro,point,average,status)
  if(status/=0)error stop 'optical full-matrix reference failure'
  call longwave_lowrank_cell([4,4,4],2,g,rec,vol,1,basis,values,m,p,a,status)
  if(status/=0)error stop 'factored optical cell rejected'
  if(maxval(abs(m-macro))>1d-11.or.maxval(abs(p-point))>1d-11.or. &
     maxval(abs(a-average))>1d-11)error stop 'factored head/wings differ from full optical cell'
  values(5)=ieee_value(0d0,ieee_quiet_nan)
  call longwave_lowrank_cell([4,4,4],2,g,rec,vol,1,basis,values,m,p,a,status)
  if(status==0)error stop 'NaN low-rank response accepted'
  values(5)=ieee_value(0d0,ieee_positive_inf)
  call longwave_lowrank_cell([4,4,4],2,g,rec,vol,1,basis,values,m,p,a,status)
  if(status==0)error stop 'infinite low-rank response accepted'
  values=0d0;values(1)=-1d0
  call longwave_lowrank_cell([4,4,4],2,g,rec,vol,1,basis,values,m,p,a,status)
  if(status==0)error stop 'negative polarization accepted'
  print *, 'GAMMA_PDEP_OPTICAL_WINGS_IEEE_PASS'
end program
