program test_dg_wpw_core_point_builder
  use dg_wpw_core_point_builder,only:evaluate_dg_wpw_core_point
  implicit none
  integer,parameter::nk=2,ng=2
  integer::core_lo(3,nk),core_hi(3,nk),fragment_ids(nk),column_ids(nk*ng),info
  integer::total_grid(3),buffer(3),width(3),grid(3)
  real(8)::hgs(3),gvec(3,ng),chi(nk),grad_chi(3,nk),omega,position(3)
  complex(8)::p(nk*ng),grad_p(3,nk*ng),phase,expected

  fragment_ids=[1,2]
  core_lo(:,1)=[1,1,1];core_hi(:,1)=[4,1,1]
  core_lo(:,2)=[5,1,1];core_hi(:,2)=[8,1,1]
  total_grid=[8,1,1];hgs=[0.5d0,1d0,1d0];buffer=[2,0,0];width=[2,0,0]
  gvec(:,1)=[0d0,0d0,0d0];gvec(:,2)=[0.4d0,-0.2d0,0.1d0]
  grid=[4,1,1];omega=16d0;position=dble(grid)*hgs
  call evaluate_dg_wpw_core_point(core_lo,core_hi,fragment_ids,total_grid,hgs,buffer,width,gvec,grid,omega,&
    column_ids,chi,grad_chi,p,grad_p,info)
  if(info/=0)error stop 1
  if(any(column_ids/=[1,2,3,4]))error stop 2
  if(abs(sum(chi**2)-1d0)>1d-13)error stop 3
  phase=exp(cmplx(0d0,sum(gvec(:,2)*position),8))/sqrt(omega)
  expected=chi(2)*phase
  if(abs(p(4)-expected)>1d-13)error stop 4
  if(maxval(abs(grad_p(:,4)-(cmplx(grad_chi(:,2),0d0,8)+cmplx(0d0,1d0,8)*gvec(:,2)*chi(2))*phase))>1d-13)&
    error stop 5
  call evaluate_dg_wpw_core_point(core_lo,core_hi,[2,1],total_grid,hgs,buffer,width,gvec,grid,omega,&
    column_ids,chi,grad_chi,p,grad_p,info)
  if(info==0)error stop 6
  print '(a)','PASS support-local stable-column WPW core point builder'
end program test_dg_wpw_core_point_builder
