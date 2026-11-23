! Regression for cancellation in the unscreened WPBEh term of Libxc HSE06.
program test_hse_low_gradient
  use iso_c_binding
  use hybrid_semilocal
  implicit none
  real(c_double) :: rho(8),sigma(8),e(8),v(8),vs(8),g(8)
  real(c_double) :: density(8),omega(3),reduced(10),rr(1),ss(1),ee(1),vv(1),sv(1),k,grad
  real(c_double) :: r2(2,1),gg(2,3,1),core(1),cg(3,1),en(1),pot(2,1),dg(2,3,1)
  real(c_double) :: rot(2,3,1),en2(1),pot2(2,1),dg2(2,3,1)
  real(c_double) :: effective_r(2,1),effective_g(2,3,1),zero_core(1),zero_gradient(3,1)
  real(c_double) :: saved,delta,fd
  character(4096) :: output
  integer :: ierr,i,j,l,m,unit
  rho=.008d0
  g=[0d0,1d-24,1d-18,1d-17,1d-16,1d-14,1d-12,1d-10]
  sigma=g*g
  call hse_semilocal_evaluate(rho,sigma,e,v,vs,ierr)
  if(ierr/=0)stop 1
  do i=1,8
    write(*,'(4es25.16)')g(i),e(i),v(i),vs(i)
  enddo
  ! Analytic WPBEh quadratic limit plus PBE correlation at rho=.008.
  if(maxval(abs(vs-.905679716979567d0))>2d-8)stop 2
  if(maxval(abs(2d0*g*vs))>2d-10)stop 3
  ! Unequal polarized densities, with one finite-gradient spin channel.
  r2(:,1)=[.004d0,.006d0];gg=0d0;gg(1,1,1)=1d-17;gg(2,2,1)=.001d0
  core=0d0;cg=0d0
  call hybrid_density_evaluate(r2,gg,core,cg,en,pot,dg,ierr)
  if(ierr/=0.or.abs(dg(1,1,1))>1d-12)stop 4
  rot(:,1,:)=gg(:,2,:);rot(:,2,:)=gg(:,3,:);rot(:,3,:)=gg(:,1,:)
  call hybrid_density_evaluate(r2,rot,core,cg,en2,pot2,dg2,ierr)
  if(ierr/=0.or.any(en/=en2).or.any(pot/=pot2))stop 5
  if(any(dg(:,1,:)/=dg2(:,3,:)).or.any(dg(:,2,:)/=dg2(:,1,:)))stop 6
  ! Polarized density derivatives, both with and without a fixed NLCC core.
  do i=1,2
    core=.002d0*real(i-1,c_double)
    cg=0d0;cg(1,1)=1d-6*real(i-1,c_double)
    call hybrid_density_evaluate(r2,gg,core,cg,en,pot,dg,ierr)
    if(ierr/=0)stop 13
    effective_r(:,1)=r2(:,1)+core(1)/2d0
    effective_g(:,:,1)=gg(:,:,1)+spread(cg(:,1)/2d0,1,2)
    zero_core=0d0;zero_gradient=0d0
    call hybrid_density_evaluate(effective_r,effective_g,zero_core,zero_gradient,en2,pot2,dg2,ierr)
    if(ierr/=0.or.maxval(abs(en-en2))>1d-15.or.maxval(abs(pot-pot2))>1d-14)stop 14
    if(maxval(abs(dg-dg2))>1d-14)stop 15
    do j=1,2
      call hybrid_density_evaluate(r2,gg,core,cg,en,pot,dg,ierr)
      if(ierr/=0)stop 8
      saved=r2(j,1);delta=saved*1d-4
      r2(j,1)=saved+delta
      call hybrid_density_evaluate(r2,gg,core,cg,en2,pot2,dg2,ierr)
      if(ierr/=0)stop 9
      r2(j,1)=saved-delta
      call hybrid_density_evaluate(r2,gg,core,cg,ee,pot2,dg2,ierr)
      r2(j,1)=saved
      fd=(en2(1)-ee(1))/(2d0*delta)
      if(ierr/=0.or.abs(fd-pot(j,1))>2d-8)stop 10
    enddo
    ! Differentiate energy w.r.t. the finite-gradient spin component.
    delta=gg(2,2,1)*1d-4
    call hybrid_density_evaluate(r2,gg,core,cg,en,pot,dg,ierr)
    gg(2,2,1)=gg(2,2,1)+delta
    call hybrid_density_evaluate(r2,gg,core,cg,en2,pot2,dg2,ierr)
    gg(2,2,1)=gg(2,2,1)-2d0*delta
    call hybrid_density_evaluate(r2,gg,core,cg,ee,pot2,dg2,ierr)
    gg(2,2,1)=gg(2,2,1)+delta
    fd=(en2(1)-ee(1))/(2d0*delta)
    if(ierr/=0.or.abs(fd-dg(2,2,1))>2d-8)stop 11
  enddo
  ! Spin-screening boundary: one channel vanishes or crosses the threshold.
  core=0d0;cg=0d0;gg=0d0
  do i=0,4
    r2(1,1)=real(i,c_double)*.5d-15;r2(2,1)=.004d0
    call hybrid_density_evaluate(r2,gg,core,cg,en,pot,dg,ierr)
    if(ierr/=0.or.maxval(abs(dg))>1d-20)stop 12
  enddo
  call get_command_argument(1,output)
  if(len_trim(output)==0)stop
  open(newunit=unit,file=trim(output),status='new',action='write')
  density=[1d-15,1.999d-15,2.001d-15,1d-12,.001d0,.008d0,.1d0,1d0];omega=[.03d0,.11d0,.3d0]
  reduced=[0d0,1d-18,1d-15,1d-10,1d-6,5d-5,9.99d-5,1.001d-4,1d-3,1d-2]
  do j=1,size(density)
    rr=density(j);k=(3d0*acos(-1d0)**2*rr(1))**(1d0/3d0)
    do l=1,3
      do m=1,size(reduced)
        grad=2d0*k*rr(1)*reduced(m);ss=grad*grad
        call hse_semilocal_evaluate(rr,ss,ee,vv,sv,ierr,screening=omega(l))
        if(ierr/=0)stop 7
        write(unit,'(7es26.17)')rr(1),omega(l),reduced(m),ss(1),ee(1),vv(1),sv(1)
      enddo
    enddo
  enddo
  close(unit)
end program
