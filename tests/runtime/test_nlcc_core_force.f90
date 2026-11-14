! Frozen-density derivative of the actual sampled core interpolation.
program test_nlcc_core_force
  use structures, only: s_pp_info,s_dft_system,s_rgrid,s_pp_nlcc,s_scalar
  use salmon_pp, only: calc_nlcc,calc_nlcc_force
  use salmon_global, only: kion,iperiodic
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_finite
  implicit none
  integer,parameter :: n=6,nr=500,na=2
  real(8),parameter :: dx=.71d0,pi=acos(-1d0),eps=1d-6
  type(s_pp_info) :: pp
  type(s_dft_system) :: sys
  type(s_rgrid) :: rg
  type(s_pp_nlcc) :: ppn
  type(s_scalar) :: vxc(2)
  real(8) :: force(3,na),spinforce(3,na),ep,em,fd,old,err,maxerr,r
  integer :: i,j,k,a,d,geom,status
  pp%flag_nlcc=.true.;pp%nrmax=nr
  allocate(pp%rad(nr,1),pp%rho_nlcc_tbl(nr,1),pp%tau_nlcc_tbl(nr,1))
  do i=1,nr
    r=(i-1)*.012d0
    pp%rad(i,1)=r
    pp%rho_nlcc_tbl(i,1)=exp(-(r/.65d0)**2)
    pp%tau_nlcc_tbl(i,1)=0d0
  end do
  allocate(kion(na));kion=1
  sys%nion=na;sys%Nspin=1;sys%hgs=dx;sys%Hvol=dx**3
  allocate(sys%Rion(3,na))
  rg%is=1;rg%ie=n;rg%num=n
  do a=1,2
    allocate(vxc(a)%f(n,n,n))
    do k=1,n
    do j=1,n
    do i=1,n
      vxc(a)%f(i,j,k)=.3d0+.07d0*sin(2*pi*(i-1)/n)+.05d0*cos(2*pi*(j+k-2)/n)
    end do
    end do
    end do
  end do
  do geom=1,3
    iperiodic=3
    if(geom==1) iperiodic=0
    sys%rmatrix_a=0d0
    do d=1,3
      sys%rmatrix_a(d,d)=1d0
    end do
    if(geom==3)sys%rmatrix_a(2,1)=.19d0 ! Lower-triangular skew, not a cuboid.
    sys%primitive_a=sys%rmatrix_a*(n*dx)
    sys%Rion(:,1)=[.314d0,.217d0,.173d0]
    sys%Rion(:,2)=[1.213d0,.873d0,2.197d0]
    call calc_nlcc_force(pp,sys,rg,vxc,force,status)
    if(status/=0)stop 10
    maxerr=0d0
    do a=1,na
    do d=1,3
      old=sys%Rion(d,a)
      sys%Rion(d,a)=old+eps
      call calc_nlcc(pp,sys,rg,ppn)
      ep=sum(ppn%rho_nlcc*vxc(1)%f)*sys%Hvol
      sys%Rion(d,a)=old-eps
      call calc_nlcc(pp,sys,rg,ppn)
      em=sum(ppn%rho_nlcc*vxc(1)%f)*sys%Hvol
      sys%Rion(d,a)=old
      fd=-(ep-em)/(2*eps)
      err=abs(fd-force(d,a));maxerr=max(maxerr,err)
      if(.not.ieee_is_finite(fd) .or. .not.ieee_is_finite(force(d,a)))stop 11
      if(err>1d-7)then
        write(*,*) 'CORE_FORCE_FAIL',geom,a,d,fd,force(d,a),err
        stop 12
      end if
    end do
    end do
    write(*,*) 'CORE_FORCE_VARIATION',geom,maxerr
    vxc(2)%f=2d0*vxc(1)%f
    sys%Nspin=2
    call calc_nlcc_force(pp,sys,rg,vxc,spinforce,status)
    if(status/=0 .or. maxval(abs(1.5d0*force-spinforce))>1d-14)stop 13
    sys%Nspin=1
  end do
  old=vxc(1)%f(1,1,1)
  vxc(1)%f(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=3 .or. maxval(abs(force))/=0d0)stop 14
  vxc(1)%f(1,1,1)=ieee_value(0d0,ieee_positive_inf)
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=3)stop 15
  vxc(1)%f(1,1,1)=old
  pp%rad(2,1)=pp%rad(1,1)
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=2)stop 16
  pp%rad(2,1)=.012d0
  old=pp%rho_nlcc_tbl(5,1)
  pp%rho_nlcc_tbl(5,1)=ieee_value(0d0,ieee_quiet_nan)
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=2)stop 17
  pp%rho_nlcc_tbl(5,1)=old
  vxc(1)%f=huge(0d0)
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=4 .or. maxval(abs(force))/=0d0)stop 20
  vxc(1)%f=.5d0*vxc(2)%f
  sys%Rion(:,1)=0d0
  call calc_nlcc_force(pp,sys,rg,vxc,force,status)
  if(status/=0)stop 18
  do a=1,na
  do d=1,3
    if(.not.ieee_is_finite(force(d,a)))stop 19
  end do
  end do
  write(*,*) 'NLCC_CORE_FORCE_PASS'
end program
