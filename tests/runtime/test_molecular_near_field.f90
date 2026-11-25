! GS-free oracle for the existing isolated dipole-source potential.
program test_molecular_near_field
  use structures, only: s_rgrid,s_rt,s_scalar,s_dft_system
  use salmon_global, only: iperiodic,num_dipole_source,cood_dipole_source, &
    vec_dipole_source,rad_dipole_source,yn_md,yn_out_rvf_rt,tw1,tw2,omega1,phi_cep1,ae_shape1,t1_t2
  use em_field, only: set_vonf,calcVbox
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(s_rgrid) :: mg,lg
  type(s_rt) :: rt
  type(s_scalar) :: potential
  type(s_dft_system) :: system
  real(8) :: d(3),expected,profile,rr,envelope,t,pi,theta,carrier,alpha
  integer :: ix,iy,iz,i,kind,n
  mg%is=1;mg%ie=3;lg%is=1;lg%ie=3;lg%num=3
  allocate(lg%coordinate(1:3,3),rt%vonf%f(1:3,1:3,1:3),potential%f(1:3,1:3,1:3))
  allocate(rt%E_ext(3,0:4),rt%Ac_ext(3,0:4))
  do i=1,3
    lg%coordinate(:,i)=[-1d0,0d0,1d0]
  enddo
  iperiodic=0;yn_md='n';yn_out_rvf_rt='n'
  num_dipole_source=2;rad_dipole_source=.5d0
  cood_dipole_source(:,1)=0d0;cood_dipole_source(:,2)=[1d0,.8d0,-.9d0]
  vec_dipole_source(:,1)=[.3d0,-.2d0,.1d0];vec_dipole_source(:,2)=[-.1d0,.4d0,.2d0]
  rt%time_step=.2d0;rt%E_ext=0d0;rt%Ac_ext=0d0
  do n=0,4
    rt%E_ext(:,n)=[.02d0,-.03d0,.04d0]
  enddo
  tw1=2d0;tw2=0d0;omega1=.7d0;phi_cep1=.13d0;t1_t2=0d0;pi=acos(-1d0)
  call set_vonf(mg,lg,[1d0,1d0,1d0],rt)
  do kind=1,2
    ae_shape1='Ecos2'
    if(kind==2)ae_shape1='Acos2'
    do n=0,4
      t=n*rt%time_step;theta=pi/tw1*(t-.5d0*tw1)
      carrier=omega1*(t-.5d0*tw1)+2d0*pi*phi_cep1;alpha=pi/tw1
      envelope=cos(theta)**2*cos(carrier)
      if(kind==2)envelope=(alpha*sin(2d0*theta)*cos(carrier)+omega1*cos(theta)**2*sin(carrier))/omega1
      call calcVbox(mg,lg,n,system,rt,potential)
      do iz=1,3
      do iy=1,3
      do ix=1,3
        profile=0d0
        do i=1,2
          d=[lg%coordinate(ix,1),lg%coordinate(iy,2),lg%coordinate(iz,3)]-cood_dipole_source(:,i)
          rr=sqrt(sum(d*d))
          profile=profile-dot_product(vec_dipole_source(:,i),d)/max(rr,rad_dipole_source)**3
        enddo
        expected=dot_product(rt%E_ext(:,n),[lg%coordinate(ix,1),lg%coordinate(iy,2),lg%coordinate(iz,3)]) &
                 +profile*envelope
        if(.not.ieee_is_finite(potential%f(ix,iy,iz)))error stop 'nonfinite potential'
        if(abs(rt%vonf%f(ix,iy,iz)-profile)>2d-14)error stop 'source profile mismatch'
        if(abs(potential%f(ix,iy,iz)-expected)>2d-14)error stop 'source pulse or uniform E mismatch'
      enddo
      enddo
      enddo
    enddo
  enddo
  print *, 'PASS: two sources, regularized core, Acos2/Ecos2 and uniform E superposition'
end program
