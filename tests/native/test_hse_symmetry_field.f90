! Copyright 2026 SALMONpy contributors. Apache-2.0.
program owned_field
 use structures,only:s_dft_system
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear
 use exx_symmetry,only:s_exx_symmetry_map,symmetry_init,symmetry_check_field
 use exx_native,only:exx_bind_system,exx_unbind_system,exx_validate_vector_potential
 use sym_sub,only:use_symmetry,SymMatA,SymMatB,Bmat,Binv
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 type(s_dft_system)::a,b,c,unbound
 type(s_exx_rt_state),target::oa,ob,oc
 type(s_exx_symmetry_map)::before,bad,empty
 real(8)::rotation(3,3),tilt(3,3),s,pi,nan,inf
 integer::i,status
 integer,target::wrong_owner
 pi=acos(-1d0);s=sqrt(.5d0)
 rotation=0d0
 do i=1,3
  rotation(i,i)=1d0
 enddo
 tilt=0d0;tilt(:,1)=[s,0d0,-s];tilt(:,2)=[0d0,1d0,0d0];tilt(:,3)=[s,0d0,s]
 call setup(a,oa,rotation,.false.)
 call setup(b,ob,rotation,.true.)
 call setup(c,oc,tilt,.false.)
 ! Wrong shared symmetry must not choose any owned result.
 use_symmetry=.false.;Bmat=99d0;Binv=-99d0
 allocate(SymMatA(3,4,2),SymMatB(3,4,2));SymMatA=99d0;SymMatB=-99d0
 before=oa%native_symmetry
 do i=1,3
  call check(a,[0d0,0d0,1d0],0)
  call check(a,[1d0,0d0,0d0],2)
  call check(a,[0d0,1d0,0d0],2)
  call check(b,[0d0,1d0,0d0],0)
  call check(b,[1d0,0d0,0d0],2)
  call check(c,[1d0,0d0,1d0],0)
  call check(c,[0d0,0d0,1d0],2)
  call check(a,[0d0,0d0,0d0],0)
  call check(a,[0d0,0d0,1d-250],0)
  call check(a,[1d-250,0d0,0d0],2)
  call check(a,[0d0,0d0,1d300],0)
  call check(a,[1d300,0d0,0d0],2)
  call check(a,[4d-11,0d0,1d0],0)
  call check(a,[6d-11,0d0,1d0],2)
 enddo
 nan=ieee_value(0d0,ieee_quiet_nan);inf=ieee_value(0d0,ieee_positive_inf)
 call check(a,[nan,0d0,0d0],1);call check(a,[0d0,inf,0d0],1)
 call symmetry_check_field(empty,[0d0,0d0,0d0],status)
 if(status/=1)stop 4
 bad=before;bad%max_little=0;call bad_map(bad)
 bad=before;deallocate(bad%cart_rotation);call bad_map(bad)
 bad=before;deallocate(bad%cart_rotation);allocate(bad%cart_rotation(3,2,2));bad%cart_rotation=0d0
 call bad_map(bad)
 bad=before;deallocate(bad%cart_rotation);allocate(bad%cart_rotation(0:2,0:2,0:1));bad%cart_rotation=0d0
 call bad_map(bad)
 bad=before;bad%cart_rotation(1,1,1)=nan;call bad_map(bad)
 bad=before;bad%cart_rotation=0d0;call bad_map(bad)
 if(any(oa%native_symmetry%cart_rotation/=before%cart_rotation))stop 5
 if(any(oa%native_symmetry%first/=before%first))stop 6
 ! A mutable model must not silently bypass the initialized map.
 a%kpoint_input%use_symmetry=.false.;call check(a,[0d0,0d0,0d0],1);a%kpoint_input%use_symmetry=.true.
 a%primitive_a(1,1)=5d0;call check(a,[0d0,0d0,1d0],1);a%primitive_a(1,1)=4d0
 a%kpoint_input%symmat_b(1,1,2)=1d0;call check(a,[0d0,0d0,1d0],1)
 a%kpoint_input%symmat_b(1,1,2)=-1d0
 unbound%kpoint_input%use_symmetry=.true.;call check(unbound,[0d0,0d0,0d0],1)
 wrong_owner=7;unbound%salmonpy_hse_owner=>wrong_owner
 call check(unbound,[0d0,0d0,0d0],1);nullify(unbound%salmonpy_hse_owner)
 call exx_unbind_system(a);call rt_state_clear(oa)
 call check(a,[0d0,0d0,0d0],1)
 call check(b,[0d0,1d0,0d0],0)
 call exx_unbind_system(b);call exx_unbind_system(c);call rt_state_clear(ob);call rt_state_clear(oc)
 print *,'Owned HSE field little group, scale/tolerance, rotated cell, invalid/closed state PASS'
contains
 subroutine setup(system,owner,basis,mirror)
 type(s_dft_system),intent(out)::system
 type(s_exx_rt_state),target,intent(inout)::owner
 real(8),intent(in)::basis(3,3)
 logical,intent(in)::mirror
 real(8)::fractional(3),weights(4)
 integer::d,x,y,z,k
 system%ngrid=64;system%nspin=1;system%no=2;system%nk=4;system%kmesh=2
 system%primitive_a=4d0*basis;system%primitive_b=(2*pi/4d0)*basis;system%hgs=1d0
 system%kpoint_input%use_symmetry=.true.;system%kpoint_input%num_kgrid=2
 allocate(system%kpoint_input%symmat_a(3,4,2),system%kpoint_input%symmat_b(3,4,2))
 system%kpoint_input%symmat_a=0d0
 do d=1,3
  system%kpoint_input%symmat_a(d,d,1)=1d0;system%kpoint_input%symmat_a(d,d,2)=1d0
 enddo
 system%kpoint_input%symmat_a(1,1,2)=-1d0
 if(.not.mirror)system%kpoint_input%symmat_a(2,2,2)=-1d0
 system%kpoint_input%symmat_b=system%kpoint_input%symmat_a
 allocate(system%vec_k(3,4));k=0
 do z=0,1;do y=0,1
  k=k+1;fractional=[-.25d0,-.25d0+.5d0*y,-.25d0+.5d0*z]
  system%vec_k(:,k)=matmul(system%primitive_b,fractional)
 enddo;enddo
 weights=.25d0
 call exx_bind_system(system,owner,status);if(status/=0)stop 1
 call symmetry_init(owner%native_symmetry,[4,4,4],system%hgs,system%vec_k,weights, &
                    system%kpoint_input%symmat_a,system%kpoint_input%symmat_b,2,status,cell=system%primitive_a)
 if(status/=0)stop 2
 end subroutine
 subroutine check(system,ac,expected)
 type(s_dft_system),intent(in)::system
 real(8),intent(in)::ac(3)
 integer,intent(in)::expected
 call exx_validate_vector_potential(system,ac,status)
 if(status/=expected)then
  print *,'field status mismatch',status,expected,ac;stop 3
 endif
 end subroutine
 subroutine bad_map(map)
 type(s_exx_symmetry_map),intent(in)::map
 call symmetry_check_field(map,[0d0,0d0,0d0],status)
 if(status/=1)stop 7
 end subroutine
end program
