program orbital_initialization
 use initialization_sub,only:s_orbital_initial_config,s_orbital_initial_layout, &
                            orbital_input_config,resolve_orbital_layout
 use salmon_global,only:calc_mode,iperiodic,spin,yn_spinorbit,nstate,nelec,nelec_spin, &
                       temperature,num_kgrid,file_kw
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(s_orbital_initial_config)::a,b,bad
 type(s_orbital_initial_layout)::x,y,saved
 integer::status,i
 calc_mode='GS';iperiodic=3;spin='unpolarized';yn_spinorbit='n'
 nstate=8;nelec=5;nelec_spin=0;temperature=.2d0;num_kgrid=1;file_kw='none'
 call orbital_input_config(a)
 b=a;b%mode='RT';b%temperature=-1d0
 ! Changes to live input cannot choose either captured layout.
 calc_mode='XX';iperiodic=99;spin='polarized';yn_spinorbit='y'
 nstate=1;nelec=99;nelec_spin=99;temperature=99d0;num_kgrid=99;file_kw='must-not-open'
 do i=1,3
  call resolve_orbital_layout(a,x,status)
  if(status/=0.or.x%nspin/=1.or.x%no/=8.or..not.x%real_orbital)stop 1
  saved=x
  call resolve_orbital_layout(b,y,status)
  if(status/=0.or.y%nspin/=1.or.y%no/=3.or.y%real_orbital)stop 2
  if(x%nspin/=saved%nspin.or.x%no/=saved%no.or.(x%real_orbital.neqv.saved%real_orbital))stop 3
 enddo
 b=a;b%kmesh=[2,1,1];call check(b,1,8,.false.)
 b=a;b%kmesh=-1;call check(b,1,8,.false.)
 b=a;b%orthogonal=.false.;call check(b,1,8,.false.)
 b=a;b%periodicity=0;b%orthogonal=.false.;call check(b,1,8,.true.)
 b=a;b%spinorbit='y';call check(b,1,8,.false.)
 b=a;b%hybrid=.true.;call check(b,1,8,.false.)
 b=a;b%mode='RT';b%temperature=-1d0;b%nelec=6;call check(b,1,3,.false.)
 b%nelec=5;b%spinorbit='y';call check(b,1,5,.false.)
 b%spinorbit='n';b%spin='polarized';b%nelec=0;b%nelec_spin=[3,5];call check(b,2,5,.false.)
 b%nelec=5;b%nelec_spin=0;call check(b,2,3,.false.)
 b%temperature=-1d-12;call check(b,2,8,.false.)
 b%temperature=-1.1d-12;call check(b,2,3,.false.)
 b%temperature=0d0;call check(b,2,8,.false.)
 bad=a;bad%mode='XX';call reject(bad)
 bad=a;bad%temperature=ieee_value(0d0,ieee_quiet_nan);call reject(bad)
 bad=a;bad%nstate=0;call reject(bad)
 bad=a;bad%kmesh=[-1,1,1];call reject(bad)
 ! File-derived mesh uses the legacy sentinel without opening a file.
 calc_mode='GS';iperiodic=3;spin='unpolarized';yn_spinorbit='n';nstate=8
 nelec=5;nelec_spin=0;temperature=.2d0;num_kgrid=1;file_kw='not-opened.dat'
 call orbital_input_config(b)
 if(any(b%kmesh/=-1))stop 4
 call check(b,1,8,.false.)
 print *,'Owned orbital initialization PASS'
contains
 subroutine check(value,ns,no,real_orbital)
 type(s_orbital_initial_config),intent(in)::value
 integer,intent(in)::ns,no
 logical,intent(in)::real_orbital
 call resolve_orbital_layout(value,y,status)
 if(status/=0.or.y%nspin/=ns.or.y%no/=no.or.(y%real_orbital.neqv.real_orbital))stop 5
 end subroutine
 subroutine reject(value)
 type(s_orbital_initial_config),intent(in)::value
 call resolve_orbital_layout(value,x,status)
 if(status==0.or.x%nspin/=saved%nspin.or.x%no/=saved%no.or. &
    (x%real_orbital.neqv.saved%real_orbital))stop 6
 end subroutine
end program
