! Exercise the optimized native refresh/action path with independent owners.
program probe
 use iso_fortran_env,only:int64
 use structures
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear,rt_state_bytes
 use exx_native,only:exx_refresh,exx_add_action,exx_bind_system,exx_unbind_system
 use salmon_global
 use sym_sub,only:use_symmetry
 use communication,only:COMM_GROUP_SELF
 implicit none
 type(s_exx_rt_state),target::a,b,reference,reference_b
 type(s_dft_system)::sa,sb,sref,srefb
 type(s_rgrid)::mg
 type(s_parallel_info)::info
 type(s_orbital)::psi,hpsi
 complex(8)::u(8,2,2),v(8,2,2),target(8,2,2),expected(8,2,2),actual(8,2,2),saved(8,2,2)
 integer::i,j,k,status
 integer(int64)::state_bytes
 character(8)::mode,isolated
 xc='hse06';theory='tddft_response';propagator='hse_taylor4';num_kgrid=[2,1,1]
  sa%kmesh=[2,1,1];sa%kpoint_input%num_kgrid=[2,1,1];num_rgrid=2
 yn_periodic='y';yn_spinorbit='n';yn_hse_wannier='n';yn_dc='n';yn_md='n';yn_jm='n'
 yn_symmetrized_stencil='n';yn_conventional_from_dcdft='n';temperature=-1d0;nstate=2;nelec=4
 call get_command_argument(1,mode)
 if(trim(mode)=='full')propagator='hse_taylor4_full'
 exx_factor_mode='off';exx_pre_scf_active=.false.;use_symmetry=.false.
 exx_mlwf_radius=0d0;exx_mlwf_norm_fraction=0d0;exx_kpoint_backend='cpu';hse_omega=.11d0
 hse_q_reduction=1;hse_block_rows=0;hse_fft_layout='strided';pbeh_coulomb_radius=0d0
 yn_exx_k_pair_symmetry='y';yn_exx_fft_orbits='n';yn_hse_sr_pair_symmetry='n'
 mg%is=1;mg%ie=2;mg%num=2
 info%io_s=1;info%io_e=2;info%ik_s=1;info%ik_e=2;info%numo=2;info%numk=2
 info%isize_o=1;info%isize_r=1;info%isize_k=1;info%isize_ro=1;info%numm=1
 info%id_k=0;info%id_r=0;info%id_ro=0;info%icomm_k=COMM_GROUP_SELF;info%icomm_rko=COMM_GROUP_SELF
 info%icomm_r=COMM_GROUP_SELF;info%icomm_o=COMM_GROUP_SELF;info%icomm_ro=COMM_GROUP_SELF
 sa%nspin=1;sa%no=2;sa%nk=2;sa%hvol=1d0;sa%hgs=1d0;sa%primitive_a=0d0;sa%primitive_b=0d0
 do i=1,3
  sa%primitive_a(i,i)=2d0;sa%primitive_b(i,i)=acos(-1d0)
 enddo
 allocate(sa%vec_k(3,2),sa%wtk(2),sa%rocc(2,2,1))
 sa%vec_k=0d0;sa%vec_k(1,2)=acos(-1d0)/2d0;sa%wtk=.5d0;sa%rocc=2d0
 sa%xc_input_ready=.true.;sa%xc_input%xcname='hse06'
 sa%xc_input%exchange_fraction=.25d0;sa%xc_input%screened_exchange=.true.
 sa%xc_input%hse_omega=.11d0
 sb=sa;sref=sa
 sb%xc_input%hse_omega=.17d0
 sb%hgs=1.2d0;sb%hvol=1.2d0**3
 sb%primitive_a=1.2d0*sa%primitive_a;sb%primitive_b=sa%primitive_b/1.2d0
 sb%vec_k=sa%vec_k/1.2d0;srefb=sb
 allocate(psi%zwf(2,2,2,1,2,2,1),hpsi%zwf(2,2,2,1,2,2,1))
 u=0d0;v=0d0
 do k=1,2
  u(1,1,k)=1d0;u(2,2,k)=1d0;v(3,1,k)=1d0;v(4,2,k)=1d0
 enddo
 do k=1,2;do j=1,2;do i=1,8
  target(i,j,k)=cmplx(.01d0*i,.02d0*j*k,8)
 enddo;enddo;enddo
 call get_command_argument(2,isolated)
 call exx_bind_system(sa,a,status);if(status/=0)stop 1
 call exx_bind_system(sb,b,status);if(status/=0)stop 2
 call exx_bind_system(sref,reference,status);if(status/=0)stop 3
 call exx_bind_system(srefb,reference_b,status);if(status/=0)stop 12
 call rt_state_bytes(a,state_bytes,status)
 if(status/=2.or.state_bytes/=-1)stop 14
 call exx_bind_system(sa,b,status);if(status==0)stop 15
 info%isize_k=2;hpsi%zwf=cmplx(.37d0,.19d0,8)
 call exx_add_action(psi,hpsi,sa,mg,info,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 16
 info%isize_k=1
 if(trim(isolated)=='solo_b')then
  call refresh(sb,v/sqrt(sb%hvol));call action(sb,actual)
  call save_action('b_native.dat',actual)
  call exx_unbind_system(sb);call rt_state_clear(b)
  stop
 endif
 call refresh(sa,u);call action(sa,expected)
 call save_action('a_native.dat',expected)
 if(trim(isolated)=='solo_a')then
  call exx_unbind_system(sa);call rt_state_clear(a)
  stop
 endif
 ! Owned inputs must survive unrelated legacy input changes between calls.
 num_kgrid=[9,9,9];use_symmetry=.true.;xc='pz';hse_omega=9d0;pbeh_coulomb_radius=7d0
 call refresh(sb,v/sqrt(sb%hvol));call refresh(srefb,v/sqrt(sb%hvol));call action(sa,actual)
 if(any(actual/=expected))stop 4
 call action(sb,saved)
 call save_action('b_native.dat',saved)
 call action(srefb,actual);if(any(actual/=saved))stop 13
 call refresh(sa,v);call refresh(sref,v)
 call action(sa,actual);call action(sref,expected)
 if(any(actual/=expected))stop 5
 call action(sb,actual);if(any(actual/=saved))stop 6
 if(abs(a%exx_exchange_energy-reference%exx_exchange_energy)>1d-14)stop 7
 call exx_unbind_system(sa);call rt_state_clear(a)
 call action(sb,actual);if(any(actual/=saved))stop 8
 call exx_bind_system(sa,a,status);if(status/=0)stop 9
 call refresh(sa,u);call action(sa,actual)
 call refresh(sref,u);call action(sref,expected)
 if(any(actual/=expected))stop 10
 call exx_unbind_system(sb);call rt_state_clear(b)
 call action(sa,actual);if(any(actual/=expected))stop 11
 call exx_unbind_system(sa);call exx_unbind_system(sref)
 call rt_state_clear(a);call rt_state_clear(reference)
 call exx_unbind_system(srefb);call rt_state_clear(reference_b)
 print *,'Native optimized HSE resources: refresh/action, energy, interleave, close/reopen PASS'
contains
 subroutine save_action(name,x)
 character(*),intent(in)::name
 complex(8),intent(in)::x(:,:,:)
 integer::unit,i,j,k
 open(newunit=unit,file=name,status='new')
 do k=1,size(x,3);do j=1,size(x,2);do i=1,size(x,1)
  write(unit,'(2es26.17)')real(x(i,j,k),8),aimag(x(i,j,k))
 enddo;enddo;enddo
 close(unit)
 end subroutine
 subroutine refresh(s,source)
 type(s_dft_system),intent(in)::s
 complex(8),intent(in)::source(8,2,2)
 psi%zwf=reshape(source,shape(psi%zwf))
 call exx_refresh(s,mg,info,psi,rt_status=status)
 if(status/=0)stop 20
 end subroutine
 subroutine action(s,out)
 type(s_dft_system),intent(in)::s
 complex(8),intent(out)::out(8,2,2)
 psi%zwf=reshape(target,shape(psi%zwf));hpsi%zwf=0d0
 call exx_add_action(psi,hpsi,s,mg,info,rt_status=status)
 if(status/=0)stop 21
 out=reshape(hpsi%zwf,shape(out))
 end subroutine
end program
