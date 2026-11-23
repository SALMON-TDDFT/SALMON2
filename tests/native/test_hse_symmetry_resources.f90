! Fast-PR exchange: independent reduced-k owners versus the same full-k kernel.
! Copyright 2026 SALMONpy contributors. Apache-2.0.
program symmetry_resources
 use structures
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use occupation,only:occupation_input_config
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear
 use exx_native,only:exx_bind_system,exx_bind_gs_system,exx_unbind_system,exx_set_phase,exx_refresh,exx_add_action
 use salmon_global
 use sym_sub,only:use_symmetry
 use communication,only:COMM_GROUP_SELF
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global
 implicit none
 type(s_exx_rt_state),target::oa,ob,ofa,ofb
 type(s_dft_system)::a,b,fa,fb
 type(s_rgrid)::mg
 type(s_parallel_info)::info
 complex(8),allocatable::ar(:,:,:),br(:,:,:),af(:,:,:),bf(:,:,:),saved(:,:,:)
 real(8)::ea,eb,efa,efb,err,energy_err,pi
 integer::status,i,j,k,t,rep(4),unit
 character(12)::mode,poison,algorithm
 call setup_parallel
 call get_command_argument(1,mode);call get_command_argument(2,poison);call get_command_argument(3,algorithm)
 pi=acos(-1d0);err=0d0;energy_err=0d0
 xc='hse06';theory='tddft_response';propagator='hse_taylor4'
 if(mode=='full')propagator='hse_taylor4_full'
 yn_periodic='y';yn_spinorbit='n';yn_hse_wannier='n';yn_dc='n';yn_md='n';yn_jm='n'
 yn_symmetrized_stencil='n';yn_conventional_from_dcdft='n';temperature=-1d0;nstate=2;nelec=4
 exx_factor_mode='off';exx_pre_scf_active=.false.;use_symmetry=.false.
 exx_mlwf_radius=0d0;exx_mlwf_norm_fraction=0d0;exx_kpoint_backend='cpu';hse_omega=.11d0
 hse_q_reduction=1;hse_block_rows=0;hse_fft_layout='strided';pbeh_coulomb_radius=0d0
 yn_exx_k_pair_symmetry='y';yn_exx_fft_orbits='n';yn_hse_sr_pair_symmetry='n'
 if(algorithm=='orbits')yn_exx_fft_orbits='y'
 num_kgrid=2;num_rgrid=4
 mg%is=1;mg%ie=4;mg%num=4
 info%io_s=1;info%io_e=2;info%numo=2
 info%isize_o=1;info%isize_r=1;info%isize_k=1;info%isize_ro=1;info%numm=1
 info%id_k=0;info%id_r=0;info%id_ro=0
 info%icomm_k=COMM_GROUP_SELF;info%icomm_rko=COMM_GROUP_SELF
 info%icomm_r=COMM_GROUP_SELF;info%icomm_o=COMM_GROUP_SELF;info%icomm_ro=COMM_GROUP_SELF
 call configure(a,.true.,1d0+.05d0*nproc_id_global,.11d0+.01d0*nproc_id_global,.false.)
 call configure(b,.true.,1.2d0+.07d0*nproc_id_global,.17d0+.01d0*nproc_id_global,.true.)
 call configure(fa,.false.,1d0+.05d0*nproc_id_global,.11d0+.01d0*nproc_id_global,.false.)
 call configure(fb,.false.,1.2d0+.07d0*nproc_id_global,.17d0+.01d0*nproc_id_global,.true.)
 theory='dft';yn_opt='n'
 call exx_bind_gs_system(a,info,oa,status)
 if(status/=2.or.associated(a%salmonpy_hse_owner).or.oa%native_resources_enabled)stop 25
 theory='tddft_response'
 call exx_bind_system(a,oa,status);if(status/=0)stop 1
 call exx_bind_system(b,ob,status);if(status/=0)stop 2
 call exx_bind_system(fa,ofa,status);if(status/=0)stop 3
 call exx_bind_system(fb,ofb,status);if(status/=0)stop 4
 ! Full mesh is independently enumerated with x fastest; reps have negative x.
 rep=[1,3,5,7]
 call evaluate(fa,.false.,0d0,af,efa)
 call evaluate(fb,.true.,0d0,bf,efb)
 if(poison=='poison')then
  num_rgrid=99;num_kgrid=99;use_symmetry=.true.;xc='pz';hse_omega=9d0
  nelec=99;temperature=.5d0;nstate=0
 endif
 if(nproc_id_global==0)then
  open(newunit=unit,file='full-reference.bin',access='stream',form='unformatted',status='new')
  write(unit)efa,efb,af,bf
  close(unit)
 endif
 call evaluate(a,.false.,0d0,ar,ea)
 call evaluate(b,.true.,0d0,br,eb)
 call compare(ar,af,ea,efa);call compare(br,bf,eb,efb)
 call reject_grid(a)
 call action(a,.false.,ar);call compare(ar,af,ea,efa)
 if(oa%native_symmetry%nfull/=8.or.oa%native_symmetry%nrep/=4)stop 5
 if(ob%native_symmetry%nfull/=8.or.ob%native_symmetry%nrep/=4)stop 6
 if(all(oa%native_symmetry%rotation==ob%native_symmetry%rotation))stop 7
 ! Alter the second source; first owner/action must remain unchanged.
 saved=ar
 call evaluate(fb,.true.,.23d0,bf,efb)
 call evaluate(b,.true.,.23d0,br,eb)
 call compare(br,bf,eb,efb)
 call action(a,.false.,ar)
 if(any(ar/=saved))stop 8
 ! Retiring A must not retire B's map/phase/FFT resources.
 saved=br
 call exx_unbind_system(a);call rt_state_clear(oa)
 call action(b,.true.,br)
 if(any(br/=saved))stop 9
 call exx_bind_system(a,oa,status);if(status/=0)stop 10
 call evaluate(a,.false.,0d0,ar,ea);call compare(ar,af,ea,efa)
 call exx_unbind_system(b);call rt_state_clear(ob)
 call action(a,.false.,ar);call compare(ar,af,ea,efa)
 ! A different per-process schedule detects unintended COMM_WORLD collectives.
 if(nproc_id_global/=0)then
  call evaluate(a,.false.,.07d0,ar,ea)
  call evaluate(fa,.false.,.07d0,af,efa)
  call compare(ar,af,ea,efa)
 endif
 call exx_unbind_system(a);call exx_unbind_system(fa);call exx_unbind_system(fb)
 call rt_state_clear(oa);call rt_state_clear(ofa);call rt_state_clear(ofb)
 print *,'max_action_error=',err,' max_energy_error=',energy_err
 print *,'Fast-PR HSE symmetry/full-k, independent owners, refresh/close/reopen PASS'
 call end_parallel
contains
 subroutine configure(s,reduced,spacing,omega,mirror)
 type(s_dft_system),intent(out)::s
 logical,intent(in)::reduced,mirror
 real(8),intent(in)::spacing,omega
 integer::x,y,z,p,d
 s%nspin=1;s%no=2;s%nk=8;s%nion=1;s%ngrid=64
 if(reduced)s%nk=4
 s%hgs=spacing;s%hvol=spacing**3;s%primitive_a=0d0;s%primitive_b=0d0
 s%vec_Ac=0d0
 do d=1,3
  s%primitive_a(d,d)=4d0*spacing;s%primitive_b(d,d)=2d0*pi/(4d0*spacing)
 enddo
 allocate(s%Rion(3,1),s%kion(1));s%Rion=0d0;s%kion=1
 s%kmesh=2;s%kpoint_input%num_kgrid=2;s%kpoint_input%use_symmetry=reduced
 if(reduced)then
  allocate(s%kpoint_input%symmat_a(3,4,2),s%kpoint_input%symmat_b(3,4,2))
  s%kpoint_input%symmat_a=0d0
  do d=1,3
   s%kpoint_input%symmat_a(d,d,1)=1d0;s%kpoint_input%symmat_a(d,d,2)=1d0
  enddo
  s%kpoint_input%symmat_a(1,1,2)=-1d0
  if(.not.mirror)s%kpoint_input%symmat_a(2,2,2)=-1d0
  s%kpoint_input%symmat_b=s%kpoint_input%symmat_a
 endif
 allocate(s%vec_k(3,s%nk),s%wtk(s%nk),s%rocc(2,s%nk,1))
 p=0
 do z=0,1;do y=0,1;do x=0,1
  if(reduced.and.x==1)cycle
  p=p+1;s%vec_k(:,p)=[(-.25d0+.5d0*x)*s%primitive_b(1,1), &
                    (-.25d0+.5d0*y)*s%primitive_b(2,2),(-.25d0+.5d0*z)*s%primitive_b(3,3)]
 enddo;enddo;enddo
 s%wtk=1d0/s%nk;s%rocc=2d0
 s%xc_input_ready=.true.;s%xc_input%xcname='hse06'
 s%xc_input%exchange_fraction=.25d0;s%xc_input%screened_exchange=.true.;s%xc_input%hse_omega=omega
 call occupation_input_config(s%occupation_input);s%occupation_input_ready=.true.
 call exx_set_phase(s,.false.)
 end subroutine
 subroutine layout(s)
 type(s_dft_system),intent(in)::s
 info%ik_s=1;info%ik_e=s%nk;info%numk=s%nk
 end subroutine
 subroutine evaluate(s,mirror,angle,out,energy)
 type(s_dft_system),intent(in)::s
 logical,intent(in)::mirror
 real(8),intent(in)::angle
 complex(8),allocatable,intent(out)::out(:,:,:)
 real(8),intent(out)::energy
 type(s_orbital)::psi
 integer::x,y,z,k,q
 real(8)::v,normalization
 call layout(s)
 allocate(psi%zwf(4,4,4,1,2,s%nk,1));normalization=1d0/sqrt(64d0*s%hvol)
 do k=1,s%nk;do z=0,3;do y=0,3;do x=0,3
  q=z;if(mirror)q=y
  v=2*pi*q/4
  psi%zwf(x+1,y+1,z+1,1,1,k,1)=normalization
  psi%zwf(x+1,y+1,z+1,1,2,k,1)=normalization*(cos(angle)*exp(cmplx(0d0,v,8)) &
                                                    +sin(angle)*exp(cmplx(0d0,2*v,8)))
 enddo;enddo;enddo;enddo
 call exx_refresh(s,mg,info,psi,rt_status=status)
 if(status/=0)then
  print *,'refresh rejection',status,s%nk,s%hgs
  stop 20
 endif
 select type(owner=>s%salmonpy_hse_owner)
 type is(s_exx_rt_state)
 energy=owner%exx_exchange_energy
 class default
 stop 21
 end select
 call action(s,mirror,out)
 end subroutine
 subroutine action(s,mirror,out)
 type(s_dft_system),intent(in)::s
 logical,intent(in)::mirror
 complex(8),allocatable,intent(out)::out(:,:,:)
 type(s_orbital)::psi,hpsi
 integer::x,y,z,k,q
 call layout(s)
 allocate(psi%zwf(4,4,4,1,2,s%nk,1),hpsi%zwf(4,4,4,1,2,s%nk,1))
 do k=1,s%nk;do z=0,3;do y=0,3;do x=0,3
  q=z;if(mirror)q=y
  psi%zwf(x+1,y+1,z+1,1,1,k,1)=cmplx(.01d0*(1+x+4*y+16*z),.02d0*q,8)
  psi%zwf(x+1,y+1,z+1,1,2,k,1)=exp(cmplx(0d0,2*pi*q/4d0,8))+.13d0
 enddo;enddo;enddo;enddo
 hpsi%zwf=0d0
 call exx_add_action(psi,hpsi,s,mg,info,rt_status=status)
 if(status/=0)stop 22
 out=reshape(hpsi%zwf,[64,2,s%nk])
 end subroutine
 subroutine reject_grid(s)
 type(s_dft_system),intent(in)::s
 type(s_rgrid)::wrong
 type(s_orbital)::psi,hpsi
 integer::i
 character(256)::saved_propagator
 call layout(s)
 allocate(psi%zwf(4,4,4,1,2,s%nk,1),hpsi%zwf(4,4,4,1,2,s%nk,1))
 psi%zwf=cmplx(.13d0,.27d0,8)
 ! Equal cell count, wrong per-axis dimensions must still reject.
 saved_propagator=propagator
 do i=1,4
  wrong=mg;info%isize_k=1;exx_factor_mode='off';propagator=saved_propagator
  hpsi%zwf=cmplx(.37d0,.19d0,8)
  select case(i)
  case(1)
   wrong%num=[4,2,8];wrong%ie=wrong%num
  case(2)
   info%isize_k=3
  case(3)
   exx_factor_mode='invalid'
  case(4)
   propagator='unsupported'
  end select
  call exx_refresh(s,wrong,info,psi,rt_status=status)
  if(status/=2)stop 23
  call exx_add_action(psi,hpsi,s,wrong,info,rt_status=status)
  if(status/=2.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 24
 enddo
 info%isize_k=1;exx_factor_mode='off';propagator=saved_propagator
 end subroutine
 subroutine compare(reduced,full,energy,full_energy)
 complex(8),intent(in)::reduced(:,:,:),full(:,:,:)
 real(8),intent(in)::energy,full_energy
 real(8)::difference
 integer::x,y,z
 if(.not.ieee_is_finite(energy).or..not.ieee_is_finite(full_energy))stop 32
 do z=1,size(full,3);do y=1,size(full,2);do x=1,size(full,1)
  if(.not.ieee_is_finite(real(full(x,y,z),8)).or..not.ieee_is_finite(aimag(full(x,y,z))))stop 33
 enddo;enddo;enddo
 do z=1,size(reduced,3);do y=1,size(reduced,2);do x=1,size(reduced,1)
  if(.not.ieee_is_finite(real(reduced(x,y,z),8)).or..not.ieee_is_finite(aimag(reduced(x,y,z))))stop 34
 enddo;enddo;enddo
 difference=maxval(abs(reduced-full(:,:,rep)))
 err=max(err,difference);energy_err=max(energy_err,abs(energy-full_energy))
 if(difference>2d-11*max(1d0,maxval(abs(full))))then
  print *,'action mismatch',difference;stop 30
 endif
 if(abs(energy-full_energy)>2d-11*max(1d0,abs(full_energy)))then
  print *,'energy mismatch',energy,full_energy;stop 31
 endif
 end subroutine
end program
