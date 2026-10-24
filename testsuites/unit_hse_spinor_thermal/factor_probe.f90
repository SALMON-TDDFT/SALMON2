program factor_thermal_probe
 use mpi
 use structures
 use salmon_global
 use exx_native
 use exx_k_exchange,only:exx_k_kernel,exx_k_kernel_init
 use exx_spinor,only:exx_spinor_weight_sources,exx_spinor_apply_distributed
 use sym_sub,only:use_symmetry
 use plusU_global,only:PLUS_U_ON
 implicit none
 type(s_dft_system)::sys
 type(s_rgrid)::mg
 type(s_parallel_info)::info
 type(s_orbital)::psi,hp
 type(exx_k_kernel)::op
 complex(8)::u(16,3,2),q(16,3,2),ref(16,3,2),w(16,3,2),old(16,3,2)
 real(8)::expected,oldenergy,pi,residual
 integer::err,status,i,j,k,route,step
 character(16)::mode
 call MPI_Init(status)
 call get_command_argument(1,mode)
 xc='hse06';yn_periodic='y';yn_spinorbit='y';yn_dc='n';yn_md='n';yn_jm='n'
 yn_symmetrized_stencil='n';yn_hse_wannier='n';yn_hse_wannier_snapshot='n';yn_conventional_from_dcdft='n'
 yn_exx_spinor_pair_fft='n';yn_exx_k_pair_symmetry='n';yn_exx_fft_orbits='n';yn_hse_profile='n'
 exx_pre_scf_active=.false.;PLUS_U_ON=.false.;use_symmetry=.false.;temperature=.001d0
 nelec=1;nstate=3;num_kgrid=[2,1,1];num_rgrid=2;hse_q_reduction=1;hse_omega=.11d0;pbeh_coulomb_radius=0d0
 exx_mlwf_radius=0d0;exx_mlwf_norm_fraction=0d0;hse_sr_tolerance=0d0;hse_block_rows=0
 exx_kpoint_backend='cpu';exx_k_algorithm='reference';hse_fft_layout='strided';theory='dft'
 propagator='hse_taylor4';if(trim(mode)=='direct')propagator='hse_taylor4_full'
 mg%is=1;mg%ie=2;mg%num=2
 info%isize_r=1;info%isize_o=1;info%isize_k=1;info%numm=1;info%numk=2;info%numo=3
 info%id_k=0;info%ik_s=1;info%ik_e=2;info%io_s=1;info%io_e=3
 info%icomm_k=MPI_COMM_SELF;info%icomm_rko=MPI_COMM_SELF
 sys%nk=2;sys%no=3;sys%nspin=2;sys%hvol=.343d0
 sys%primitive_a=reshape([0d0,1.4d0,1.4d0,1.4d0,0d0,1.4d0,1.4d0,1.4d0,0d0],[3,3])
 pi=acos(-1d0);sys%primitive_b=pi/1.4d0
 do i=1,3;sys%primitive_b(i,i)=-pi/1.4d0;enddo
 sys%hgs=sqrt(2d0)*.7d0;sys%hvol=.7d0**3*2d0
 allocate(sys%vec_k(3,2),sys%wtk(2),sys%rocc(3,2,2))
 sys%vec_k=0d0;sys%vec_k(:,2)=.5d0*sys%primitive_b(:,1);sys%wtk=.5d0
 sys%rocc(:,1,1)=[.65d0,.35d0,0d0];sys%rocc(:,2,1)=[.8d0,.2d0,0d0];sys%rocc(:,:,2)=sys%rocc(:,:,1)
 allocate(psi%zwf(2,2,2,2,3,2,1),hp%zwf(2,2,2,2,3,2,1))
 do k=1,2;do j=1,3;do i=1,16
  u(i,j,k)=cmplx(sin(real(i+3*j+k,8)),cos(real(2*i-j+k,8)),8)/5d0
 enddo;enddo;enddo
 call exx_unpack(u,psi,mg,info)
 call exx_k_kernel_init(op,mg%num,num_kgrid,sys%hgs,sys%vec_k,hse_omega,3,err, &
  cell=sys%primitive_a,reciprocal=sys%primitive_b,fft_layout='strided')
 if(err/=0)error stop 'reference initialization'
 call reference_action(ref)
 call exx_refresh(sys,mg,info,psi)
 call check_action(ref,'initial native')
 expected=exchange_energy(ref)
 if(abs(exx_exchange_energy-expected)>1d-11)error stop 'native energy'
 theory='tddft_pulse';yn_predictor_corrector='y';exx_factor_mode='learned'
 exx_ace_support='occupied';exx_pair_screening='off'
 exx_factor_exact_interval=2;exx_factor_n_teach=2;exx_factor_warmup_steps=4
 info%icomm_r=MPI_COMM_SELF;info%icomm_o=MPI_COMM_SELF;info%icomm_ro=MPI_COMM_SELF
 call exx_factor_begin(1,sys,mg,info,psi)
 call check_action(ref,'fractional spinor factor startup')
 ! Give the regression a nonsingular teaching signal, then retain three
 ! identical endpoints. The ensuing prediction must equal physical exchange.
 exx_factor_warmup_steps=20
 do step=2,20
  call exx_factor_begin(step,sys,mg,info,psi)
  if(mod(step,2)==0)then
   if(step<=6)then
    call exx_unpack(u*(1d0+.01d0*step),psi,mg,info)
   else
    call exx_unpack(u,psi,mg,info)
   endif
   call exx_refresh(sys,mg,info,psi)
   call exx_factor_endpoint(sys,info,step)
  endif
 enddo
 exx_factor_warmup_steps=0
 call exx_factor_begin(21,sys,mg,info,psi)
 if(.not.exx_freeze)error stop 'fixture did not enter learned prediction'
 if(trim(mode)=='changed')then
  sys%rocc(1,1,1)=sys%rocc(1,1,1)-.01d0
  sys%rocc(2,1,1)=sys%rocc(2,1,1)+.01d0
  call exx_refresh(sys,mg,info,psi)
  error stop 'changed occupation accepted'
 endif
 call exx_refresh(sys,mg,info,psi)
 call check_action(ref,'predicted spinor all-target action')
 if(abs(exx_exchange_energy-expected)>1d-10)error stop 'predicted spinor energy'
 print *,'PASS thermal spinor actual learned prediction, action and energy'
 call MPI_Finalize(status)

contains
 subroutine reference_action(action)
 complex(8),intent(out)::action(:,:,:)
 call exx_spinor_weight_sources(u,sys%rocc(:,:,1),q,err)
 if(err/=0)error stop 'source weighting'
 call exx_spinor_apply_distributed(op,q,u,action,[1],[2],0,0,transpose_local,err)
 if(err/=0)error stop 'reference Fock'
 end subroutine
 real(8) function exchange_energy(action) result(en)
 complex(8),intent(in)::action(:,:,:)
 en=0d0
 do k=1,2
 en=en+.5d0*.25d0*sys%hvol*sys%wtk(k)*sum(sum(real(conjg(u(:,:,k))*action(:,:,k),8),dim=1)*sys%rocc(:,k,1))
 enddo
 end function
 subroutine check_action(expected,label)
 complex(8),intent(in)::expected(:,:,:)
 character(*),intent(in)::label
 hp%zwf=0d0
 call exx_add_action(psi,hp,sys,mg,info)
 call exx_pack(hp,mg,info,w)
 if(maxval(abs(w-.25d0*expected))>1d-10)then
 print *,label,maxval(abs(w-.25d0*expected))
 error stop 'native action'
 endif
 end subroutine
 subroutine transpose_local(send,recv,count,comm)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count,comm
 recv=send
 end subroutine
end program
