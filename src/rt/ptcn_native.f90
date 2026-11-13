#include "config.h"
! Native HSE PT-CN adapter. Full local real-space domain.
module ptcn_native
 implicit none
 private
 public::native_ptcn_step
contains
 subroutine native_ptcn_step(dt,ac0,ac1,lg,mg,system,info,stencil,xc_func,srg,srg_scalar,pp,ppg,ppn, &
   input,output,rho,rho_s,vlocal,vh,vxc,vpsl,fg,poisson,energy)
 use structures
 use ptcn_core,only:exx_ptcn_solve,pt_orthogonalize,pt_uniform_occupation
#ifdef USE_SCALAPACK
 use ptcn_native_matrix
#endif
#ifdef USE_HSE
 use exx_native,only:exx_freeze,exx_ptcn_begin,exx_ptcn_predict_inner,exx_ptcn_accept
#endif
 use salmon_global,only:xc,yn_spinorbit
 use ptcn_rectangular_fft
 use ptcn_kinetic_preconditioner
 use omp_lib,only:omp_get_wtime
 use communication,only:comm_get_max,comm_summation,comm_is_root
 use density_matrix,only:calc_density
 use hartree_sub,only:hartree
 use salmon_xc,only:exchange_correlation
 use hamiltonian,only:hpsi,update_vlocal,update_kvector_nonlocalpt
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(8),intent(in)::dt,ac0(3),ac1(3)
 type(s_rgrid),intent(in)::lg,mg
 type(s_dft_system),intent(inout)::system
 type(s_parallel_info),intent(in)::info
 type(s_stencil),intent(in)::stencil
 type(s_xc_functional),intent(in)::xc_func
 type(s_sendrecv_grid),intent(inout)::srg,srg_scalar
 type(s_pp_info),intent(in)::pp
 type(s_pp_grid),intent(inout)::ppg
 type(s_pp_nlcc),intent(in)::ppn
 type(s_orbital),intent(in)::input
 type(s_orbital),intent(inout)::output
 type(s_scalar),intent(inout)::rho,rho_s(:),vlocal(:),vh,vxc(:)
 type(s_scalar),intent(in)::vpsl
 type(s_reciprocal_grid),intent(inout)::fg
 type(s_poisson),intent(inout)::poisson
 type(s_dft_energy),intent(inout)::energy
 type(s_orbital)::trial,htrial
#ifdef USE_SCALAPACK
 type(ptcn_owner_matrix)::owner
#endif
 type(ptcn_fft_plan)::fft
 complex(8),allocatable::u(:,:,:),x(:,:,:),den(:,:),pre_in(:,:),pre_out(:,:)
 complex(8),allocatable::initial_hu(:,:,:),rhs_cn(:,:,:),work_hu(:,:,:)
 real(8)::residual,gram,gram_k,ne_local,ne,expected,qvec(3),kvec(3),symbol,fmetric(6),pre_seconds,tick
 real(8)::gram_before,ne_before,post_residual,correction_norm,lowdin_seconds
 integer::ng,no,nk,builds,apps,err,k,j,ik,io,is,g,ix,iy,iz,ns,status
 logical::first_action,hybrid
 character(12)::label
 if(info%isize_r/=1.or.any(mg%num/=lg%num))error stop 'Native PT-CN requires full local spatial domain'
 hybrid=xc=='hse06'
 if(.not.hybrid)error stop 'Native PT-CN: unsupported functional'
 label='HSE_PT_CN'
#ifndef USE_HSE
 if(hybrid)error stop 'Native PT-CN: HSE not built'
#endif
 if(info%numm/=1)error stop 'Native PT-CN requires one system'
 if(hybrid.and.system%nspin/=1)error stop 'Restored HSE Lowdin requires the certified scalar occupation route'
 ! RT initialization truncates the zero-temperature scalar system to occupied
 ! orbitals (nelec/2), independently of the GS nstate. Full-column Lowdin is
 ! valid here only when every propagated orbital has the same occupation 2.
 if(hybrid)then
  status=0
  if(.not.pt_uniform_occupation(system%rocc,2d0))status=1
  call comm_get_max(status,info%icomm_rko)
  if(status/=0)error stop 'Standard HSE PT-CN requires occupied-only scalar RT orbitals'
  if(comm_is_root(info%id_rko))write(*,'(a,i0,a)') &
   'HSE_PT_CN_OCCUPIED_ONLY states=',system%no,' occupation=2'
 endif
 ns=system%nspin
 if(ns/=1.and.(ns/=2.or.yn_spinorbit/='y'))error stop 'Native PT-CN requires scalar or SOC spinor orbitals'
#ifndef USE_SCALAPACK
 if(info%isize_o/=1.or.info%numo/=system%no.or.info%io_s/=1) &
   error stop 'Native PT-CN orbital distribution requires ScaLAPACK'
#endif
 ng=product(mg%num);no=info%numo;nk=info%numk
 call allocate_orbital_complex(ns,mg,info,trial)
 call allocate_orbital_complex(ns,mg,info,htrial)
 allocate(u(ns*ng,no,nk),x(ns*ng,no,nk),den(ng,nk),pre_in(ng,ns*no),pre_out(ng,ns*no))
 call ptcn_fft_create(fft,mg%num,ns*no,err)
 status=err;call comm_get_max(status,info%icomm_rko)
 if(status/=0)error stop 'Native PT-CN FFT setup failed'
 fmetric=stencil%coef_F
 if(stencil%if_orthogonal)fmetric=[1d0,1d0,1d0,0d0,0d0,0d0]
 do k=1,nk
  ik=k+info%ik_s-1;kvec=system%vec_k(:,ik)+.5d0*(ac0+ac1)
  g=0
  do iz=0,mg%num(3)-1;do iy=0,mg%num(2)-1;do ix=0,mg%num(1)-1
   g=g+1;qvec=2*acos(-1d0)*real([ix,iy,iz],8)/mg%num
   call ptcn_kinetic_symbol(qvec,stencil%coef_lap0,stencil%coef_lap,stencil%coef_nab, &
     fmetric,system%rmatrix_B,kvec,symbol,err)
   if(err/=0)error stop 'Native PT-CN kinetic symbol rejected'
   den(g,k)=1+(0d0,.5d0)*dt*symbol
  enddo;enddo;enddo
 enddo
 pre_seconds=0d0
#ifdef USE_SCALAPACK
 call ptcn_owner_create(owner,ns*ng,system%no,info%io_s,info%io_e,info%icomm_o,err)
#else
 err=0
#endif
 status=err;call comm_get_max(status,info%icomm_rko)
 if(status/=0)error stop 'Native PT-CN owner setup failed'
 call pack(input,u)
#ifdef USE_HSE
 if(hybrid)call exx_ptcn_begin(dt)
#endif
 first_action=.true.
 if(hybrid)then
  allocate(initial_hu(ns*ng,no,nk),rhs_cn(ns*ng,no,nk),work_hu(ns*ng,no,nk))
  ! Keep the exact initial RHS for an independent post-retraction residual.
  ! This is the same starting callback the solver otherwise executes, not a duplicate.
  call action(u,initial_hu,.true.)
  call project(u,initial_hu,system%hvol,work_hu,system%no)
  rhs_cn=u-(0d0,.5d0)*dt*work_hu
  call exx_ptcn_solve(u,dt,system%hvol,action,precondition,global_norm,x,residual,builds,apps,err, &
   noccupied=system%no,use_anderson=.true.,initial_hu=initial_hu, &
   projection=project,consensus=consensus,mixing_reduce=reduce_gram)
  builds=builds+1 ! Account for the exact RHS callback performed above.
 else
  call exx_ptcn_solve(u,dt,system%hvol,action,precondition,global_norm,x,residual,builds,apps,err, &
   noccupied=system%no,use_anderson=.true.,projection=project,consensus=consensus,mixing_reduce=reduce_gram)
 endif
#ifdef USE_HSE
 if(hybrid)exx_freeze=.false.
#endif
 status=err;call comm_get_max(status,info%icomm_rko)
 if(status/=0)then
  if(comm_is_root(info%id_rko))write(*,*)trim(label)//'_REJECT residual/builds/apps',residual,builds,apps
  error stop 'Native PT-CN nonlinear convergence failed'
 endif
 ! Preserve the original raw Gram/electron gate: do not hide large drift.
 call state_metrics(x,gram,ne,status)
 if(status/=0)error stop 'Native PT-CN raw Gram/electron gate failed'
 if(hybrid)then
  gram_before=gram;ne_before=ne;initial_hu=x
  tick=omp_get_wtime()
#ifdef USE_SCALAPACK
  do k=1,nk
   call ptcn_owner_lowdin(owner,x(:,:,k),system%hvol,gram_k,err)
   if(err/=0)exit
  enddo
#else
  call pt_orthogonalize(x,system%hvol,err)
#endif
  lowdin_seconds=omp_get_wtime()-tick
  status=abs(err);call comm_get_max(status,info%icomm_rko)
  if(status/=0)error stop 'Native HSE PT-CN Lowdin failed'
  correction_norm=global_norm(x-initial_hu)/global_norm(u)
  call state_metrics(x,gram,ne,status)
  if(status/=0)error stop 'Native HSE PT-CN post-Lowdin Gram/electron gate failed'
  ! Retraction is after the CN nonlinear solve, as in the original HSE route.
  ! Rebuild H/ACE at the ACTUAL output state, before committing an exact teacher.
  call action(x,work_hu,.true.)
  call project(x,work_hu,system%hvol,initial_hu,system%no)
  work_hu=x+(0d0,.5d0)*dt*initial_hu-rhs_cn
  post_residual=global_norm(work_hu)/global_norm(u)
  status=0
  if(.not.ieee_is_finite(post_residual).or..not.ieee_is_finite(correction_norm))status=1
  call comm_get_max(status,info%icomm_rko)
  if(status/=0)error stop 'Native HSE post-Lowdin residual is nonfinite'
  if(comm_is_root(info%id_rko))write(*,'(5(a,es22.13))') &
   'HSE_PT_CN_ORTH pre_gram=',gram_before,' post_gram=',gram, &
   ' pre_residual=',residual,' post_residual=',post_residual,' correction=',correction_norm
  if(comm_is_root(info%id_rko))write(*,'(a,es22.13,a)') &
   'HSE_PT_CN_ORTH seconds=',lowdin_seconds,' post_fresh_actions=1'
 endif
#ifdef USE_HSE
 if(hybrid)call exx_ptcn_accept(system,info)
#endif
 call unpack(x,output)
 if(comm_is_root(info%id_rko))write(*,'(a,es22.13,a,i0,a,i0,a,es22.13,a,es22.13)') &
  trim(label)//' residual=',residual,' builds=',builds,' applications=',apps,' Gram=',gram,' electrons=',ne
 if(comm_is_root(info%id_rko))write(*,'(a,es22.13)')trim(label)//'_PRECONDITION seconds=',pre_seconds
 call ptcn_fft_destroy(fft)
#ifdef USE_SCALAPACK
 call ptcn_owner_destroy(owner)
#endif
 call deallocate_orbital(trial);call deallocate_orbital(htrial)
 contains
 subroutine state_metrics(a,gr,number,code)
 complex(8),intent(in)::a(:,:,:)
 real(8),intent(out)::gr,number
 integer,intent(out)::code
 real(8)::gk,local_number,target,gmax(1)
 integer::kk,jj,global_k,global_o,ec
 gr=0d0;local_number=0d0;code=0
 do kk=1,size(a,3)
  call check_gram(a(:,:,kk),gk,ec)
  if(.not.ieee_is_finite(gk))ec=1
  code=max(code,ec);gr=max(gr,gk);global_k=kk+info%ik_s-1
  do jj=1,no
   global_o=jj+info%io_s-1
   local_number=local_number+sum(abs(a(:,jj,kk))**2)*system%hvol* &
    system%rocc(global_o,global_k,1)*system%wtk(global_k)
  enddo
 enddo
 call comm_summation(local_number,number,info%icomm_ko)
 call comm_get_max([gr],gmax,1,info%icomm_rko)
 gr=gmax(1)
 target=sum(system%rocc(:,:,1)*spread(system%wtk,1,system%no))
 if(.not.ieee_is_finite(number).or..not.ieee_is_finite(gr))code=1
 if(gr>1d-8.or.abs(number-target)>1d-7)code=1
 call comm_get_max(code,info%icomm_rko)
 end subroutine
 subroutine check_gram(a,error,code)
 complex(8),intent(in)::a(:,:)
 real(8),intent(out)::error
 integer,intent(out)::code
#ifndef USE_SCALAPACK
 complex(8)::matrix(size(a,2),size(a,2))
 integer::j
 matrix=system%hvol*matmul(conjg(transpose(a)),a)
 do j=1,size(a,2);matrix(j,j)=matrix(j,j)-1;enddo
 error=maxval(abs(matrix));code=0
 if(.not.ieee_is_finite(error))code=1
#else
 call ptcn_owner_gram(owner,a,system%hvol,error,code)
#endif
 end subroutine
 subroutine reduce_gram(matrix)
 complex(8),intent(inout)::matrix(:,:)
 complex(8)::total(size(matrix,1),size(matrix,2))
 call comm_summation(matrix,total,size(matrix),info%icomm_rko)
 matrix=total
 end subroutine
 subroutine consensus(code)
 integer,intent(inout)::code
 call comm_get_max(code,info%icomm_rko)
 end subroutine
 subroutine project(a,b,dv,r,noccupied)
 complex(8),intent(in)::a(:,:,:),b(:,:,:)
 complex(8),intent(out)::r(:,:,:)
 real(8),intent(in)::dv
 integer,intent(in)::noccupied
 integer::kk,ii,code
#ifndef USE_SCALAPACK
 complex(8)::matrix(size(a,2),size(a,2))
 integer::i,j
#endif
 do kk=1,size(a,3)
  ii=kk+info%ik_s-1
#ifdef USE_SCALAPACK
  call ptcn_owner_project(owner,a(:,:,kk),b(:,:,kk),dv,noccupied,r(:,:,kk),code,system%rocc(:,ii,1))
#else
  matrix=dv*matmul(conjg(transpose(a(:,:,kk))),b(:,:,kk))
  do j=1,size(a,2);do i=1,size(a,2)
   if(system%rocc(i,ii,1)/=system%rocc(j,ii,1))matrix(i,j)=0
  enddo;enddo
  r(:,:,kk)=b(:,:,kk)-matmul(a(:,:,kk),matrix);code=0
#endif
  if(code/=0)error stop 'Native PT-CN occupation projection failed'
 enddo
 end subroutine
 subroutine action(a,b,refresh)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 logical,intent(in)::refresh
 if(first_action)then
  system%vec_Ac=ac0;first_action=.false.
 else
  system%vec_Ac=ac1
 endif
#ifdef USE_HSE
 ! Install at most once. Failed fresh checks retain their new exact ACE for correction.
 if(hybrid.and..not.refresh)call exx_ptcn_predict_inner(system,info)
#endif
 call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
 call unpack(a,trial)
 call calc_density(system,rho_s,trial,info,mg)
 rho%f=rho_s(1)%f
 if(system%nspin==2)rho%f=rho%f+rho_s(2)%f
 call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,vh)
 ! Local XC is rebuilt on every action. HSE fixes only exchange internally;
 ! every outer/final callback reconstructs the full physical exchange operator.
#ifdef USE_HSE
 if(hybrid)exx_freeze=.not.refresh
#endif
 call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,trial,stencil,vxc,energy%E_xc)
 call update_vlocal(mg,system%nspin,vh,vpsl,vxc,vlocal)
 call hpsi(trial,htrial,info,mg,vlocal,system,stencil,srg,ppg)
 call pack(htrial,b)
 end subroutine
 subroutine precondition(a,b)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 integer::kk,code
 real(8)::start
 start=omp_get_wtime()
 do kk=1,size(a,3)
  pre_in=reshape(a(:,:,kk),[ng,ns*no])
  call ptcn_fft_apply(fft,pre_in,den(:,kk),pre_out,code)
  if(code/=0)error stop 'Native PT-CN FFT precondition failed'
  b(:,:,kk)=reshape(pre_out,[ns*ng,no])
 enddo
 pre_seconds=pre_seconds+omp_get_wtime()-start
 end subroutine
 real(8) function global_norm(a)
 complex(8),intent(in)::a(:,:,:)
 real(8)::local,total
 local=sum(abs(a)**2)
 call comm_summation(local,total,info%icomm_rko)
 global_norm=sqrt(total)
 end function
 subroutine pack(psi,a)
 type(s_orbital),intent(in)::psi
 complex(8),intent(out)::a(:,:,:)
 integer::kk,jj,ss,xx,yy,zz,gg
 do kk=info%ik_s,info%ik_e;do jj=info%io_s,info%io_e;do ss=1,ns
 do zz=mg%is(3),mg%ie(3);do yy=mg%is(2),mg%ie(2);do xx=mg%is(1),mg%ie(1)
  gg=1+xx-mg%is(1)+mg%num(1)*(yy-mg%is(2)+mg%num(2)*(zz-mg%is(3)))+(ss-1)*ng
  a(gg,jj-info%io_s+1,kk-info%ik_s+1)=psi%zwf(xx,yy,zz,ss,jj,kk,1)
 enddo;enddo;enddo;enddo;enddo;enddo
 end subroutine
 subroutine unpack(a,psi)
 complex(8),intent(in)::a(:,:,:)
 type(s_orbital),intent(inout)::psi
 integer::kk,jj,ss,xx,yy,zz,gg
 do kk=info%ik_s,info%ik_e;do jj=info%io_s,info%io_e;do ss=1,ns
 do zz=mg%is(3),mg%ie(3);do yy=mg%is(2),mg%ie(2);do xx=mg%is(1),mg%ie(1)
  gg=1+xx-mg%is(1)+mg%num(1)*(yy-mg%is(2)+mg%num(2)*(zz-mg%is(3)))+(ss-1)*ng
  psi%zwf(xx,yy,zz,ss,jj,kk,1)=a(gg,jj-info%io_s+1,kk-info%ik_s+1)
 enddo;enddo;enddo;enddo;enddo;enddo
 psi%update_zwf_overlap=.false.
 end subroutine
 end subroutine
end module
