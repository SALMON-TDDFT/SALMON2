! Native RT adapter for core-partitioned LCFO screened exchange.
! Each rank owns one core. Hartree, XC, propagation and current stay native.
! Full density-factor reference or opt-in initial-MLWF U reuse and source masks.
module hse_lcfo_rt
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
  use communication, only: comm_bcast,comm_summation,comm_get_max
  use salmon_global, only: hse_omega,ae_shape1
  use lcfo_rt_basis
  use lcfo_rt_wannier, only: lcfo_mlwf_enabled,lcfo_mlwf_configure,lcfo_mlwf_source, &
    lcfo_mlwf_stage,lcfo_mlwf_accept_cached,lcfo_mlwf_track
  use hse_wannier, only: s_hse_wannier,wannier_init,wannier_apply,wannier_forward
  use lcfo_ace_local, only: lcfo_ace_local_action,lcfo_ace_half_trace
  use hse_ace, only: hse_ace_state,hse_ace_average
  use lcfo_dist_rows, only: s_lcfo_halo,lcfo_halo_init,lcfo_halo_get,lcfo_halo_sum
  use lcfo_projection, only: s_lcfo_projection,lcfo_projection_init,lcfo_projection_apply
  use lcfo_dist_dense, only: lcfo_distributed_ace_build
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_hse_refresh,lcfo_hse_add_action,lcfo_hse_stage
  complex(8),allocatable,save :: fragment_basis(:,:),core_basis(:,:),hx(:,:),initial_hx(:,:),midpoint_hx(:,:)
  integer,allocatable,save :: selected(:),fragment_global_index(:)
  real(8),allocatable,save :: core_weight(:)
  logical,save :: measure_continuity=.false.
  type(s_hse_wannier),save :: fragment_operator
  type(s_lcfo_halo),save :: exchange_plan
  type(s_lcfo_projection),save :: projection_plan
  type(hse_ace_state),save :: ace,initial_ace,midpoint_ace
  logical,save :: ace_valid=.false.,initial_ace_valid=.false.,midpoint_ace_valid=.false.,use_midpoint=.false.
  integer,save :: refresh_count=0,ace_interval=1,rt_step=0,step_start_builds=0,refresh_origin=0
  complex(8),allocatable,save :: cached_coeff(:,:)
  real(8),save :: cached_energy=0d0
  real(8),allocatable,save :: cached_occupation(:)
contains
  subroutine initialize_fragment()
    complex(8),allocatable :: block(:,:)
    integer,allocatable :: mapping(:,:),first(:)
    integer :: nf,ns(3),ng,f,g,x,y,z,p(3),global_point(3),rel(3),j,nsel,ierr,env_status,ios,bad,fft_batch
    character(16) :: value
    nf=size(lcfo_counts);ns=lcfo_core+2*lcfo_buffer;ng=product(ns)
    if(any(ns>lcfo_grid))error stop 'LCFO HSE: fragment exceeds global periodic grid'
    allocate(mapping(ng,nf),first(nf),fragment_global_index(ng),core_weight(ng));mapping=0;first=0
    core_weight=0d0
    bad=0
    if(lcfo_rank==0)then
      call get_environment_variable('SALMON_LCFO_RT_CONTINUITY',value,status=env_status)
      measure_continuity=env_status==0.and.trim(value)=='1'
      call get_environment_variable('SALMON_LCFO_RT_ACE_INTERVAL',value,status=env_status)
      if(env_status==0.and.len_trim(value)>0)then
        read(value,*,iostat=ios)ace_interval
        if(ios/=0)bad=1
      else if(env_status/=1.and.env_status/=0)then
        bad=1
      endif
      if(ace_interval<1)bad=1
      if(ae_shape1=='impulse')refresh_origin=1
    endif
    call comm_bcast(bad,lcfo_comm,0)
    if(bad/=0)error stop 'LCFO HSE: ACE interval must be a positive integer'
    call comm_bcast(ace_interval,lcfo_comm,0)
    call comm_bcast(refresh_origin,lcfo_comm,0)
    if(lcfo_rank==0)write(*,'(a,2i8)')'LCFO HSE ACE interval/refresh origin: ',ace_interval,refresh_origin
    call comm_bcast(measure_continuity,lcfo_comm,0)
    fft_batch=1;bad=0
    if(lcfo_rank==0)then
      call get_environment_variable('SALMON_LCFO_RT_FFT_BATCH',value,status=env_status)
      if(env_status==0.and.len_trim(value)>0)then
        read(value,*,iostat=ios)fft_batch
        if(ios/=0)bad=1
      else if(env_status/=1.and.env_status/=0)then
        bad=1
      endif
      if(fft_batch<1.or.fft_batch>32)bad=1
    endif
    call comm_bcast(bad,lcfo_comm,0)
    if(bad/=0)error stop 'LCFO HSE: FFT batch must be an integer from 1 to 32'
    call comm_bcast(fft_batch,lcfo_comm,0)
    ! FFT order is core/right-buffer then the periodic left buffer.
    g=0
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      g=g+1;p=[x,y,z]
      where(p>=lcfo_core+lcfo_buffer)p=p-ns
      global_point=modulo(lcfo_origins(:,lcfo_rank+1)+p,lcfo_grid)
      fragment_global_index(g)=1+global_point(1)+lcfo_grid(1)*(global_point(2)+lcfo_grid(2)*global_point(3))
      if(all([x,y,z]<lcfo_core))core_weight(g)=1d0
      do f=1,nf
        rel=modulo(global_point-lcfo_origins(:,f),lcfo_grid)
        if(all(rel<lcfo_core))mapping(g,f)=1+rel(1)+lcfo_core(1)*(rel(2)+lcfo_core(2)*rel(3))
      enddo
      if(count(mapping(g,:)>0)/=1)error stop 'LCFO HSE: core partition is not unique'
    enddo;enddo;enddo
    nsel=0
    do f=1,nf
      if(any(mapping(:,f)>0))then
        first(f)=nsel+1;nsel=nsel+lcfo_counts(f)
      endif
    enddo
    allocate(fragment_basis(ng,nsel),core_basis(ng,nsel),selected(nsel));fragment_basis=0d0
    do f=1,nf
      allocate(block(product(lcfo_core),lcfo_counts(f)));block=0d0
      if(f==lcfo_rank+1)block=lcfo_basis
      call comm_bcast(block,lcfo_comm,f-1)
      if(first(f)>0)then
        do j=1,lcfo_counts(f)
          selected(first(f)+j-1)=lcfo_offsets(f)+j
        enddo
        do g=1,ng
          if(mapping(g,f)>0)fragment_basis(g,first(f):first(f)+lcfo_counts(f)-1)=block(mapping(g,f),:)
        enddo
      endif
      deallocate(block)
    enddo
    core_basis=fragment_basis;g=0
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      g=g+1
      if(any([x,y,z]>=lcfo_core))core_basis(g,:)=0d0
    enddo;enddo;enddo
    call lcfo_projection_init(projection_plan,core_basis,.25d0*lcfo_dv)
    call wannier_init(fragment_operator,ns,[1,1,1],lcfo_h,reshape([0d0,0d0,0d0],[3,1]),hse_omega,ierr)
    if(ierr/=0)error stop 'LCFO HSE: fragment periodic exchange initialization failed'
    fragment_operator%fft_batch_size=fft_batch
    write(*,'(a,5i10)')'LCFO compact projection rank/rows/columns/fullrows/fullcolumns:', &
      lcfo_rank,size(projection_plan%rows),size(projection_plan%columns),ng,nsel
    call lcfo_mlwf_configure()
    if(lcfo_rank==0)write(*,'(a,3i6)')'LCFO HSE fragment grid:',ns
  end subroutine

  subroutine pack_coefficients(psi,system,mg,info,coeff)
    type(s_orbital),intent(in) :: psi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable,intent(out) :: coeff(:,:)
    complex(8),allocatable :: grid(:,:),local(:,:),full_local(:,:)
    integer :: io,is(3),ie(3),n
    if(.not.lcfo_rt_active)error stop 'LCFO HSE: inactive LCFO basis'
    if(system%nk/=1.or.system%nspin/=1.or.info%numo<1.or. &
       info%ik_s/=1.or.info%ik_e/=1.or.info%numm/=1)error stop 'LCFO HSE: requires Gamma/local orbital block'
    if(any(mg%num/=lcfo_core))error stop 'LCFO HSE: native grid does not match LCFO core'
    if(.not.allocated(psi%zwf))error stop 'LCFO HSE: complex wavefunctions required'
    is=mg%is;ie=mg%ie;n=size(lcfo_basis,2)
    allocate(grid(product(mg%num),info%numo),local(n,system%no),full_local(n,system%no));local=0d0
    do io=info%io_s,info%io_e
      grid(:,io-info%io_s+1)=reshape(psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[product(mg%num)])
    enddo
    local(:,info%io_s:info%io_e)=matmul(conjg(transpose(lcfo_basis)),grid)*lcfo_dv
    ! Refresh masters alone receive all orbital columns. Real-space WFs stay local.
    call comm_summation(local,full_local,size(local),lcfo_orb_comm,0)
    if(lcfo_orb_rank==0)then
      call move_alloc(full_local,coeff)
    else
      allocate(coeff(0,0))
    endif
  end subroutine

  subroutine lcfo_hse_refresh(system,mg,info,psi,exchange_energy)
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    real(8),intent(out) :: exchange_energy
    complex(8),allocatable :: coeff(:,:)
    integer :: old_count,nb,factor_shape(3),nsel
    real(8) :: started,pack_seconds
    started=wall_seconds();old_count=refresh_count;nb=sum(lcfo_counts)
    call pack_coefficients(psi,system,mg,info,coeff)
    pack_seconds=wall_seconds()-started
    if(.not.exchange_plan%ready)then
      if(lcfo_orb_rank==0)call initialize_fragment()
      nsel=0
      if(lcfo_orb_rank==0)nsel=size(selected)
      call comm_bcast(nsel,lcfo_orb_comm,0)
      if(lcfo_orb_rank/=0)allocate(selected(nsel))
      call comm_bcast(selected,lcfo_orb_comm,0)
      call lcfo_halo_init(exchange_plan,lcfo_counts,selected,lcfo_comm)
      if(lcfo_orb_rank==0)write(*,'(a,4i8)')'LCFO distributed storage rank/local/global/halo rows:', &
        lcfo_rank,size(coeff,1),nb,nsel
    endif
    if(lcfo_orb_rank==0)call refresh_master(system,mg,info,coeff,exchange_energy,pack_seconds)
    call comm_bcast(exchange_energy,lcfo_orb_comm,0)
    call comm_bcast(refresh_count,lcfo_orb_comm,0)
    call comm_bcast(refresh_origin,lcfo_orb_comm,0)
    call comm_bcast(ace_interval,lcfo_orb_comm,0)
    call comm_bcast(ace_valid,lcfo_orb_comm,0)
    if(refresh_count/=old_count)then
      if(.not.allocated(hx))allocate(hx(size(selected),size(selected)))
      call comm_bcast(hx,lcfo_orb_comm,0)
      if(ace_valid)then
        if(lcfo_orb_rank==0)factor_shape=shape(ace%factors)
        call comm_bcast(factor_shape,lcfo_orb_comm,0)
        if(lcfo_orb_rank/=0)then
          if(allocated(ace%factors))deallocate(ace%factors)
          allocate(ace%factors(factor_shape(1),factor_shape(2),factor_shape(3)))
        endif
        call comm_bcast(ace%factors,lcfo_orb_comm,0)
        call comm_bcast(ace%dv,lcfo_orb_comm,0)
        call comm_bcast(ace%condition,lcfo_orb_comm,0)
      else
        if(allocated(ace%factors))deallocate(ace%factors)
      endif
    endif
  end subroutine

  subroutine refresh_master(system,mg,info,coeff,exchange_energy,pack_seconds)
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),intent(in) :: coeff(:,:)
    real(8),intent(out) :: exchange_energy
    complex(8),allocatable :: weighted(:,:),density(:,:),work(:),source(:,:,:),action(:,:,:)
    complex(8),allocatable :: projected(:,:),near_coeff(:,:),w(:,:),contribution(:,:)
    real(8),allocatable :: eigenvalues(:),rwork(:)
    real(8),intent(in) :: pack_seconds
    real(8) :: threshold,discarded,local_energy,started,source_seconds,exchange_seconds
    real(8) :: action_seconds,projection_seconds,projection_started
    integer :: nsel,ng,no,j,ierr,nrank,first,changed,total_changed
    external :: zheev
    started=wall_seconds()
    no=size(coeff,2);nsel=size(selected);ng=size(fragment_basis,1)
    if(any(system%rocc(:,1,1)<0d0).or.any(system%rocc(:,1,1)>2d0).or. &
       .not.all(ieee_is_finite(system%rocc(:,1,1))))error stop 'LCFO HSE: invalid occupations'
    ! On impulse step1 rebuild both endpoints. Smooth fields can retain the
    ! zero-field initial operator; subsequent refreshes are physical-step based.
    if(ace_interval>1.and.rt_step>0.and.ace_valid.and.allocated(cached_occupation))then
      if(mod(rt_step-refresh_origin,ace_interval)/=0.and.all(system%rocc(:,1,1)==cached_occupation))then
        if(lcfo_mlwf_enabled)call lcfo_mlwf_track(coeff)
        call lcfo_ace_half_trace(coeff,ace%factors(:,:,1),system%rocc(:,1,1), &
          ace%dv,lcfo_comm,exchange_energy)
        ! A transported current_frame no longer belongs to the exact cache key.
        if(allocated(cached_coeff))deallocate(cached_coeff)
        if(lcfo_rank==0)write(*,'(a,i8,a,es20.10)')'LCFO HSE ACE retained at step ',rt_step, &
          ' frozen-operator trace energy ',exchange_energy
        call report_timings('retained',pack_seconds,wall_seconds()-started,0d0,0d0)
        return
      endif
    endif
    changed=1
    if(allocated(cached_coeff).and..not.(rt_step==1.and.refresh_origin==1))then
      if(all(coeff==cached_coeff).and.all(system%rocc(:,1,1)==cached_occupation))changed=0
    endif
    ! A local cache hit is insufficient: every core must agree before returning.
    call comm_summation(changed,total_changed,lcfo_comm)
    if(total_changed==0)then
      exchange_energy=cached_energy
      call lcfo_mlwf_accept_cached()
      return
    endif
    call lcfo_halo_get(exchange_plan,coeff,near_coeff)
    if(lcfo_mlwf_enabled)then
      if(any(system%rocc(:,1,1)/=2d0))error stop 'LCFO MLWF reuse: fixed fully occupied states required'
      call lcfo_mlwf_source(coeff,fragment_basis,selected,source,exchange_plan)
      nrank=size(source,2);threshold=0d0;discarded=0d0
    else
    allocate(weighted(nsel,no),density(nsel,nsel),eigenvalues(nsel),work(max(1,2*nsel)),rwork(max(1,3*nsel-2)))
    do j=1,no
      weighted(:,j)=near_coeff(:,j)*sqrt(system%rocc(j,1,1)/2d0)
    enddo
    density=matmul(weighted,transpose(conjg(weighted)))
    density=.5d0*(density+transpose(conjg(density)))
    call zheev('V','U',nsel,density,nsel,eigenvalues,work,size(work),rwork,ierr)
    if(ierr/=0)error stop 'LCFO HSE: fragment density diagonalization failed'
    threshold=1d-14*max(0d0,eigenvalues(nsel))
    if(eigenvalues(1)<-max(threshold,1d-12*maxval(abs(eigenvalues)))) &
      error stop 'LCFO HSE: fragment density is not positive semidefinite'
    nrank=count(eigenvalues>threshold);first=nsel-nrank+1
    discarded=sum(max(0d0,eigenvalues(:first-1)))
    allocate(source(ng,nrank,1))
    do j=1,nrank
      source(:,j,1)=matmul(fragment_basis,density(:,first+j-1))*sqrt(eigenvalues(first+j-1))
    enddo
    endif
    source_seconds=wall_seconds()-started;started=wall_seconds()
    if(allocated(fragment_operator%source))deallocate(fragment_operator%source)
    allocate(fragment_operator%source(ng,nrank))
    call wannier_forward(fragment_operator,source,fragment_operator%source)
    allocate(action(ng,nsel,1))
    call wannier_apply(fragment_operator,reshape(fragment_basis,[ng,nsel,1]),action,ierr)
    if(ierr/=0)error stop 'LCFO HSE: fragment exchange action failed'
    action_seconds=wall_seconds()-started
    write(*,'(a,i6,2i14)')'LCFO exchange FFT pairs rank/executed/possible:', &
      lcfo_rank,fragment_operator%fft_pairs_executed,fragment_operator%fft_pairs_total
    write(*,'(a,i6,i6,i14)')'LCFO exchange FFT batches rank/width/calls:', &
      lcfo_rank,fragment_operator%worker_batch,fragment_operator%fft_batches_executed
    if(measure_continuity)call exchange_continuity(near_coeff,system%rocc(:,1,1),action)
    projection_started=wall_seconds()
    call lcfo_projection_apply(projection_plan,action(:,:,1),projected)
    projection_seconds=wall_seconds()-projection_started
    ! Retain this fragment's Hermitian contribution, not a replicated global Hx.
    hx=projected
    if(.not.all(ieee_is_finite(real(hx))).or..not.all(ieee_is_finite(aimag(hx)))) &
      error stop 'LCFO HSE: nonfinite projected exchange'
    contribution=matmul(hx,near_coeff)
    call lcfo_halo_sum(exchange_plan,contribution,w)
    local_energy=0d0
    do j=1,no
      local_energy=local_energy+.5d0*system%rocc(j,1,1)*real(sum(conjg(coeff(:,j))*w(:,j)),8)
    enddo
    call comm_summation(local_energy,exchange_energy,lcfo_comm)
    exchange_seconds=wall_seconds()-started;started=wall_seconds()
    call lcfo_distributed_ace_build(ace,coeff,w,1d0,lcfo_comm,ierr);ace_valid=ierr==0
    call report_timings('build',pack_seconds,source_seconds,exchange_seconds,wall_seconds()-started, &
      action_seconds,projection_seconds)
    cached_coeff=coeff;cached_occupation=system%rocc(:,1,1);cached_energy=exchange_energy
    refresh_count=refresh_count+1
    if(lcfo_rank==0)write(*,'(a,i8)')'LCFO HSE exchange rebuilt at step ',rt_step
    ! Every rank reports its actual source rank and discarded density weight.
    write(*,'(a,i6,a,i6,a,i6,a,es12.4,a,es12.4)')'LCFO HSE rank ',lcfo_rank,' refresh ',refresh_count, &
      ' density factors ',nrank,' eig threshold ',threshold,' discarded trace ',discarded
    if(lcfo_rank==0)then
      if(ace_valid)then
        write(*,'(a,i8,a,es12.4)')'LCFO HSE ACE build ',refresh_count,' valid, condition ',ace%condition
      else
        write(*,'(a,i8,a)')'LCFO HSE ACE build ',refresh_count, &
          ' rejected (indefinite/singular metric); using full projected exchange, no eigenvalue clipping'
      endif
    endif
  end subroutine

  subroutine exchange_continuity(coeff,occupation,action)
    ! Density source 2 Im sum_i f_i psi_i^* (K psi_i), evaluated BEFORE
    ! LCFO output projection. For full Fock it cancels pointwise; masks may
    ! break cancellation even though the assembled K remains Hermitian.
    ! Include BOTH core-weighted adjoint halves, then sum overlapping buffers.
    complex(8),intent(in) :: coeff(:,:),action(:,:,:)
    real(8),intent(in) :: occupation(:)
    complex(8),allocatable :: core_action(:,:,:),psi(:,:),left(:,:),right(:,:)
    real(8),allocatable :: local(:),global(:)
    integer :: ng,ns,g,j,ierr
    ng=size(fragment_basis,1);ns=size(fragment_basis,2)
    allocate(core_action(ng,ns,1))
    call wannier_apply(fragment_operator,reshape(core_basis,[ng,ns,1]),core_action,ierr)
    if(ierr/=0)error stop 'LCFO HSE: continuity diagnostic action failed'
    psi=matmul(fragment_basis,coeff)
    left=matmul(action(:,:,1),coeff)
    right=matmul(core_action(:,:,1),coeff)
    allocate(local(product(lcfo_grid)),global(product(lcfo_grid)));local=0d0
    do j=1,size(coeff,2);do g=1,ng
      local(fragment_global_index(g))=local(fragment_global_index(g))+ &
        .25d0*occupation(j)*aimag(conjg(psi(g,j))*(core_weight(g)*left(g,j)+right(g,j)))
    enddo;enddo
    call comm_summation(local,global,size(local),lcfo_comm)
    if(lcfo_rank==0)write(*,'(a,i8,3es19.10)')'LCFO EXX continuity build/signed/L1/max: ', &
      refresh_count+1,sum(global)*lcfo_dv,sum(abs(global))*lcfo_dv,maxval(abs(global))
  end subroutine

  subroutine lcfo_hse_stage(stage,system,mg,info,psi)
    integer,intent(in) :: stage
    type(s_dft_system),intent(in),optional :: system
    type(s_rgrid),intent(in),optional :: mg
    type(s_parallel_info),intent(in),optional :: info
    type(s_orbital),intent(in),optional :: psi
    real(8) :: rebuilt_energy
    integer :: ierr
    if(stage==0)then
      rt_step=rt_step+1
      if(rt_step==1.and.refresh_origin==1)then
        if(.not.(present(system).and.present(mg).and.present(info).and.present(psi))) &
          error stop 'LCFO HSE: impulse initial refresh requires current orbitals'
        call lcfo_hse_refresh(system,mg,info,psi,rebuilt_energy)
        if(lcfo_rank==0.and.lcfo_orb_rank==0)write(*,'(a)')'LCFO HSE impulse ACE rebuilt before first predictor'
      endif
    endif
    if(lcfo_orb_rank==0)call lcfo_mlwf_stage(stage)
    select case(stage)
    case(0)
      step_start_builds=refresh_count
      if(.not.allocated(hx))error stop 'LCFO HSE: initial operator missing'
      initial_hx=hx;initial_ace=ace;initial_ace_valid=ace_valid;use_midpoint=.false.
    case(1)
      if(.not.allocated(initial_hx).or..not.allocated(hx))error stop 'LCFO HSE: midpoint endpoints missing'
      midpoint_hx=.5d0*(initial_hx+hx)
      midpoint_ace_valid=initial_ace_valid.and.ace_valid
      if(midpoint_ace_valid)then
        if(step_start_builds==refresh_count)then
          midpoint_ace=ace
        else
          call hse_ace_average(initial_ace,ace,midpoint_ace,ierr)
          midpoint_ace_valid=ierr==0
        endif
      endif
      use_midpoint=.true.
    case(2)
      use_midpoint=.false.;initial_ace_valid=.false.;midpoint_ace_valid=.false.
      if(allocated(initial_hx))deallocate(initial_hx)
      if(allocated(midpoint_hx))deallocate(midpoint_hx)
      if(allocated(initial_ace%factors))deallocate(initial_ace%factors)
      if(allocated(midpoint_ace%factors))deallocate(midpoint_ace%factors)
    case default
      error stop 'LCFO HSE: unknown Taylor stage'
    end select
  end subroutine

  subroutine lcfo_hse_add_action(psi,hpsi,system,mg,info)
    type(s_orbital),intent(in) :: psi
    type(s_orbital),intent(inout) :: hpsi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable :: coeff(:,:),grid(:,:),hgrid(:,:),local_action(:,:),near_coeff(:,:),contribution(:,:)
    integer :: ng,no,io,j,is(3),ie(3)
    if(.not.allocated(hx))error stop 'LCFO HSE: refresh required before action'
    no=info%numo;ng=product(mg%num);is=mg%is;ie=mg%ie
    allocate(grid(ng,no),hgrid(ng,no))
    do j=1,no
      io=info%io_s+j-1
      grid(:,j)=reshape(psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[ng])
      hgrid(:,j)=reshape(hpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[ng])
    enddo
    if(use_midpoint.and.midpoint_ace_valid)then
      call lcfo_ace_local_action(lcfo_basis,grid,hgrid,midpoint_ace%factors(:,:,1), &
        lcfo_dv,midpoint_ace%dv,lcfo_comm)
    else if(.not.use_midpoint.and.ace_valid)then
      call lcfo_ace_local_action(lcfo_basis,grid,hgrid,ace%factors(:,:,1),lcfo_dv,ace%dv,lcfo_comm)
    else
      coeff=matmul(conjg(transpose(lcfo_basis)),grid)*lcfo_dv
      call lcfo_halo_get(exchange_plan,coeff,near_coeff)
      if(use_midpoint)then
        contribution=matmul(midpoint_hx,near_coeff)
      else
        contribution=matmul(hx,near_coeff)
      endif
      call lcfo_halo_sum(exchange_plan,contribution,local_action)
      hgrid=matmul(lcfo_basis,matmul(conjg(transpose(lcfo_basis)),hgrid)*lcfo_dv+local_action)
    endif
    do j=1,no
      io=info%io_s+j-1
      hpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)=reshape(hgrid(:,j),mg%num)
    enddo
    hpsi%update_zwf_overlap=.false.
  end subroutine
  real(8) function wall_seconds() result(seconds)
    integer(8) :: tick,rate
    call system_clock(tick,rate)
    seconds=dble(tick)/dble(rate)
  end function

  subroutine report_timings(label,pack,source,exchange,build,action,projection)
    character(*),intent(in) :: label
    real(8),intent(in) :: pack,source,exchange,build
    real(8),intent(in),optional :: action,projection
    real(8) :: local(6),maximum(6)
    local=0d0;local(1:4)=[pack,source,exchange,build]
    if(present(action))local(5)=action
    if(present(projection))local(6)=projection
    call comm_get_max(local,maximum,6,lcfo_comm)
    if(lcfo_rank==0)then
      write(*,'(a,a,4es16.7)')'LCFO timing pack/source/exchange/ACE ',trim(label),maximum(1:4)
      if(present(action))write(*,'(a,2es16.7)')'LCFO exchange detail action/projection ',maximum(5:6)
    endif
  end subroutine
end module
