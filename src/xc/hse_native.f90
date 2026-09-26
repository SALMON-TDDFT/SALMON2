#include "config.h"
! SALMON adapter: initial certified layout is complete grid/orbitals per rank,
! distributed k points. The legacy backend transfers density tiles; the Wannier
! baseline gathers a fragment on its k root. ACE applications stay local.
module hse_native
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use structures
  use plusU_global, only: PLUS_U_ON
  use hse_exchange
  use hse_ace
  use hse_wannier
  use hse_symmetry
  use sym_sub, only: use_symmetry,SymMatA,SymMatB
  use communication, only: comm_summation,comm_alltoall
  use salmon_global, only: xc,yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_md,yn_symmetrized_stencil,propagator,num_kgrid,hse_omega, &
    yn_hse_wannier,hse_mlwf_interval,hse_mlwf_maxiter,hse_mlwf_tolerance
  implicit none
  private
  public :: hse_export_snapshot,hse_eigen_diagnostic_enabled,hse_export_eigen_pair
  public :: hse_enabled,hse_refresh,hse_add_action,hse_exchange_energy,hse_freeze
  public :: hse_pack,hse_unpack,hse_timings,hse_walltime
  public :: hse_taylor_stage,hse_core_exchange,hse_force_full_action
  type(hse_symmetry_map),save :: symmetry_map
  type(hse_kernel),save :: kernel
  type(s_hse_wannier),save :: wannier
  real(8),allocatable,save :: cached_occupation(:,:)
  complex(8),allocatable,save :: cached_action(:,:,:)
  logical,save :: hse_force_full_action=.false.
  type(hse_ace_state),save :: ace
  type(hse_ace_state),save :: initial_ace,midpoint_ace
  complex(8),allocatable,save :: full_source(:,:,:),initial_source(:,:,:),midpoint_source(:,:,:)
  logical,save :: taylor_active=.false.,taylor_midpoint=.false.
  complex(8),allocatable,save :: cached_source(:,:,:),target_work(:,:,:),action_work(:,:,:),output_work(:,:,:)
  real(8),save :: hse_exchange_energy=0d0
  real(8),save :: hse_timings(4)=0d0 ! full EXX, ACE build, ACE apply, EXX collectives
  logical,save :: hse_freeze=.false.,reported_team=.false.,timing_enabled=.false.
contains
  logical function hse_eigen_diagnostic_enabled(info,variable) result(enabled)
    use communication, only: comm_bcast
    type(s_parallel_info),intent(in) :: info
    character(*),optional,intent(in) :: variable
    integer :: flag,status
    character(8) :: setting
    enabled=.false.
    if(.not.hse_enabled())return
    flag=0
    if(info%id_rko==0)then
      if(present(variable))then
        call get_environment_variable(variable,setting,status=status)
      else
        call get_environment_variable('SALMON_HSE_EIGEN_DIAGNOSTIC',setting,status=status)
      endif
      if(status==0.and.trim(setting)=='1')flag=1
    endif
    call comm_bcast(flag,info%icomm_rko,0)
    enabled=flag==1
  end function

  subroutine hse_export_eigen_pair(system,mg,info,psi,hpsi,tag)
    use iso_fortran_env, only: int32
    use salmon_global, only: base_directory
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi,hpsi
    character(*),optional,intent(in) :: tag
    character(:),allocatable :: filename
    complex(8),allocatable :: p(:,:,:),hp(:,:,:)
    integer :: iu,status,closed,total,ng,no
    ! Diagnostic format deliberately supports only a full Gamma fragment.
    if(system%nk/=1.or.info%isize_r/=1.or.info%isize_o/=1) &
      error stop 'HSE eigen diagnostic requires Gamma and full grid/orbitals'
    ng=product(mg%num);no=system%no
    allocate(p(ng,no,1),hp(ng,no,1))
    call hse_pack(psi,mg,info,p);call hse_pack(hpsi,mg,info,hp)
    filename=trim(base_directory)//'hse_eigen_pair.bin'
    if(present(tag))filename=trim(base_directory)//'hse_eigen_'//tag//'.bin'
    open(newunit=iu,file=filename, &
      status='replace',access='stream',form='unformatted',iostat=status)
    if(status==0)then
      write(iu,iostat=status)int([16909060,1,ng,no,1],int32)
      if(status==0)write(iu,iostat=status)system%hvol,system%rocc(:,:,1),p,hp
      close(iu,iostat=closed)
      if(status==0)status=closed
    endif
    call comm_summation(abs(status),total,info%icomm_rko)
    if(total/=0)error stop 'HSE eigen diagnostic write failed'
  end subroutine

  subroutine hse_export_snapshot(system,mg,info,psi,iteration,residual,converged)
    use salmon_global, only: base_directory
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    integer,intent(in) :: iteration
    real(8),intent(in) :: residual
    logical,intent(in) :: converged
    character(8) :: setting
    integer :: status,enabled
    if(.not.hse_enabled().or..not.use_wannier_exchange())return
    enabled=0
    if(info%id_k==0)then
      call get_environment_variable('SALMON_HSE_WANNIER_SNAPSHOT',setting,status=status)
      if(status==0.and.trim(setting)=='1')enabled=1
    endif
    call comm_bcast(enabled,info%icomm_k,0)
    if(enabled==0)return
    call hse_refresh(system,mg,info,psi)
    status=0
    if(info%id_k==0)call wannier_snapshot(wannier,wannier%source_occupation,hse_omega,hse_exchange_energy, &
      residual,iteration,converged,trim(base_directory)//'hse_wannier_snapshot.bin',status)
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier snapshot: write failed'
  end subroutine
  subroutine hse_taylor_stage(stage)
    integer,intent(in) :: stage
    integer :: ierr,no
    select case(stage)
    case(0)
      taylor_active=.true.;taylor_midpoint=.false.
      if(propagator=='hse_taylor4_full')then
        initial_source=full_source
      else
        initial_ace=ace
      endif
    case(1)
      if(propagator=='hse_taylor4_full')then
        no=size(full_source,2)
        if(allocated(midpoint_source))deallocate(midpoint_source)
        allocate(midpoint_source(size(full_source,1),2*no,size(full_source,3)))
        midpoint_source(:,:no,:)=initial_source/sqrt(2d0)
        midpoint_source(:,no+1:,:)=full_source/sqrt(2d0)
      else
        call hse_ace_average(initial_ace,ace,midpoint_ace,ierr)
        if(ierr/=0)error stop 'HSE Taylor midpoint ACE failed'
      endif
      taylor_midpoint=.true.
    case(2)
      taylor_active=.false.;taylor_midpoint=.false.
      if(allocated(initial_ace%factors))deallocate(initial_ace%factors)
      if(allocated(midpoint_ace%factors))deallocate(midpoint_ace%factors)
      if(allocated(initial_source))deallocate(initial_source)
      if(allocated(midpoint_source))deallocate(midpoint_source)
    case default
      error stop 'Invalid HSE Taylor stage'
    end select
  end subroutine

  real(8) function hse_walltime()
    integer(int64) :: count,rate
    call system_clock(count,rate)
    hse_walltime=real(count,8)/real(rate,8)
  end function

  logical function hse_enabled()
    hse_enabled=trim(xc)=='hse06'
  end function

  subroutine hse_pack(psi,mg,info,a)
    type(s_orbital),intent(in) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),intent(out) :: a(:,:,:)
    integer :: ik,io,is(3),ie(3)
    is=mg%is;ie=mg%ie
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e
      a(:,io-info%io_s+1,ik-info%ik_s+1)=reshape( &
        psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,ik,1),[product(mg%num)])
    enddo;enddo
  end subroutine

  subroutine hse_unpack(a,psi,mg,info)
    complex(8),intent(in) :: a(:,:,:)
    type(s_orbital),intent(inout) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ik,io,is(3),ie(3)
    is=mg%is;ie=mg%ie
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e
      psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,ik,1)= &
        reshape(a(:,io-info%io_s+1,ik-info%ik_s+1),mg%num)
    enddo;enddo
    psi%update_zwf_overlap=.false.
  end subroutine

  subroutine hse_refresh(system,mg,info,psi)
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: w(:,:,:),local(:,:,:)
    real(8) :: ex,offdiag(3,3),tick,communication_before
    integer :: ierr,total_error,ng,nk,no,n,mesh,j,first_full,count_full
    if(.not.hse_enabled().or.hse_freeze)return
    if(info%isize_r/=1.or.info%isize_o/=1.or.info%numm/=1) &
      error stop 'HSE06: initial native support requires k-only MPI distribution'
    if(yn_periodic/='y'.or.system%nspin/=1.or..not.allocated(psi%zwf)) &
      error stop 'HSE06: periodic complex unpolarized orbitals required'
    if(yn_spinorbit/='n'.or.yn_jm/='n'.or.yn_md/='n'.or.yn_symmetrized_stencil=='y') &
      error stop 'HSE06: unsupported Hamiltonian/ionic extension'
    if(PLUS_U_ON)error stop 'HSE06: DFT+U combination unsupported'
    if(allocated(system%Ac_micro%v))error stop 'HSE06: microscopic vector potential unsupported'
    if(use_wannier_exchange())then
      call refresh_wannier(system,mg,info,psi)
      return
    endif
    ng=product(mg%num);nk=system%nk;no=system%no;n=mg%num(1);mesh=nint(real(nk,8)**(1d0/3d0))
    if(use_symmetry)mesh=num_kgrid(1)
    if(any(mg%num/=n).or.(.not.use_symmetry.and.mesh**3/=nk).or. &
      maxval(abs(system%hgs-system%hgs(1)))>1d-12) &
      error stop 'HSE06: cubic grid and full cubic k mesh required'
    offdiag=system%primitive_a
    do j=1,3;offdiag(j,j)=0;enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE06: orthogonal cell required'
    if((.not.use_symmetry.and.maxval(abs(system%wtk-1d0/nk))>1d-12).or. &
      maxval(abs(system%rocc-2d0))>1d-12) &
      error stop 'HSE06: uniform k weights and fully occupied spin pairs required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE06: unsupported orbital layout'
    if(kernel%n==0)then
      ! Keep representatives persistent; expand only rank-local stars for EXX.
      if(use_symmetry)then
        if(any(num_kgrid/=mesh))error stop 'HSE symmetry: cubic full mesh required'
        call symmetry_init(symmetry_map,mg%num,system%hgs,system%vec_k(:,:nk),system%wtk(:nk), &
          SymMatA,SymMatB,mesh,ierr)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: invalid grid, stars or weights'
        call symmetry_validate_atoms(symmetry_map,system%Rion,system%kion,ierr)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: operation does not preserve atoms'
        first_full=symmetry_map%first(info%ik_s)
        count_full=symmetry_map%first(info%ik_e+1)-first_full
        call hse_kernel_init(kernel,n,mesh,system%hgs(1),symmetry_map%full_k,hse_omega, &
          max(1,min(16,64/info%isize_k)),ierr,first_full,count_full)
      else
        call hse_kernel_init(kernel,n,mesh,system%hgs(1),system%vec_k,hse_omega, &
          max(1,min(16,64/info%isize_k)),ierr,info%ik_s,info%numk)
      endif
      call comm_summation(ierr,total_error,info%icomm_rko)
      if(total_error/=0)error stop 'HSE06: kernel initialization failed'
      timing_enabled=kernel%profile.or.propagator=='hse_ptcn'
    endif
    allocate(local(ng,no,info%numk))
    call hse_pack(psi,mg,info,local)
    ierr=1
    if(allocated(cached_source))then
      if(all(shape(cached_source)==shape(local)))then
        if(all(cached_source==local))ierr=0
      endif
    endif
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error==0)return
    allocate(w(ng,no,info%numk))
    if(propagator=='hse_taylor4_full')full_source=local
    if(timing_enabled)then
      tick=hse_walltime();communication_before=hse_timings(4)
    endif
    call apply_distributed(local,local,w,info,ierr)
    if(timing_enabled)hse_timings(1)=hse_timings(1)+hse_walltime()-tick-(hse_timings(4)-communication_before)
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: distributed exchange action failed'
    if(kernel%profile.and.info%id_k==0)write(*,'(a,i0,a,6es14.5)') &
      'HSE_PROFILE rank0 block=',kernel%block,' density comm fft kernel packing action=',kernel%seconds
    if(.not.reported_team.and.info%id_k==0)write(*,'(a,i0)')'HSE_OPENMP threads=',kernel%threads_used
    if(.not.reported_team.and.info%id_k==0)write(*,'(a,l1)')'HSE_FFT_CONTIGUOUS=',kernel%contiguous_fft
    if(.not.reported_team.and.info%id_k==0.and.kernel%auto_fft) &
      write(*,'(a,2es14.5)')'HSE_FFT_AUTO trial strided contiguous seconds=',kernel%fft_trial_seconds
    reported_team=.true.
    if(timing_enabled)tick=hse_walltime()
    call hse_ace_build(ace,local,w,system%hvol,ierr)
    if(timing_enabled)hse_timings(2)=hse_timings(2)+hse_walltime()-tick
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: ACE metric failed'
    cached_source=local
    ex=0d0
    do j=1,info%numk
      ex=ex+.25d0*real(sum(conjg(local(:,:,j))*w(:,:,j)),8)*system%hvol*system%wtk(info%ik_s+j-1)
    enddo
    call comm_summation(ex,hse_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_distributed(source,target,action,info,ierr)
    complex(8),intent(in) :: source(:,:,:),target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    integer,intent(out) :: ierr
    integer :: local_layout(2*info%isize_k),layout(2*info%isize_k),np,j,index,total_error
    complex(8),allocatable :: expanded_target(:,:,:),expanded_action(:,:,:),rotated(:,:)
    np=info%isize_k;local_layout=0
    local_layout(info%id_k+1)=info%ik_s
    local_layout(np+info%id_k+1)=info%numk
    if(use_symmetry)then
      local_layout(info%id_k+1)=symmetry_map%first(info%ik_s)
      local_layout(np+info%id_k+1)=symmetry_map%first(info%ik_e+1)-symmetry_map%first(info%ik_s)
    endif
    call comm_summation(local_layout,layout,size(layout),info%icomm_k)
    if(use_symmetry)then
      ierr=0
      if(.not.all(ieee_is_finite(real(source))).or..not.all(ieee_is_finite(aimag(source))))ierr=1
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)then
        ierr=1;return
      endif
      allocate(expanded_target(size(target,1),size(target,2),layout(np+info%id_k+1)))
      allocate(rotated(size(source,1),size(source,2)))
      do j=1,size(expanded_target,3)
        index=layout(info%id_k+1)+j-1
        call symmetry_transform(symmetry_map,index,1,target(:,:,symmetry_map%owner(index)-info%ik_s+1), &
          expanded_target(:,:,j),ierr)
        if(ierr/=0)error stop 'HSE symmetry target transformation failed'
      enddo
      allocate(expanded_action(size(expanded_target,1),size(expanded_target,2),size(expanded_target,3)))
      ! Source argument is only a valid layout placeholder when the density callback is supplied.
      call hse_kernel_apply_distributed(kernel,expanded_target,expanded_target,expanded_action,layout(:np), &
        layout(np+1:),info%id_k,transpose_tiles,ierr,fill_density)
      if(ierr==0)then
        do j=1,info%numk
          index=symmetry_map%first(info%ik_s+j-1)-symmetry_map%first(info%ik_s)+1
          action(:,:,j)=expanded_action(:,:,index)
        enddo
      endif
      return
    endif
    call hse_kernel_apply_distributed(kernel,source,target,action,layout(:np),layout(np+1:), &
      info%id_k,transpose_tiles,ierr)
  contains
    subroutine fill_density(j,lo,rows,density)
      integer,intent(in) :: j,lo,rows
      complex(8),intent(out) :: density(:,:)
      integer :: full_index,rep,op,g,stat,ng,no
      complex(8) :: coefficient
      external :: zgemm
      ng=size(source,1);no=size(source,2)
      full_index=layout(info%id_k+1)+j-1
      rep=symmetry_map%owner(full_index)-info%ik_s+1
      coefficient=cmplx(1d0/symmetry_map%multiplicity(full_index),0d0,8)
      density=0d0
      ! Average little-group projectors one orbital block at a time, not copies of all orbitals.
      do op=1,symmetry_map%multiplicity(full_index)
        call symmetry_transform(symmetry_map,full_index,op,source(:,:,rep),rotated,stat)
        if(stat/=0)error stop 'HSE symmetry source transformation failed'
        do g=1,no
          rotated(:,g)=rotated(:,g)*kernel%phase(:,j)
        enddo
        call zgemm('N','C',rows,ng,no,coefficient,rotated(lo,1),ng,rotated(1,1),ng, &
          (1d0,0d0),density(1,1),size(density,1))
      enddo
    end subroutine
    subroutine transpose_tiles(send,recv,count)
      complex(8),intent(in) :: send(:)
      complex(8),intent(out) :: recv(:)
      integer,intent(in) :: count
      real(8) :: start
      if(timing_enabled)start=hse_walltime()
      call comm_alltoall(send,recv,info%icomm_k,count)
      if(timing_enabled)hse_timings(4)=hse_timings(4)+hse_walltime()-start
    end subroutine
  end subroutine

  subroutine hse_add_action(psi,hpsi,system,mg,info)
    type(s_orbital),intent(in) :: psi
    type(s_orbital),intent(inout) :: hpsi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ierr,ng,total_error
    real(8) :: tick,communication_before
    if(.not.hse_enabled())return
    if(.not.allocated(ace%factors))error stop 'HSE06: occupied exchange source is not initialized'
    ng=product(mg%num)
    if(allocated(target_work))then
      if(any(shape(target_work)/=[ng,info%numo,info%numk]))deallocate(target_work,action_work,output_work)
    endif
    if(.not.allocated(target_work))allocate(target_work(ng,info%numo,info%numk), &
      action_work(ng,info%numo,info%numk),output_work(ng,info%numo,info%numk))
    call hse_pack(psi,mg,info,target_work)
    if(timing_enabled)tick=hse_walltime()
    if(use_wannier_exchange().and.hse_force_full_action)then
      call apply_wannier_collective(target_work,action_work,info)
      ierr=0
    else if(taylor_active.and.propagator=='hse_taylor4_full')then
      communication_before=hse_timings(4)
      if(taylor_midpoint)then
        call apply_distributed(midpoint_source,target_work,action_work,info,ierr)
      else
        call apply_distributed(initial_source,target_work,action_work,info,ierr)
      endif
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'HSE Taylor full target action failed'
      if(timing_enabled)hse_timings(1)=hse_timings(1)+hse_walltime()-tick-(hse_timings(4)-communication_before)
    else
      if(taylor_active.and.taylor_midpoint)then
        call hse_ace_apply(midpoint_ace,target_work,action_work,ierr)
      else
        call hse_ace_apply(ace,target_work,action_work,ierr)
      endif
      if(timing_enabled)hse_timings(3)=hse_timings(3)+hse_walltime()-tick
    endif
    if(ierr/=0)error stop 'HSE06: ACE application failed'
    call hse_pack(hpsi,mg,info,output_work)
    output_work=output_work+.25d0*action_work
    call hse_unpack(output_work,hpsi,mg,info)
  end subroutine
  logical function use_wannier_exchange()
    implicit none
    use_wannier_exchange=yn_dc=='y'.or.yn_hse_wannier=='y'
  end function

  subroutine refresh_wannier(system,mg,info,psi)
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: local(:,:,:),send(:,:,:),allpsi(:,:,:),allw(:,:,:),w(:,:,:)
    real(8) :: ex,offdiag(3,3)
    integer :: ng,no,nk,ik,j,status,changed,total_changed,maxiter
    ng=product(mg%num);no=system%no;nk=system%nk
    if(use_symmetry.or.nk/=product(num_kgrid))error stop 'HSE Wannier: full k mesh required'
    if(maxval(abs(system%wtk-1d0/nk))>1d-12)error stop 'HSE Wannier: uniform k weights required'
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE Wannier: orthogonal cell required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE Wannier: invalid orbital layout'
    allocate(local(ng,no,info%numk))
    call hse_pack(psi,mg,info,local)
    changed=1
    if(allocated(cached_source).and.allocated(cached_occupation))then
      if(all(shape(cached_source)==shape(local)).and.all(shape(cached_occupation)==shape(system%rocc(:,:,1))))then
        if(all(cached_source==local).and.all(cached_occupation==system%rocc(:,:,1)))changed=0
      endif
    endif
    call comm_summation(changed,total_changed,info%icomm_k)
    if(total_changed==0)return
    allocate(send(ng,no,nk),allpsi(ng,no,nk),allw(ng,no,nk),w(ng,no,info%numk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=local
    call comm_summation(send,allpsi,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)then
      if(wannier%ng==0)call wannier_init(wannier,mg%num,num_kgrid,system%hgs,system%vec_k,hse_omega,status)
      if(status==0)then
        maxiter=0
        if(mod(wannier%updates,hse_mlwf_interval)==0)maxiter=hse_mlwf_maxiter
        call wannier_refresh_source(wannier,allpsi,system%rocc(:,:,1),maxiter,hse_mlwf_tolerance,status)
      endif
      if(status==0)call wannier_apply(wannier,allpsi,allw,status)
      if(status==0.and.(wannier%updates==1.or.maxiter>0))then
        write(*,'(a,3i7,3es16.7)')'HSE_WANNIER refresh/iterations/status/spread/gradient/overlap: ', &
        wannier%updates,wannier%localization_iterations,wannier%localization_status, &
        wannier%spread,wannier%gradient,wannier%min_singular
        if(wannier%localization_status/=0) &
          write(*,'(a)')'HSE_WANNIER: localization not converged; retaining full-support exact exchange.'
      endif
    endif
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: collective exchange refresh failed'
    call comm_bcast(allw,info%icomm_k,0)
    w=allw(:,:,info%ik_s:info%ik_e)
    call hse_ace_build(ace,local,w,system%hvol,status)
    call comm_summation(status,total_changed,info%icomm_k)
    if(total_changed/=0)error stop 'HSE Wannier: ACE construction metric failed'
    cached_source=local;cached_occupation=system%rocc(:,:,1);cached_action=w
    ex=0d0
    do ik=1,info%numk
      do j=1,no
        ex=ex+.125d0*system%rocc(j,info%ik_s+ik-1,1)*system%wtk(info%ik_s+ik-1)*system%hvol &
          *real(sum(conjg(local(:,j,ik))*w(:,j,ik)),8)
      enddo
    enddo
    call comm_summation(ex,hse_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_wannier_collective(target,action,info)
    use communication, only: comm_bcast
    implicit none
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable :: send(:,:,:),alltarget(:,:,:),allaction(:,:,:)
    integer :: ng,nt,nk,status
    ng=size(target,1);nt=size(target,2);nk=product(num_kgrid)
    allocate(send(ng,nt,nk),alltarget(ng,nt,nk),allaction(ng,nt,nk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=target
    call comm_summation(send,alltarget,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)call wannier_apply(wannier,alltarget,allaction,status)
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: full trial action failed'
    call comm_bcast(allaction,info%icomm_k,0)
    action=allaction(:,:,info%ik_s:info%ik_e)
  end subroutine

  subroutine hse_core_exchange(system,mg,info,psi,core,energy)
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    integer,intent(in) :: core(3)
    real(8),intent(out) :: energy
    real(8) :: local
    integer :: ix,iy,iz,io,ik,g
    if(.not.allocated(cached_action))error stop 'DC HSE: refreshed exchange action missing'
    local=0d0
    do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do iz=mg%is(3),min(mg%ie(3),core(3))
          do iy=mg%is(2),min(mg%ie(2),core(2))
            do ix=mg%is(1),min(mg%ie(1),core(1))
              g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))
              local=local+.125d0*system%rocc(io,ik,1)*system%wtk(ik)*system%hvol &
                *real(conjg(psi%zwf(ix,iy,iz,1,io,ik,1))*cached_action(g,io,ik-info%ik_s+1),8)
            enddo
          enddo
        enddo
      enddo
    enddo
    call comm_summation(local,energy,info%icomm_k)
  end subroutine
end module
