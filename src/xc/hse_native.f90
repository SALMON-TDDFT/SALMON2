#include "config.h"
! SALMON adapter: initial certified layout is complete grid/orbitals per rank,
! distributed k points. Exchange transfers density tiles; ACE applications stay local.
module hse_native
  use iso_fortran_env, only: int64
  use structures
  use plusU_global, only: PLUS_U_ON
  use hse_exchange
  use hse_ace
  use communication, only: comm_summation,comm_alltoall
  use salmon_global, only: xc,yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_md,yn_symmetrized_stencil,propagator
  implicit none
  private
  public :: hse_enabled,hse_refresh,hse_add_action,hse_exchange_energy,hse_freeze
  public :: hse_pack,hse_unpack,hse_timings,hse_walltime
  public :: hse_taylor_stage
  type(hse_kernel),save :: kernel
  type(hse_ace_state),save :: ace
  type(hse_ace_state),save :: initial_ace,midpoint_ace
  complex(8),allocatable,save :: full_source(:,:,:),initial_source(:,:,:),midpoint_source(:,:,:)
  logical,save :: taylor_active=.false.,taylor_midpoint=.false.
  complex(8),allocatable,save :: cached_source(:,:,:),target_work(:,:,:),action_work(:,:,:),output_work(:,:,:)
  real(8),save :: hse_exchange_energy=0d0
  real(8),save :: hse_timings(4)=0d0 ! full EXX, ACE build, ACE apply, EXX collectives
  logical,save :: hse_freeze=.false.,reported_team=.false.
contains
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
    integer :: ierr,total_error,ng,nk,no,n,mesh,j
    if(.not.hse_enabled().or.hse_freeze)return
    if(info%isize_r/=1.or.info%isize_o/=1.or.info%numm/=1) &
      error stop 'HSE06: initial native support requires k-only MPI distribution'
    if(yn_periodic/='y'.or.system%nspin/=1.or..not.allocated(psi%zwf)) &
      error stop 'HSE06: periodic complex unpolarized orbitals required'
    if(yn_spinorbit/='n'.or.yn_jm/='n'.or.yn_dc/='n'.or.yn_md/='n'.or.yn_symmetrized_stencil=='y') &
      error stop 'HSE06: unsupported Hamiltonian/ionic extension'
    if(PLUS_U_ON)error stop 'HSE06: DFT+U combination unsupported'
    if(allocated(system%Ac_micro%v))error stop 'HSE06: microscopic vector potential unsupported'
    ng=product(mg%num);nk=system%nk;no=system%no;n=mg%num(1);mesh=nint(real(nk,8)**(1d0/3d0))
    if(any(mg%num/=n).or.mesh**3/=nk.or.maxval(abs(system%hgs-system%hgs(1)))>1d-12) &
      error stop 'HSE06: cubic grid and full cubic k mesh required'
    offdiag=system%primitive_a
    do j=1,3;offdiag(j,j)=0;enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE06: orthogonal cell required'
    if(maxval(abs(system%wtk-1d0/nk))>1d-12.or.maxval(abs(system%rocc-2d0))>1d-12) &
      error stop 'HSE06: uniform k weights and fully occupied spin pairs required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE06: unsupported orbital layout'
    if(kernel%n==0)then
      ! Up to 64 rows per collective round, at most 16 per rank; bound tile memory.
      call hse_kernel_init(kernel,n,mesh,system%hgs(1),system%vec_k,.11d0, &
        max(1,min(16,64/info%isize_k)),ierr,info%ik_s,info%numk)
      call comm_summation(ierr,total_error,info%icomm_rko)
      if(total_error/=0)error stop 'HSE06: kernel initialization failed'
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
    tick=hse_walltime();communication_before=hse_timings(4)
    call apply_distributed(local,local,w,info,ierr)
    hse_timings(1)=hse_timings(1)+hse_walltime()-tick-(hse_timings(4)-communication_before)
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: distributed exchange action failed'
    if(.not.reported_team.and.info%id_k==0)write(*,'(a,i0)')'HSE_OPENMP threads=',kernel%threads_used
    reported_team=.true.
    tick=hse_walltime()
    call hse_ace_build(ace,local,w,system%hvol,ierr)
    hse_timings(2)=hse_timings(2)+hse_walltime()-tick
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: ACE metric failed'
    cached_source=local
    ex=.25d0*real(sum(conjg(local)*w),8)*system%hvol/nk
    call comm_summation(ex,hse_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_distributed(source,target,action,info,ierr)
    complex(8),intent(in) :: source(:,:,:),target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    integer,intent(out) :: ierr
    integer :: local_layout(2*info%isize_k),layout(2*info%isize_k),np
    np=info%isize_k;local_layout=0
    local_layout(info%id_k+1)=info%ik_s
    local_layout(np+info%id_k+1)=info%numk
    call comm_summation(local_layout,layout,size(layout),info%icomm_k)
    call hse_kernel_apply_distributed(kernel,source,target,action,layout(:np),layout(np+1:), &
      info%id_k,transpose_tiles,ierr)
  contains
    subroutine transpose_tiles(send,recv,count)
      complex(8),intent(in) :: send(:)
      complex(8),intent(out) :: recv(:)
      integer,intent(in) :: count
      real(8) :: start
      start=hse_walltime()
      call comm_alltoall(send,recv,info%icomm_k,count)
      hse_timings(4)=hse_timings(4)+hse_walltime()-start
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
    tick=hse_walltime()
    if(taylor_active.and.propagator=='hse_taylor4_full')then
      communication_before=hse_timings(4)
      if(taylor_midpoint)then
        call apply_distributed(midpoint_source,target_work,action_work,info,ierr)
      else
        call apply_distributed(initial_source,target_work,action_work,info,ierr)
      endif
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'HSE Taylor full target action failed'
      hse_timings(1)=hse_timings(1)+hse_walltime()-tick-(hse_timings(4)-communication_before)
    else
      if(taylor_active.and.taylor_midpoint)then
        call hse_ace_apply(midpoint_ace,target_work,action_work,ierr)
      else
        call hse_ace_apply(ace,target_work,action_work,ierr)
      endif
      hse_timings(3)=hse_timings(3)+hse_walltime()-tick
    endif
    if(ierr/=0)error stop 'HSE06: ACE application failed'
    call hse_pack(hpsi,mg,info,output_work)
    output_work=output_work+.25d0*action_work
    call hse_unpack(output_work,hpsi,mg,info)
  end subroutine
end module
