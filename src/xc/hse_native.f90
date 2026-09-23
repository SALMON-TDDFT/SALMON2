! SALMON adapter: initial certified layout is complete grid/orbitals per rank,
! distributed k points. Exchange builds share sources; ACE applications stay local.
module hse_native
  use iso_fortran_env, only: int64
  use structures
  use plusU_global, only: PLUS_U_ON
  use hse_exchange
  use hse_ace
  use communication, only: comm_summation
  use salmon_global, only: xc,yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_md,yn_symmetrized_stencil
  implicit none
  private
  public :: hse_enabled,hse_refresh,hse_add_action,hse_exchange_energy,hse_freeze
  public :: hse_pack,hse_unpack,hse_timings,hse_walltime
  type(hse_kernel),save :: kernel
  type(hse_ace_state),save :: ace
  complex(8),allocatable,save :: cached_source(:,:,:),target_work(:,:,:),action_work(:,:,:),output_work(:,:,:)
  real(8),save :: hse_exchange_energy=0d0
  real(8),save :: hse_timings(4)=0d0 ! full EXX, ACE build, ACE apply, EXX collectives
  logical,save :: hse_freeze=.false.
contains
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
    complex(8),allocatable :: part(:,:,:),source(:,:,:),w(:,:,:),local(:,:,:)
    real(8) :: ex,offdiag(3,3),tick
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
      call hse_kernel_init(kernel,n,mesh,system%hgs(1),system%vec_k,.11d0,16,ierr)
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
    allocate(part(ng,no,nk),source(ng,no,nk),w(ng,no,nk))
    part=0;part(:,:,info%ik_s:info%ik_e)=local
    tick=hse_walltime()
    call comm_summation(part,source,size(source),info%icomm_k)
    hse_timings(4)=hse_timings(4)+hse_walltime()-tick
    tick=hse_walltime()
    call hse_kernel_apply(kernel,source,source,part,info%id_k,info%isize_k,ierr)
    hse_timings(1)=hse_timings(1)+hse_walltime()-tick
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: full exchange action failed'
    tick=hse_walltime()
    call comm_summation(part,w,size(w),info%icomm_k)
    hse_timings(4)=hse_timings(4)+hse_walltime()-tick
    tick=hse_walltime()
    call hse_ace_build(ace,local,w(:,:,info%ik_s:info%ik_e),system%hvol,ierr)
    hse_timings(2)=hse_timings(2)+hse_walltime()-tick
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: ACE metric failed'
    cached_source=local
    ex=.25d0*real(sum(conjg(local)*w(:,:,info%ik_s:info%ik_e)),8)*system%hvol/nk
    call comm_summation(ex,hse_exchange_energy,info%icomm_k)
  end subroutine

  subroutine hse_add_action(psi,hpsi,system,mg,info)
    type(s_orbital),intent(in) :: psi
    type(s_orbital),intent(inout) :: hpsi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ierr,ng
    real(8) :: tick
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
    call hse_ace_apply(ace,target_work,action_work,ierr)
    hse_timings(3)=hse_timings(3)+hse_walltime()-tick
    if(ierr/=0)error stop 'HSE06: ACE application failed'
    call hse_pack(hpsi,mg,info,output_work)
    output_work=output_work+.25d0*action_work
    call hse_unpack(output_work,hpsi,mg,info)
  end subroutine
end module
