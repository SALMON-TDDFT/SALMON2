! Native RT adapter for core-partitioned LCFO screened exchange.
! Each rank owns one core. Hartree, XC, propagation and current stay native.
! Full fragment support only: density eigfactors are not MLWFs or pair pruning.
module hse_lcfo_rt
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
  use communication, only: comm_bcast,comm_summation
  use salmon_global, only: hse_omega
  use lcfo_rt_basis
  use hse_wannier, only: s_hse_wannier,wannier_init,wannier_apply,wannier_forward
  use hse_ace, only: hse_ace_state,hse_ace_build,hse_ace_apply,hse_ace_average
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_hse_refresh,lcfo_hse_add_action,lcfo_hse_stage
  complex(8),allocatable,save :: fragment_basis(:,:),core_basis(:,:),hx(:,:),initial_hx(:,:),midpoint_hx(:,:)
  integer,allocatable,save :: selected(:)
  type(s_hse_wannier),save :: fragment_operator
  type(hse_ace_state),save :: ace,initial_ace,midpoint_ace
  logical,save :: ace_valid=.false.,initial_ace_valid=.false.,midpoint_ace_valid=.false.,use_midpoint=.false.
  integer,save :: refresh_count=0
contains
  subroutine initialize_fragment()
    complex(8),allocatable :: block(:,:)
    integer,allocatable :: mapping(:,:),first(:)
    integer :: nf,ns(3),ng,f,g,x,y,z,p(3),global_point(3),rel(3),j,nsel,ierr
    nf=size(lcfo_counts);ns=lcfo_core+2*lcfo_buffer;ng=product(ns)
    if(any(ns>lcfo_grid))error stop 'LCFO HSE: fragment exceeds global periodic grid'
    allocate(mapping(ng,nf),first(nf));mapping=0;first=0
    ! FFT order is core/right-buffer then the periodic left buffer.
    g=0
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      g=g+1;p=[x,y,z]
      where(p>=lcfo_core+lcfo_buffer)p=p-ns
      global_point=modulo(lcfo_origins(:,lcfo_rank+1)+p,lcfo_grid)
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
    call wannier_init(fragment_operator,ns,[1,1,1],lcfo_h,reshape([0d0,0d0,0d0],[3,1]),hse_omega,ierr)
    if(ierr/=0)error stop 'LCFO HSE: fragment periodic exchange initialization failed'
    if(lcfo_rank==0)write(*,'(a,3i6,a)')'LCFO HSE fragment grid:',ns,'; full support; no MLWF cutoff'
  end subroutine

  subroutine pack_coefficients(psi,system,mg,info,coeff)
    type(s_orbital),intent(in) :: psi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable,intent(out) :: coeff(:,:)
    complex(8),allocatable :: grid(:,:)
    integer :: io,is(3),ie(3)
    if(.not.lcfo_rt_active)error stop 'LCFO HSE: inactive LCFO basis'
    if(system%nk/=1.or.system%nspin/=1.or.info%io_s/=1.or.info%io_e/=system%no.or. &
       info%ik_s/=1.or.info%ik_e/=1.or.info%numm/=1)error stop 'LCFO HSE: requires Gamma/all orbitals per core'
    if(any(mg%num/=lcfo_core))error stop 'LCFO HSE: native grid does not match LCFO core'
    if(.not.allocated(psi%zwf))error stop 'LCFO HSE: complex wavefunctions required'
    is=mg%is;ie=mg%ie
    allocate(grid(product(mg%num),system%no),coeff(sum(lcfo_counts),system%no))
    do io=1,system%no
      grid(:,io)=reshape(psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[product(mg%num)])
    enddo
    call lcfo_collect_coefficients(grid,coeff)
  end subroutine

  subroutine lcfo_hse_refresh(system,mg,info,psi,exchange_energy)
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    real(8),intent(out) :: exchange_energy
    complex(8),allocatable :: coeff(:,:),weighted(:,:),density(:,:),work(:),source(:,:,:),action(:,:,:)
    complex(8),allocatable :: projected(:,:),local_hx(:,:),u(:,:,:),w(:,:,:)
    real(8),allocatable :: eigenvalues(:),rwork(:)
    real(8) :: threshold,discarded
    integer :: nb,nsel,ng,no,j,k,ierr,nrank,first
    external :: zheev
    call pack_coefficients(psi,system,mg,info,coeff)
    if(.not.allocated(fragment_basis))call initialize_fragment()
    nb=size(coeff,1);no=size(coeff,2);nsel=size(selected);ng=size(fragment_basis,1)
    if(any(system%rocc(:,1,1)<0d0).or.any(system%rocc(:,1,1)>2d0).or. &
       .not.all(ieee_is_finite(system%rocc(:,1,1))))error stop 'LCFO HSE: invalid occupations'
    allocate(weighted(nsel,no),density(nsel,nsel),eigenvalues(nsel),work(max(1,2*nsel)),rwork(max(1,3*nsel-2)))
    do j=1,no
      weighted(:,j)=coeff(selected,j)*sqrt(system%rocc(j,1,1)/2d0)
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
    if(allocated(fragment_operator%source))deallocate(fragment_operator%source)
    allocate(fragment_operator%source(ng,nrank))
    call wannier_forward(fragment_operator,source,fragment_operator%source)
    allocate(action(ng,nsel,1))
    call wannier_apply(fragment_operator,reshape(fragment_basis,[ng,nsel,1]),action,ierr)
    if(ierr/=0)error stop 'LCFO HSE: fragment exchange action failed'
    projected=.25d0*lcfo_dv*matmul(transpose(conjg(core_basis)),action(:,:,1))
    projected=.5d0*(projected+transpose(conjg(projected)))
    allocate(local_hx(nb,nb));local_hx=0d0
    do k=1,nsel;do j=1,nsel
      local_hx(selected(j),selected(k))=projected(j,k)
    enddo;enddo
    if(.not.allocated(hx))allocate(hx(nb,nb))
    call comm_summation(local_hx,hx,size(hx),lcfo_comm)
    if(.not.all(ieee_is_finite(real(hx))).or..not.all(ieee_is_finite(aimag(hx)))) &
      error stop 'LCFO HSE: nonfinite projected exchange'
    allocate(u(nb,no,1),w(nb,no,1));u(:,:,1)=coeff;w(:,:,1)=matmul(hx,coeff)
    exchange_energy=0d0
    do j=1,no
      exchange_energy=exchange_energy+.5d0*system%rocc(j,1,1)*real(sum(conjg(coeff(:,j))*w(:,j,1)),8)
    enddo
    call hse_ace_build(ace,u,w,1d0,ierr);ace_valid=ierr==0
    refresh_count=refresh_count+1
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

  subroutine lcfo_hse_stage(stage)
    integer,intent(in) :: stage
    integer :: ierr
    select case(stage)
    case(0)
      if(.not.allocated(hx))error stop 'LCFO HSE: initial operator missing'
      initial_hx=hx;initial_ace=ace;initial_ace_valid=ace_valid;use_midpoint=.false.
    case(1)
      if(.not.allocated(initial_hx).or..not.allocated(hx))error stop 'LCFO HSE: midpoint endpoints missing'
      midpoint_hx=.5d0*(initial_hx+hx)
      midpoint_ace_valid=initial_ace_valid.and.ace_valid
      if(midpoint_ace_valid)then
        call hse_ace_average(initial_ace,ace,midpoint_ace,ierr)
        midpoint_ace_valid=ierr==0
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
    complex(8),allocatable :: coeff(:,:),u(:,:,:),w(:,:,:),grid(:,:)
    integer :: nb,no,ierr,lo,hi,io,is(3),ie(3)
    if(.not.allocated(hx))error stop 'LCFO HSE: refresh required before action'
    call pack_coefficients(psi,system,mg,info,coeff)
    nb=size(coeff,1);no=size(coeff,2);allocate(u(nb,no,1),w(nb,no,1));u(:,:,1)=coeff
    ierr=0
    if(use_midpoint)then
      if(midpoint_ace_valid)then
        call hse_ace_apply(midpoint_ace,u,w,ierr)
      else
        w(:,:,1)=matmul(midpoint_hx,coeff)
      endif
    else
      if(ace_valid)then
        call hse_ace_apply(ace,u,w,ierr)
      else
        w(:,:,1)=matmul(hx,coeff)
      endif
    endif
    if(ierr/=0)error stop 'LCFO HSE: ACE application failed'
    lo=lcfo_offsets(lcfo_rank+1)+1;hi=lcfo_offsets(lcfo_rank+2)
    grid=matmul(lcfo_basis,w(lo:hi,:,1));is=mg%is;ie=mg%ie
    do io=1,no
      hpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)= &
        hpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)+reshape(grid(:,io),mg%num)
    enddo
    hpsi%update_zwf_overlap=.false.
  end subroutine
end module
