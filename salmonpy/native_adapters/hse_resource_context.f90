! Copyright 2026 SALMON developers. Apache-2.0.
! Opt-in ownership adapter for pinned fast-PR kernels, not a full TDDFT API.
! Callers serialize calls, never assign/copy contexts, and explicitly close them.
! No process-global settings, active pointer, MPI communicator or output unit.
module hse_resource_context
 use iso_fortran_env,only:int64
 use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
 use exx_k_exchange,only:exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply,exx_k_kernel_destroy
 use exx_ace,only:s_exx_ace,exx_ace_build,exx_ace_apply
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear
 implicit none
 private
 public :: hse_context,hse_open,hse_refresh,hse_apply,hse_close
 type hse_context
  private
  ! Allocatable ownership is moved/allocated, never copied: kernel contains FFTW pointers.
  type(exx_k_kernel),allocatable :: kernel
  type(s_exx_rt_state) :: rt
  integer :: ng=0,nk=0
  real(8) :: dv=0d0,fraction=0d0
 end type
contains
 subroutine hse_open(ctx,n,mesh,h,k,omega,fraction,status)
  type(hse_context),target,intent(inout)::ctx
  integer,intent(in)::n(3),mesh(3)
  real(8),intent(in)::h(3),k(:,:),omega,fraction
  integer,intent(out)::status
  integer::axis,j,err
  integer(int64)::ng,nk
  real(8)::dv
  status=1
  if(allocated(ctx%kernel))return
  if(any(n<1).or.any(mesh<1))return
  ng=1;nk=1
  do axis=1,3
   if(ng>int(huge(1),int64)/n(axis).or.nk>int(huge(1),int64)/mesh(axis))return
   ng=ng*n(axis);nk=nk*mesh(axis)
   if(.not.ieee_is_finite(h(axis)))return
   if(h(axis)<=0d0)return
  enddo
  if(ng>int(huge(1),int64)/nk)return
  if(size(k,1)/=3.or.size(k,2)/=nk)return
  do j=1,size(k,2)
   do axis=1,3
    if(.not.ieee_is_finite(k(axis,j)))return
   enddo
  enddo
  if(.not.ieee_is_finite(omega).or..not.ieee_is_finite(fraction))return
  if(omega<0d0.or.fraction<0d0.or.fraction>1d0)return
  dv=product(h)
  if(.not.ieee_is_finite(dv).or.dv<=0d0)return
  allocate(ctx%kernel,stat=err)
  if(err/=0)then
   status=2;return
  endif
  ! Deterministic layout; no shared automatic-layout tuning or accelerators.
  call exx_k_kernel_init(ctx%kernel,n,mesh,h,k,omega,min(16,int(ng)),err,fft_layout='strided')
  if(err/=0)then
   call hse_close(ctx);status=3;return
  endif
  ctx%ng=int(ng);ctx%nk=int(nk);ctx%dv=dv;ctx%fraction=fraction
  status=0
 end subroutine

 subroutine hse_refresh(ctx,source,status)
  type(hse_context),target,intent(inout)::ctx
  complex(8),intent(in)::source(:,:,:)
  integer,intent(out)::status
  type(s_exx_ace)::candidate
  complex(8),allocatable::action(:,:,:),accepted(:,:,:)
  integer::err,i,j,k
  complex(8)::overlap
  status=1
  if(.not.allocated(ctx%kernel))return
  if(size(source,1)/=ctx%ng.or.size(source,3)/=ctx%nk.or.size(source,2)<1)return
  if(.not.finite_orbitals(source))return
  ! This primitive accepts an occupied orthonormal frame, not arbitrary density factors.
  if(size(source,2)>ctx%ng)return
  do k=1,ctx%nk
   do j=1,size(source,2)
    do i=1,j
     overlap=dot_product(source(:,i,k),source(:,j,k))*ctx%dv
     if(.not.ieee_is_finite(real(overlap,8)).or..not.ieee_is_finite(aimag(overlap)))return
     if(i==j)overlap=overlap-1d0
     if(abs(overlap)>1d-8)return
    enddo
   enddo
  enddo
  allocate(action(ctx%ng,size(source,2),ctx%nk),accepted(ctx%ng,size(source,2),ctx%nk),stat=err)
  if(err/=0)then
   status=2;return
  endif
  accepted=source
  call exx_k_kernel_apply(ctx%kernel,source,source,action,0,1,err)
  if(err/=0)then
   status=3;return
  endif
  call exx_ace_build(candidate,source,action,ctx%dv,err)
  if(err/=0)then
   status=3;return
  endif
  ! Only the value-owned ACE/source commit; never assign an FFT-owning kernel.
  ctx%rt%ace=candidate
  call move_alloc(accepted,ctx%rt%cached_source)
  status=0
 end subroutine

 subroutine hse_apply(ctx,target,action,status)
  type(hse_context),intent(in)::ctx
  complex(8),intent(in)::target(:,:,:)
  complex(8),intent(inout)::action(:,:,:)
  integer,intent(out)::status
  complex(8),allocatable::candidate(:,:,:)
  integer::err
  status=1
  if(.not.allocated(ctx%kernel).or..not.allocated(ctx%rt%cached_source))return
  if(size(target,1)/=ctx%ng.or.size(target,3)/=ctx%nk.or.size(target,2)<1)return
  if(any(shape(action)/=shape(target)).or..not.finite_orbitals(target))return
  allocate(candidate(size(target,1),size(target,2),size(target,3)),stat=err)
  if(err/=0)then
   status=2;return
  endif
  call exx_ace_apply(ctx%rt%ace,target,candidate,err)
  if(err/=0)then
   status=3;return
  endif
  candidate=ctx%fraction*candidate
  if(.not.finite_orbitals(candidate))then
   status=3;return
  endif
  action=candidate;status=0
 end subroutine

 subroutine hse_close(ctx)
  type(hse_context),intent(inout)::ctx
  call rt_state_clear(ctx%rt)
  if(allocated(ctx%kernel))then
   call exx_k_kernel_destroy(ctx%kernel)
   deallocate(ctx%kernel)
  endif
  ctx%ng=0;ctx%nk=0;ctx%dv=0;ctx%fraction=0
 end subroutine

 logical function finite_orbitals(u) result(ok)
  complex(8),intent(in)::u(:,:,:)
  integer::i,j,k
  ok=.false.
  do k=1,size(u,3)
   do j=1,size(u,2)
    do i=1,size(u,1)
     if(.not.ieee_is_finite(real(u(i,j,k),8)))return
     if(.not.ieee_is_finite(aimag(u(i,j,k))))return
    enddo
   enddo
  enddo
  ok=.true.
 end function
end module
