! Full spin density matrix exchange. Rows pack (real-space point, spin),
! sources may be weighted by sqrt(f), with no closed-shell factor.
module exx_spinor
 implicit none
 private
 public :: exx_spinor_apply_distributed,exx_spinor_weight_sources
contains
 subroutine exx_spinor_weight_sources(orbitals,occupation,source,ierr)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  complex(8),intent(in)::orbitals(:,:,:)
  real(8),intent(in)::occupation(:,:)
  complex(8),intent(out)::source(:,:,:)
  integer,intent(out)::ierr
  integer::i,j,k
  ierr=1;source=0d0
  if(any(shape(source)/=shape(orbitals)))return
  if(size(occupation,1)/=size(orbitals,2).or.size(occupation,2)/=size(orbitals,3))return
  do k=1,size(occupation,2);do j=1,size(occupation,1)
   if(.not.ieee_is_finite(occupation(j,k)))return
   if(occupation(j,k)<0d0.or.occupation(j,k)>1d0)return
   do i=1,size(orbitals,1)
    if(.not.ieee_is_finite(real(orbitals(i,j,k),8)).or..not.ieee_is_finite(aimag(orbitals(i,j,k))))return
   enddo
  enddo;enddo
  do k=1,size(occupation,2);do j=1,size(occupation,1)
   source(:,j,k)=orbitals(:,j,k)*sqrt(occupation(j,k))
  enddo;enddo
  ierr=0
 end subroutine
 subroutine exx_spinor_apply_distributed(op,source,target,action,starts,counts,rank,comm,transpose_tiles,ierr,pair_tiles)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use exx_k_exchange, only: exx_k_kernel,exx_k_kernel_apply_distributed,exx_k_kernel_apply_pair_distributed
  implicit none
  type(exx_k_kernel),intent(inout)::op
  complex(8),intent(in)::source(:,:,:),target(:,:,:)
  complex(8),intent(out)::action(:,:,:)
  integer,intent(in)::starts(:),counts(:),rank,comm
  integer,intent(out)::ierr
  logical,optional,intent(in)::pair_tiles
  logical::use_pairs
  interface
   subroutine transpose_tiles(send,recv,count,comm)
    implicit none
    complex(8),intent(in)::send(:)
    complex(8),intent(out)::recv(:)
    integer,intent(in)::count,comm
   end subroutine
  end interface
  complex(8),allocatable::left(:,:,:),right(:,:,:),part(:,:,:)
  complex(8)::sent(size(counts)),received(size(counts))
  logical::valid
  integer::ng,no,nk,s,t,j,b,i
  interface
   subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(1),intent(in)::ta,tb
    integer,intent(in)::m,n,k,lda,ldb,ldc
    complex(8),intent(in)::alpha,beta,a(lda,*),b(ldb,*)
    complex(8),intent(inout)::c(ldc,*)
   end subroutine
  end interface
  use_pairs=.false.
  if(present(pair_tiles))use_pairs=pair_tiles
  ierr=1;action=0d0
  ng=op%ng;no=size(source,2);nk=size(source,3)
  valid=size(source,1)==2*ng.and.size(target,1)==2*ng.and.all(shape(action)==shape(target))
  valid=valid.and.allocated(op%phase)
  if(allocated(op%phase))valid=valid.and.size(op%phase,2)==nk
  do j=1,nk;do b=1,no;do i=1,size(source,1)
   if(.not.ieee_is_finite(real(source(i,b,j),8)).or..not.ieee_is_finite(aimag(source(i,b,j))))valid=.false.
  enddo;enddo;enddo
  do j=1,size(target,3);do b=1,size(target,2);do i=1,size(target,1)
   if(.not.ieee_is_finite(real(target(i,b,j),8)).or..not.ieee_is_finite(aimag(target(i,b,j))))valid=.false.
  enddo;enddo;enddo
  ! All k peers agree before any callback can consume cross-spin density.
  sent=cmplx(merge(0,1,valid),merge(1,0,use_pairs),8)
  call transpose_tiles(sent,received,1,comm)
  if(any(real(received)/=0d0).or.any(aimag(received)/=real(merge(1,0,use_pairs),8)))return
  allocate(left(ng,no,nk),right(ng,no,nk),part(ng,size(target,2),nk))
  do s=1,2;do t=1,2
   do j=1,nk;do b=1,no
    left(:,b,j)=source((s-1)*ng+1:s*ng,b,j)*op%phase(:,j)
    right(:,b,j)=source((t-1)*ng+1:t*ng,b,j)*op%phase(:,j)
   enddo;enddo
   ! Callback supplies cross-spin density; scalar backend handles k transposes,
   ! the screened interaction, target phases and uniform Brillouin-zone weights.
   if(use_pairs)then
    call exx_k_kernel_apply_pair_distributed(op,target((t-1)*ng+1:t*ng,:,:),part, &
     starts,counts,rank,comm,transpose_tiles,fill_density,ierr)
   else
   call exx_k_kernel_apply_distributed(op,source((s-1)*ng+1:s*ng,:,:), &
    target((t-1)*ng+1:t*ng,:,:),part,starts,counts,rank,comm,transpose_tiles,ierr,fill_density)
   endif
   if(ierr/=0)return
   action((s-1)*ng+1:s*ng,:,:)=action((s-1)*ng+1:s*ng,:,:)+part
  enddo;enddo
 contains
  subroutine fill_density(j,lo,rows,density)
   implicit none
   integer,intent(in)::j,lo,rows
   complex(8),intent(out)::density(:,:)
   density=0d0
   call zgemm('N','C',rows,ng,no,(1d0,0d0),left(lo,1,j),ng,right(1,1,j),ng, &
    (0d0,0d0),density,size(density,1))
  end subroutine
 end subroutine
end module
