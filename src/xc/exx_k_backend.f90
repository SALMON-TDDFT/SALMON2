! Owning accelerator interface for k-mesh density-tile convolution.
module exx_k_backend
 implicit none
 private
 public :: s_exx_k_backend,k_backend_factory
 type,abstract :: s_exx_k_backend
 contains
  procedure(prepare_k),deferred :: prepare
  procedure(apply_k),deferred :: apply
  procedure(release_k),deferred :: release
 end type
 abstract interface
  subroutine prepare_k(self,n,mesh,block,kernel,point,shift,slot,nslots,status)
   import s_exx_k_backend
   implicit none
   class(s_exx_k_backend),target,intent(inout) :: self
   integer,intent(in) :: n,mesh,block,point(:,:),shift(:,:),slot(:),nslots
   real(8),intent(in) :: kernel(0:,0:,0:)
   integer,intent(out) :: status
  end subroutine
  subroutine apply_k(self,lo,rows,buffer,action,status)
   import s_exx_k_backend
   implicit none
   class(s_exx_k_backend),target,intent(inout) :: self
   integer,intent(in) :: lo,rows
   complex(8),intent(in) :: buffer(:,:,:)
   complex(8),intent(out) :: action(:,:,:)
   integer,intent(out) :: status
  end subroutine
  subroutine release_k(self,status)
   import s_exx_k_backend
   implicit none
   class(s_exx_k_backend),target,intent(inout) :: self
   integer,intent(out) :: status
  end subroutine
  subroutine k_backend_factory(backend)
   import s_exx_k_backend
   implicit none
   class(s_exx_k_backend),allocatable,intent(out) :: backend
  end subroutine
 end interface
end module
