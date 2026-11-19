! Stateful compact-pair backend. Ownership belongs to one exchange action.
module exx_batch_backend
 implicit none
 private
 public :: s_exx_batch_backend,local_backend_factory
 type,abstract :: s_exx_batch_backend
 contains
  procedure(prepare_batch),deferred :: prepare
  procedure(apply_batch),deferred :: apply
  procedure(release_batch),deferred :: release
 end type
 abstract interface
  subroutine prepare_batch(self,padded,indices,filter,source,capacity,status)
   import s_exx_batch_backend
   implicit none
   class(s_exx_batch_backend),target,intent(inout) :: self
   integer,intent(in) :: padded(3),indices(:),capacity
   complex(8),intent(in) :: filter(:,:,:),source(:)
   integer,intent(out) :: status
  end subroutine
  subroutine apply_batch(self,targets,action,status)
   import s_exx_batch_backend
   implicit none
   class(s_exx_batch_backend),target,intent(inout) :: self
   complex(8),intent(in) :: targets(:,:)
   complex(8),intent(out) :: action(:,:)
   integer,intent(out) :: status
  end subroutine
  subroutine release_batch(self,status)
   import s_exx_batch_backend
   implicit none
   class(s_exx_batch_backend),target,intent(inout) :: self
   integer,intent(out) :: status
  end subroutine
  subroutine local_backend_factory(backend)
   import s_exx_batch_backend
   implicit none
   class(s_exx_batch_backend),allocatable,intent(out) :: backend
  end subroutine
 end interface
end module
