program cufft_stub_probe
  use exx_batch_backend, only: s_exx_batch_backend,local_backend_factory
  use exx_cufft, only: exx_cufft_apply,s_exx_cufft,exx_cufft_create
  implicit none
  class(s_exx_batch_backend),allocatable :: backend
  procedure(local_backend_factory),pointer :: factory=>exx_cufft_create
  integer :: padded(3)=[2,3,4],indices(2)=[1,24],status
  complex(8) :: filter(2,3,4),source(2),targets(2,2),action(2,2)
  complex(8) :: no_targets(2,0),no_action(2,0),no_source(0),empty_targets(0,2),empty_action(0,2)
  integer :: no_indices(0)
  complex(8) :: empty_batch(0,0),empty_batch_action(0,0)
  filter=1d0;source=cmplx(.3d0,.7d0,8);targets=cmplx(.2d0,-.4d0,8)
  call exx_cufft_apply(padded,indices,filter,source,targets,action,status)
  if(status/=-1)error stop 'Disabled cuFFT must reject nonempty work explicitly'
  if(any(action/=(0d0,0d0)))error stop 'Disabled cuFFT output must be zero'
  call exx_cufft_apply(padded,indices,filter,source,no_targets,no_action,status)
  if(status/=0)error stop 'Empty batch is a successful no-op'
  call exx_cufft_apply(padded,no_indices,filter,no_source,empty_targets,empty_action,status)
  if(status/=0)error stop 'Empty support is a successful no-op'
  call exx_cufft_apply(padded,[1,25],filter,source,targets,action,status)
  if(status/=-2)error stop 'Out-of-range index must fail validation'
  call exx_cufft_apply(padded,[1,1],filter,source,targets,action,status)
  if(status/=-2)error stop 'Duplicate index must fail validation'
  call exx_cufft_apply(padded,indices,filter(:,1:2,:),source,targets,action,status)
  if(status/=-2)error stop 'Filter shape mismatch must fail validation'
  ! Exercise dispatch through the neutral abstract type/factory, not a copied owner.
  call factory(backend)
  if(.not.allocated(backend))error stop 'Factory failed to allocate backend'
  call check_zero_counters()
  call backend%apply(targets,action,status)
  if(status/=-2.or.any(action/=(0d0,0d0)))error stop 'Apply before prepare must fail'
  call backend%prepare(padded,no_indices,filter,no_source,2,status)
  if(status/=0)error stop 'Empty stateful prepare failed'
  call backend%apply(empty_targets,empty_action,status)
  if(status/=0)error stop 'Empty stateful support apply failed'
  call backend%apply(empty_batch,empty_batch_action,status)
  if(status/=0)error stop 'Empty stateful batch failed'
  call backend%prepare(padded,no_indices,filter,no_source,2,status)
  if(status/=0)error stop 'Repeated empty prepare failed'
  call check_zero_counters()
  call backend%release(status)
  if(status/=0)error stop 'First stateful release failed'
  call backend%release(status)
  if(status/=0)error stop 'Repeated stateful release failed'
  call backend%prepare(padded,indices,filter,source,0,status)
  if(status/=-2)error stop 'Zero capacity accepted'
  call backend%prepare(padded,indices,filter,source,-1,status)
  if(status/=-2)error stop 'Negative capacity accepted'
  call backend%prepare(padded,indices,filter,source,2,status)
  if(status/=-1)error stop 'Disabled stateful backend accepted nonempty prepare'
  call backend%apply(targets,action,status)
  if(status/=-2.or.any(action/=(0d0,0d0)))error stop 'Failed prepare left a usable backend'
  call check_zero_counters()
  call backend%prepare(padded,no_indices,filter,no_source,2,status)
  if(status/=0)error stop 'Prepare before implicit finalization failed'
  deallocate(backend) ! Finalization must release a prepared object without an explicit release.
  call factory(backend)
  call check_zero_counters()
  call backend%release(status)
  if(status/=0)error stop 'Fresh unprepared release failed'
  deallocate(backend)
  print *, 'PASS CPU cuFFT stub: disabled backend, empty work, invalid inputs, factory and stateful lifecycle'
contains
  subroutine check_zero_counters()
    implicit none
    select type(backend)
    type is(s_exx_cufft)
      if(any([backend%plan_builds,backend%filter_uploads,backend%index_uploads, &
        backend%source_uploads,backend%batch_uploads]/=0))error stop 'CPU stub counters must remain zero'
    class default
      error stop 'Factory returned the wrong dynamic type'
    end select
  end subroutine
end program
