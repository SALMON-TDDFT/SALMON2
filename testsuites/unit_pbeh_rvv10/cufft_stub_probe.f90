program cufft_stub_probe
  use exx_cufft, only: exx_cufft_apply
  implicit none
  integer :: padded(3)=[2,3,4],indices(2)=[1,24],status
  complex(8) :: filter(2,3,4),source(2),targets(2,2),action(2,2)
  complex(8) :: no_targets(2,0),no_action(2,0),no_source(0),empty_targets(0,2),empty_action(0,2)
  integer :: no_indices(0)
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
  print *, 'PASS CPU cuFFT stub: disabled backend, empty work, invalid indices and shape'
end program
