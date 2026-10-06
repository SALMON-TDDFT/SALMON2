program cufft_gpu_probe
  use exx_local_fft, only: s_exx_local_fft,exx_local_init,exx_local_prepare,exx_local_apply,exx_local_destroy
  use exx_batch_backend, only: s_exx_batch_backend
  use exx_cufft, only: exx_cufft_apply,s_exx_cufft,exx_cufft_create
  implicit none
  integer,parameter :: ns=7,nt=7
  type(s_exx_local_fft) :: plan
  class(s_exx_batch_backend),allocatable :: backend
  integer :: expected(5),saved_index
  integer :: points(3,ns),n(3)=[12,16,20],status,x,y,z,j,i,batch,first,last,pass
  integer :: batch_sizes(4)=[1,2,3,8]
  real(8) :: multiplier(12,16,20),error,scale
  complex(8) :: source(ns),targets(ns,nt),action(ns,nt),reference(ns,nt),density(ns),potential(ns)
  complex(8) :: no_targets(ns,0),no_action(ns,0)
  logical :: used
  points(:,1)=[11,15,19];points(:,2)=[0,0,0];points(:,3)=[11,1,2]
  points(:,4)=[0,15,2];points(:,5)=[11,0,1];points(:,6)=[0,1,19];points(:,7)=[0,0,2]
  do z=1,n(3);do y=1,n(2);do x=1,n(1)
    multiplier(x,y,z)=.7d0+.1d0*sin(.31d0*x+.17d0*y+.23d0*z)
  enddo;enddo;enddo
  call exx_local_init(plan,multiplier,status)
  if(status/=0)error stop 'FFTW kernel initialization'
  call exx_local_prepare(plan,points,used,status)
  if(status/=0.or..not.used)error stop 'FFTW compact preparation'
  if(any(plan%padded/=[3,5,8]))error stop 'Expected unequal noncubic padded dimensions'
  do i=1,ns
    source(i)=cmplx(.2d0+.03d0*i,.11d0*cos(real(i,8)),8)
    do j=1,nt
      targets(i,j)=cmplx(sin(.27d0*i+.31d0*j),cos(.43d0*i-.13d0*j),8)
    enddo
  enddo
  targets(:,4)=0d0
  do pass=1,2
    ! Repeated calls change source phases, so stale device data cannot pass.
    source=source*exp(cmplx(0d0,.37d0*pass,8))
    targets=targets*cmplx(.9d0,.1d0*pass,8)
    do j=1,nt
      density=conjg(source)*targets(:,j)
      call exx_local_apply(plan,density,potential,status)
      if(status/=0)error stop 'FFTW reference application'
      reference(:,j)=-source*potential
    enddo
    scale=max(1d0,maxval(abs(reference)))
    call exx_cufft_apply(plan%padded,plan%indices,plan%filter,source,targets,action,status)
    if(status/=0)error stop 'cuFFT application failed (GPU and CUDA runtime required)'
    error=maxval(abs(action-reference))/scale
    if(error>2d-12)error stop 'cuFFT full-batch differs from FFTW'
    if(maxval(abs(action(:,4)))>2d-14)error stop 'Zero target produced nonzero action'
    do batch=1,size(batch_sizes)
      action=cmplx(999d0,-999d0,8)
      do first=1,nt,batch_sizes(batch)
        last=min(nt,first+batch_sizes(batch)-1)
        call exx_cufft_apply(plan%padded,plan%indices,plan%filter,source, &
          targets(:,first:last),action(:,first:last),status)
        if(status/=0)error stop 'cuFFT partial batch failed'
      enddo
      error=maxval(abs(action-reference))/scale
      if(error>2d-12)error stop 'cuFFT partial batches differ from FFTW'
    enddo
  enddo
  call exx_cufft_apply(plan%padded,plan%indices,plan%filter,source,no_targets,no_action,status)
  if(status/=0)error stop 'Empty GPU batch failed'
  ! Keep one resident owner over several tail batches, then update each input independently.
  call exx_cufft_create(backend)
  expected=0;call check_counters()
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,3,status)
  if(status/=0)error stop 'Initial stateful GPU prepare failed'
  expected=[1,1,1,1,0];call check_counters()
  call apply_resident(3)
  call backend%apply(no_targets,no_action,status)
  if(status/=0)error stop 'Resident empty batch failed'
  call check_counters() ! Empty apply must not upload a target batch.
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,3,status)
  if(status/=0)error stop 'Unchanged stateful GPU prepare failed'
  call check_counters() ! No plan construction or constant uploads for unchanged data.
  targets=targets*cmplx(.8d0,-.2d0,8)
  call apply_resident(3)

  ! Nonuniform amplitude/phase change: a global phase alone cancels from exchange.
  source(2)=source(2)*cmplx(1.3d0,.4d0,8)
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,3,status)
  if(status/=0)error stop 'Changed source prepare failed'
  expected(4)=expected(4)+1;call check_counters()
  call apply_resident(3)
  plan%filter=plan%filter*cmplx(.7d0,.15d0,8)
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,3,status)
  if(status/=0)error stop 'Changed filter prepare failed'
  expected(2)=expected(2)+1;call check_counters()
  call apply_resident(3)
  saved_index=plan%indices(1);plan%indices(1)=plan%indices(2);plan%indices(2)=saved_index
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,3,status)
  if(status/=0)error stop 'Changed indices prepare failed'
  expected(3)=expected(3)+1;call check_counters()
  call apply_resident(3)

  call backend%prepare(plan%padded,plan%indices,plan%filter,source,4,status)
  if(status/=0)error stop 'Changed capacity prepare failed'
  expected(1:4)=expected(1:4)+1;call check_counters()
  call apply_resident(4)
  points(1,1)=10 ! Change the support box and FFT dimensions, keeping seven source points.
  call exx_local_prepare(plan,points,used,status)
  if(status/=0.or..not.used.or.any(plan%padded/=[5,5,8]))error stop 'Changed geometry FFTW prepare failed'
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,4,status)
  if(status/=0)error stop 'Changed geometry GPU prepare failed'
  expected(1:4)=expected(1:4)+1;call check_counters()
  call apply_resident(4)
  call backend%release(status)
  if(status/=0)error stop 'Explicit resident release failed'
  call check_counters() ! Counters are lifetime totals, not current allocation counts.
  call backend%release(status)
  if(status/=0)error stop 'Repeated resident release failed'
  call check_counters()
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,4,status)
  if(status/=0)error stop 'Prepare after release failed'
  expected(1:4)=expected(1:4)+1;call check_counters()
  call apply_resident(4)
  deallocate(backend) ! Prepared owner: finalization must release its device data and plan.
  call exx_cufft_create(backend)
  expected=0;call check_counters()
  call backend%prepare(plan%padded,plan%indices,plan%filter,source,4,status)
  if(status/=0)error stop 'Fresh owner after finalization failed'
  expected=[1,1,1,1,0];call check_counters()
  call apply_resident(4)
  deallocate(backend)
  call exx_local_destroy(plan)
  print *, 'PASS cuFFT/FFTW complex128 parity: one-shot, resident reuse, independent uploads, rebuilds and lifecycle'
contains
  subroutine check_counters()
    implicit none
    integer :: actual(5)
    select type(backend)
    type is(s_exx_cufft)
      actual=[backend%plan_builds,backend%filter_uploads,backend%index_uploads, &
        backend%source_uploads,backend%batch_uploads]
    class default
      error stop 'Factory returned the wrong GPU backend type'
    end select
    if(any(actual/=expected))then
      print *, 'Operation counters actual/expected: ',actual,expected
      error stop 'Resident operation counters differ'
    endif
  end subroutine

  subroutine apply_resident(capacity)
    implicit none
    integer,intent(in) :: capacity
    integer :: column,begin_column,end_column
    do column=1,nt
      density=conjg(source)*targets(:,column)
      call exx_local_apply(plan,density,potential,status)
      if(status/=0)error stop 'Resident FFTW oracle failed'
      reference(:,column)=-source*potential
    enddo
    action=cmplx(999d0,-999d0,8)
    do begin_column=1,nt,capacity
      end_column=min(nt,begin_column+capacity-1)
      call backend%apply(targets(:,begin_column:end_column),action(:,begin_column:end_column),status)
      if(status/=0)error stop 'Resident GPU batch failed'
      expected(5)=expected(5)+1
      call check_counters() ! apply may upload targets, never rebuild/reupload constant data.
    enddo
    scale=max(1d0,maxval(abs(reference)))
    if(maxval(abs(action-reference))/scale>2d-12)error stop 'Resident GPU/FFTW parity failed'
    if(maxval(abs(action(:,4)))>2d-14)error stop 'Resident zero target was contaminated'
  end subroutine
end program
