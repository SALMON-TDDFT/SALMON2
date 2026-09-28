program cufft_gpu_probe
  use exx_local_fft, only: s_exx_local_fft,exx_local_init,exx_local_prepare,exx_local_apply,exx_local_destroy
  use exx_cufft, only: exx_cufft_apply
  implicit none
  integer,parameter :: ns=7,nt=7
  type(s_exx_local_fft) :: plan
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
  call exx_local_destroy(plan)
  print *, 'PASS cuFFT/FFTW complex128 parity: 3x5x8, complex source, zero column, tails, repeated and empty batches'
end program
