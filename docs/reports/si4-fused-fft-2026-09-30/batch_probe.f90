program batch_probe
  use iso_c_binding
  use omp_lib
  implicit none
  include 'fftw3.f03'
  complex(c_double_complex),allocatable :: data(:),initial(:),reference(:)
  type(c_ptr) :: plan,refplan
  integer :: n,batch,threads,mode,flag,dims(1),i,k,r,base,nlines
  integer :: sizes(3)=[12,16,64],batches(4)=[1,4,16,64],teams(3)=[1,2,4],ii,jj,tt
  real(8) :: tick,setup,elapsed,error
  nlines=2048
  do ii=1,size(sizes)
    n=sizes(ii);dims=n
    allocate(data(n*nlines),initial(n*nlines),reference(n*nlines))
    do i=1,size(data)
      initial(i)=cmplx(sin(.017d0*i),cos(.031d0*i),8)
    enddo
    refplan=fftw_plan_dft_1d(n,data,data,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    data=initial
    do k=0,nlines-1
      base=1+k*n
      call fftw_execute_dft(refplan,data(base:base+n-1),data(base:base+n-1))
    enddo
    reference=data
    call fftw_destroy_plan(refplan)
    do mode=1,2
      flag=FFTW_ESTIMATE
      if(mode==2)flag=FFTW_MEASURE
      do jj=1,size(batches)
        batch=batches(jj)
        tick=omp_get_wtime()
        plan=fftw_plan_many_dft(1,dims,batch,data,dims,1,n,data,dims,1,n, &
          FFTW_FORWARD,ior(flag,FFTW_UNALIGNED))
        if(.not.c_associated(plan))error stop 'plan failure'
        setup=omp_get_wtime()-tick
        do tt=1,size(teams)
          threads=teams(tt)
          call omp_set_num_threads(threads)
          ! Warm up the plan and team before timing repeated line batches.
          data=initial
!$omp parallel do default(none) private(k,base) shared(nlines,batch,n,plan,data)
          do k=0,nlines/batch-1
            base=1+k*n*batch
            call fftw_execute_dft(plan,data(base:base+n*batch-1),data(base:base+n*batch-1))
          enddo
!$omp end parallel do
          tick=omp_get_wtime()
          do r=1,5000
!$omp parallel do default(none) private(k,base) shared(nlines,batch,n,plan,data,initial)
            do k=0,nlines/batch-1
              base=1+k*n*batch
              data(base:base+n*batch-1)=initial(base:base+n*batch-1)
              call fftw_execute_dft(plan,data(base:base+n*batch-1),data(base:base+n*batch-1))
            enddo
!$omp end parallel do
          enddo
          elapsed=omp_get_wtime()-tick
          error=maxval(abs(data-reference))
          if(error>1d-11)error stop 'FFT reference mismatch'
          write(*,*) 'BATCH',n,batch,threads,mode,setup,elapsed,error
        enddo
        call fftw_destroy_plan(plan)
      enddo
    enddo
    deallocate(data,initial,reference)
  enddo
end program
