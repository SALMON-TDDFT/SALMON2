! Synthetic cached pair FFTs; no MPI, density generation or contraction timed.
program benchmark
 use iso_c_binding
 use iso_fortran_env,only:int64
 use exx_collinear_fft
 implicit none
 integer,parameter::batch=4096,nk=64
 complex(c_double_complex),allocatable,target::data(:)
 type(collinear_fft_plans)::plans
 integer::pass,w,i,err,threads
 integer(int64)::t0,t1,rate
 real(8)::times(5)
 character(16)::layout
 logical::contiguous
 call get_command_argument(1,layout)
 contiguous=trim(layout)=='contiguous'
 allocate(data(batch*nk));data=cmplx(.25d0,.125d0,8)
 call collinear_fft_create(plans,[4,4,4],nk,batch,batch,contiguous,data,err)
 if(err/=0)stop 1
 call step()
 call system_clock(count_rate=rate)
 do pass=1,5
 call system_clock(t0)
 do w=1,10
 call step()
 enddo
 call system_clock(t1)
 times(pass)=real(t1-t0,8)/real(rate,8)/10
 enddo
 if(maxval(abs(data-cmplx(.25d0,.125d0,8)))>1d-12)stop 2
 print *,'COLLINEAR_FFT_BENCH',threads,trim(layout),' BYTES ',size(data)*16,' SECONDS ',times
 call collinear_fft_destroy(plans)
contains
 subroutine step()
 call collinear_fft_execute(plans,data,FFTW_FORWARD_DIRECTION,threads)
 !$omp parallel do schedule(static)
 do i=1,size(data)
 data(i)=data(i)/nk
 enddo
 !$omp end parallel do
 call collinear_fft_execute(plans,data,FFTW_BACKWARD_DIRECTION,threads)
 end subroutine
end program
