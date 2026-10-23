! Synthetic representative tile: timings exclude setup, communication and contraction.
program benchmark
 use iso_c_binding
 use iso_fortran_env,only:int64
 use exx_spinor_wave_fft
 implicit none
 integer,parameter::n(3)=[12,12,84],mesh(3)=[4,2,1],b=16,cols=756,nk=8,ng=12096
 complex(c_double_complex),allocatable,target::rho(:,:,:,:)
 real(c_double),allocatable::kernel(:,:,:)
 integer::point(3,ng),shift(3,nk),first(b),i,j,k,w,err,threads,pass
 integer(int64)::t0,t1,rate
 real(8)::seconds(5),checksum
 type(wave_fft_plans)::plans
 allocate(rho(b,cols,nk,4),kernel(0:n(1)*mesh(1)-1,0:n(2)*mesh(2)-1,0:n(3)*mesh(3)-1))
 do i=1,ng
 point(:,i)=[mod(i-1,n(1)),mod((i-1)/n(1),n(2)),(i-1)/(n(1)*n(2))]
 enddo
 do i=1,b
 first(i)=i
 enddo
 do k=1,nk
 shift(:,k)=[mod(k-1,mesh(1)),mod((k-1)/mesh(1),mesh(2)),0]
 enddo
 kernel=1d0/real(nk,8)
 rho=cmplx(.25d0,.125d0,8)
 call wave_fft_create(plans,mesh,rho,err)
 if(err/=0)stop 1
 call wave_fft_execute(plans,rho,n,point,shift,kernel,first,1,1,b,cols,threads)
 call system_clock(count_rate=rate)
 do pass=1,5
 call system_clock(t0)
 do w=1,10
 call wave_fft_execute(plans,rho,n,point,shift,kernel,first,1,1,b,cols,threads)
 enddo
 call system_clock(t1)
 seconds(pass)=real(t1-t0,8)/real(rate,8)/10
 enddo
 checksum=sum(abs(rho))
 if(abs(checksum-real(size(rho),8)*abs(cmplx(.25d0,.125d0,8)))>1d-6)stop 2
 print *, 'WAVE_FFT_TILE_BENCH',threads,' RHO_BYTES ',size(rho)*16,' SECONDS ',seconds
 call wave_fft_destroy(plans)
end program
