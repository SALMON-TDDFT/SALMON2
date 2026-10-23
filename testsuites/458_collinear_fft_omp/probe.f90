program probe
 use iso_c_binding
 use exx_collinear_fft
 !$ use omp_lib,only:omp_get_max_threads
 implicit none
 integer,parameter::batch=7,nk=12,nq=6
 complex(c_double_complex),target::data(batch*nk)
 complex(c_double_complex)::original(batch*nk),expected(batch*nk),v(nq)
 type(collinear_fft_plans)::plans
 integer::layout,active,pass,err,threads,want,p,k,j,g,ki(3),ji(3),stride,distance,offset,where
 real(8)::pi,angle,error
 logical::contiguous
 pi=acos(-1d0);error=0
 want=1
 !$ want=omp_get_max_threads()
 do layout=0,1
 contiguous=layout==1
 stride=batch;distance=1
 if(contiguous)then
 stride=1;distance=nk
 endif
 do p=1,size(data)
 original(p)=cmplx(sin(real(p,8)),cos(real(3*p+1,8)),8)
 enddo
 do pass=1,3
 active=batch
 if(pass==2)active=5
 if(pass==3)active=0
 expected=original
 do g=0,nk/nq-1;do p=1,active
 offset=(p-1)*distance+g*nq*stride
 v=0
 do k=0,nq-1;do j=0,nq-1
 ki=[mod(k,2),k/2,0];ji=[mod(j,2),j/2,0]
 angle=-2*pi*sum(real(ki*ji,8)/[2,3,1])
 v(k+1)=v(k+1)+original(offset+1+j*stride)*exp(cmplx(0d0,angle,8))
 enddo;enddo
 do k=0,nq-1
 expected(offset+1+k*stride)=v(k+1)
 enddo
 enddo;enddo
 data=original
 call collinear_fft_create(plans,[2,3,1],nk,batch,active,contiguous,data,err)
 if(err/=0)stop 1
 call collinear_fft_execute(plans,data,FFTW_FORWARD_DIRECTION,threads)
 if(active>0.and.threads/=want)stop 2
 error=max(error,maxval(abs(data-expected)))
 call collinear_fft_execute(plans,data,FFTW_BACKWARD_DIRECTION,threads)
 do g=0,nk/nq-1;do p=1,active;do k=0,nq-1
 where=1+(p-1)*distance+(g*nq+k)*stride
 expected(where)=original(where)*nq
 enddo;enddo;enddo
 if(active==0)expected=original
 error=max(error,maxval(abs(data-expected)))
 call collinear_fft_destroy(plans)
 ! Failed creation must leave no cached plan and allow a fresh valid creation.
 call collinear_fft_create(plans,[5,3,1],nk,batch,active,contiguous,data,err)
 if(err==0.or.allocated(plans%forward))stop 4
 enddo
 enddo
 if(error>1d-11)stop 3
 print *,'COLLINEAR_FFT_OMP_PASS',want,error
end program
