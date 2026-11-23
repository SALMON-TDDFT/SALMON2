! Independent collinear pair FFT batches over existing storage, with serial plans.
module exx_collinear_fft
 use iso_c_binding
 !$ use omp_lib,only:omp_get_max_threads,omp_get_num_threads
 implicit none
 private
 include 'fftw3.f03'
 integer,parameter,public::FFTW_FORWARD_DIRECTION=FFTW_FORWARD,FFTW_BACKWARD_DIRECTION=FFTW_BACKWARD
 public::collinear_fft_plans,collinear_fft_create,collinear_fft_execute,collinear_fft_destroy
 type collinear_fft_plans
  type(c_ptr),allocatable::forward(:),backward(:)
  integer,allocatable::offset(:)
  integer::groups=0,group_stride=0
 end type
contains
 subroutine collinear_fft_create(plans,mesh,nk,batch,active,contiguous,buffer,ierr)
 type(collinear_fft_plans),intent(inout)::plans
 integer,intent(in)::mesh(3),nk,batch,active
 logical,intent(in)::contiguous
 complex(c_double_complex),contiguous,target,intent(inout)::buffer(:)
 integer,intent(out)::ierr
 integer::workers,w,first,last,start,dims(3),flags,stride,distance,nq
 call collinear_fft_destroy(plans)
 ierr=1
 if(any(mesh<1).or.nk<1.or.batch<1.or.active<0.or.active>batch)return
 nq=product(mesh)
 if(modulo(nk,nq)/=0.or.size(buffer)/nk<batch)return
 ierr=0
 if(active==0)return
 workers=1
 !$ workers=omp_get_max_threads()
 workers=max(1,min(workers,active))
 stride=batch;distance=1
 if(contiguous)then
 stride=1;distance=nk
 endif
 plans%groups=nk/nq;plans%group_stride=nq*stride
 allocate(plans%forward(workers),plans%backward(workers),plans%offset(workers))
 plans%forward=c_null_ptr;plans%backward=c_null_ptr
 dims=mesh(3:1:-1);flags=ior(FFTW_ESTIMATE,FFTW_UNALIGNED)
 ! Planning is deliberately outside application OpenMP. Only FFTW's serial
 ! library is needed; threads execute disjoint new-array transform batches.
 do w=1,workers
 first=((w-1)*active)/workers;last=(w*active)/workers
 start=1+first*distance;plans%offset(w)=start
 plans%forward(w)=fftw_plan_many_dft(3,dims,last-first,buffer(start:),dims,stride,distance, &
 buffer(start:),dims,stride,distance,FFTW_FORWARD,flags)
 plans%backward(w)=fftw_plan_many_dft(3,dims,last-first,buffer(start:),dims,stride,distance, &
 buffer(start:),dims,stride,distance,FFTW_BACKWARD,flags)
 if(.not.c_associated(plans%forward(w)).or..not.c_associated(plans%backward(w)))then
 call collinear_fft_destroy(plans)
 ierr=1
 return
 endif
 enddo
 end subroutine

 subroutine collinear_fft_execute(plans,buffer,direction,threads)
 type(collinear_fft_plans),intent(in)::plans
 complex(c_double_complex),contiguous,intent(inout)::buffer(:)
 integer,intent(in)::direction
 integer,intent(out)::threads
 integer::w,g,start
 threads=1
 if(.not.allocated(plans%forward))return
 !$omp parallel default(none) shared(plans,buffer,direction,threads) private(w,g,start)
 !$omp single
 !$ threads=omp_get_num_threads()
 !$omp end single
 !$omp do collapse(2) schedule(static)
 do g=0,plans%groups-1
 do w=1,size(plans%forward)
 start=plans%offset(w)+g*plans%group_stride
 if(direction==FFTW_FORWARD)then
 call fftw_execute_dft(plans%forward(w),buffer(start:),buffer(start:))
 else
 call fftw_execute_dft(plans%backward(w),buffer(start:),buffer(start:))
 endif
 enddo
 enddo
 !$omp end do
 !$omp end parallel
 end subroutine

 subroutine collinear_fft_destroy(plans)
 type(collinear_fft_plans),intent(inout)::plans
 integer::w
 if(allocated(plans%forward))then
 do w=1,size(plans%forward)
 if(c_associated(plans%forward(w)))call fftw_destroy_plan(plans%forward(w))
 if(c_associated(plans%backward(w)))call fftw_destroy_plan(plans%backward(w))
 enddo
 deallocate(plans%forward,plans%backward,plans%offset)
 endif
 plans%groups=0;plans%group_stride=0
 end subroutine
end module
