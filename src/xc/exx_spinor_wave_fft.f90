! Disjoint pair lanes of the existing SOC density tile; no thread-local tile copies.
module exx_spinor_wave_fft
 use iso_c_binding
 !$ use omp_lib,only:omp_get_max_threads,omp_get_num_threads
 implicit none
 private
 include 'fftw3.f03'
 public::wave_fft_plans,wave_fft_create,wave_fft_execute,wave_fft_destroy
 type wave_fft_plans
  type(c_ptr),allocatable::forward(:),backward(:)
  integer,allocatable::offset(:),count(:)
  integer::batch=0,nk=0,mesh(3)=0
 end type
contains
 subroutine wave_fft_create(plans,mesh,rho,ierr)
 type(wave_fft_plans),intent(inout)::plans
 integer,intent(in)::mesh(3)
 complex(c_double_complex),contiguous,target,intent(inout)::rho(:,:,:,:)
 integer,intent(out)::ierr
 complex(c_double_complex),pointer,contiguous::flat(:)
 integer::workers,w,dims(3),start,flags
 call wave_fft_destroy(plans)
 ierr=1
 if(any(mesh<1).or.size(rho,1)<1.or.size(rho,2)<1.or.(size(rho,4)/=1.and.size(rho,4)/=4))return
 if(product(mesh)/=size(rho,3))return
 plans%batch=size(rho,1)*size(rho,2);plans%nk=size(rho,3);plans%mesh=mesh
 workers=1
 !$ workers=omp_get_max_threads()
 workers=max(1,min(workers,plans%batch))
 allocate(plans%forward(workers),plans%backward(workers),plans%offset(workers),plans%count(workers))
 plans%forward=c_null_ptr;plans%backward=c_null_ptr
 call c_f_pointer(c_loc(rho),flat,[size(rho)])
 dims=mesh(3:1:-1);flags=ior(FFTW_ESTIMATE,FFTW_UNALIGNED)
 ! FFTW planning/destruction stays serial. Link the serial FFTW library.
 do w=1,workers
 plans%offset(w)=((w-1)*plans%batch)/workers+1
 plans%count(w)=(w*plans%batch)/workers-plans%offset(w)+1
 start=plans%offset(w)
 plans%forward(w)=fftw_plan_many_dft(3,dims,plans%count(w),flat(start:),dims,plans%batch,1, &
  flat(start:),dims,plans%batch,1,FFTW_FORWARD,flags)
 plans%backward(w)=fftw_plan_many_dft(3,dims,plans%count(w),flat(start:),dims,plans%batch,1, &
  flat(start:),dims,plans%batch,1,FFTW_BACKWARD,flags)
 if(.not.c_associated(plans%forward(w)).or..not.c_associated(plans%backward(w)))then
 call wave_fft_destroy(plans)
 return
 endif
 enddo
 ierr=0
 end subroutine

 subroutine wave_fft_execute(plans,rho,n,point,shift,kernel,first,lo,col,rows,ncols,threads)
 type(wave_fft_plans),intent(in)::plans
 complex(c_double_complex),contiguous,target,intent(inout)::rho(:,:,:,:)
 integer,intent(in)::n(3),point(:,:),shift(:,:),first(:),lo,col,rows,ncols
 real(c_double),intent(in)::kernel(0:,0:,0:)
 integer,intent(out)::threads
 complex(c_double_complex),pointer,contiguous::flat(:)
 integer::w,c,start,p,r,g,k,idx(3),b
 b=size(rho,1)
 call c_f_pointer(c_loc(rho),flat,[size(rho)])
 threads=1
 ! Each task owns distinct (spin,pair) lanes across every k, including padded
 ! lanes. The original tile and reduction order remain unchanged. No MPI here.
 !$omp parallel default(none) &
 !$omp shared(plans,rho,n,point,shift,kernel,first,lo,col,rows,ncols,threads,flat,b) &
 !$omp private(w,c,start,p,r,g,k,idx)
 !$omp single
 !$ threads=omp_get_num_threads()
 !$omp end single
 !$omp do collapse(2) schedule(static)
 do c=1,size(rho,4)
 do w=1,size(plans%forward)
 start=plans%offset(w)+(c-1)*plans%batch*plans%nk
 call fftw_execute_dft(plans%forward(w),flat(start:),flat(start:))
 do k=1,plans%nk
 do p=plans%offset(w),plans%offset(w)+plans%count(w)-1
 r=1+modulo(p-1,b);g=1+(p-1)/b
 if(r>rows.or.g>ncols)cycle
 idx=modulo(point(:,first(lo+r-1))-point(:,col+g-1)-shift(:,k),n*plans%mesh)
 rho(r,g,k,c)=rho(r,g,k,c)*kernel(idx(1),idx(2),idx(3))
 enddo
 enddo
 call fftw_execute_dft(plans%backward(w),flat(start:),flat(start:))
 enddo
 enddo
 !$omp end do
 !$omp end parallel
 end subroutine

 subroutine wave_fft_destroy(plans)
 type(wave_fft_plans),intent(inout)::plans
 integer::w
 if(allocated(plans%forward))then
 do w=1,size(plans%forward)
 if(c_associated(plans%forward(w)))call fftw_destroy_plan(plans%forward(w))
 if(c_associated(plans%backward(w)))call fftw_destroy_plan(plans%backward(w))
 enddo
 deallocate(plans%forward,plans%backward,plans%offset,plans%count)
 endif
 plans%batch=0;plans%nk=0;plans%mesh=0
 end subroutine
end module
