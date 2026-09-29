! API test double only: vendor ABI and runtime still require MIYABI validation.
module nvpl_test_state
 use iso_c_binding, only: c_int
 implicit none
 integer(c_int) :: blas=0,lapack=5
 integer :: calls=0
!$omp threadprivate(blas,lapack,calls)
contains
 function blas_local(n) bind(C,name='nvpl_blas_set_num_threads_local') result(previous)
  implicit none
  integer(c_int),value :: n
  integer(c_int) :: previous
  previous=blas;blas=n;calls=calls+1
 end function
 function lapack_local(n) bind(C,name='nvpl_lapack_set_num_threads_local') result(previous)
  implicit none
  integer(c_int),value :: n
  integer(c_int) :: previous
  previous=lapack;lapack=n
 end function
end module
program nvpl_probe
 use omp_lib, only: omp_set_dynamic,omp_set_num_threads
 use nvpl_test_state, only: blas,lapack,calls
 use exx_ace, only: s_exx_ace,exx_ace_build,exx_ace_apply
 use exx_blas_threads, only: exx_blas_thread_control
 implicit none
 type(s_exx_ace) :: ace
 complex(8) :: u(8,2,6),w(8,2,6),action(8,2,6)
 integer :: status,k,nk,errors,state(3),outer(3)
 call omp_set_dynamic(.false.)
 call omp_set_num_threads(3)
 u=0d0;w=0d0
 do k=1,6
  u(1,1,k)=1d0;u(2,2,k)=1d0
  w(1,1,k)=-2d0;w(2,2,k)=-3d0
 enddo
 state=0;outer=0
 call exx_blas_thread_control(4,outer)
 call exx_blas_thread_control(1,state)
 if(blas/=1.or.lapack/=1)error stop 'NVPL nested enter'
 call exx_blas_thread_control(0,state)
 if(blas/=4.or.lapack/=4)error stop 'NVPL nested restore'
 call exx_blas_thread_control(0,outer)
 if(blas/=0.or.lapack/=5)error stop 'NVPL independent defaults'
 do nk=1,6,5
!$omp parallel
 calls=0
!$omp end parallel
  call exx_ace_build(ace,u(:,:,:nk),w(:,:,:nk),1d0,status,thread_control=exx_blas_thread_control)
  if(status/=0.or.blas/=0.or.lapack/=5)error stop 'NVPL parent restoration'
  errors=0
!$omp parallel reduction(+:errors)
  if(blas/=0.or.lapack/=5)errors=errors+1
  if(nk>1.and.calls<2)errors=errors+1
  calls=0
!$omp end parallel
  if(errors/=0)error stop 'NVPL build worker enter/restore'
  call exx_ace_apply(ace,u(:,:,:nk),action(:,:,:nk),status,thread_control=exx_blas_thread_control)
  if(status/=0.or.maxval(abs(action(:,:,:nk)-w(:,:,:nk)))>1d-12)error stop 'NVPL action'
  errors=0
!$omp parallel reduction(+:errors)
  if(blas/=0.or.lapack/=5)errors=errors+1
  if(nk>1.and.calls<2)errors=errors+1
!$omp end parallel
  if(errors/=0)error stop 'NVPL worker enter/restore'
 enddo
 w=-w
 call exx_ace_build(ace,u,w,1d0,status,thread_control=exx_blas_thread_control)
 if(status==0.or.allocated(ace%factors))error stop 'NVPL failure cleanup'
 errors=0
!$omp parallel reduction(+:errors)
 if(blas/=0.or.lapack/=5)errors=errors+1
!$omp end parallel
 if(errors/=0)error stop 'NVPL failure restoration'
 print *, 'PASS NVPL API double: per-worker scopes, independent defaults, action and failure restoration'
end program
