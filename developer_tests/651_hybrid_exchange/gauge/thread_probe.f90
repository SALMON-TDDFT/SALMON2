#include "config.h"
! Test-only adapter records scope entry; the production backend still sets BLAS.
module exx_blas_threads
 use gauge_test_backend, only: backend_control=>exx_blas_thread_control
 implicit none
 private
 public :: exx_blas_thread_control,scope_entries,scope_active
 integer :: scope_entries=0
 logical :: scope_active=.false.
contains
 subroutine exx_blas_thread_control(requested,state)
  use iso_c_binding, only: c_int
  use omp_lib, only: omp_get_max_threads,omp_in_parallel
  implicit none
  integer,intent(in) :: requested
  integer,intent(inout) :: state(3)
#ifdef HAVE_EXX_OPENBLAS_THREADS
  integer :: previous
  interface
   function get_threads() bind(C,name='openblas_get_num_threads') result(n)
    import c_int
    implicit none
    integer(c_int) :: n
   end function
  end interface
  previous=int(get_threads())
  if(requested==0.and.state(1)/=0)previous=state(2)
#endif
  if(omp_in_parallel())error stop 'gauge scope within OMP'
  if(requested>0)then
   if(scope_active.or.requested/=omp_get_max_threads())error stop 'gauge scope entry'
   scope_active=.true.;scope_entries=scope_entries+1
  endif
  call backend_control(requested,state)
#ifdef HAVE_EXX_OPENBLAS_THREADS
  if(requested>0.and.state(1)/=0)then
   if(get_threads()/=requested)error stop 'BLAS count not applied'
  elseif(requested==0)then
   if(get_threads()/=previous)error stop 'BLAS count not restored'
  endif
#endif
  if(requested==0)scope_active=.false.
 end subroutine
end module
