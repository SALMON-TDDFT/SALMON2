#include "config.h"
! Scope state: active flag, previous BLAS/task count, previous LAPACK count.
! Request zero restores the scope on the same thread; no module state is kept.
module exx_blas_threads
  implicit none
  private
  public :: exx_blas_thread_control
contains
  subroutine exx_blas_thread_control(requested,state)
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS)
    use iso_c_binding, only: c_int
#endif
!$  use omp_lib, only: omp_in_parallel,omp_get_max_threads,omp_set_num_threads
    implicit none
    integer,intent(in) :: requested
    integer,intent(inout) :: state(3)
#ifdef HAVE_EXX_NVPL_THREADS
    integer(c_int) :: ignored
    interface
      function blas_local(n) bind(C,name='nvpl_blas_set_num_threads_local') result(previous)
        import c_int
        implicit none
        integer(c_int),value :: n
        integer(c_int) :: previous
      end function
      function lapack_local(n) bind(C,name='nvpl_lapack_set_num_threads_local') result(previous)
        import c_int
        implicit none
        integer(c_int),value :: n
        integer(c_int) :: previous
      end function
    end interface
#elif defined(HAVE_EXX_OPENBLAS_THREADS)
    interface
      function get_threads() bind(C,name='openblas_get_num_threads') result(n)
        import c_int
        implicit none
        integer(c_int) :: n
      end function
      function get_parallel() bind(C,name='openblas_get_parallel') result(mode)
        import c_int
        implicit none
        integer(c_int) :: mode
      end function
      subroutine set_threads(n) bind(C,name='openblas_set_num_threads')
        import c_int
        implicit none
        integer(c_int),value :: n
      end subroutine
    end interface
#endif
    if(requested==0)then
      if(state(1)==0)return
#ifdef HAVE_EXX_NVPL_THREADS
      ignored=blas_local(int(state(2),c_int))
      ignored=lapack_local(int(state(3),c_int))
#elif defined(HAVE_EXX_OPENBLAS_THREADS)
      call set_threads(int(state(2),c_int))
#elif defined(HAVE_EXX_FUJITSU_THREADS)
!$    call omp_set_num_threads(state(2))
#endif
      state=0
      return
    endif
    state=0
    if(requested<1)return
#ifdef HAVE_EXX_NVPL_THREADS
    ! Zero is a valid prior value: follow the library's global/default setting.
    state(2)=int(blas_local(int(requested,c_int)))
    state(3)=int(lapack_local(int(requested,c_int)))
    state(1)=1
#elif defined(HAVE_EXX_OPENBLAS_THREADS)
    ! Process-wide control belongs to the parent only, never to concurrent workers.
!$  if(omp_in_parallel())return
    ! A sequential OpenBLAS build may lack concurrent-call locking.
    if(get_parallel()==0)return
    state(2)=int(get_threads());state(1)=1
    call set_threads(int(requested,c_int))
#elif defined(HAVE_EXX_FUJITSU_THREADS)
    ! SSL2BLAMP follows the calling OpenMP task's thread-count setting.
!$  state(2)=omp_get_max_threads();state(1)=1
!$  call omp_set_num_threads(requested)
#endif
  end subroutine
end module
