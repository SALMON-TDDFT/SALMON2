#ifdef SALMONPY_TIMER_STUBS
module misc_routines
 implicit none
 real(8)::test_clock=0d0
contains
 function get_wtime() result(t)
  real(8)::t
  t=test_clock
 end function
end module
module backup_routines
 implicit none
contains
 subroutine load_value(fd,value)
  integer,intent(in)::fd
  real(8),allocatable,intent(out)::value(:,:)
  stop 91
 end subroutine
 subroutine save_value(fd,value)
  integer,intent(in)::fd
  real(8),intent(in)::value(:,:)
  stop 92
 end subroutine
end module
#else
! Deterministic clock supplied by the driver; no sleep or wall-time tolerance.
program test_timer_context
 use timer
 use misc_routines,only:test_clock
 implicit none
 type(s_timer_context)::host,a,b
 call timer_initialize
 call timer_set(LOG_TOTAL,7d0)
 call timer_begin(LOG_CALC_RHO)
 call timer_thread_begin(LOG_CALC_HARTREE)
 test_clock=2d0
 call timer_context_swap(host)
 call timer_context_swap(a)
 call timer_initialize
 call timer_enable_sub
 call timer_set(LOG_TOTAL,11d0)
 call timer_begin(LOG_CALC_RHO)
 test_clock=5d0
 call timer_end(LOG_CALC_RHO)
 call timer_context_swap(a)
 call timer_context_swap(host)
 if(timer_get(LOG_TOTAL)/=7d0)stop 1
 if(timer_now(LOG_CALC_RHO)/=2d0)stop 2
 if(timer_context_get(a,LOG_TOTAL)/=11d0)stop 3
 if(timer_context_get(a,LOG_CALC_RHO)/=3d0)stop 4
 test_clock=6d0
 call timer_end(LOG_CALC_RHO)
 call timer_thread_end(LOG_CALC_HARTREE)
 if(timer_get(LOG_CALC_RHO)/=3d0)stop 5
 if(timer_thread_get(LOG_CALC_HARTREE,0)/=3d0)stop 6
 call timer_context_swap(host)
 call timer_context_swap(b)
 call timer_initialize
 call timer_set(LOG_TOTAL,19d0)
 call timer_context_swap(b)
 call timer_context_swap(host)
 if(timer_get(LOG_TOTAL)/=7d0)stop 7
 if(timer_context_get(a,LOG_TOTAL)/=11d0)stop 8
 if(timer_context_get(b,LOG_TOTAL)/=19d0)stop 9
 call timer_context_swap(host)
 call timer_context_swap(a)
 call timer_initialize
 if(timer_get_sub(LOG_TOTAL)/=11d0)stop 10
 call timer_context_swap(a)
 call timer_context_swap(host)
 print *, 'timer contexts: host/running/thread/subrange/sibling PASS'
end program
#endif
