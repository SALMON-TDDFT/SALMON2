program scope_probe
 use omp_lib, only: omp_get_max_threads,omp_set_num_threads,omp_get_thread_num
 use exx_blas_threads, only: exx_blas_thread_control
 implicit none
 integer :: state(3),worker_state(3),before,errors,t
 call omp_set_num_threads(4)
 before=omp_get_max_threads()
 state=0
 call exx_blas_thread_control(1,state)
 if(state(1)==0.or.omp_get_max_threads()/=1)error stop 'scope entry'
 errors=0
!$omp parallel num_threads(3) private(worker_state,t) reduction(+:errors)
 t=2+omp_get_thread_num()
 call omp_set_num_threads(t)
 worker_state=0
 call exx_blas_thread_control(1,worker_state)
 if(omp_get_max_threads()/=1)errors=errors+1
 call exx_blas_thread_control(0,worker_state)
 if(omp_get_max_threads()/=t)errors=errors+1
!$omp end parallel
 call exx_blas_thread_control(0,state)
 if(errors/=0.or.omp_get_max_threads()/=before)error stop 'scope restoration'
 print *, 'PASS OpenMP worker and parent setting restoration'
end program
