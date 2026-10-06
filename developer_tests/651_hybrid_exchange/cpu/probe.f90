program probe
!$ use omp_lib, only: omp_get_thread_num
 use exx_local_fft, only: s_exx_local_fft,exx_local_prepare_compact,exx_local_cpu_prepare, &
                         exx_local_cpu_pair,exx_local_apply,exx_local_destroy
 implicit none
 type(s_exx_local_fft) :: plan
 complex(8) :: kernel(5,8,8),source(192),targets(192,7),actual(192,7),reference(192,7),potential(192)
 integer :: workers(4)=[1,4,2,1],iteration,w,j,i,status,worker,box(3),n(3),count
 logical :: used,failed
 n=[16,8,8];box=[3,8,8]
 kernel=(.01d0,.003d0)
 source=(.2d0,.01d0)
 do j=1,7
  do i=1,192
   targets(i,j)=cmplx(sin(.03d0*i*j),cos(.07d0*i-j),8)
  enddo
 enddo
 targets(:,3)=0d0
 do iteration=1,2
  if(iteration==2)box=[2,8,8]
  count=product(box)
  call exx_local_prepare_compact(plan,n,box,kernel(:2*box(1)-1,:,:),status)
  if(status/=0)error stop 'compact prepare'
  do j=1,7
   call exx_local_apply(plan,conjg(source(:count))*targets(:count,j),potential(:count),status)
   if(status/=0)error stop 'serial action'
   reference(:count,j)=-source(:count)*potential(:count)
  enddo
  do w=1,4
   call exx_local_cpu_prepare(plan,workers(w),status)
   if(status/=0)error stop 'worker prepare'
   failed=.false.
!$omp parallel do default(none) num_threads(workers(w)) &
!$omp shared(plan,source,targets,actual,count) private(j,worker,used) reduction(.or.:failed)
   do j=1,7
    worker=1
!$  worker=omp_get_thread_num()+1
    call exx_local_cpu_pair(plan,worker,source(:count),targets(:count,j),actual(:count,j),used)
    if(used.neqv.(j/=3))failed=.true.
   enddo
!$omp end parallel do
   if(failed)error stop 'zero density detection'
   if(maxval(abs(actual(:count,:)-reference(:count,:)))>1d-12)error stop 'worker action mismatch'
  enddo
 enddo
 call exx_local_destroy(plan)
 print *, 'PASS CPU pair workspaces, worker count and geometry changes, zero pairs'
end program
