program probe
!$ use omp_lib, only: omp_in_parallel,omp_set_num_threads
 use exx_ace, only: s_exx_ace,exx_ace_build,exx_ace_apply
 implicit none
 type(s_exx_ace) :: ace,empty_grid_ace
 complex(8) :: u(8,3,5),w(8,3,5),action(8,3,5),target(8,2,5),applied(8,2,5),reference(8,2,5)
 complex(8) :: empty_target(0,2,5),empty_action(0,2,5)
 integer :: nk,i,k,status,threads,previous,calls,requested_first,expected_workers
 threads=7;calls=0;requested_first=0;expected_workers=1
!$ expected_workers=3
!$ call omp_set_num_threads(3)
 u=0d0;w=0d0
 do k=1,5
  do i=1,3
   u(i,i,k)=1d0;w(i,i,k)=-real(i+k,8)
  enddo
 enddo
 w(:,:,5)=0d0
 do nk=1,5,4
  requested_first=0
  call exx_ace_build(ace,u(:,:,:nk),w(:,:,:nk),1d0,status,thread_control=control)
  if(status/=0.or.threads/=7)error stop 'build/restore'
  if(expected_workers>1)then
   if(nk==1.and.requested_first/=3)error stop 'single k BLAS policy'
   if(nk>1.and.requested_first/=1)error stop 'multiple k BLAS policy'
  else
   if(requested_first/=0)error stop 'serial BLAS setting changed'
  endif
  requested_first=0
  call exx_ace_apply(ace,u(:,:,:nk),action(:,:,:nk),status,thread_control=control)
  if(threads/=7)error stop 'apply restore'
  if(expected_workers>1)then
   if(nk==1.and.requested_first/=3)error stop 'apply single k policy'
   if(nk>1.and.requested_first/=1)error stop 'apply multiple k policy'
  endif
  if(status/=0.or.maxval(abs(action(:,:,:nk)-w(:,:,:nk)))>1d-12)error stop 'action'
 enddo
 target(:,1,:)=u(:,1,:)+(0d0,0.3d0)*u(:,2,:)
 target(:,2,:)=2d0*u(:,3,:)-u(:,2,:)
 call exx_ace_apply(ace,target,reference,status)
 if(status/=0)error stop 'apply reference'
 call exx_ace_apply(ace,target,applied,status,thread_control=control)
 if(status/=0.or.maxval(abs(applied-reference))>1d-12)error stop 'rectangular apply'
 calls=0
 call exx_ace_apply(ace,target,applied,status,reduce_grid,control)
 if(status/=0.or.calls/=7.or.threads/=7)error stop 'apply collective policy'
 if(maxval(abs(applied-reference))>1d-12)error stop 'apply collective result'
 call exx_ace_apply(ace,target(:,:0,:),applied(:,:0,:),status,thread_control=control)
 if(status==0.or.threads/=7)error stop 'invalid target state'
 allocate(empty_grid_ace%factors(0,3,5))
 empty_grid_ace%dv=1d0;calls=0
 call exx_ace_apply(empty_grid_ace,empty_target,empty_action,status,reduce_grid,control)
 if(status/=0.or.calls/=7.or.threads/=7)error stop 'empty grid collective path'
 calls=0
 call exx_ace_build(ace,u,w,1d0,status,reduce_grid,control)
 if(status/=0.or.calls/=10.or.threads/=7)error stop 'callback build'
 w=-w
 call exx_ace_build(ace,u,w,1d0,status,thread_control=control)
 if(status==0.or.allocated(ace%factors).or.threads/=7)error stop 'failure cleanup'
 print *, 'PASS k ACE, zero exchange, serial callback, failure and BLAS setting restore'
contains
 subroutine control(requested,saved)
  implicit none
  integer,intent(in) :: requested
  integer,intent(inout) :: saved(3)
  if(requested==0)then
   if(saved(1)/=0)threads=saved(2)
   saved=0
   return
  endif
  saved=0
!$ if(omp_in_parallel())return
  if(requested_first==0)requested_first=requested
  saved=[1,threads,0];threads=requested
 end subroutine
 subroutine reduce_grid(matrix)
  implicit none
  complex(8),intent(inout) :: matrix(:,:)
!$ if(omp_in_parallel())error stop 'collective in parallel region'
  calls=calls+1
 end subroutine
end program
