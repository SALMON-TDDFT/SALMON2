program backend_probe
 use iso_c_binding, only: c_int
 use exx_ace, only: s_exx_ace,exx_ace_build,exx_ace_apply
 use exx_blas_threads, only: exx_blas_thread_control
 implicit none
 interface
  function get_threads() bind(C,name='openblas_get_num_threads') result(n)
   import c_int
   implicit none
   integer(c_int) :: n
  end function
 end interface
 type(s_exx_ace) :: ace
 complex(8) :: u(8,2,3),w(8,2,3),action(8,2,3)
 integer :: old,status,k,nk
 old=get_threads();u=0d0;w=0d0
 do k=1,3
  u(1,1,k)=1d0;u(2,2,k)=1d0
  w(1,1,k)=-2d0;w(2,2,k)=-3d0
 enddo
 do nk=1,3,2
  call exx_ace_build(ace,u(:,:,:nk),w(:,:,:nk),1d0,status,thread_control=exx_blas_thread_control)
  if(status/=0.or.get_threads()/=old)error stop 'real backend restoration'
  call exx_ace_apply(ace,u(:,:,:nk),action(:,:,:nk),status,thread_control=exx_blas_thread_control)
  if(get_threads()/=old)error stop 'apply real backend restoration'
  if(status/=0.or.maxval(abs(action(:,:,:nk)-w(:,:,:nk)))>1d-12)error stop 'backend action'
 enddo
 w=-w
 call exx_ace_build(ace,u,w,1d0,status,thread_control=exx_blas_thread_control)
 if(status==0.or.get_threads()/=old)error stop 'backend error restoration'
 print *, 'PASS real BLAS restoration and single/multiple-k action'
end program
