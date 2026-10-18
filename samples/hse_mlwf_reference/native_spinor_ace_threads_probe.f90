program probe
 use exx_ace
 implicit none
 type(s_exx_ace)::ace
 complex(8)::u(1024,8,3),w(1024,8,3),action(1024,8,3)
 integer::i,j,k,e
 ! SPD metric with a diagonal K and linearly independent Fourier vectors.
 do k=1,3;do j=1,8;do i=1,1024
 u(i,j,k)=exp(cmplx(0d0,real(i*j,8)*.004d0,8))/32
 w(i,j,k)=-u(i,j,k)*(.01d0*i+1)
 enddo;enddo;enddo
 call exx_ace_build(ace,u,w,1d0,e,thread_control=threads)
 if(e/=0)stop 1
 call exx_ace_apply(ace,u,action,e,thread_control=threads)
 if(e/=0.or.maxval(abs(action-w))>1d-8)stop 2
 print *, "ACE_THREADED_K_METRIC_APPLY_PASS"
contains
 subroutine threads(n,state)
 integer,intent(in)::n
 integer,intent(inout)::state(3)
 state=0
 if(n>0)state(1)=1
 end subroutine
end program
