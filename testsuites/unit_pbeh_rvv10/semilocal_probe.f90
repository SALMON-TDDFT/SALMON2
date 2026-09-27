program probe
 use hse_semilocal
 implicit none
 real(8) :: r(4),s(4),e(4),v(4),w(4)
 integer :: i,status
 r=[1d-7,.001d0,.1d0,1d0];s=[1d-12,.0001d0,.03d0,.2d0]
 call pbeh_semilocal_evaluate(r,s,e,v,w,status)
 if(status/=0)error stop 'semilocal'
 do i=1,4
  write(*,*)e(i),v(i),w(i)
 enddo
end program
