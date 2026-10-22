module crps_output
 implicit none
 private
 public :: crps_write_v1
contains

! Versioned CRPS interchange. Density is rho(r), NOT 4*pi*r^2*rho(r).
subroutine crps_write_v1(path,n,h,lmax,zval,nc,ipc,lloc,ls,js,rc,occ,eps,rho,vl,v,u)
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 character(*),intent(in)::path
 integer,intent(in)::n,lmax,nc,ipc,lloc,ls(nc),js(nc),rc(nc)
 real(8),intent(in)::h,zval,occ(nc),eps(nc),rho(0:n),vl(0:n),v(0:n,nc),u(0:n,nc)
 real(8)::d(0:n),dd(0:n)
 integer::i,k,unit
 if(n<4.or.h<=0.or.ipc<0.or.ipc>=n.or.nc/=2*lmax+1) error stop 'invalid CRPS dimensions'
 if(.not.all(ieee_is_finite(rho)).or.any(rho<0d0)) error stop 'invalid CRPS core density'
 if(.not.all(ieee_is_finite(v)).or..not.all(ieee_is_finite(u)).or.&
    .not.all(ieee_is_finite(vl))) error stop 'invalid CRPS radial data'
 if(ipc==0.and.any(rho/=0d0)) error stop 'core-off density must be zero'
 ! Five point derivatives; origin follows even radial density symmetry.
 d=0d0; dd=0d0
 do i=2,n-2
  d(i)=(rho(i-2)-8*rho(i-1)+8*rho(i+1)-rho(i+2))/(12*h)
  dd(i)=(-rho(i-2)+16*rho(i-1)-30*rho(i)+16*rho(i+1)-rho(i+2))/(12*h*h)
 end do
 dd(0)=(-30*rho(0)+32*rho(1)-2*rho(2))/(12*h*h)
 d(1)=(rho(1)-8*rho(0)+8*rho(2)-rho(3))/(12*h)
 dd(1)=(-rho(1)+16*rho(0)-30*rho(1)+16*rho(2)-rho(3))/(12*h*h)
 d(n-1)=(rho(n)-rho(n-2))/(2*h); d(n)=(3*rho(n)-4*rho(n-1)+rho(n-2))/(2*h)
 dd(n-1)=(rho(n)-2*rho(n-1)+rho(n-2))/(h*h);dd(n)=dd(n-1)
 open(newunit=unit,file=path,status='new',action='write')
 write(unit,'(i8,es25.16,i5,es25.16)')n,h,lmax,zval
 write(unit,'(a)')'KY_SOC_NLCC_V1'
 write(unit,'(a)')'ANGSTROM EV U_R_INV_SQRT_ANGSTROM RHO_INV_ANGSTROM3'
 write(unit,'(3i8)')ipc,lloc,nc
 do k=1,nc
  write(unit,'(2i5,3es25.16)')ls(k),js(k),rc(k)*h,occ(k),eps(k)
 end do
 do i=0,n
  write(unit,'(*(es25.16,1x))')i*h,vl(i),rho(i),d(i),dd(i),(v(i,k),u(i,k),k=1,nc)
 end do
 write(unit,'(a)')'END_KY_SOC_NLCC_V1'
 close(unit)
end subroutine crps_write_v1

end module crps_output
