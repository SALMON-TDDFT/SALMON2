! Stateless serial ABI. Every invocation owns its work arrays; no saved state.
module sbe_c_api
 use iso_c_binding
 use sbe_types
 use sbe_velocity_kernel
 implicit none
contains
 subroutine salmonpy_sbe(nb,nk,ne,order,vnl,volume,dt,eigen,weight,p,rvnl,rho,aold,anew,advance,j) bind(C)
  integer(c_int), value :: nb,nk,ne,order,vnl,advance
  real(c_double), value :: volume,dt
  real(c_double), intent(in) :: eigen(nb,nk),weight(nk),aold(3),anew(3)
  complex(c_double_complex), intent(in) :: p(nb,nb,3,nk),rvnl(nb,nb,3,nk)
  complex(c_double_complex), intent(inout) :: rho(nb,nb,nk)
  real(c_double), intent(out) :: j(3)
  type(s_sbe_gs_info) :: gs
  type(s_sbe_bloch_solver) :: sbe
  complex(8) :: tmp(3)
  real(8) :: trace
  integer :: ik,ib,jb
  gs%nb=nb; gs%nk=nk; gs%ne=ne; gs%volume=volume
  gs%eigen=eigen; gs%kweight=weight; gs%p_tm_matrix=p; gs%rvnl_tm_matrix=rvnl
  allocate(gs%delta_omega(nb,nb,nk),gs%occup(nb,nk))
  gs%occup=0d0
  gs%occup(1:ne/2,:)=2d0
  do ik=1,nk
   do jb=1,nb
    do ib=1,nb
     gs%delta_omega(ib,jb,ik)=eigen(ib,ik)-eigen(jb,ik)
    end do
   end do
  end do
  sbe%nb=nb; sbe%nk=nk; sbe%ik_min=1; sbe%ik_max=nk
  sbe%flag_vnl_correction=vnl/=0; sbe%rho=rho
  if(advance/=0) call velocity_evolve(sbe,gs,0.5d0*(aold+anew),dt)
  call velocity_current_local(sbe,gs,anew,order,tmp)
  j=real(tmp)/sum(weight)/volume
  if(order<1) then
   trace=0d0
   do ik=1,nk
    do ib=1,nb
     trace=trace+real(sbe%rho(ib,ib,ik))*weight(ik)
    end do
   end do
   j=j+anew*(trace/sum(weight))/volume
  end if
  rho=sbe%rho
 end subroutine
end module
