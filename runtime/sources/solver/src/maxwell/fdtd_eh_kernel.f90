! Copyright 2026 SALMON developers. Apache-2.0.
module fdtd_eh_kernel
  implicit none
  private
  public :: eh_fd
contains
  subroutine eh_fd(ista,iend,ng_is,ng_ie,Nd,c1,c2,f1,f2,f3,var,dir,&
                   cb_x,cb_y,cb_z)
    implicit none
    integer,     intent(in)          :: ista(3),iend(3),ng_is(3),ng_ie(3)
    integer,     intent(in)          :: Nd
    real(8),     intent(in)          :: c1(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd),&
                                        c2(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    real(8),     intent(inout)       :: f1(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    real(8),     intent(in)          :: f2(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd),&
                                        f3(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    character(1),intent(in)          :: var,dir
    real(8),     intent(in),optional :: cb_x(ng_is(1)-Nd:ng_ie(1)+Nd), &
                                        cb_y(ng_is(2)-Nd:ng_ie(2)+Nd), &
                                        cb_z(ng_is(3)-Nd:ng_ie(3)+Nd)
    integer :: ix,iy,iz
    
    if(var=='e') then
      if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_x(ix-1)*(f2(ix-1,iy,iz)+f3(ix-1,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_y(iy-1)*(f2(ix,iy-1,iz)+f3(ix,iy-1,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_z(iz-1)*(f2(ix,iy,iz-1)+f3(ix,iy,iz-1)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(var=='h') then
      if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_x(ix+1)*(f2(ix+1,iy,iz)+f3(ix+1,iy,iz)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_y(iy+1)*(f2(ix,iy+1,iz)+f3(ix,iy+1,iz)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_z(iz+1)*(f2(ix,iy,iz+1)+f3(ix,iy,iz+1)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    end if
    
    return
  end subroutine eh_fd
end module fdtd_eh_kernel
