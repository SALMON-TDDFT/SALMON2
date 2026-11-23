! Copyright 2026 SALMON developers. Apache-2.0.
! Extracted from fdtd_eh; keep arithmetic and loop order for reproducibility.
module lorentz_drude
  implicit none
  private
  public :: s_ld_state, ld_allocate, ld_coefficients, ld_store_old, ld_update, ld_finalize
  type s_ld_state
    logical :: enabled=.false.
    integer             :: num_ld=0                                !LD: number of LD media
    integer             :: max_pole_num_ld=0                       !LD: maximum of pole_num_ld
    integer,allocatable :: media_ld(:)                             !LD: imedia number for LD model(num_ld)
    real(8),allocatable ::   rjx_ld(:,:,:,:),rjx_old_ld(:,:,:,:),& !LD: poparization current density J
                             rjy_ld(:,:,:,:),rjy_old_ld(:,:,:,:),& !    and its stock
                             rjz_ld(:,:,:,:),rjz_old_ld(:,:,:,:)   !    (x,y,z,max_pole_num_ld)
    real(8),allocatable ::    px_ld(:,:,:,:),&                     !LD: poparization vector P
                              py_ld(:,:,:,:),&                     !    (x,y,z,max_pole_num_ld)
                              pz_ld(:,:,:,:)                       !
    real(8),allocatable :: c1_jx_ld(:,:,:,:),&                     !LD: coefficient for Jx
                           c2_jx_ld(:,:,:,:),&                     !    (x,y,z,max_pole_num_ld)
                           c3_jx_ld(:,:,:,:)                       !
    real(8),allocatable :: c1_jy_ld(:,:,:,:),&                     !LD: coefficient for Jy
                           c2_jy_ld(:,:,:,:),&                     !    (x,y,z,max_pole_num_ld)
                           c3_jy_ld(:,:,:,:)                       !
    real(8),allocatable :: c1_jz_ld(:,:,:,:),&                     !LD: coefficient for Jz
                           c2_jz_ld(:,:,:,:),&                     !    (x,y,z,max_pole_num_ld)
                           c3_jz_ld(:,:,:,:)                       !
    real(8),allocatable :: rjx_fdtd_ld(:,:,:),&                    !LD: J substituted into FDTD
                           rjy_fdtd_ld(:,:,:),&                    !    (x,y,z)
                           rjz_fdtd_ld(:,:,:)                      !
    real(8),allocatable :: ex_old_ld(:,:,:),&                      !LD: old E
                           ey_old_ld(:,:,:),&                      !    (x,y,z)
                           ez_old_ld(:,:,:)                        !
  end type s_ld_state
contains
  subroutine ld_allocate(ld,lo,hi,status &
#ifdef SALMONPY_RECOVERY_TESTING
                         ,failure_slot &
#endif
                         )
    type(s_ld_state), intent(inout) :: ld
    integer, intent(in) :: lo(3),hi(3)
    integer,optional,intent(out) :: status
    logical :: failed
#ifdef SALMONPY_RECOVERY_TESTING
    integer,optional,intent(in) :: failure_slot
    integer :: slot
    slot=0
#endif
    failed=.false.
    if(present(status))status=0
    call allocate4(ld%rjx_ld)
    if(failed)return
    ld%rjx_ld=0.0d0
    call allocate4(ld%rjx_old_ld)
    if(failed)return
    ld%rjx_old_ld=0.0d0
    call allocate4(ld%rjy_ld)
    if(failed)return
    ld%rjy_ld=0.0d0
    call allocate4(ld%rjy_old_ld)
    if(failed)return
    ld%rjy_old_ld=0.0d0
    call allocate4(ld%rjz_ld)
    if(failed)return
    ld%rjz_ld=0.0d0
    call allocate4(ld%rjz_old_ld)
    if(failed)return
    ld%rjz_old_ld=0.0d0
    call allocate4(ld%px_ld)
    if(failed)return
    ld%px_ld=0.0d0
    call allocate4(ld%py_ld)
    if(failed)return
    ld%py_ld=0.0d0
    call allocate4(ld%pz_ld)
    if(failed)return
    ld%pz_ld=0.0d0
    call allocate4(ld%c1_jx_ld)
    if(failed)return
    ld%c1_jx_ld=0.0d0
    call allocate4(ld%c2_jx_ld)
    if(failed)return
    ld%c2_jx_ld=0.0d0
    call allocate4(ld%c3_jx_ld)
    if(failed)return
    ld%c3_jx_ld=0.0d0
    call allocate4(ld%c1_jy_ld)
    if(failed)return
    ld%c1_jy_ld=0.0d0
    call allocate4(ld%c2_jy_ld)
    if(failed)return
    ld%c2_jy_ld=0.0d0
    call allocate4(ld%c3_jy_ld)
    if(failed)return
    ld%c3_jy_ld=0.0d0
    call allocate4(ld%c1_jz_ld)
    if(failed)return
    ld%c1_jz_ld=0.0d0
    call allocate4(ld%c2_jz_ld)
    if(failed)return
    ld%c2_jz_ld=0.0d0
    call allocate4(ld%c3_jz_ld)
    if(failed)return
    ld%c3_jz_ld=0.0d0
    call allocate3(ld%rjx_fdtd_ld)
    if(failed)return
    ld%rjx_fdtd_ld=0.0d0
    call allocate3(ld%rjy_fdtd_ld)
    if(failed)return
    ld%rjy_fdtd_ld=0.0d0
    call allocate3(ld%rjz_fdtd_ld)
    if(failed)return
    ld%rjz_fdtd_ld=0.0d0
    call allocate3(ld%ex_old_ld)
    if(failed)return
    ld%ex_old_ld=0.0d0
    call allocate3(ld%ey_old_ld)
    if(failed)return
    ld%ey_old_ld=0.0d0
    call allocate3(ld%ez_old_ld)
    if(failed)return
    ld%ez_old_ld=0.0d0
  contains
    ! Keep allocation and zero-initialization order; an optional status gives
    ! native owners a recoverable local return, without an MPI dependency here.
    subroutine allocate4(array)
      real(8),allocatable,intent(out) :: array(:,:,:,:)
      integer :: ios
      if(present(status))then
        ios=0
#ifdef SALMONPY_RECOVERY_TESTING
        slot=slot+1
        if(present(failure_slot))then
          if(slot==failure_slot)ios=1
        endif
#endif
        if(ios==0)allocate(array(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld),stat=ios)
        if(ios/=0)then
          status=1
          failed=.true.
        endif
      else
        allocate(array(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
      endif
    end subroutine allocate4
    subroutine allocate3(array)
      real(8),allocatable,intent(out) :: array(:,:,:)
      integer :: ios
      if(present(status))then
        ios=0
#ifdef SALMONPY_RECOVERY_TESTING
        slot=slot+1
        if(present(failure_slot))then
          if(slot==failure_slot)ios=1
        endif
#endif
        if(ios==0)allocate(array(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=ios)
        if(ios/=0)then
          status=1
          failed=.true.
        endif
      else
        allocate(array(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      endif
    end subroutine allocate3
  end subroutine ld_allocate
  subroutine ld_coefficients(media_ld,pole_num_ld,gamma_ld,omega_ld,f_ld,omega_p_ld,dt_em,pi, &
                             tmp_c0_ld,tmp_c1_ld,tmp_c2_ld,tmp_c3_ld,tmp_c4_ld,tmp_c5_ld,epsilon_infinity)
    integer, intent(in) :: media_ld(:),pole_num_ld(0:)
    real(8), intent(in) :: gamma_ld(0:,:),omega_ld(0:,:),f_ld(0:,:),omega_p_ld(0:),dt_em,pi
    real(8), intent(inout) :: tmp_c0_ld(:,:),tmp_c1_ld(:,:),tmp_c2_ld(:,:),tmp_c3_ld(:,:)
    real(8), intent(inout) :: tmp_c4_ld(:),tmp_c5_ld(:)
    real(8),optional,intent(in) :: epsilon_infinity(0:)
    real(8) :: background
    integer :: ii,ij
      do ii = 1,size(media_ld)
        do ij = 1,pole_num_ld(media_ld(ii))
          tmp_c0_ld(ij,ii) = 1.0d0                                          &
                             + ( gamma_ld(media_ld(ii),ij)*dt_em/2.0d0 ) &
                             + ((omega_ld(media_ld(ii),ij)*dt_em)**2.0d0)/4.0d0
          tmp_c1_ld(ij,ii) = 1.0d0                                          &
                             - ( gamma_ld(media_ld(ii),ij)*dt_em/2.0d0 ) &
                             - ((omega_ld(media_ld(ii),ij)*dt_em)**2.0d0)/4.0d0
          tmp_c1_ld(ij,ii) = tmp_c1_ld(ij,ii) / tmp_c0_ld(ij,ii)
          tmp_c2_ld(ij,ii) = f_ld(media_ld(ii),ij)*(omega_p_ld(media_ld(ii))**2.0d0)*dt_em &
                             /(8.0d0*pi*tmp_c0_ld(ij,ii))
          tmp_c3_ld(ij,ii) = (omega_ld(media_ld(ii),ij)**2.0d0)*dt_em/tmp_c0_ld(ij,ii)
          tmp_c4_ld(ii) = tmp_c4_ld(ii) + tmp_c2_ld(ij,ii)
          tmp_c5_ld(ii) = tmp_c5_ld(ii) + tmp_c2_ld(ij,ii)
        end do
        background=1.0d0
        if(present(epsilon_infinity))background=epsilon_infinity(media_ld(ii))
        tmp_c4_ld(ii) = background + 2.0d0*pi*dt_em*tmp_c4_ld(ii)
        tmp_c5_ld(ii) = background - 2.0d0*pi*dt_em*tmp_c5_ld(ii)
      end do
      
  end subroutine ld_coefficients
    subroutine ld_store_old(ld,lo,hi,ista,iend,dt_em,ex_y,ex_z,ey_z,ey_x,ez_x,ez_y)
      implicit none
      integer, intent(in) :: lo(3),hi(3),ista(3),iend(3)
      type(s_ld_state), intent(inout) :: ld
      real(8), intent(in) :: dt_em
      real(8), intent(in) :: ex_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ex_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      integer :: ii,ix,iy,iz
      
      if(ld%enabled) then
!$omp parallel
!$omp do private(ii,ix,iy,iz) collapse(3)
        do ii=1,ld%max_pole_num_ld
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          ld%rjx_old_ld(ix,iy,iz,ii) = ld%rjx_ld(ix,iy,iz,ii)
          ld%rjy_old_ld(ix,iy,iz,ii) = ld%rjy_ld(ix,iy,iz,ii)
          ld%rjz_old_ld(ix,iy,iz,ii) = ld%rjz_ld(ix,iy,iz,ii)
        end do
        end do
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          ld%ex_old_ld(ix,iy,iz) = ex_y(ix,iy,iz) + ex_z(ix,iy,iz);
          ld%ey_old_ld(ix,iy,iz) = ey_z(ix,iy,iz) + ey_x(ix,iy,iz);
          ld%ez_old_ld(ix,iy,iz) = ez_x(ix,iy,iz) + ez_y(ix,iy,iz);
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
      
      return
    end subroutine ld_store_old
    subroutine ld_update(ld,lo,hi,ista,iend,dt_em,ex_y,ex_z,ey_z,ey_x,ez_x,ez_y,cross_drive)
      implicit none
      integer, intent(in) :: lo(3),hi(3),ista(3),iend(3)
      type(s_ld_state), intent(inout) :: ld
      real(8), intent(in) :: dt_em
      real(8), intent(in) :: ex_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ex_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      real(8),optional,intent(in)::cross_drive(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld)
      integer :: ii,ix,iy,iz
      
      !initialize
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        ld%rjx_fdtd_ld(ix,iy,iz)=0.0d0; ld%rjy_fdtd_ld(ix,iy,iz)=0.0d0; ld%rjz_fdtd_ld(ix,iy,iz)=0.0d0;
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !update LD polarization current & vector and FDTD current
      do ii=1,ld%max_pole_num_ld
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          ld%rjx_ld(ix,iy,iz,ii) = ld%c1_jx_ld(ix,iy,iz,ii)*ld%rjx_old_ld(ix,iy,iz,ii)            &
                                  +ld%c2_jx_ld(ix,iy,iz,ii)*( ex_y(ix,iy,iz)+ex_z(ix,iy,iz) &
                                                             +ld%ex_old_ld(ix,iy,iz) )            &
                                  -ld%c3_jx_ld(ix,iy,iz,ii)*ld%px_ld(ix,iy,iz,ii)
          ld%rjy_ld(ix,iy,iz,ii) = ld%c1_jy_ld(ix,iy,iz,ii)*ld%rjy_old_ld(ix,iy,iz,ii)            &
                                  +ld%c2_jy_ld(ix,iy,iz,ii)*( ey_z(ix,iy,iz)+ey_x(ix,iy,iz) &
                                                             +ld%ey_old_ld(ix,iy,iz) )            &
                                  -ld%c3_jy_ld(ix,iy,iz,ii)*ld%py_ld(ix,iy,iz,ii)
          ld%rjz_ld(ix,iy,iz,ii) = ld%c1_jz_ld(ix,iy,iz,ii)*ld%rjz_old_ld(ix,iy,iz,ii)            &
                                  +ld%c2_jz_ld(ix,iy,iz,ii)*( ez_x(ix,iy,iz)+ez_y(ix,iy,iz) &
                                                             +ld%ez_old_ld(ix,iy,iz) )            &
                                  -ld%c3_jz_ld(ix,iy,iz,ii)*ld%pz_ld(ix,iy,iz,ii)
          if(present(cross_drive))then
            ld%rjx_ld(ix,iy,iz,ii)=ld%rjx_ld(ix,iy,iz,ii)+cross_drive(1,ix,iy,iz,ii)
            ld%rjy_ld(ix,iy,iz,ii)=ld%rjy_ld(ix,iy,iz,ii)+cross_drive(2,ix,iy,iz,ii)
            ld%rjz_ld(ix,iy,iz,ii)=ld%rjz_ld(ix,iy,iz,ii)+cross_drive(3,ix,iy,iz,ii)
          endif
          ld%px_ld(ix,iy,iz,ii)  = ld%px_ld(ix,iy,iz,ii) &
                                  +0.5d0*dt_em*( ld%rjx_ld(ix,iy,iz,ii) + ld%rjx_old_ld(ix,iy,iz,ii) )
          ld%py_ld(ix,iy,iz,ii)  = ld%py_ld(ix,iy,iz,ii) &
                                  +0.5d0*dt_em*( ld%rjy_ld(ix,iy,iz,ii) + ld%rjy_old_ld(ix,iy,iz,ii) )
          ld%pz_ld(ix,iy,iz,ii)  = ld%pz_ld(ix,iy,iz,ii) &
                                  +0.5d0*dt_em*( ld%rjz_ld(ix,iy,iz,ii) + ld%rjz_old_ld(ix,iy,iz,ii) )
          ld%rjx_fdtd_ld(ix,iy,iz)  = ld%rjx_fdtd_ld(ix,iy,iz)                                          &
                                     +0.5d0*( (1.0d0+ld%c1_jx_ld(ix,iy,iz,ii)) * ld%rjx_ld(ix,iy,iz,ii) &
                                                    -ld%c3_jx_ld(ix,iy,iz,ii)  * ld%px_ld(ix,iy,iz,ii) )
          ld%rjy_fdtd_ld(ix,iy,iz)  = ld%rjy_fdtd_ld(ix,iy,iz)                                          &
                                     +0.5d0*( (1.0d0+ld%c1_jy_ld(ix,iy,iz,ii)) * ld%rjy_ld(ix,iy,iz,ii) &
                                                    -ld%c3_jy_ld(ix,iy,iz,ii)  * ld%py_ld(ix,iy,iz,ii) )
          ld%rjz_fdtd_ld(ix,iy,iz)  = ld%rjz_fdtd_ld(ix,iy,iz)                                          &
                                     +0.5d0*( (1.0d0+ld%c1_jz_ld(ix,iy,iz,ii)) * ld%rjz_ld(ix,iy,iz,ii) &
                                                    -ld%c3_jz_ld(ix,iy,iz,ii)  * ld%pz_ld(ix,iy,iz,ii) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end do
      
      return
    end subroutine ld_update
  subroutine ld_finalize(ld)
    type(s_ld_state), intent(inout) :: ld
    type(s_ld_state) :: empty
    ld=empty
  end subroutine ld_finalize
end module lorentz_drude
