! Copyright 2026 SALMON developers. Apache-2.0.
! Extracted from fdtd_eh; keep arithmetic and loop order for reproducibility.
module legacy_ld
  implicit none
  private
  public :: s_legacy_ld, legacy_allocate, legacy_coefficients, legacy_store_old, legacy_update, legacy_finalize
  type s_legacy_ld
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
  end type s_legacy_ld
contains
  subroutine legacy_allocate(ld,lo,hi)
    type(s_legacy_ld), intent(inout) :: ld
    integer, intent(in) :: lo(3),hi(3)
    allocate(ld%rjx_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjx_ld=0.0d0
    allocate(ld%rjx_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjx_old_ld=0.0d0
    allocate(ld%rjy_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjy_ld=0.0d0
    allocate(ld%rjy_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjy_old_ld=0.0d0
    allocate(ld%rjz_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjz_ld=0.0d0
    allocate(ld%rjz_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%rjz_old_ld=0.0d0
    allocate(ld%px_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%px_ld=0.0d0
    allocate(ld%py_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%py_ld=0.0d0
    allocate(ld%pz_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%pz_ld=0.0d0
    allocate(ld%c1_jx_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c1_jx_ld=0.0d0
    allocate(ld%c2_jx_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c2_jx_ld=0.0d0
    allocate(ld%c3_jx_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c3_jx_ld=0.0d0
    allocate(ld%c1_jy_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c1_jy_ld=0.0d0
    allocate(ld%c2_jy_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c2_jy_ld=0.0d0
    allocate(ld%c3_jy_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c3_jy_ld=0.0d0
    allocate(ld%c1_jz_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c1_jz_ld=0.0d0
    allocate(ld%c2_jz_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c2_jz_ld=0.0d0
    allocate(ld%c3_jz_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ld%max_pole_num_ld))
    ld%c3_jz_ld=0.0d0
    allocate(ld%rjx_fdtd_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%rjx_fdtd_ld=0.0d0
    allocate(ld%rjy_fdtd_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%rjy_fdtd_ld=0.0d0
    allocate(ld%rjz_fdtd_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%rjz_fdtd_ld=0.0d0
    allocate(ld%ex_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%ex_old_ld=0.0d0
    allocate(ld%ey_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%ey_old_ld=0.0d0
    allocate(ld%ez_old_ld(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    ld%ez_old_ld=0.0d0
  end subroutine legacy_allocate
  subroutine legacy_coefficients(media_ld,pole_num_ld,gamma_ld,omega_ld,f_ld,omega_p_ld,dt_em,pi, &
                             tmp_c0_ld,tmp_c1_ld,tmp_c2_ld,tmp_c3_ld,tmp_c4_ld,tmp_c5_ld)
    integer, intent(in) :: media_ld(:),pole_num_ld(0:)
    real(8), intent(in) :: gamma_ld(0:,:),omega_ld(0:,:),f_ld(0:,:),omega_p_ld(0:),dt_em,pi
    real(8), intent(inout) :: tmp_c0_ld(:,:),tmp_c1_ld(:,:),tmp_c2_ld(:,:),tmp_c3_ld(:,:)
    real(8), intent(inout) :: tmp_c4_ld(:),tmp_c5_ld(:)
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
        tmp_c4_ld(ii) = 1.0d0 + 2.0d0*pi*dt_em*tmp_c4_ld(ii)
        tmp_c5_ld(ii) = 1.0d0 - 2.0d0*pi*dt_em*tmp_c5_ld(ii)
      end do
      
  end subroutine legacy_coefficients
    subroutine legacy_store_old(ld,lo,hi,ista,iend,dt_em,ex_y,ex_z,ey_z,ey_x,ez_x,ez_y)
      implicit none
      integer, intent(in) :: lo(3),hi(3),ista(3),iend(3)
      type(s_legacy_ld), intent(inout) :: ld
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
    end subroutine legacy_store_old
    subroutine legacy_update(ld,lo,hi,ista,iend,dt_em,ex_y,ex_z,ey_z,ey_x,ez_x,ez_y)
      implicit none
      integer, intent(in) :: lo(3),hi(3),ista(3),iend(3)
      type(s_legacy_ld), intent(inout) :: ld
      real(8), intent(in) :: dt_em
      real(8), intent(in) :: ex_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ex_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ey_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
                            ez_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
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
    end subroutine legacy_update
  subroutine legacy_finalize(ld)
    type(s_legacy_ld), intent(inout) :: ld
    type(s_legacy_ld) :: empty
    ld=empty
  end subroutine legacy_finalize
end module legacy_ld
