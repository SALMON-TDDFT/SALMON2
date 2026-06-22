!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Momentum (kinetic velocity) matrix element for the GW module (spec-b1-5a).
!
!   vmat(1:3, n, m) = < n,k | -i nabla | m,k >    (Cartesian, a.u.)
!
! This is the kinetic part of the velocity operator v = -i nabla (+ k + [V_nl,r]).
! The +k term vanishes for the off-diagonal (n /= m) pairs the absorption head
! needs (<n,k|k|m,k> = k <u_nk|u_mk> = k delta_nm), so only the gradient is taken.
! The nonlocal pseudopotential commutator [V_nl, r] is a separate refinement
! (b1-5b); it is NOT included here.
!
! Mirrors the transition-matrix assembly in io/write.f90 (halo exchange ->
! calc_gradient_psi -> grid inner product * Hvol), kept independent (no ported
! external code; the gradient stencil and overlap exchange are SALMON's own).
!
! No proper nouns appear in this file per the project constraint.
!
module gw_velocity_sub
  implicit none
  private
  public :: calc_velocity_mtxel

contains

  subroutine calc_velocity_mtxel(system, info, mg, stencil, srg, spsi, ik, ispin, vmat)
    use structures,    only: s_dft_system, s_parallel_info, s_rgrid, s_stencil, &
                             s_sendrecv_grid, s_orbital
    use communication, only: comm_summation
    use sendrecv_grid, only: update_overlap_complex8
    use stencil_sub,   only: calc_gradient_psi
    implicit none
    type(s_dft_system),    intent(in)    :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_stencil),       intent(in)    :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_orbital),       intent(inout) :: spsi   ! inout: update_overlap fills the halo
    integer,               intent(in)    :: ik, ispin
    complex(8),            intent(out)   :: vmat(3, system%no, system%no)

    complex(8), allocatable :: gtpsi_l(:,:,:,:), gtpsi(:,:,:,:)
    complex(8), allocatable :: vmat_l(:,:,:)
    complex(8), parameter   :: zI = (0.0d0, 1.0d0)
    complex(8) :: wrk(3)
    integer    :: no, m, n, i, narray, im
    integer    :: is(3), ie(3)

    no = system%no
    im = 1
    is = mg%is
    ie = mg%ie

    allocate(gtpsi_l(3, mg%is_array(1):mg%ie_array(1), &
                        mg%is_array(2):mg%ie_array(2), &
                        mg%is_array(3):mg%ie_array(3)))
    allocate(gtpsi  (3, mg%is_array(1):mg%ie_array(1), &
                        mg%is_array(2):mg%ie_array(2), &
                        mg%is_array(3):mg%ie_array(3)))
    allocate(vmat_l(3, no, no))
    vmat_l(:,:,:) = (0.0d0, 0.0d0)

    narray = 3 * (mg%ie_array(1) - mg%is_array(1) + 1) &
               * (mg%ie_array(2) - mg%is_array(2) + 1) &
               * (mg%ie_array(3) - mg%is_array(3) + 1)

    ! halo exchange only when the real-space grid is distributed; on a single
    ! rank the idx/idy/idz index maps wrap periodically (no halo needed).
    if (info%if_divide_rspace) call update_overlap_complex8(srg, mg, spsi%zwf)

    do m = 1, no
      gtpsi_l(:,:,:,:) = (0.0d0, 0.0d0)
      call calc_gradient_psi(spsi%zwf(:,:,:,ispin,m,ik,im), gtpsi_l, &
           mg%is_array, mg%ie_array, mg%is, mg%ie, mg%idx, mg%idy, mg%idz, &
           stencil%coef_nab, system%rmatrix_B)
      call comm_summation(gtpsi_l, gtpsi, narray, info%icomm_o)

      do n = 1, no
        do i = 1, 3
          wrk(i) = sum( conjg(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,n,ik,im)) &
                            * gtpsi(i, is(1):ie(1), is(2):ie(2), is(3):ie(3)) )
        end do
        vmat_l(:,n,m) = - zI * wrk * system%Hvol
      end do
    end do

    call comm_summation(vmat_l, vmat, 3*no*no, info%icomm_r)

    deallocate(gtpsi_l, gtpsi, vmat_l)
  end subroutine calc_velocity_mtxel

end module gw_velocity_sub
