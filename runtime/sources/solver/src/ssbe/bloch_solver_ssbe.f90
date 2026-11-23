module bloch_solver_ssbe
    use sbe_types, only: s_sbe_bloch_solver
    use sbe_velocity_kernel, only: velocity_evolve, velocity_current_local
    use math_constants, only: pi, zi
    use communication, only: comm_get_groupinfo, comm_summation
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none







contains



subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))

    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do

    sbe%flag_vnl_correction = .false.
end subroutine


subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    use salmon_global, only: norder_correction
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(3)
    real(8), intent(out) :: jmat(3)
    integer, intent(in) :: icomm
    complex(8) :: tmp1(3), tmp(3)
    call velocity_current_local(sbe, gs, Ac, norder_correction, tmp1)
    call comm_summation(tmp1, tmp, 3, icomm)

    jmat(:) = real(tmp(1:3)) / sum(gs%kweight(:)) / gs%volume
    if (.not. (norder_correction >= 1)) then
        ! Gauge current cancellation
        jmat(:) = jmat(:) + Ac * calc_trace(sbe, gs, sbe%nb, icomm) / gs%volume
    end if

    return
end subroutine calc_current_bloch


subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: Ac(3), dt
    call velocity_evolve(sbe, gs, Ac, dt)
end subroutine


function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication
    use salmon_global
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    select case(trim(gauge_sbe))
    case ("velocity_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    case ("length_gauge")
        !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
        do ik = sbe%ik_min, sbe%ik_max
            do ib = 1, nb_max
                tmp1 = tmp1 + real(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
            end do
        end do
        !$omp end parallel do
    end select
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)

    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    ! real(8) :: kvec(1:3)
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                !& + dot_product(kvec(:), Ac(:))
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return
end function calc_energy

subroutine prepare_qnm(sbe, gs, icomm)
  use salmon_global, only: epdir_re1, am_s
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(in) :: gs
  integer, intent(in) :: icomm
  complex(8) :: dnm(sbe%nb, sbe%nb, sbe%nk)
  integer :: nb,nk
  integer :: ik,ib,jb,ii,jj

  nk = sbe%nk
  nb = sbe%nb

  allocate(sbe%qnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%grad_qnm(sbe%nb, sbe%nb, 1:3, sbe%nk))
  allocate(sbe%qnm_new(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%dqnm_stock(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max, am_s))
  allocate(sbe%dnm_i(sbe%nb, sbe%nb, 1:3, sbe%ik_min:sbe%ik_max))
  allocate(sbe%abs_dnm(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))
  allocate(sbe%exp_iphi(sbe%nb, sbe%nb, sbe%ik_min:sbe%ik_max))

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
    sbe%qnm(ib, jb, ik) = 0.d0
    sbe%qnm_new(ib, jb, ik) = 0.d0
    dnm(ib, jb, ik)=0.d0
    sbe%abs_dnm(ib, jb, ik)=0.d0
    sbe%exp_iphi(ib, jb, ik)=0.d0
  end do
  end do
  end do

  do ii=1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib=1,nb
  do jb=1,nb
      sbe%dqnm_stock(ib, jb, ik, ii) = 0.d0
  end do
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb, jj) 
  do ik = sbe%ik_min, sbe%ik_max
    do jj=1,3
      do ib=1,nb
        do jb=1,nb
          dnm(ib, jb, ik) = dnm(ib, jb, ik) + epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
        end do
      end do
    end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb, jj) collapse(4)
  do ik = sbe%ik_min, sbe%ik_max
    do jj=1,3
      do ib=1,nb
        do jb=1,nb
          sbe%dnm_i(ib, jb, jj, ik) = epdir_re1(jj) * gs%d_matrix(ib, jb, jj, ik)
        end do
      end do
    end do
  end do

  sbe%abs_dnm=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        sbe%abs_dnm(ib, jb, ik) = abs(dnm(ib, jb, ik))
      end do
    end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
    do ib=1,nb
      do jb=1,nb
        if(sbe%abs_dnm(ib, jb, ik) >= 1.d-13) then
          sbe%exp_iphi(ib, jb, ik) = dnm(ib, jb, ik)/sbe%abs_dnm(ib, jb, ik)
        end if
        if(ib==jb)then
          sbe%qnm_new(ib, jb, ik) = sbe%rho(ib, jb, ik)
        else
          sbe%qnm_new(ib, jb, ik) = conjg(sbe%exp_iphi(ib, jb, ik)) *  &
                                 &  sbe%rho(ib, jb, ik)
        end if
      end do
    end do
  end do

end subroutine prepare_qnm

subroutine dt_evolve_bloch_lg(sbe, gs, E, bj_am, dt, icomm)
  use salmon_global, only: am_s, num_kgrid, t_2, epdir_re1
  use common_ssbe, only: grad_k_array_nb2d_dcomplex
  implicit none
  type(s_sbe_bloch_solver), intent(inout) :: sbe
  type(s_sbe_gs_info), intent(inout) :: gs
  real(8), intent(in) :: E(1:3)
  real(8), intent(in) :: bj_am(8,8)
  real(8), intent(in) :: dt
  integer, intent(in) :: icomm
  integer :: nb,nk
  integer :: ik,ib,jb
  integer :: ii,lb
  complex(8) :: shift_vector(3)
  real(8) :: abs_E, proj_E
  complex(8) :: qnm_tmp1(sbe%nb, sbe%nb, sbe%nk)
  complex(8) :: qnm_tmp(sbe%nb, sbe%nb, sbe%nk)

  nb = sbe%nb 
  nk = sbe%nk

  shift_vector(:) = 0.d0 ! shift_vector is set to be 0

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm(ib, jb, ik) = sbe%qnm_new(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s-1
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, ii) = sbe%dqnm_stock(ib, jb, ik, ii+1)
  end do
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%dqnm_stock(ib, jb, ik, am_s) = 0.d0
  end do
  end do
  end do

  qnm_tmp1=0.d0
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    qnm_tmp1(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  call comm_summation(qnm_tmp1, qnm_tmp, nb*nb*nk, icomm)
  call grad_k_array_nb2d_dcomplex(nb,nk,gs%b_matrix,qnm_tmp,sbe%grad_qnm)

  abs_E = sqrt(E(1)**2+E(2)**2+E(3)**2)
  proj_E = E(1)*epdir_re1(1) + E(2)*epdir_re1(2) + E(3)*epdir_re1(3) 

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    if(ib == jb) then
    ! qnn (diagonal part)
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, ib, ik, am_s) = E(1) * sbe%grad_qnm(ib, ib, 1, ik) + &
                                           E(2) * sbe%grad_qnm(ib, ib, 2, ik) + &
                                           E(3) * sbe%grad_qnm(ib, ib, 3, ik)
        do lb = 1,nb
          if(ib /= lb)then
            sbe%dqnm_stock(ib, ib, ik, am_s) = &
              & sbe%dqnm_stock(ib, ib, ik, am_s) + &
              & 2.d0 * proj_E * sbe%abs_dnm(ib, lb, ik) * &
              &        aimag(sbe%qnm(lb, ib, ik))
          end if
        end do
      end if
    else
    ! qnm (off-diagonal part)
      sbe%dqnm_stock(ib, jb, ik, am_s) = - sbe%qnm(ib, jb, ik)/t_2
      if(abs_E >= 1.d-13) then
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & E(1) * sbe%grad_qnm(ib, jb, 1, ik) + &
          & E(2) * sbe%grad_qnm(ib, jb, 2, ik) + &
          & E(3) * sbe%grad_qnm(ib, jb, 3, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * gs%delta_omega(ib, jb, ik) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) - &
          & zi * (E(1) * shift_vector(1) + &
                  E(2) * shift_vector(2) + &
                  E(3) * shift_vector(3)) * sbe%qnm(ib, jb, ik)
        sbe%dqnm_stock(ib, jb, ik, am_s) = &
          & sbe%dqnm_stock(ib, jb, ik, am_s) + &
          & zi * proj_E * sbe%abs_dnm(ib, jb, ik) * &
          &     (sbe%qnm(ib, ib, ik) - sbe%qnm(jb, jb, ik))
        do lb = 1,nb
          if(lb /= ib .and. lb /= jb) then
            sbe%dqnm_stock(ib, jb, ik, am_s) = &
              & sbe%dqnm_stock(ib, jb, ik, am_s) - &
              & zi * proj_E * &
              &     (sbe%abs_dnm(ib, lb, ik) * &
              &      sbe%qnm(lb, jb, ik) - &
              &      sbe%abs_dnm(lb, jb, ik) * &
              &      sbe%qnm(ib, lb, ik)) * &
              &      sbe%exp_iphi(ib, lb, ik) * &
              &      sbe%exp_iphi(lb, jb, ik) * &
              &      conjg(sbe%exp_iphi(ib, jb, ik))
          end if
        end do
      end if
    end if
  end do
  end do
  end do

  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
    sbe%qnm_new(ib, jb, ik) = sbe%qnm(ib, jb, ik)
  end do
  end do
  end do

  do ii = 1,am_s
  !$omp parallel do default(shared) private(ik, ib, jb) collapse(3)
  do ik = sbe%ik_min, sbe%ik_max
  do ib = 1,nb
  do jb = 1,nb
      sbe%qnm_new(ib, jb, ik) = sbe%qnm_new(ib, jb, ik) + &
       &  bj_am(am_s+1-ii,am_s) * sbe%dqnm_stock(ib, jb, ik, ii) * dt
  end do
  end do
  end do
  end do

end subroutine dt_evolve_bloch_lg

subroutine calc_current_bloch_lg(sbe, gs, jmat, icomm)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(out) :: jmat(3)
    integer, intent(in) :: icomm
    integer :: ik, ib, jb, jj
    integer :: nk, nb
    complex(8) :: tmp1(3),tmp(3)

    nk = sbe%nk
    nb = sbe%nb

    tmp1(:) = 0.d0

    !$omp parallel do default(shared) private(ik,ib,jb,jj) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
    do jj = 1,3
    do ib = 1,nb
    do jb = 1,nb
        if(ib == jb) then
            tmp1(jj) = tmp1(jj) + gs%grad_k_eigen(ib, jj, ik) *  &
                               & dble(sbe%qnm_new(ib, ib, ik)) * gs%kweight(ik)
        else
            tmp1(jj) = tmp1(jj) + gs%delta_omega(ib, jb, ik) * &
                               & abs(sbe%dnm_i(ib, jb, jj, ik)) * &
                               & aimag(sbe%qnm_new(jb, ib, ik)) * gs%kweight(ik)
        end if
    end do
    end do
    end do
    end do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = -dble(tmp(:)) / sum(gs%kweight(:)) / gs%volume

end subroutine calc_current_bloch_lg

subroutine adams_moulton_coefs(bj_am)
  implicit none
  real(8) :: bj_am(8,8)

  bj_am(1,1) = 1.d0

  !bj_am(1,2) = 3.d0/2.d0
  !bj_am(2,2) = -1.d0/2.d0

  !bj_am(1,3) = 23.d0/12.d0
  !bj_am(2,3) = -4.d0/3.d0
  !bj_am(3,3) = 5.d0/12.d0

  bj_am(1,4) = 55.d0/24.d0
  bj_am(2,4) = -59.d0/24.d0
  bj_am(3,4) = 37.d0/24.d0
  bj_am(4,4) = -3.d0/8.d0

  !bj_am(1,5) = 1901.d0/720.d0
  !bj_am(2,5) = -1387.d0/360.d0
  !bj_am(3,5) = 109.d0/30.d0
  !bj_am(4,5) = -637.d0/360.d0
  !bj_am(5,5) = 251.d0/720.d0

  !bj_am(1,6) = 4277.d0/1440.d0
  !bj_am(2,6) = -2641.d0/480.d0
  !bj_am(3,6) = 4991.d0/720.d0
  !bj_am(4,6) = -3649.d0/720.d0
  !bj_am(5,6) = 959.d0/480.d0
  !bj_am(6,6) = -95.d0/288.d0

  !bj_am(1,7) = 198721.d0/60480.d0
  !bj_am(2,7) = -18637.d0/2520.d0
  !bj_am(3,7) = 235183.d0/20160.d0
  !bj_am(4,7) = -10754.d0/945.d0
  !bj_am(5,7) = 135713.d0/20160.d0
  !bj_am(6,7) = -5603.d0/2520.d0
  !bj_am(7,7) = 19087.d0/60480.d0

  bj_am(1,8) = 16083.d0/4480.d0
  bj_am(2,8) = -1152169.d0/120960.d0
  bj_am(3,8) = 242653.d0/13440.d0
  bj_am(4,8) = -296053.d0/13440.d0
  bj_am(5,8) = 2102243.d0/120960.d0
  bj_am(6,8) = -115747.d0/13440.d0
  bj_am(7,8) = 32863.d0/13440.d0
  bj_am(8,8) = -5257.d0/17280.d0

end subroutine adams_moulton_coefs


end module

