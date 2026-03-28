!
!  SALMON2 stress tensor implementation
!  Nielsen-Martin formula: sigma_{ab} = -(1/V) dE/d(eps_{ab})
!
module stress_sub
  implicit none
  private
  public :: calc_stress, prepare_stress_field_state, refresh_stress_output_state, s_stress_field_state

  type s_stress_field_state
    logical :: use_micro_ac = .false.
    real(8) :: Ac_uniform(3) = 0d0
  end type s_stress_field_state

contains

  subroutine prepare_stress_field_state(system, theory, field_state)
    use structures, only: s_dft_system
    implicit none
    type(s_dft_system),         intent(in)  :: system
    character(*),               intent(in)  :: theory
    type(s_stress_field_state), intent(out) :: field_state

    if(trim(theory) == 'single_scale_maxwell_tddft') then
      field_state%use_micro_ac = .true.
      field_state%Ac_uniform = 0d0
    else
      field_state%use_micro_ac = .false.
      field_state%Ac_uniform = system%vec_Ac
    end if
  end subroutine prepare_stress_field_state

  subroutine refresh_stress_output_state(system, theory, info, mg, stencil, srg, ppg, &
                                         energy, V_local, live_psi, work_psi1, work_psi2, field_state)
    use structures
    use hamiltonian,  only: update_kvector_nonlocalpt, update_kvector_nonlocalpt_microAc
    use Total_Energy, only: calc_eigen_energy
    implicit none
    type(s_dft_system),         intent(inout) :: system
    character(*),               intent(in)    :: theory
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_sendrecv_grid),      intent(inout) :: srg
    type(s_pp_grid),            intent(inout) :: ppg
    type(s_dft_energy),         intent(inout) :: energy
    type(s_scalar),             intent(in)    :: V_local(system%nspin)
    type(s_orbital),            intent(inout) :: live_psi
    type(s_orbital),            intent(inout) :: work_psi1
    type(s_orbital),            intent(inout) :: work_psi2
    type(s_stress_field_state), intent(out)   :: field_state

    if(trim(theory) == 'single_scale_maxwell_tddft') then
      call update_kvector_nonlocalpt_microAc(info%ik_s, info%ik_e, system, ppg)
    else
      call update_kvector_nonlocalpt(info%ik_s, info%ik_e, system, ppg)
    end if
    call calc_eigen_energy(energy, live_psi, work_psi1, work_psi2, system, info, mg, V_local, stencil, srg, ppg)
    call prepare_stress_field_state(system, theory, field_state)
  end subroutine refresh_stress_output_state

  subroutine calc_stress(system, pp, fg, info, mg, stencil, poisson, &
                         srg, ppg, ppn, tpsi, ewald, energy, xc_func, rho_s, Vxc, field_state)
    use structures
    use plusU_global,  only: PLUS_U_ON
    use salmon_global, only: xc, yn_stress_loc_fd
    use sendrecv_grid, only: update_overlap_complex8
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_pp_info),            intent(in)    :: pp
    type(s_reciprocal_grid),    intent(in)    :: fg
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_poisson),            intent(in)    :: poisson
    type(s_sendrecv_grid),      intent(inout) :: srg
    type(s_pp_grid),            intent(in)    :: ppg
    type(s_pp_nlcc),            intent(in)    :: ppn
    type(s_orbital),            intent(inout) :: tpsi
    type(s_ewald_ion_ion),      intent(in)    :: ewald
    type(s_dft_energy),         intent(in)    :: energy
    type(s_xc_functional),      intent(in)    :: xc_func
    type(s_scalar),             intent(in)    :: rho_s(:)
    type(s_scalar),             intent(in)    :: Vxc(:)
    type(s_stress_field_state), intent(in)    :: field_state

    if(PLUS_U_ON) stop "stress tensor is not supported with PLUS_U_ON"
    if(system%nspin /= 1) stop "stress tensor supports only unpolarized calculations"
    if(trim(xc) /= 'pz') stop "stress tensor supports only built-in xc='PZ'"

    if(info%if_divide_rspace) then
      if(.not. tpsi%update_zwf_overlap) then
        call update_overlap_complex8(srg, mg, tpsi%zwf)
        tpsi%update_zwf_overlap = .true.
      end if
    end if

    system%stress_kin = 0d0
    system%stress_har = 0d0
    system%stress_xc = 0d0
    system%stress_loc = 0d0
    system%stress_loc_fd = 0d0
    system%stress_nl = 0d0
    system%stress_ewa = 0d0
    system%stress_tensor = 0d0

    call calc_stress_kin(system, info, mg, stencil, ppg, tpsi, field_state)
    call calc_stress_har(system, info, mg, fg, poisson, energy)
    call calc_stress_xc(system, info, mg, ppn, rho_s, Vxc, energy, xc_func)
    call calc_stress_loc(system, pp, fg, info, mg, ppg, poisson, energy)
    call calc_stress_nl(system, info, mg, stencil, ppg, tpsi, energy, field_state)
    call calc_stress_ewa(system, pp, fg, info, mg, ewald)

    if(yn_stress_loc_fd == 'y') then
      call calc_stress_loc_fd(system, pp, fg, info, mg, ppg, poisson, energy)
    end if

    call symmetrize_stress_term(system%stress_kin)
    call symmetrize_stress_term(system%stress_har)
    call symmetrize_stress_term(system%stress_xc)
    call symmetrize_stress_term(system%stress_loc)
    if(yn_stress_loc_fd == 'y') call symmetrize_stress_term(system%stress_loc_fd)
    call symmetrize_stress_term(system%stress_nl)
    call symmetrize_stress_term(system%stress_ewa)

    system%stress_tensor = system%stress_kin + system%stress_har + system%stress_xc &
                         + system%stress_loc + system%stress_nl + system%stress_ewa
    call symmetrize_stress_term(system%stress_tensor)
  end subroutine calc_stress

  pure subroutine symmetrize_stress_term(strs)
    real(8), intent(inout) :: strs(3,3)

    strs = 0.5d0 * (strs + transpose(strs))
  end subroutine symmetrize_stress_term

  subroutine get_g_bounds(fg, ig_s, ig_e)
    use structures, only: s_reciprocal_grid
    implicit none
    type(s_reciprocal_grid), intent(in)  :: fg
    integer,                 intent(out) :: ig_s(3), ig_e(3)

    ig_s(1) = lbound(fg%if_Gzero, 1)
    ig_s(2) = lbound(fg%if_Gzero, 2)
    ig_s(3) = lbound(fg%if_Gzero, 3)
    ig_e(1) = ubound(fg%if_Gzero, 1)
    ig_e(2) = ubound(fg%if_Gzero, 2)
    ig_e(3) = ubound(fg%if_Gzero, 3)
  end subroutine get_g_bounds

  subroutine sum_stress_tensor(icomm, strs, strs_sum)
    use communication, only: comm_summation
    implicit none
    integer, intent(in)  :: icomm
    real(8), intent(in)  :: strs(3,3)
    real(8), intent(out) :: strs_sum(3,3)

    call comm_summation(strs, strs_sum, 9, icomm)
  end subroutine sum_stress_tensor

  subroutine calc_stress_har(system, info, mg, fg, poisson, energy)
    use structures
    use math_constants, only: pi
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_poisson),         intent(in)    :: poisson
    type(s_dft_energy),      intent(in)    :: energy
    integer :: ix, iy, iz, a, b, ig_s(3), ig_e(3)
    real(8) :: g(3), G2, coeff, strs(3,3), strs_sum(3,3), V

    V = system%det_a
    strs = 0d0
    call get_g_bounds(fg, ig_s, ig_e)

    do iz = ig_s(3), ig_e(3)
    do iy = ig_s(2), ig_e(2)
    do ix = ig_s(1), ig_e(1)
      if(fg%if_Gzero(ix,iy,iz)) cycle
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      G2 = g(1)**2 + g(2)**2 + g(3)**2
      coeff = abs(poisson%zrhoG_ele(ix,iy,iz))**2 * (4d0*pi / G2**2)
      do b = 1, 3
      do a = 1, 3
        strs(a,b) = strs(a,b) - coeff * g(a) * g(b)
      end do
      end do
    end do
    end do
    end do

    call sum_stress_tensor(info%icomm_r, strs, strs_sum)
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_h / V
    end do
    system%stress_har = strs_sum
  end subroutine calc_stress_har

  subroutine calc_stress_xc(system, info, mg, ppn, rho_s, Vxc, energy, xc_func)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_nlcc),       intent(in)    :: ppn
    type(s_scalar),        intent(in)    :: rho_s(:)
    type(s_scalar),        intent(in)    :: Vxc(:)
    type(s_dft_energy),    intent(in)    :: energy
    type(s_xc_functional), intent(in)    :: xc_func
    integer :: ix, iy, iz, ispin, a
    real(8) :: rho_xc, E_vxc_loc, E_vxc, V

    if(xc_func%use_gradient) stop "stress tensor supports only built-in PZ LDA XC"

    V = system%det_a
    E_vxc_loc = 0d0
    do ispin = 1, system%nspin
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_xc = rho_s(ispin)%f(ix,iy,iz)
        if(allocated(ppn%rho_nlcc)) rho_xc = rho_xc + 0.5d0 * ppn%rho_nlcc(ix,iy,iz)
        E_vxc_loc = E_vxc_loc + rho_xc * Vxc(ispin)%f(ix,iy,iz)
      end do
      end do
      end do
    end do
    E_vxc_loc = E_vxc_loc * system%Hvol
    call comm_summation(E_vxc_loc, E_vxc, info%icomm_r)

    system%stress_xc = 0d0
    do a = 1, 3
      system%stress_xc(a,a) = (E_vxc - energy%E_xc) / V
    end do
  end subroutine calc_stress_xc

  subroutine calc_stress_kin(system, info, mg, stencil, ppg, tpsi, field_state)
    use structures
    use math_constants, only: zi
    use stencil_sub,    only: calc_gradient_psi
    use communication,  only: comm_summation
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_pp_grid),            intent(in)    :: ppg
    type(s_orbital),            intent(inout) :: tpsi
    type(s_stress_field_state), intent(in)    :: field_state
    integer :: ix, iy, iz, ik, io, ispin, im, a, b
    real(8) :: rtmp, kAc(3), strs(3,3), strs_sum(3,3), V
    complex(8) :: w(3), psi_r
    complex(8), allocatable :: gtpsi(:,:,:,:)

    V = system%det_a
    im = 1
    strs = 0d0
    allocate(gtpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))

    do ik = info%ik_s, info%ik_e
    do io = info%io_s, info%io_e
    do ispin = 1, system%nspin
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im), gtpsi, &
           mg%is_array, mg%ie_array, mg%is, mg%ie, mg%idx, mg%idy, mg%idz, stencil%coef_nab, system%rmatrix_B)
      rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol
      !$omp parallel do collapse(2) private(ix,iy,iz,kAc,psi_r,w,a,b) reduction(+:strs)
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        if(field_state%use_micro_ac) then
          kAc(1) = system%vec_k(1,ik) + system%Ac_micro%v(1,ix,iy,iz)
          kAc(2) = system%vec_k(2,ik) + system%Ac_micro%v(2,ix,iy,iz)
          kAc(3) = system%vec_k(3,ik) + system%Ac_micro%v(3,ix,iy,iz)
        else
          kAc(1) = system%vec_k(1,ik) + field_state%Ac_uniform(1)
          kAc(2) = system%vec_k(2,ik) + field_state%Ac_uniform(2)
          kAc(3) = system%vec_k(3,ik) + field_state%Ac_uniform(3)
        end if
        psi_r = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        w(1) = gtpsi(1,ix,iy,iz) + zi * kAc(1) * psi_r
        w(2) = gtpsi(2,ix,iy,iz) + zi * kAc(2) * psi_r
        w(3) = gtpsi(3,ix,iy,iz) + zi * kAc(3) * psi_r
        do b = 1, 3
        do a = 1, 3
          strs(a,b) = strs(a,b) - rtmp * dble(conjg(w(a)) * w(b))
        end do
        end do
      end do
      end do
      end do
      !$omp end parallel do
    end do
    end do
    end do

    deallocate(gtpsi)
    call comm_summation(strs, strs_sum, 9, info%icomm_rko)
    system%stress_kin = strs_sum / V
  end subroutine calc_stress_kin

  subroutine calc_stress_loc(system, pp, fg, info, mg, ppg, poisson, energy)
    use structures
    use math_constants, only: pi, zi
    use communication,  only: comm_summation
    use salmon_global,  only: cutoff_g, kion
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_pp_grid),         intent(in)    :: ppg
    type(s_poisson),         intent(in)    :: poisson
    type(s_dft_energy),      intent(in)    :: energy
    integer :: ix, iy, iz, ia, a, b, ig_s(3), ig_e(3)
    real(8) :: g(3), G2, Gd, coeff_lr, strs(3,3), strs_sum(3,3), E_sr, E_lr, E_sr_loc, E_lr_loc, V
    complex(8) :: rho_e, V_sr_sum, V_lr_sum, dVsr_dG2_sum, phase

    V = system%det_a
    strs = 0d0
    E_sr_loc = 0d0
    E_lr_loc = 0d0
    call get_g_bounds(fg, ig_s, ig_e)

    do iz = ig_s(3), ig_e(3)
    do iy = ig_s(2), ig_e(2)
    do ix = ig_s(1), ig_e(1)
      if(fg%if_Gzero(ix,iy,iz)) cycle
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      G2 = g(1)**2 + g(2)**2 + g(3)**2
      if(G2 > cutoff_g**2) cycle

      rho_e = poisson%zrhoG_ele(ix,iy,iz)
      V_sr_sum = (0d0, 0d0)
      V_lr_sum = (0d0, 0d0)
      dVsr_dG2_sum = (0d0, 0d0)

      do ia = 1, system%nion
        Gd = g(1)*system%Rion(1,ia) + g(2)*system%Rion(2,ia) + g(3)*system%Rion(3,ia)
        phase = exp(-zi*Gd)
        V_sr_sum = V_sr_sum + ppg%zVG_ion(ix,iy,iz,kion(ia)) * phase
        V_lr_sum = V_lr_sum - (4d0*pi / G2) * pp%zps(kion(ia)) * phase
        dVsr_dG2_sum = dVsr_dG2_sum + ppg%dVG_ion_dG2(ix,iy,iz,kion(ia)) * phase
      end do

      E_sr_loc = E_sr_loc + dble(conjg(rho_e) * V_sr_sum)
      E_lr_loc = E_lr_loc + dble(conjg(rho_e) * V_lr_sum)
      coeff_lr = 2d0 * dble(conjg(rho_e) * V_lr_sum) / G2

      do b = 1, 3
      do a = 1, 3
        strs(a,b) = strs(a,b) - 2d0 * dble(conjg(rho_e) * dVsr_dG2_sum) * g(a) * g(b)
        strs(a,b) = strs(a,b) + coeff_lr * g(a) * g(b)
      end do
      end do
    end do
    end do
    end do

    call comm_summation(E_sr_loc, E_sr, info%icomm_r)
    call comm_summation(E_lr_loc, E_lr, info%icomm_r)
    call sum_stress_tensor(info%icomm_r, strs, strs_sum)

    strs_sum = strs_sum / V**2
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + (E_sr + E_lr) / V
    end do
    system%stress_loc = strs_sum
  end subroutine calc_stress_loc

  subroutine calc_stress_nl(system, info, mg, stencil, ppg, tpsi, energy, field_state)
    use structures
    use math_constants,     only: zi
    use stencil_sub,        only: calc_gradient_psi
    use nonlocal_potential, only: calc_uVpsi, calc_uVpsi_rdivided
    use communication,      only: comm_summation
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_pp_grid),            intent(in)    :: ppg
    type(s_orbital),            intent(inout) :: tpsi
    type(s_dft_energy),         intent(in)    :: energy
    type(s_stress_field_state), intent(in)    :: field_state
    integer :: ix, iy, iz, ik, io, ispin, im, ilma, ia, j, a, b
    real(8) :: rtmp, kAc(3), strs(3,3), strs_sum(3,3), V
    complex(8) :: psi_r, w(3), r_uVpsi_b(3), uVpsi_ilma
    complex(8), allocatable :: gtpsi(:,:,:,:), uVpsibox(:,:,:,:,:), uVpsibox2(:,:,:,:,:)

    V = system%det_a
    im = 1
    strs = 0d0
    allocate(gtpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))

    if(info%if_divide_rspace) then
      call calc_uVpsi_rdivided(system%nspin, info, ppg, tpsi, uVpsibox, uVpsibox2)
    else
      call calc_uVpsi(system%nspin, info, ppg, tpsi, uVpsibox2)
    end if

    do ik = info%ik_s, info%ik_e
    do io = info%io_s, info%io_e
    do ispin = 1, system%nspin
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im), gtpsi, &
           mg%is_array, mg%ie_array, mg%is, mg%ie, mg%idx, mg%idy, mg%idz, stencil%coef_nab, system%rmatrix_B)
      rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol

      do ilma = 1, ppg%nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi_ilma = uVpsibox2(ispin, io, ik, im, ilma)
        do a = 1, 3
          r_uVpsi_b = (0d0, 0d0)
          ! Follow-up micro-optimization: kAc(1:3) is rebuilt inside the tensor-row loop
          ! even though only kAc(a) is consumed on each pass.
          do j = 1, ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            if(field_state%use_micro_ac) then
              kAc(1) = system%vec_k(1,ik) + system%Ac_micro%v(1,ix,iy,iz)
              kAc(2) = system%vec_k(2,ik) + system%Ac_micro%v(2,ix,iy,iz)
              kAc(3) = system%vec_k(3,ik) + system%Ac_micro%v(3,ix,iy,iz)
            else
              kAc(1) = system%vec_k(1,ik) + field_state%Ac_uniform(1)
              kAc(2) = system%vec_k(2,ik) + field_state%Ac_uniform(2)
              kAc(3) = system%vec_k(3,ik) + field_state%Ac_uniform(3)
            end if
            psi_r = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
            w(a) = gtpsi(a,ix,iy,iz) + zi * kAc(a) * psi_r
            do b = 1, 3
              r_uVpsi_b(b) = r_uVpsi_b(b) + ppg%rxyz(b,j,ia) * conjg(ppg%zekr_uV(j,ilma,ik)) * w(a)
            end do
          end do
          do b = 1, 3
            strs(a,b) = strs(a,b) + rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))
          end do
        end do
      end do
    end do
    end do
    end do

    deallocate(gtpsi)
    if(allocated(uVpsibox)) deallocate(uVpsibox)
    if(allocated(uVpsibox2)) deallocate(uVpsibox2)

    call comm_summation(strs, strs_sum, 9, info%icomm_rko)
    strs_sum = strs_sum / V
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_ion_nloc / V
    end do
    system%stress_nl = strs_sum
  end subroutine calc_stress_nl

  subroutine calc_stress_ewa(system, pp, fg, info, mg, ewald)
    use structures
    use math_constants, only: pi
    use salmon_math
    use communication, only: comm_summation
    use salmon_global, only: aEwald, cutoff_r, kion
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_ewald_ion_ion),   intent(in)    :: ewald
    integer :: ix, iy, iz, ia, ib, a, b, iia, ipair, ig_s(3), ig_e(3)
    real(8) :: g(3), G2, fact, rr, r_abs, rab(3), r(3), Qtot, strs_G(3,3), strs_R(3,3), strs_G_sum(3,3), strs_R_sum(3,3), V
    complex(8) :: SG

    V = system%det_a
    strs_G = 0d0
    strs_R = 0d0
    Qtot = 0d0
    do ia = 1, system%nion
      Qtot = Qtot + pp%zps(kion(ia))
    end do

    call get_g_bounds(fg, ig_s, ig_e)
    do iz = ig_s(3), ig_e(3)
    do iy = ig_s(2), ig_e(2)
    do ix = ig_s(1), ig_e(1)
      if(fg%if_Gzero(ix,iy,iz)) cycle
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      G2 = g(1)**2 + g(2)**2 + g(3)**2

      SG = (0d0, 0d0)
      do ia = 1, system%nion
        SG = SG + pp%zps(kion(ia)) * exp((0d0,1d0) * (g(1)*system%Rion(1,ia) + g(2)*system%Rion(2,ia) + g(3)*system%Rion(3,ia)))
      end do
      fact = (2d0*pi/G2) * exp(-G2/(4d0*aEwald)) / V**2 * abs(SG)**2
      do b = 1, 3
      do a = 1, 3
        strs_G(a,b) = strs_G(a,b) - fact * 2d0 * g(a) * g(b) / G2 * (1d0 + G2/(4d0*aEwald))
      end do
      end do
      do a = 1, 3
        strs_G(a,a) = strs_G(a,a) + fact
      end do
    end do
    end do
    end do

    call sum_stress_tensor(info%icomm_r, strs_G, strs_G_sum)
    do a = 1, 3
      strs_G_sum(a,a) = strs_G_sum(a,a) - pi * Qtot**2 / (2d0*aEwald*V**2)
    end do

    do iia = 1, info%nion_mg
      ia = info%ia_mg(iia)
      do ipair = 1, ewald%npair_bk(iia)
        r(1) = ewald%bk(1,ipair,iia)*system%primitive_a(1,1) &
             + ewald%bk(2,ipair,iia)*system%primitive_a(1,2) &
             + ewald%bk(3,ipair,iia)*system%primitive_a(1,3)
        r(2) = ewald%bk(1,ipair,iia)*system%primitive_a(2,1) &
             + ewald%bk(2,ipair,iia)*system%primitive_a(2,2) &
             + ewald%bk(3,ipair,iia)*system%primitive_a(2,3)
        r(3) = ewald%bk(1,ipair,iia)*system%primitive_a(3,1) &
             + ewald%bk(2,ipair,iia)*system%primitive_a(3,2) &
             + ewald%bk(3,ipair,iia)*system%primitive_a(3,3)
        ib = ewald%bk(4,ipair,iia)
        rab = system%Rion(:,ia) - r(:) - system%Rion(:,ib)
        rr = sum(rab(:)**2)
        if(rr > cutoff_r**2) cycle
        r_abs = sqrt(rr)
        fact = 0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib)) / (V * r_abs**3) &
             * (erfc_salmon(sqrt(aEwald)*r_abs) + 2d0*r_abs*sqrt(aEwald/pi)*exp(-aEwald*rr))
        do b = 1, 3
        do a = 1, 3
          strs_R(a,b) = strs_R(a,b) + fact * rab(a) * rab(b)
        end do
        end do
      end do
    end do

    call comm_summation(strs_R, strs_R_sum, 9, info%icomm_r)
    system%stress_ewa = strs_G_sum + strs_R_sum
  end subroutine calc_stress_ewa

  subroutine calc_stress_loc_fd(system, pp, fg, info, mg, ppg, poisson, energy)
    use structures
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_pp_grid),         intent(in)    :: ppg
    type(s_poisson),         intent(in)    :: poisson
    type(s_dft_energy),      intent(in)    :: energy

    stop "calc_stress_loc_fd pre-Task-10 guard reached; finish Task 10 or keep yn_stress_loc_fd='n'"
  end subroutine calc_stress_loc_fd

end module stress_sub
