!
!  SALMON2 stress tensor implementation
!  Nielsen-Martin formula: sigma_{ab} = +(1/V) dE/d(eps_{ab})
!
module stress_sub
  implicit none
  private
  public :: calc_stress, calc_stress_har_shadow, calc_stress_loc_sr_rs, &
    prepare_stress_field_state, refresh_stress_output_state, s_stress_field_state

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

    if(allocated(system%stress_nl_l)) deallocate(system%stress_nl_l)

    system%stress_kin = 0d0
    system%stress_har = 0d0
    system%stress_har_shadow = 0d0
    system%stress_xc = 0d0
    system%stress_loc = 0d0
    system%stress_loc_fd = 0d0
    system%stress_loc_grad = 0d0
    system%stress_loc_diag = 0d0
    system%stress_loc_fullobj_grad = 0d0
    system%stress_loc_fullobj_diag = 0d0
    system%stress_loc_sr_grad = 0d0
    system%stress_loc_lr_grad = 0d0
    system%stress_loc_sr_diag = 0d0
    system%stress_loc_lr_diag = 0d0
    system%stress_loc_sr_scr_grad = 0d0
    system%stress_loc_lr_scr_grad = 0d0
    system%stress_loc_sr_scr_diag = 0d0
    system%stress_loc_lr_scr_diag = 0d0
    system%stress_nl = 0d0
    system%stress_ewa = 0d0
    system%stress_loc_sr_energy = 0d0
    system%stress_loc_lr_energy = 0d0
    system%stress_loc_sr_scr_energy = 0d0
    system%stress_loc_lr_scr_energy = 0d0
    system%stress_xc_e_vxc = 0d0
    system%stress_ewa_g = 0d0
    system%stress_ewa_r = 0d0
    system%stress_ewa_g_grad = 0d0
    system%stress_ewa_g_diag = 0d0
    system%stress_ewa_g_self = 0d0
    system%stress_ewa_energy_G = 0d0
    system%stress_ewa_energy_R = 0d0
    system%stress_tensor = 0d0
    system%stress_kin_dbg_grad2 = 0d0
    system%stress_kin_dbg_cross = 0d0
    system%stress_kin_dbg_k2 = 0d0

    call calc_stress_kin(system, info, mg, stencil, ppg, tpsi, field_state)
    call calc_stress_har(system, info, mg, fg, poisson, energy)
    call calc_stress_xc(system, info, mg, ppn, rho_s, Vxc, energy, xc_func)
    call calc_stress_loc(system, pp, fg, info, mg, ppg, poisson, energy)
    call calc_stress_nl(system, pp, info, mg, stencil, ppg, tpsi, energy, field_state)
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
    system%stress_har = -strs_sum
  end subroutine calc_stress_har

  subroutine calc_stress_har_shadow(system, info, mg, stencil, srg_scalar, Vh, energy)
    use structures
    use math_constants, only: pi
    use sendrecv_grid, only: update_overlap_real8
    use stencil_sub, only: calc_gradient_field
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_stencil),       intent(in)    :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg_scalar
    type(s_scalar),        intent(in)    :: Vh
    type(s_dft_energy),    intent(in)    :: energy
    integer :: ix, iy, iz, a, b
    real(8) :: box(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
    real(8) :: grad_vh(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8) :: strs(3,3), strs_sum(3,3), coeff, V

    V = system%det_a
    coeff = system%Hvol / (4d0 * pi)
    strs = 0d0
    box = 0d0

    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      box(ix,iy,iz) = Vh%f(ix,iy,iz)
    end do
    end do
    end do

    if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, box)
    call calc_gradient_field(mg, stencil%coef_nab, system%rmatrix_B, box, grad_vh)

    ! Keep the diagnostic accumulation order deterministic across runs.
    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      do b = 1, 3
      do a = 1, 3
        strs(a,b) = strs(a,b) - coeff * grad_vh(a,ix,iy,iz) * grad_vh(b,ix,iy,iz)
      end do
      end do
    end do
    end do
    end do

    call sum_stress_tensor(info%icomm_r, strs, strs_sum)
    strs_sum = strs_sum / V
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_h / V
    end do
    call symmetrize_stress_term(strs_sum)
    system%stress_har_shadow = -strs_sum
  end subroutine calc_stress_har_shadow

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
    system%stress_xc_e_vxc = E_vxc
    do a = 1, 3
      system%stress_xc(a,a) = -(E_vxc - energy%E_xc) / V
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
    real(8) :: grad2_loc, cross_loc, k2_loc, grad2_sum, cross_sum, k2_sum, psi_abs2
    complex(8) :: w(3), psi_r
    complex(8), allocatable :: gtpsi(:,:,:,:)

    V = system%det_a
    im = 1
    strs = 0d0
    grad2_loc = 0d0
    cross_loc = 0d0
    k2_loc = 0d0
    allocate(gtpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))

    do ik = info%ik_s, info%ik_e
    do io = info%io_s, info%io_e
    do ispin = 1, system%nspin
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im), gtpsi, &
           mg%is_array, mg%ie_array, mg%is, mg%ie, mg%idx, mg%idy, mg%idz, stencil%coef_nab, system%rmatrix_B)
      rtmp = system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol
      !$omp parallel do collapse(2) private(ix,iy,iz,kAc,psi_r,w,a,b,psi_abs2) &
      !$omp& reduction(+:strs,grad2_loc,cross_loc,k2_loc)
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
        psi_abs2 = abs(psi_r)**2
        w(1) = gtpsi(1,ix,iy,iz) + zi * kAc(1) * psi_r
        w(2) = gtpsi(2,ix,iy,iz) + zi * kAc(2) * psi_r
        w(3) = gtpsi(3,ix,iy,iz) + zi * kAc(3) * psi_r
        do a = 1, 3
          grad2_loc = grad2_loc + rtmp * abs(gtpsi(a,ix,iy,iz))**2
          cross_loc = cross_loc + rtmp * dble(conjg(gtpsi(a,ix,iy,iz)) * (zi * kAc(a) * psi_r))
          k2_loc = k2_loc + rtmp * kAc(a)**2 * psi_abs2
        end do
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
    call comm_summation(grad2_loc, grad2_sum, info%icomm_rko)
    call comm_summation(cross_loc, cross_sum, info%icomm_rko)
    call comm_summation(k2_loc, k2_sum, info%icomm_rko)
    system%stress_kin = strs_sum / V
    system%stress_kin_dbg_grad2 = grad2_sum
    system%stress_kin_dbg_cross = cross_sum
    system%stress_kin_dbg_k2 = k2_sum
  end subroutine calc_stress_kin

  subroutine calc_stress_loc(system, pp, fg, info, mg, ppg, poisson, energy)
    use structures
    use math_constants, only: pi, zi
    use communication,  only: comm_summation
    use salmon_global,  only: aEwald, cutoff_g, kion
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_pp_grid),         intent(in)    :: ppg
    type(s_poisson),         intent(in)    :: poisson
    type(s_dft_energy),      intent(in)    :: energy
    integer :: ix, iy, iz, ia, ik, a, b, ig_s(3), ig_e(3)
    real(8) :: g(3), r(3), G2, Gd, coeff_lr, coeff_lr_scr, scr_fac, strs(3,3), strs_sum(3,3), &
      & strs_grad_sum(3,3), strs_diag_sum(3,3), strs_sr(3,3), strs_lr(3,3), strs_lr_scr(3,3), &
      & strs_sr_sum(3,3), strs_lr_sum(3,3), strs_lr_scr_sum(3,3), E_sr, E_lr, E_sr_scr, E_lr_scr, &
      & E_sr_loc, E_lr_loc, E_lr_scr_loc, V
    complex(8) :: rho_e, V_sr_sum, V_lr_sum, V_lr_scr_sum, dVsr_dG2_sum, phase

    V = system%det_a
    strs = 0d0
    strs_sr = 0d0
    strs_lr = 0d0
    strs_lr_scr = 0d0
    E_sr_loc = 0d0
    E_lr_loc = 0d0
    E_lr_scr_loc = 0d0
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
      V_lr_scr_sum = (0d0, 0d0)
      dVsr_dG2_sum = (0d0, 0d0)

      do ia = 1, system%nion
        ik = kion(ia)
        if(ik < 1 .or. ik > size(pp%zps)) stop "calc_stress_loc: invalid species index"
        if(.not. allocated(ppg%dVG_ion_dG2)) stop "calc_stress_loc: dVG_ion_dG2 is not allocated"
        r(:) = system%Rion(1:3,ia)
        Gd = sum(g(:) * r(:))
        phase = dcmplx(cos(Gd), -sin(Gd))
        V_sr_sum = V_sr_sum + ppg%zVG_ion(ix,iy,iz,ik) * phase
        V_lr_sum = V_lr_sum - (4d0*pi / G2) * pp%zps(ik) * phase
        dVsr_dG2_sum = dVsr_dG2_sum + ppg%dVG_ion_dG2(ix,iy,iz,ik) * phase
      end do

      E_sr_loc = E_sr_loc + dble(conjg(rho_e) * V_sr_sum)
      E_lr_loc = E_lr_loc + dble(conjg(rho_e) * V_lr_sum)
      scr_fac = fg%exp_ewald(ix,iy,iz)
      V_lr_scr_sum = scr_fac * V_lr_sum
      E_lr_scr_loc = E_lr_scr_loc + dble(conjg(rho_e) * V_lr_scr_sum)
      coeff_lr = -2d0 * dble(conjg(rho_e) * V_lr_sum) / G2
      coeff_lr_scr = -2d0 * dble(conjg(rho_e) * V_lr_scr_sum) / G2 * (1d0 + G2 / (4d0 * aEwald))

      do b = 1, 3
      do a = 1, 3
        strs_sr(a,b) = strs_sr(a,b) + 2d0 * dble(conjg(rho_e) * dVsr_dG2_sum) * g(a) * g(b)
        strs_lr(a,b) = strs_lr(a,b) + coeff_lr * g(a) * g(b)
        strs_lr_scr(a,b) = strs_lr_scr(a,b) + coeff_lr_scr * g(a) * g(b)
      end do
      end do
    end do
    end do
    end do

    call comm_summation(E_sr_loc, E_sr, info%icomm_r)
    call comm_summation(E_lr_loc, E_lr, info%icomm_r)
    call comm_summation(E_lr_scr_loc, E_lr_scr, info%icomm_r)
    strs = strs_sr + strs_lr
    call sum_stress_tensor(info%icomm_r, strs, strs_sum)
    call sum_stress_tensor(info%icomm_r, strs_sr, strs_sr_sum)
    call sum_stress_tensor(info%icomm_r, strs_lr, strs_lr_sum)
    call sum_stress_tensor(info%icomm_r, strs_lr_scr, strs_lr_scr_sum)

    ! poisson%zrhoG_ele already carries the reciprocal-space 1/V normalization,
    ! while dVG_ion_dG2 is stored as a bare pseudopotential derivative.
    strs_sum = strs_sum / V
    strs_sr_sum = strs_sr_sum / V
    strs_lr_sum = strs_lr_sum / V
    strs_lr_scr_sum = strs_lr_scr_sum / V
    strs_grad_sum = strs_sum
    strs_diag_sum = 0d0
    E_sr_scr = (E_sr + E_lr) - E_lr_scr
    system%stress_loc_grad = -strs_grad_sum
    system%stress_loc_diag = 0d0
    system%stress_loc_fullobj_grad = -strs_grad_sum
    system%stress_loc_fullobj_diag = 0d0
    system%stress_loc_sr_grad = -strs_sr_sum
    system%stress_loc_lr_grad = -strs_lr_sum
    system%stress_loc_sr_scr_grad = -(strs_grad_sum - strs_lr_scr_sum)
    system%stress_loc_lr_scr_grad = -strs_lr_scr_sum
    system%stress_loc_sr_diag = 0d0
    system%stress_loc_lr_diag = 0d0
    system%stress_loc_sr_scr_diag = 0d0
    system%stress_loc_lr_scr_diag = 0d0
    do a = 1, 3
      ! Keep the diagonal consistent with the periodic local-energy decomposition:
      ! E_ion_loc = E_sr + E_lr.
      strs_diag_sum(a,a) = (E_sr + E_lr) / V
      strs_sum(a,a) = strs_grad_sum(a,a) + strs_diag_sum(a,a)
      system%stress_loc_sr_diag(a,a) = -E_sr / V
      system%stress_loc_lr_diag(a,a) = -E_lr / V
      system%stress_loc_sr_scr_diag(a,a) = -E_sr_scr / V
      system%stress_loc_lr_scr_diag(a,a) = -E_lr_scr / V
      system%stress_loc_diag(a,a) = -strs_diag_sum(a,a)
      system%stress_loc_fullobj_diag(a,a) = -energy%E_ion_loc / V
    end do
    system%stress_loc_sr_energy = E_sr
    system%stress_loc_lr_energy = E_lr
    system%stress_loc_sr_scr_energy = E_sr_scr
    system%stress_loc_lr_scr_energy = E_lr_scr
    system%stress_loc = -strs_sum
  end subroutine calc_stress_loc

  subroutine calc_stress_nl(system, pp, info, mg, stencil, ppg, tpsi, energy, field_state)
    use structures
    use math_constants,     only: zi
    use stencil_sub,        only: calc_gradient_psi
    use nonlocal_potential, only: calc_uVpsi, calc_uVpsi_rdivided
    use communication,      only: comm_summation
    use salmon_global,      only: stress_fd_detail, kion
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_pp_info),            intent(in)    :: pp
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_pp_grid),            intent(in)    :: ppg
    type(s_orbital),            intent(inout) :: tpsi
    type(s_dft_energy),         intent(in)    :: energy
    type(s_stress_field_state), intent(in)    :: field_state
    integer :: ix, iy, iz, ik, io, ispin, im, ilma, ia, j, a, b, ll, lmax_nl
    real(8) :: rtmp, kAc(3), strs(3,3), strs_sum(3,3), V, nl_energy_part
    complex(8) :: psi_r, w(3), r_uVpsi_b(3), uVpsi_ilma
    complex(8), allocatable :: gtpsi(:,:,:,:), uVpsibox(:,:,:,:,:), uVpsibox2(:,:,:,:,:)
    integer, allocatable :: ll_of_ilma(:)
    real(8), allocatable :: strs_l(:,:,:), strs_l_sum(:,:,:), e_nl_l(:), e_nl_l_sum(:)
    logical :: want_l_detail

    V = system%det_a
    im = 1
    strs = 0d0
    want_l_detail = (trim(stress_fd_detail) == 'high')

    if(want_l_detail) then
      call build_nl_l_channel_map(pp, kion, ppg%nlma, ll_of_ilma, lmax_nl)
      allocate(strs_l(0:lmax_nl,3,3), strs_l_sum(0:lmax_nl,3,3))
      allocate(e_nl_l(0:lmax_nl), e_nl_l_sum(0:lmax_nl))
      strs_l = 0d0
      e_nl_l = 0d0
    end if

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
      rtmp = system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol

      do ilma = 1, ppg%nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi_ilma = uVpsibox2(ispin, io, ik, im, ilma)
        if(want_l_detail) then
          ll = ll_of_ilma(ilma)
          nl_energy_part = rtmp * dble(conjg(uVpsi_ilma) * uVpsi_ilma) / ppg%rinv_uvu(ilma)
          e_nl_l(ll) = e_nl_l(ll) + nl_energy_part
        end if
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
            strs(a,b) = strs(a,b) + 2d0 * rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))
            if(want_l_detail) strs_l(ll,a,b) = strs_l(ll,a,b) + 2d0 * rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))
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
    if(want_l_detail) then
      call comm_summation(strs_l, strs_l_sum, size(strs_l), info%icomm_rko)
      call comm_summation(e_nl_l, e_nl_l_sum, size(e_nl_l), info%icomm_rko)
    end if
    strs_sum = strs_sum / V
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_ion_nloc / V
    end do
    system%stress_nl = -strs_sum

    if(want_l_detail) then
      allocate(system%stress_nl_l(0:lmax_nl,3,3))
      system%stress_nl_l = -strs_l_sum / V
      do ll = 0, lmax_nl
        do a = 1, 3
          system%stress_nl_l(ll,a,a) = system%stress_nl_l(ll,a,a) - e_nl_l_sum(ll) / V
        end do
      end do
      deallocate(strs_l, strs_l_sum, e_nl_l, e_nl_l_sum, ll_of_ilma)
    end if
  end subroutine calc_stress_nl

  subroutine build_nl_l_channel_map(pp, kion, nlma, ll_of_ilma, lmax_nl)
    use structures
    implicit none
    type(s_pp_info),      intent(in)  :: pp
    integer,              intent(in)  :: kion(:)
    integer,              intent(in)  :: nlma
    integer, allocatable, intent(out) :: ll_of_ilma(:)
    integer,              intent(out) :: lmax_nl
    integer :: ia, ik, ll, l, l0, m, ilma

    lmax_nl = 0
    do ia = 1, size(kion)
      ik = kion(ia)
      lmax_nl = max(lmax_nl, pp%mlps(ik))
    end do

    allocate(ll_of_ilma(nlma))
    ll_of_ilma = -1

    ilma = 0
    do ia = 1, size(kion)
      ik = kion(ia)
      l0 = 0
      do ll = 0, pp%mlps(ik)
        do l = l0, l0 + pp%nproj(ll,ik) - 1
          if(pp%inorm(l,ik) == 0) cycle
          do m = -ll, ll
            ilma = ilma + 1
            ll_of_ilma(ilma) = ll
          end do
        end do
        l0 = l
      end do
    end do

    if(ilma /= nlma) stop "error: build_nl_l_channel_map size mismatch"
  end subroutine build_nl_l_channel_map

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
    real(8) :: g(3), G2, fact, rr, r_abs, rab(3), r(3), Qtot, zps2, strs_G(3,3), strs_R(3,3), strs_G_sum(3,3), strs_R_sum(3,3), V
    real(8) :: strs_G_grad(3,3), strs_G_diag(3,3), strs_G_self(3,3), strs_G_grad_sum(3,3), strs_G_diag_sum(3,3)
    real(8) :: E_ewa_G_loc, E_ewa_R_loc, E_ewa_G_sum, E_ewa_R_sum
    complex(8) :: SG

    V = system%det_a
    strs_G = 0d0
    strs_R = 0d0
    strs_G_grad = 0d0
    strs_G_diag = 0d0
    strs_G_self = 0d0
    E_ewa_G_loc = 0d0
    E_ewa_R_loc = 0d0
    Qtot = 0d0
    zps2 = 0d0
    do ia = 1, system%nion
      Qtot = Qtot + pp%zps(kion(ia))
      zps2 = zps2 + pp%zps(kion(ia))**2
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
      E_ewa_G_loc = E_ewa_G_loc + (2d0*pi/G2) * exp(-G2/(4d0*aEwald)) / V * abs(SG)**2
      do b = 1, 3
      do a = 1, 3
        strs_G_grad(a,b) = strs_G_grad(a,b) + fact * 2d0 * g(a) * g(b) / G2 * (1d0 + G2/(4d0*aEwald))
      end do
      end do
      do a = 1, 3
        strs_G_diag(a,a) = strs_G_diag(a,a) - fact
      end do
    end do
    end do
    end do

    strs_G = strs_G_grad + strs_G_diag
    call sum_stress_tensor(info%icomm_r, strs_G, strs_G_sum)
    call sum_stress_tensor(info%icomm_r, strs_G_grad, strs_G_grad_sum)
    call sum_stress_tensor(info%icomm_r, strs_G_diag, strs_G_diag_sum)
    do a = 1, 3
      strs_G_sum(a,a) = strs_G_sum(a,a) + pi * Qtot**2 / (2d0*aEwald*V**2)
      strs_G_self(a,a) = pi * Qtot**2 / (2d0*aEwald*V**2)
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
        E_ewa_R_loc = E_ewa_R_loc + 0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib)) &
                   * erfc_salmon(sqrt(aEwald)*r_abs) / r_abs
        do b = 1, 3
        do a = 1, 3
          strs_R(a,b) = strs_R(a,b) + fact * rab(a) * rab(b)
        end do
        end do
      end do
    end do

    call comm_summation(strs_R, strs_R_sum, 9, info%icomm_r)
    call comm_summation(E_ewa_G_loc, E_ewa_G_sum, info%icomm_r)
    call comm_summation(E_ewa_R_loc, E_ewa_R_sum, info%icomm_r)
    E_ewa_G_sum = E_ewa_G_sum - pi * Qtot**2 / (2d0*aEwald*V) - sqrt(aEwald/pi) * zps2
    system%stress_ewa_energy_G = E_ewa_G_sum
    system%stress_ewa_energy_R = E_ewa_R_sum
    system%stress_ewa_g_grad = strs_G_grad_sum
    system%stress_ewa_g_diag = strs_G_diag_sum
    system%stress_ewa_g_self = strs_G_self
    system%stress_ewa_g = strs_G_sum
    system%stress_ewa_r = strs_R_sum
    system%stress_ewa = strs_G_sum + strs_R_sum
  end subroutine calc_stress_ewa

  !----------------------------------------------------------------------------
  ! Real-space local-SR stress shadow (Nielsen-Martin Paper I, Eq. 30b)
  !
  ! Computes the short-range part of the ion-electron stress directly in real
  ! space, completely bypassing the G-space SR/LR decomposition and its
  ! associated numerical integration / cancellation issues.
  !
  ! Formula:
  !   sigma_{ab}^{sr} = -(1/V) sum_I sum_j rho(r_j) V'_sr(|s_j|)
  !                     * s_{j,a} s_{j,b} / |s_j| * hvol
  !
  ! where s_j = r_j - R_I,  V_sr(r) = V_loc(r) + Z/r,
  !       V'_sr(r) = dV_loc/dr - Z/r^2.
  !
  ! The LR Coulomb gradient and diagonal terms must be added from the existing
  ! G-space computation to form the total local stress.
  !----------------------------------------------------------------------------
  subroutine calc_stress_loc_sr_rs(system, pp, info, mg, ppg, rho_s)
    use structures
    use math_constants, only: pi
    use communication,  only: comm_summation
    use salmon_global,  only: kion
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_grid),       intent(in)    :: ppg
    type(s_scalar),        intent(in)    :: rho_s(:)
    integer  :: ia, ik, j, a, b, ix, iy, iz, intr, ispin, nspin
    real(8)  :: s(3), r_abs, r_inv, rho_val, dvsr_dr, vloc_r, dvloc_r
    real(8)  :: ratio, V, hvol
    real(8)  :: strs(3,3), strs_sum(3,3)
    real(8)  :: r_lo, r_hi, v_lo, v_hi, dv_lo, dv_hi

    V = system%det_a
    hvol = system%Hvol
    nspin = system%nspin
    strs = 0d0

    do ia = 1, system%nion
      ik = kion(ia)
      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        s(1) = ppg%rxyz(1,j,ia)
        s(2) = ppg%rxyz(2,j,ia)
        s(3) = ppg%rxyz(3,j,ia)
        r_abs = sqrt(s(1)**2 + s(2)**2 + s(3)**2)
        if(r_abs < 1d-12) cycle  ! skip on-site point
        r_inv = 1d0 / r_abs

        ! Sum electron density over spins at this grid point
        rho_val = 0d0
        do ispin = 1, nspin
          rho_val = rho_val + rho_s(ispin)%f(ix,iy,iz)
        end do

        ! Interpolate V_loc(r) and dV_loc/dr from the PP table.
        ! Find radial interval by bisection in pp%rad.
        call find_radial_index(r_abs, pp%rad(:,ik), pp%nrloc(ik), intr)
        if(intr < 2 .or. intr >= pp%nrloc(ik)) cycle

        r_lo = pp%rad(intr, ik)
        r_hi = pp%rad(intr+1, ik)
        v_lo = pp%vloctbl(intr, ik)
        v_hi = pp%vloctbl(intr+1, ik)
        dv_lo = pp%dvloctbl(intr, ik)
        dv_hi = pp%dvloctbl(intr+1, ik)
        ratio = (r_abs - r_lo) / (r_hi - r_lo)

        ! Linear interpolation of V_loc and dV_loc/dr
        vloc_r = v_lo + ratio * (v_hi - v_lo)
        dvloc_r = dv_lo + ratio * (dv_hi - dv_lo)

        ! V'_sr(r) = dV_loc/dr - Z/r^2
        ! (V_sr = V_loc + Z/r, so V'_sr = V'_loc - Z/r^2)
        dvsr_dr = dvloc_r - dble(pp%zps(ik)) * r_inv**2

        ! Accumulate: -rho * V'_sr * s_a * s_b / |s| * hvol
        do b = 1, 3
        do a = 1, 3
          strs(a,b) = strs(a,b) - rho_val * dvsr_dr * s(a) * s(b) * r_inv * hvol
        end do
        end do
      end do
    end do

    ! Use icomm_r (not icomm_rko): this is a density-only sum with no k/orbital loops.
    call comm_summation(strs, strs_sum, 9, info%icomm_r)
    system%stress_loc_sr_rs = -strs_sum / V

  contains

    subroutine find_radial_index(r, rad, nr, idx)
      implicit none
      real(8), intent(in)  :: r
      real(8), intent(in)  :: rad(:)
      integer, intent(in)  :: nr
      integer, intent(out) :: idx
      integer :: lo, hi, mid
      lo = 1
      hi = nr
      do while(hi - lo > 1)
        mid = (lo + hi) / 2
        if(rad(mid) <= r) then
          lo = mid
        else
          hi = mid
        end if
      end do
      idx = lo
    end subroutine find_radial_index

  end subroutine calc_stress_loc_sr_rs

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
