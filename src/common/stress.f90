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

  real(8), allocatable, target, save :: stress_rho_box_work(:,:,:)
  real(8), allocatable, target, save :: stress_grho_work(:,:,:,:)
  complex(8), allocatable, target, save :: stress_gtpsi_work(:,:,:,:)
  integer, save :: stress_work_is_array(3) = 0
  integer, save :: stress_work_ie_array(3) = -1
  integer, save :: stress_work_is(3) = 0
  integer, save :: stress_work_ie(3) = -1

contains

  subroutine fail_stress(message)
    use communication, only: comm_is_root
    use parallelization, only: end_parallel, nproc_id_global
    implicit none
    character(*), intent(in) :: message

    if(comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'Error in stress: '//trim(message)
    end if
    call end_parallel
    stop 1
  end subroutine fail_stress

  subroutine ensure_r2scan_stress_workspace(mg)
    use structures, only: s_rgrid
    implicit none
    type(s_rgrid), intent(in) :: mg
    logical :: resize_needed

    resize_needed = any(stress_work_is_array /= mg%is_array) .or. any(stress_work_ie_array /= mg%ie_array) .or. &
         any(stress_work_is /= mg%is) .or. any(stress_work_ie /= mg%ie)
    if (resize_needed) then
      if (allocated(stress_rho_box_work)) deallocate(stress_rho_box_work)
      if (allocated(stress_grho_work)) deallocate(stress_grho_work)
      if (allocated(stress_gtpsi_work)) deallocate(stress_gtpsi_work)
      stress_work_is_array = mg%is_array
      stress_work_ie_array = mg%ie_array
      stress_work_is = mg%is
      stress_work_ie = mg%ie
    end if

    if (.not. allocated(stress_rho_box_work)) then
      allocate(stress_rho_box_work(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))
    end if
    if (.not. allocated(stress_grho_work)) then
      allocate(stress_grho_work(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    end if
    if (.not. allocated(stress_gtpsi_work)) then
      allocate(stress_gtpsi_work(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3)))
    end if
  end subroutine ensure_r2scan_stress_workspace

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
                         srg, ppg, ppn, tpsi, ewald, energy, xc_func, rho_s, Vxc, field_state, srg_scalar)
    use structures
    use plusU_global,  only: PLUS_U_ON
    use salmon_global, only: xc
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
    type(s_sendrecv_grid),      intent(inout), optional :: srg_scalar

    if(PLUS_U_ON) call fail_stress("stress tensor is not supported with PLUS_U_ON")
    if(system%nspin /= 1) call fail_stress("stress tensor supports only unpolarized calculations")
    if(trim(xc) /= 'pz' .and. trim(xc) /= 'r2scan') then
      call fail_stress("stress tensor supports only built-in xc='PZ' or xc='r2scan'")
    end if

    if(info%if_divide_rspace) then
      if(.not. tpsi%update_zwf_overlap) then
        call update_overlap_complex8(srg, mg, tpsi%zwf)
        tpsi%update_zwf_overlap = .true.
      end if
    end if

    if(allocated(system%stress_nl_species_l)) deallocate(system%stress_nl_species_l)

    system%stress_kin = 0d0
    system%stress_har = 0d0
    system%stress_har_shadow = 0d0
    system%stress_xc = 0d0
    system%stress_xc_local = 0d0
    system%stress_xc_grad = 0d0
    system%stress_xc_grad_payload = 0d0
    system%stress_xc_grad_vsigma = 0d0
    system%stress_xc_tau = 0d0
    system%stress_x = 0d0
    system%stress_c = 0d0
    system%stress_xc_valence = 0d0
    system%stress_xc_nlcc = 0d0
    system%stress_loc = 0d0
    system%stress_loc_sr_grad = 0d0
    system%stress_loc_lr_grad = 0d0
    system%stress_loc_sr_diag = 0d0
    system%stress_loc_lr_diag = 0d0
    system%stress_loc_sr_scr_grad = 0d0
    system%stress_loc_lr_scr_grad = 0d0
    system%stress_loc_sr_scr_diag = 0d0
    system%stress_loc_lr_scr_diag = 0d0
    system%stress_loc_sr_rs = 0d0
    system%stress_loc_sr_rs_bins = 0d0
    system%stress_loc_sr_rs_legacy = 0d0
    system%stress_loc_sr_rs_legacy_bins = 0d0
    system%stress_nl = 0d0
    system%stress_ewa = 0d0
    system%stress_loc_sr_energy = 0d0
    system%stress_loc_lr_energy = 0d0
    system%stress_loc_sr_scr_energy = 0d0
    system%stress_loc_lr_scr_energy = 0d0
    system%stress_xc_e_vxc = 0d0
    system%stress_xc_energy_valence = 0d0
    system%stress_xc_energy_nlcc = 0d0
    system%stress_xc_e_vxc_valence = 0d0
    system%stress_xc_e_vxc_nlcc = 0d0
    system%stress_ewa_g = 0d0
    system%stress_ewa_r = 0d0
    system%stress_ewa_g_grad = 0d0
    system%stress_ewa_g_diag = 0d0
    system%stress_ewa_g_self = 0d0
    system%stress_ewa_energy_G = 0d0
    system%stress_ewa_energy_R = 0d0
    system%stress_tensor = 0d0
    system%stress_xc_dbg_grho_local_payload_maxdiff = 0d0
    system%stress_xc_dbg_grho_direct_payload_maxdiff = 0d0
    system%stress_xc_dbg_grho_direct_local_maxdiff = 0d0
    system%stress_xc_dbg_rdedd_dot_grho_local = 0d0
    system%stress_xc_dbg_rdedd_dot_grho_payload = 0d0
    system%stress_xc_dbg_rho_div_rdedd = 0d0

    call calc_stress_kin(system, info, mg, stencil, ppg, tpsi, field_state)
    call calc_stress_har(system, info, mg, fg, poisson, energy)
    call calc_stress_xc(system, pp, info, mg, stencil, srg, ppn, rho_s, Vxc, energy, xc_func, tpsi, field_state, srg_scalar)
    call calc_stress_loc(system, pp, fg, info, mg, ppg, poisson)
    call calc_stress_nl(system, pp, info, mg, stencil, ppg, tpsi, energy, field_state)
    call calc_stress_ewa(system, pp, fg, info, mg, ewald)

    call symmetrize_stress_term(system%stress_kin)
    call symmetrize_stress_term(system%stress_har)
    call symmetrize_stress_term(system%stress_xc)
    call symmetrize_stress_term(system%stress_xc_local)
    call symmetrize_stress_term(system%stress_xc_grad)
    call symmetrize_stress_term(system%stress_xc_grad_payload)
    call symmetrize_stress_term(system%stress_xc_grad_vsigma)
    call symmetrize_stress_term(system%stress_xc_tau)
    call symmetrize_stress_term(system%stress_xc_valence)
    call symmetrize_stress_term(system%stress_xc_nlcc)
    call symmetrize_stress_term(system%stress_loc)
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
    use communication, only: comm_summation, comm_get_max_array1d_double
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
    real(8) :: g(3), G2, coeff, strs(3,3), strs_sum(3,3), V, fourpi

    V = system%det_a
    fourpi = 4d0 * pi
    strs = 0d0
    call get_g_bounds(fg, ig_s, ig_e)

    !$omp parallel do collapse(2) default(none) &
    !$omp   private(ix,iy,iz,a,b,g,G2,coeff) &
    !$omp   shared(fg,poisson,ig_s,ig_e,fourpi) &
    !$omp   reduction(+:strs)
    do iz = ig_s(3), ig_e(3)
    do iy = ig_s(2), ig_e(2)
    do ix = ig_s(1), ig_e(1)
      if(fg%if_Gzero(ix,iy,iz)) cycle
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      G2 = g(1)**2 + g(2)**2 + g(3)**2
      coeff = abs(poisson%zrhoG_ele(ix,iy,iz))**2 * (fourpi / G2**2)
      do b = 1, 3
      do a = 1, 3
        strs(a,b) = strs(a,b) - coeff * g(a) * g(b)
      end do
      end do
    end do
    end do
    end do
    !$omp end parallel do

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

    !$omp parallel do collapse(2) default(none) private(ix,iy,iz) shared(box,Vh,mg)
    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      box(ix,iy,iz) = Vh%f(ix,iy,iz)
    end do
    end do
    end do
    !$omp end parallel do

    if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, box)
    call calc_gradient_field(mg, stencil%coef_nab, system%rmatrix_B, box, grad_vh)

    !$omp parallel do collapse(2) default(none) &
    !$omp   private(ix,iy,iz,a,b) &
    !$omp   shared(coeff,grad_vh,mg) &
    !$omp   reduction(+:strs)
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
    !$omp end parallel do

    call sum_stress_tensor(info%icomm_r, strs, strs_sum)
    strs_sum = strs_sum / V
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_h / V
    end do
    call symmetrize_stress_term(strs_sum)
    system%stress_har_shadow = -strs_sum
  end subroutine calc_stress_har_shadow

  subroutine calc_stress_xc(system, pp, info, mg, stencil, srg, ppn, rho_s, Vxc, energy, xc_func, tpsi, field_state, srg_scalar)
    use structures
    use salmon_global, only: xc
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_stencil),       intent(in)    :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_pp_nlcc),       intent(in)    :: ppn
    type(s_scalar),        intent(in)    :: rho_s(:)
    type(s_scalar),        intent(in)    :: Vxc(:)
    type(s_dft_energy),    intent(in)    :: energy
    type(s_xc_functional), intent(in)    :: xc_func
    type(s_orbital),        intent(inout) :: tpsi
    type(s_stress_field_state), intent(in) :: field_state
    type(s_sendrecv_grid), intent(inout), optional :: srg_scalar

    select case (trim(xc))
    case ('pz')
      call calc_stress_xc_builtin_pz(system, info, mg, ppn, rho_s, energy, xc_func)
      call calc_stress_xc_nlcc_cc(system, pp, info, mg, ppn, rho_s)
    case ('r2scan')
      call calc_stress_xc_builtin_r2scan(system, pp, info, mg, stencil, srg, ppn, rho_s, Vxc, energy, xc_func, tpsi, field_state, srg_scalar)
      call calc_stress_xc_nlcc_cc_r2scan(system, pp, info, mg, ppn, Vxc)
    case default
      if(xc_func%use_gradient) call fail_stress("stress tensor supports only built-in PZ or built-in r2SCAN XC")
      call fail_stress("stress tensor supports only built-in xc='PZ' or xc='r2scan'")
    end select
  end subroutine calc_stress_xc

  subroutine calc_stress_xc_builtin_pz(system, info, mg, ppn, rho_s, energy, xc_func)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_nlcc),       intent(in)    :: ppn
    type(s_scalar),        intent(in)    :: rho_s(:)
    type(s_dft_energy),    intent(in)    :: energy
    type(s_xc_functional), intent(in)    :: xc_func
    integer :: ix, iy, iz, ispin, a
    real(8) :: rho_xc, trho_xc
    real(8) :: exc_x, exc_c, vxc_x, vxc_c
    real(8) :: exc_x_val, exc_c_val, vxc_x_val, vxc_c_val
    real(8) :: E_vx_loc, E_vc_loc, E_vx, E_vc
    real(8) :: E_x_loc, E_c_loc, E_x, E_c
    real(8) :: E_vx_val_loc, E_vc_val_loc, E_vx_val, E_vc_val
    real(8) :: E_x_val_loc, E_c_val_loc, E_x_val, E_c_val
    real(8) :: E_xc_total, E_xc_valence, E_xc_nlcc
    real(8) :: E_vxc_total, E_vxc_valence, E_vxc_nlcc
    real(8) :: V

    if(xc_func%use_gradient) call fail_stress("stress tensor supports only built-in PZ LDA XC")

    V = system%det_a
    E_vx_loc = 0d0
    E_vc_loc = 0d0
    E_x_loc = 0d0
    E_c_loc = 0d0
    E_vx_val_loc = 0d0
    E_vc_val_loc = 0d0
    E_x_val_loc = 0d0
    E_c_val_loc = 0d0
    !$omp parallel do collapse(3) default(none) &
    !$omp   private(ispin,ix,iy,iz,rho_xc,trho_xc,exc_x,exc_c,vxc_x,vxc_c,exc_x_val,exc_c_val,vxc_x_val,vxc_c_val) &
    !$omp   shared(system,mg,ppn,rho_s) &
    !$omp   reduction(+:E_x_loc,E_c_loc,E_vx_loc,E_vc_loc,E_x_val_loc,E_c_val_loc,E_vx_val_loc,E_vc_val_loc)
    do ispin = 1, system%nspin
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_xc = rho_s(ispin)%f(ix,iy,iz)
        trho_xc = rho_xc
        if(allocated(ppn%rho_nlcc)) trho_xc = trho_xc + ppn%rho_nlcc(ix,iy,iz)
        call calc_builtin_pz_xc_split(trho_xc, exc_x, exc_c, vxc_x, vxc_c)
        call calc_builtin_pz_xc_split(rho_xc, exc_x_val, exc_c_val, vxc_x_val, vxc_c_val)
        E_x_loc = E_x_loc + trho_xc * exc_x
        E_c_loc = E_c_loc + trho_xc * exc_c
        E_vx_loc = E_vx_loc + rho_xc * vxc_x
        E_vc_loc = E_vc_loc + rho_xc * vxc_c
        E_x_val_loc = E_x_val_loc + rho_xc * exc_x_val
        E_c_val_loc = E_c_val_loc + rho_xc * exc_c_val
        E_vx_val_loc = E_vx_val_loc + rho_xc * vxc_x_val
        E_vc_val_loc = E_vc_val_loc + rho_xc * vxc_c_val
      end do
      end do
      end do
    end do
    !$omp end parallel do
    E_x_loc = E_x_loc * system%Hvol
    E_c_loc = E_c_loc * system%Hvol
    E_vx_loc = E_vx_loc * system%Hvol
    E_vc_loc = E_vc_loc * system%Hvol
    E_x_val_loc = E_x_val_loc * system%Hvol
    E_c_val_loc = E_c_val_loc * system%Hvol
    E_vx_val_loc = E_vx_val_loc * system%Hvol
    E_vc_val_loc = E_vc_val_loc * system%Hvol
    call comm_summation(E_x_loc, E_x, info%icomm_r)
    call comm_summation(E_c_loc, E_c, info%icomm_r)
    call comm_summation(E_vx_loc, E_vx, info%icomm_r)
    call comm_summation(E_vc_loc, E_vc, info%icomm_r)
    call comm_summation(E_x_val_loc, E_x_val, info%icomm_r)
    call comm_summation(E_c_val_loc, E_c_val, info%icomm_r)
    call comm_summation(E_vx_val_loc, E_vx_val, info%icomm_r)
    call comm_summation(E_vc_val_loc, E_vc_val, info%icomm_r)

    E_xc_total = E_x + E_c
    E_xc_valence = E_x_val + E_c_val
    E_xc_nlcc = E_xc_total - E_xc_valence
    E_vxc_total = E_vx + E_vc
    E_vxc_valence = E_vx_val + E_vc_val
    E_vxc_nlcc = E_vxc_total - E_vxc_valence

    system%stress_x = 0d0
    system%stress_c = 0d0
    system%stress_xc = 0d0
    system%stress_xc_valence = 0d0
    system%stress_xc_nlcc = 0d0
    do a = 1, 3
      system%stress_x(a,a) = -(E_vx - E_x) / V
      system%stress_c(a,a) = -(E_vc - E_c) / V
      system%stress_xc_valence(a,a) = -(E_vxc_valence - E_xc_valence) / V
    end do
    system%stress_xc = system%stress_x + system%stress_c
    system%stress_xc_nlcc = system%stress_xc - system%stress_xc_valence
    system%stress_xc_e_vxc = E_vx + E_vc
    system%stress_xc_energy_valence = E_xc_valence
    system%stress_xc_energy_nlcc = E_xc_nlcc
    system%stress_xc_e_vxc_valence = E_vxc_valence
    system%stress_xc_e_vxc_nlcc = E_vxc_nlcc
  end subroutine calc_stress_xc_builtin_pz

  subroutine calc_stress_xc_nlcc_cc(system, pp, info, mg, ppn, rho_s)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: kion
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_nlcc),       intent(in)    :: ppn
    type(s_scalar),        intent(in)    :: rho_s(:)
    integer :: j1, j2, j3, ispin
    real(8) :: trho, exc_x, exc_c, vxc_x, vxc_c, vxc_sum
    real(8) :: strs(3,3), strs_sum(3,3), vol
    real(8), allocatable :: vxc_box(:,:,:)

    ! Strain derivative of E_xc from the rigid, atom-centered NLCC core
    ! density: sigma_cc_ab = (1/V) * int vxc(r) * sum_a rhoc'(|s|) s_a s_b / |s|.
    ! Real-space equivalent of QE PW/src/stres_cc.f90 (diagonal + nondiagonal,
    ! G=0 included). The valence-dilution diagonal -(E_vxc - E_xc)/V is already
    ! handled by calc_stress_xc_builtin_pz.
    system%stress_xc_cc = 0d0
    if(.not. pp%flag_nlcc) return
    if(.not. allocated(ppn%rho_nlcc)) return

    vol = system%det_a
    strs = 0d0

    allocate(vxc_box(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    ! spin-averaged PZ vxc of the total (valence + core) density,
    ! consistent with the builtin-PZ stress energy path
    !$omp parallel do collapse(2) default(none) &
    !$omp   private(j1,j2,j3,ispin,trho,exc_x,exc_c,vxc_x,vxc_c,vxc_sum) &
    !$omp   shared(mg,system,rho_s,ppn,vxc_box)
    do j3 = mg%is(3), mg%ie(3)
    do j2 = mg%is(2), mg%ie(2)
    do j1 = mg%is(1), mg%ie(1)
      vxc_sum = 0d0
      do ispin = 1, system%nspin
        trho = rho_s(ispin)%f(j1,j2,j3) + ppn%rho_nlcc(j1,j2,j3)
        call calc_builtin_pz_xc_split(trho, exc_x, exc_c, vxc_x, vxc_c)
        vxc_sum = vxc_sum + vxc_x + vxc_c
      end do
      vxc_box(j1,j2,j3) = vxc_sum / dble(system%nspin)
    end do
    end do
    end do
    !$omp end parallel do

    call accumulate_stress_nlcc_cc(system, pp, mg, vxc_box, tbl=pp%rho_nlcc_tbl, strs=strs)
    deallocate(vxc_box)

    call comm_summation(strs, strs_sum, 9, info%icomm_r)
    system%stress_xc_cc = strs_sum / vol
    system%stress_xc = system%stress_xc + system%stress_xc_cc
  end subroutine calc_stress_xc_nlcc_cc

  subroutine accumulate_stress_nlcc_cc(system, pp, mg, pot_box, tbl, strs)
    use structures
    use salmon_global, only: kion
    implicit none
    type(s_dft_system), intent(in)    :: system
    type(s_pp_info),    intent(in)    :: pp
    type(s_rgrid),      intent(in)    :: mg
    real(8),            intent(in)    :: pot_box(mg%is(1):,mg%is(2):,mg%is(3):)
    real(8),            intent(in)    :: tbl(:,:)
    real(8),            intent(inout) :: strs(3,3)
    integer :: a, ik, i, ir, intr, i1, i2, i3, j1, j2, j3, ia, ib
    real(8) :: rc, u, v, w, r, s(3), drc_dr, contrib, hvol, Rion_repr(3)
    logical :: flag_cuboid

    hvol = system%Hvol
    flag_cuboid = .true.
    if( abs(system%primitive_a(1,2)) >= 1d-10 .or. &
        abs(system%primitive_a(1,3)) >= 1d-10 .or. &
        abs(system%primitive_a(2,3)) >= 1d-10 ) flag_cuboid = .false.

    ! geometry mirrors calc_nlcc (salmon_pp.f90): atoms x (+/-2 replicas) x local grid
    !$omp parallel do default(none) &
    !$omp   private(a,ik,rc,i,ir,intr,i1,i2,i3,j1,j2,j3,u,v,w,r,s,drc_dr,contrib,Rion_repr,ia,ib) &
    !$omp   shared(system,pp,mg,kion,pot_box,tbl,hvol,flag_cuboid) &
    !$omp   reduction(+:strs)
    do a = 1, system%nion
      ik = kion(a)
      rc = 15d0
      do i = 1, pp%nrmax
        if(pp%rho_nlcc_tbl(i,ik) + pp%tau_nlcc_tbl(i,ik) < 1d-6) then
          rc = pp%rad(i,ik)
          exit
        end if
      end do

      do i1 = -2, 2
      do i2 = -2, 2
      do i3 = -2, 2
        ! full lattice-vector translation (see calc_nlcc): the diagonal-only
        ! form misplaces NLCC periodic images for nonorthogonal cells
        Rion_repr(1) = system%Rion(1,a) + i1 * system%primitive_a(1,1) + i2 * system%primitive_a(1,2) + i3 * system%primitive_a(1,3)
        Rion_repr(2) = system%Rion(2,a) + i1 * system%primitive_a(2,1) + i2 * system%primitive_a(2,2) + i3 * system%primitive_a(2,3)
        Rion_repr(3) = system%Rion(3,a) + i1 * system%primitive_a(3,1) + i2 * system%primitive_a(3,2) + i3 * system%primitive_a(3,3)

        do j1 = mg%is(1), mg%ie(1)
          u = (j1-1) * system%hgs(1)
        do j2 = mg%is(2), mg%ie(2)
          v = (j2-1) * system%hgs(2)
        do j3 = mg%is(3), mg%ie(3)
          w = (j3-1) * system%hgs(3)
          if(flag_cuboid) then
            s(1) = u - Rion_repr(1)
            s(2) = v - Rion_repr(2)
            s(3) = w - Rion_repr(3)
          else
            s(1) = u*system%rmatrix_a(1,1) + v*system%rmatrix_a(1,2) + w*system%rmatrix_a(1,3) - Rion_repr(1)
            s(2) = u*system%rmatrix_a(2,1) + v*system%rmatrix_a(2,2) + w*system%rmatrix_a(2,3) - Rion_repr(2)
            s(3) = u*system%rmatrix_a(3,1) + v*system%rmatrix_a(3,2) + w*system%rmatrix_a(3,3) - Rion_repr(3)
          end if
          r = sqrt(s(1)**2 + s(2)**2 + s(3)**2)
          if(r > rc .or. r < 1d-12) cycle
          do ir = 1, pp%nrmax
            if(pp%rad(ir,ik) > r) exit
          end do
          intr = ir - 1
          if(intr < 1 .or. intr >= pp%nrmax) cycle
          ! segment slope of the radial NLCC table: matches the linear
          ! interpolation used to build ppn%rho_nlcc/tau_nlcc in calc_nlcc
          drc_dr = (tbl(intr+1,ik) - tbl(intr,ik)) &
                 / (pp%rad(intr+1,ik) - pp%rad(intr,ik))
          contrib = pot_box(j1,j2,j3) * drc_dr / r * hvol
          do ib = 1, 3
          do ia = 1, 3
            strs(ia,ib) = strs(ia,ib) + contrib * s(ia) * s(ib)
          end do
          end do
        end do
        end do
        end do
      end do
      end do
      end do
    end do
    !$omp end parallel do
  end subroutine accumulate_stress_nlcc_cc

  subroutine calc_stress_xc_nlcc_cc_r2scan(system, pp, info, mg, ppn, Vxc)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: yn_tau_nlcc
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_nlcc),       intent(in)    :: ppn
    type(s_scalar),        intent(in)    :: Vxc(:)
    integer :: j1, j2, j3
    real(8) :: strs(3,3), strs_sum(3,3), vol
    real(8), allocatable :: pot_box(:,:,:), vtau_box(:,:,:)

    ! NLCC cc stress for the meta-GGA path. Uses the multiplicative KS
    ! potential Vxc = vrho - div(2*vsigma*grad n), which equals dE/dn at
    ! fixed tau, so the gradient (vsigma) core response is included by
    ! integration by parts. No vtau term: tau_nlcc is not consumed by the XC energy
    ! (salmon_xc.f90 passes rho_nlcc only), so adding one would break
    ! energy/stress consistency.
    system%stress_xc_cc = 0d0
    if(.not. pp%flag_nlcc) return
    if(.not. allocated(ppn%rho_nlcc)) return

    vol = system%det_a
    strs = 0d0
    allocate(pot_box(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    !$omp parallel do collapse(2) default(none) private(j1,j2,j3) shared(mg,pot_box,Vxc)
    do j3 = mg%is(3), mg%ie(3)
    do j2 = mg%is(2), mg%ie(2)
    do j1 = mg%is(1), mg%ie(1)
      pot_box(j1,j2,j3) = Vxc(1)%f(j1,j2,j3)
    end do
    end do
    end do
    !$omp end parallel do

    call accumulate_stress_nlcc_cc(system, pp, mg, pot_box, tbl=pp%rho_nlcc_tbl, strs=strs)
    deallocate(pot_box)

    if (yn_tau_nlcc == 'y' .and. allocated(system%xc_payload%vtau%f)) then
      ! second cc term: int vtau * d(tau_core)/d(eps_ab); only meaningful
      ! when the energy itself sees tau_core (yn_tau_nlcc='y')
      allocate(vtau_box(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
      !$omp parallel do collapse(2) default(none) private(j1,j2,j3) shared(mg,vtau_box,system)
      do j3 = mg%is(3), mg%ie(3)
      do j2 = mg%is(2), mg%ie(2)
      do j1 = mg%is(1), mg%ie(1)
        vtau_box(j1,j2,j3) = system%xc_payload%vtau%f(mg%idx(j1), mg%idy(j2), mg%idz(j3))
      end do
      end do
      end do
      !$omp end parallel do
      call accumulate_stress_nlcc_cc(system, pp, mg, vtau_box, tbl=pp%tau_nlcc_tbl, strs=strs)
      deallocate(vtau_box)
    end if

    call comm_summation(strs, strs_sum, 9, info%icomm_r)
    system%stress_xc_cc = strs_sum / vol
    system%stress_xc = system%stress_xc + system%stress_xc_cc
  end subroutine calc_stress_xc_nlcc_cc_r2scan

  subroutine calc_stress_xc_builtin_r2scan(system, pp, info, mg, stencil, srg, ppn, rho_s, Vxc, energy, xc_func, tpsi, field_state, srg_scalar)
    use structures
    use stencil_sub, only: calc_gradient_field, calc_gradient_psi
    use communication, only: comm_summation, comm_get_max_array1d_double
    use sendrecv_grid, only: update_overlap_real8
    use math_constants, only: zi
    use salmon_global, only: yn_out_stress_numerics, yn_nlcc_grho
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_pp_info),            intent(in)    :: pp
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_sendrecv_grid),      intent(inout) :: srg
    type(s_pp_nlcc),            intent(in)    :: ppn
    type(s_scalar),             intent(in)    :: rho_s(:)
    type(s_dft_energy),         intent(in)    :: energy
    type(s_scalar),             intent(in)    :: Vxc(:)
    type(s_xc_functional),      intent(in)    :: xc_func
    type(s_orbital),            intent(inout) :: tpsi
    type(s_stress_field_state), intent(in)    :: field_state
    type(s_sendrecv_grid),      intent(inout), optional :: srg_scalar
    integer :: ix, iy, iz, ik, io, ispin, im, idir, a, b
    real(8) :: rho_loc, V, E_vxc_loc, E_vxc, rtmp, vtau_loc, vsigma_loc
    real(8) :: kAc(3)
    real(8) :: strs_grad(3,3), strs_grad_sum(3,3)
    real(8) :: strs_grad_payload(3,3), strs_grad_payload_sum(3,3)
    real(8) :: strs_grad_vsigma(3,3), strs_grad_vsigma_sum(3,3)
    real(8) :: strs_tau(3,3), strs_tau_sum(3,3)
    logical :: want_numerics_diag, has_grho_payload
    real(8), pointer :: rho_box(:,:,:)
    real(8), pointer :: grho_local(:,:,:,:)
    complex(8), pointer :: gtpsi(:,:,:,:)
    complex(8) :: psi_r, w(3)

    if (system%nspin /= 1) call fail_stress("r2scan stress supports only nspin=1")
    if (.not. allocated(system%xc_payload%rdedd%v)) call fail_stress("r2scan stress requires rdedd payload")
    if (.not. system%xc_payload%rdedd_has_shadow_values) call fail_stress("r2scan stress requires rdedd payload with shadow values")
    if (.not. allocated(system%xc_payload%vtau%f)) call fail_stress("r2scan stress requires vtau payload")
    if (.not. system%xc_payload%vtau_has_shadow_values) call fail_stress("r2scan stress requires vtau payload with shadow values")
    if (.not. allocated(system%xc_payload%vsigma%f)) call fail_stress("r2scan stress requires vsigma payload")
    if (.not. system%xc_payload%vsigma_has_shadow_values) call fail_stress("r2scan stress requires vsigma payload with shadow values")

    V = system%det_a
    want_numerics_diag = (trim(yn_out_stress_numerics) == 'y')
    has_grho_payload = allocated(system%xc_payload%grho%v)
    call ensure_r2scan_stress_workspace(mg)
    rho_box => stress_rho_box_work
    grho_local => stress_grho_work
    gtpsi => stress_gtpsi_work

    rho_box = 0d0
    !$omp parallel do collapse(2) default(none) private(ix,iy,iz) shared(rho_box,rho_s,mg)
    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      rho_box(ix,iy,iz) = dble(rho_s(1)%f(ix,iy,iz))
    end do
    end do
    end do
    !$omp end parallel do
    if (yn_nlcc_grho == 'y' .and. allocated(ppn%rho_nlcc)) then
      ! NLCC core-inclusive gradient (matches the energy when yn_nlcc_grho='y'):
      ! the r2scan grad stress (rdedd/vsigma) then uses |grad(rho_val+rho_core)|^2,
      ! consistent with the payload computed from the core-inclusive density.
      ! Default 'n' keeps rho_box valence-only (bit-identical, gradient unchanged).
      ! rho_box feeds only grho_local (the local E_vxc term reads rho_s directly).
!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_box(ix,iy,iz) = rho_box(ix,iy,iz) + ppn%rho_nlcc(ix,iy,iz)
      end do
      end do
      end do
!$omp end parallel do
    end if
    if (info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rho_box)
    call calc_gradient_field(mg, stencil%coef_nab, system%rmatrix_B, rho_box, grho_local)

    E_vxc_loc = 0d0
    strs_grad = 0d0
    strs_grad_payload = 0d0
    strs_grad_vsigma = 0d0
    !$omp parallel do collapse(2) default(none) &
    !$omp   private(ix,iy,iz,a,b,rho_loc,vsigma_loc) &
    !$omp   shared(system,Vxc,rho_s,grho_local,mg,has_grho_payload) &
    !$omp   reduction(+:E_vxc_loc,strs_grad,strs_grad_payload,strs_grad_vsigma)
    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      rho_loc = dble(rho_s(1)%f(ix,iy,iz))
      E_vxc_loc = E_vxc_loc + rho_loc * Vxc(1)%f(ix,iy,iz)
      vsigma_loc = system%xc_payload%vsigma%f(mg%idx(ix), mg%idy(iy), mg%idz(iz))
      do b = 1, 3
      do a = 1, 3
        strs_grad(a,b) = strs_grad(a,b) + system%Hvol * &
             system%xc_payload%rdedd%v(a, mg%idx(ix), mg%idy(iy), mg%idz(iz)) * grho_local(b,ix,iy,iz)
        if (has_grho_payload) then
          strs_grad_payload(a,b) = strs_grad_payload(a,b) + system%Hvol * &
               system%xc_payload%rdedd%v(a, mg%idx(ix), mg%idy(iy), mg%idz(iz)) * &
               system%xc_payload%grho%v(b, ix, iy, iz)
        end if
        strs_grad_vsigma(a,b) = strs_grad_vsigma(a,b) + system%Hvol * &
             2d0 * vsigma_loc * grho_local(a,ix,iy,iz) * grho_local(b,ix,iy,iz)
      end do
      end do
    end do
    end do
    end do
    !$omp end parallel do

    strs_tau = 0d0
    im = 1
    do ik = info%ik_s, info%ik_e
    do io = info%io_s, info%io_e
    do ispin = 1, system%nspin
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im), gtpsi, &
           mg%is_array, mg%ie_array, mg%is, mg%ie, mg%idx, mg%idy, mg%idz, stencil%coef_nab, system%rmatrix_B)
      rtmp = system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol
      !$omp parallel do collapse(2) default(none) &
      !$omp   private(ix,iy,iz,kAc,psi_r,w,a,b,vtau_loc) &
      !$omp   shared(mg,field_state,system,ik,tpsi,ispin,io,im,gtpsi,rtmp) &
      !$omp   reduction(+:strs_tau)
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
        vtau_loc = system%xc_payload%vtau%f(mg%idx(ix), mg%idy(iy), mg%idz(iz))
        psi_r = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        w(1) = gtpsi(1,ix,iy,iz) + zi * kAc(1) * psi_r
        w(2) = gtpsi(2,ix,iy,iz) + zi * kAc(2) * psi_r
        w(3) = gtpsi(3,ix,iy,iz) + zi * kAc(3) * psi_r
        do b = 1, 3
        do a = 1, 3
          strs_tau(a,b) = strs_tau(a,b) - rtmp * vtau_loc * dble(conjg(w(a)) * w(b))
        end do
        end do
      end do
      end do
      end do
      !$omp end parallel do
    end do
    end do
    end do

    call comm_summation(E_vxc_loc * system%Hvol, E_vxc, info%icomm_r)
    call comm_summation(strs_grad, strs_grad_sum, 9, info%icomm_r)
    if (allocated(system%xc_payload%grho%v)) then
      call comm_summation(strs_grad_payload, strs_grad_payload_sum, 9, info%icomm_r)
    else
      strs_grad_payload_sum = 0d0
    end if
    call comm_summation(strs_grad_vsigma, strs_grad_vsigma_sum, 9, info%icomm_r)
    call comm_summation(strs_tau, strs_tau_sum, 9, info%icomm_rko)

    system%stress_x = 0d0
    system%stress_c = 0d0
    system%stress_xc = 0d0
    system%stress_xc_local = 0d0
    system%stress_xc_grad = 0d0
    system%stress_xc_grad_payload = 0d0
    system%stress_xc_grad_vsigma = 0d0
    system%stress_xc_tau = 0d0
    do a = 1, 3
      system%stress_xc_local(a,a) = -(E_vxc - energy%E_xc) / V
    end do
    system%stress_xc_grad = strs_grad_sum / V
    do b = 1, 3
    do a = 1, 3
      system%stress_xc_grad_payload(a,b) = system%stress_xc_grad_payload(a,b) + strs_grad_payload_sum(a,b) / V
      system%stress_xc_grad_vsigma(a,b) = system%stress_xc_grad_vsigma(a,b) + strs_grad_vsigma_sum(a,b) / V
    end do
    end do
    system%stress_xc_tau = strs_tau_sum / V
    system%stress_xc = system%stress_xc_local + system%stress_xc_grad + system%stress_xc_tau
    system%stress_xc_e_vxc = E_vxc

    if (want_numerics_diag) then
      if (present(srg_scalar)) then
        call calc_r2scan_high_diagnostics(system, info, mg, stencil, srg_scalar, rho_s, grho_local)
      end if
    end if
  end subroutine calc_stress_xc_builtin_r2scan

  subroutine calc_r2scan_high_diagnostics(system, info, mg, stencil, srg_scalar, rho_s, grho_local)
    use structures
    use stencil_sub, only: calc_gradient_field
    use sendrecv_grid, only: update_overlap_real8
    use communication, only: comm_get_max_array1d_double, comm_summation
    implicit none
    type(s_dft_system),         intent(inout) :: system
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_stencil),            intent(in)    :: stencil
    type(s_sendrecv_grid),      intent(inout) :: srg_scalar
    type(s_scalar),             intent(in)    :: rho_s(:)
    real(8),                    intent(in)    :: grho_local(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), &
                                                             mg%is(3):mg%ie(3))
    real(8) :: rhd_direct(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
    real(8) :: rdedd_component(mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), &
                               mg%is_array(3):mg%ie_array(3))
    real(8) :: div_component(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8) :: grho_direct(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: ispin, idir, ix, iy, iz
    real(8) :: maxdiff_in(3), maxdiff_out(3)
    real(8) :: diag_in(3), diag_out(3), rho_loc, rdedd_loc

    maxdiff_in = 0d0
    diag_in = 0d0
    rhd_direct = 0d0
    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      rhd_direct(ix,iy,iz) = dble(rho_s(1)%f(ix,iy,iz))
    end do
    end do
    end do
    if (info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, rhd_direct)
    call calc_gradient_field(mg, stencil%coef_nab, system%rmatrix_B, rhd_direct, grho_direct)

    do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      rho_loc = dble(rho_s(1)%f(ix,iy,iz))
      do idir = 1, 3
        rdedd_loc = system%xc_payload%rdedd%v(idir, mg%idx(ix), mg%idy(iy), mg%idz(iz))
        maxdiff_in(1) = max(maxdiff_in(1), abs(system%xc_payload%grho%v(idir,ix,iy,iz) - &
                           grho_local(idir,ix,iy,iz)))
        maxdiff_in(2) = max(maxdiff_in(2), abs(system%xc_payload%grho%v(idir,ix,iy,iz) - &
                           grho_direct(idir,ix,iy,iz)))
        maxdiff_in(3) = max(maxdiff_in(3), abs(grho_direct(idir,ix,iy,iz) - &
                           grho_local(idir,ix,iy,iz)))
        diag_in(1) = diag_in(1) + system%Hvol * rdedd_loc * grho_local(idir,ix,iy,iz)
        diag_in(2) = diag_in(2) + system%Hvol * rdedd_loc * system%xc_payload%grho%v(idir,ix,iy,iz)
      end do
    end do
    end do
    end do

    do idir = 1, 3
      rdedd_component = system%xc_payload%rdedd%v(idir,:,:,:)
      call calc_gradient_field(mg, stencil%coef_nab, system%rmatrix_B, rdedd_component, div_component)
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_loc = dble(rho_s(1)%f(ix,iy,iz))
        diag_in(3) = diag_in(3) + system%Hvol * rho_loc * div_component(idir,ix,iy,iz)
      end do
      end do
      end do
    end do

    call comm_get_max_array1d_double(maxdiff_in, maxdiff_out, 3, info%icomm_rko)
    call comm_summation(diag_in, diag_out, 3, info%icomm_rko)
    system%stress_xc_dbg_grho_local_payload_maxdiff = maxdiff_out(1)
    system%stress_xc_dbg_grho_direct_payload_maxdiff = maxdiff_out(2)
    system%stress_xc_dbg_grho_direct_local_maxdiff = maxdiff_out(3)
    system%stress_xc_dbg_rdedd_dot_grho_local = diag_out(1)
    system%stress_xc_dbg_rdedd_dot_grho_payload = diag_out(2)
    system%stress_xc_dbg_rho_div_rdedd = diag_out(3)
  end subroutine calc_r2scan_high_diagnostics

  pure subroutine calc_builtin_pz_xc_split(trho, exc_x, exc_c, vxc_x, vxc_c)
    implicit none
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: gammaU = -0.1423d0, beta1U = 1.0529d0
    real(8), parameter :: beta2U = 0.3334d0, AU = 0.0311d0, BU = -0.048d0
    real(8), parameter :: CU = 0.002d0, DU = -0.0116d0
    real(8), parameter :: const = 0.75d0 * (3d0 / (2d0 * pi))**(2d0 / 3d0)
    real(8), intent(in)  :: trho
    real(8), intent(out) :: exc_x, exc_c, vxc_x, vxc_c
    real(8) :: ttrho, rs, rssq, rsln
    real(8) :: dexc_x_drho, dexc_c_drho

    ttrho = trho + 1d-10
    rs = (3d0 / (4d0 * pi * ttrho))**(1d0 / 3d0)

    exc_x = -const / rs
    dexc_x_drho = exc_x / (3d0 * ttrho)
    vxc_x = exc_x + ttrho * dexc_x_drho

    if(rs > 1d0) then
      rssq = sqrt(rs)
      exc_c = gammaU / (1d0 + beta1U * rssq + beta2U * rs)
      dexc_c_drho = gammaU * (0.5d0 * beta1U * rssq + beta2U * rs) / (3d0 * ttrho) &
                  / (1d0 + beta1U * rssq + beta2U * rs)**2
    else
      rsln = log(rs)
      exc_c = AU * rsln + BU + CU * rs * rsln + DU * rs
      dexc_c_drho = -rs / (3d0 * ttrho) * (AU / rs + CU * (rsln + 1d0) + DU)
    end if
    vxc_c = exc_c + ttrho * dexc_c_drho
  end subroutine calc_builtin_pz_xc_split

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
      rtmp = system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol
      !$omp parallel do collapse(2) default(none) &
      !$omp   private(ix,iy,iz,kAc,psi_r,w,a,b) &
      !$omp   shared(mg,field_state,system,ik,tpsi,ispin,io,im,gtpsi,rtmp) &
      !$omp   reduction(+:strs)
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

  subroutine calc_stress_loc(system, pp, fg, info, mg, ppg, poisson)
    use structures
    use math_constants, only: pi, zi
    use communication,  only: comm_summation
    use filesystem,     only: open_filehandle
    use inputoutput,    only: au_pressure_gpa
    use prep_pp_sub,    only: integrate_local_sr_dvg_dg2_shell
    use salmon_global,  only: aEwald, base_directory, cutoff_g, kion, &
                              npt_loc_sr_aux_2pi, sysname, yn_out_loc_sr_grad_sampled_dump, &
                              yn_out_stress_numerics
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_rgrid),           intent(in)    :: mg
    type(s_pp_grid),         intent(in)    :: ppg
    type(s_poisson),         intent(in)    :: poisson
    integer :: ix, iy, iz, ia, ik, a, b, ig_s(3), ig_e(3), fh_grad_dump
    real(8) :: g(3), r(3), G2, Gd, coeff_lr, coeff_lr_scr, scr_fac, strs(3,3), strs_sum(3,3), &
      & strs_grad_sum(3,3), strs_diag_sum(3,3), strs_sr(3,3), strs_lr(3,3), strs_lr_scr(3,3), &
      & strs_sr_sum(3,3), strs_lr_sum(3,3), strs_lr_scr_sum(3,3), E_sr, E_lr, E_sr_scr, E_lr_scr, &
      & E_sr_loc, E_lr_loc, E_lr_scr_loc, V, cutoff_g2, gmag, legacy_dvg_dg2, pressure_factor, &
      & pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa
    complex(8) :: rho_e, V_sr_sum, V_lr_sum, V_lr_scr_sum, dVsr_dG2_sum, dVsr_dG2_current, &
      & dVsr_dG2_legacy, phase
    logical :: dump_grad
    character(256) :: file_grad_dump

    V = system%det_a
    strs = 0d0
    strs_sr = 0d0
    strs_lr = 0d0
    strs_lr_scr = 0d0
    E_sr_loc = 0d0
    E_lr_loc = 0d0
    E_lr_scr_loc = 0d0
    cutoff_g2 = cutoff_g**2
    pressure_factor = au_pressure_gpa / (3d0 * V)
    call get_g_bounds(fg, ig_s, ig_e)
    if(.not. allocated(ppg%zVG_ion_stress)) call fail_stress("calc_stress_loc: zVG_ion_stress is not allocated")
    if(.not. allocated(ppg%dVG_ion_dG2)) call fail_stress("calc_stress_loc: dVG_ion_dG2 is not allocated")
    do ia = 1, system%nion
      ik = kion(ia)
      if(ik < 1 .or. ik > size(pp%zps)) call fail_stress("calc_stress_loc: invalid species index")
    end do
    dump_grad = (yn_out_loc_sr_grad_sampled_dump == 'y' .and. info%id_k == 0 .and. info%id_o == 0)
    fh_grad_dump = -1
    if(dump_grad) call open_loc_sr_grad_sampled_dump()

    !$omp parallel do collapse(2) default(none) &
    !$omp& private(ix,iy,iz,ia,ik,a,b,g,r,G2,gmag,Gd,rho_e, &
    !$omp&         V_sr_sum,V_lr_sum,V_lr_scr_sum,dVsr_dG2_sum, &
    !$omp&         dVsr_dG2_current,dVsr_dG2_legacy,phase,scr_fac, &
    !$omp&         coeff_lr,coeff_lr_scr,legacy_dvg_dg2, &
    !$omp&         pressure_current_gpa,pressure_legacy_gpa,delta_pressure_gpa) &
    !$omp& shared(fg,ig_s,ig_e,poisson,ppg,pp,system,kion,cutoff_g2, &
    !$omp&        aEwald,dump_grad,fh_grad_dump,npt_loc_sr_aux_2pi,pressure_factor) &
    !$omp& reduction(+:strs_sr,strs_lr,strs_lr_scr,E_sr_loc,E_lr_loc,E_lr_scr_loc)
    do iz = ig_s(3), ig_e(3)
    do iy = ig_s(2), ig_e(2)
    do ix = ig_s(1), ig_e(1)
      if(fg%if_Gzero(ix,iy,iz)) then
        ! G=0 short-range (alpha-Z) term: rho(0) * sum_a V_sr,a(0).
        ! Keeps E_sr consistent with the total-energy bookkeeping
        ! (total_energy.f90 includes G=0 in the electron-ion core sum).
        ! The LR term has no G=0 (1/G^2 singular, neutrality-cancelled)
        ! and the gradient term carries G_a*G_b = 0.
        rho_e = poisson%zrhoG_ele(ix,iy,iz)
        V_sr_sum = (0d0, 0d0)
        do ia = 1, system%nion
          ik = kion(ia)
          V_sr_sum = V_sr_sum + ppg%zVG_ion_stress(ix,iy,iz,ik)
        end do
        E_sr_loc = E_sr_loc + dble(conjg(rho_e) * V_sr_sum)
        cycle
      end if
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      G2 = g(1)**2 + g(2)**2 + g(3)**2
      gmag = sqrt(G2)
      if(G2 > cutoff_g2) cycle

      rho_e = poisson%zrhoG_ele(ix,iy,iz)
      V_sr_sum = (0d0, 0d0)
      V_lr_sum = (0d0, 0d0)
      V_lr_scr_sum = (0d0, 0d0)
      dVsr_dG2_sum = (0d0, 0d0)
      dVsr_dG2_legacy = (0d0, 0d0)

      do ia = 1, system%nion
        ik = kion(ia)
        r(:) = system%Rion(1:3,ia)
        Gd = sum(g(:) * r(:))
        phase = dcmplx(cos(Gd), -sin(Gd))
        V_sr_sum = V_sr_sum + ppg%zVG_ion_stress(ix,iy,iz,ik) * phase
        V_lr_sum = V_lr_sum - (4d0*pi / G2) * pp%zps(ik) * phase
        dVsr_dG2_sum = dVsr_dG2_sum + ppg%dVG_ion_dG2(ix,iy,iz,ik) * phase
        if(dump_grad) then
          legacy_dvg_dg2 = integrate_local_sr_dvg_dg2_shell(pp, ik, gmag, npt_loc_sr_aux_2pi)
          dVsr_dG2_legacy = dVsr_dG2_legacy + legacy_dvg_dg2 * phase
        end if
      end do

      E_sr_loc = E_sr_loc + dble(conjg(rho_e) * V_sr_sum)
      E_lr_loc = E_lr_loc + dble(conjg(rho_e) * V_lr_sum)
      scr_fac = fg%exp_ewald(ix,iy,iz)
      V_lr_scr_sum = scr_fac * V_lr_sum
      E_lr_scr_loc = E_lr_scr_loc + dble(conjg(rho_e) * V_lr_scr_sum)
      coeff_lr = -2d0 * dble(conjg(rho_e) * V_lr_sum) / G2
      coeff_lr_scr = -2d0 * dble(conjg(rho_e) * V_lr_scr_sum) / G2 * (1d0 + G2 / (4d0 * aEwald))
      if(dump_grad) then
        dVsr_dG2_current = dVsr_dG2_sum
        pressure_current_gpa = 2d0 * dble(conjg(rho_e) * dVsr_dG2_current) * G2 * pressure_factor
        pressure_legacy_gpa = 2d0 * dble(conjg(rho_e) * dVsr_dG2_legacy) * G2 * pressure_factor
        delta_pressure_gpa = pressure_current_gpa - pressure_legacy_gpa
        !$omp critical(loc_sr_grad_sampled_dump_io)
        call dump_loc_sr_grad_sampled_points(ix, iy, iz, g, gmag, rho_e, dVsr_dG2_current, dVsr_dG2_legacy, &
             pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa)
        !$omp end critical(loc_sr_grad_sampled_dump_io)
      end if

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
    !$omp end parallel do

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
    end do
    system%stress_loc_sr_energy = E_sr
    system%stress_loc_lr_energy = E_lr
    system%stress_loc_sr_scr_energy = E_sr_scr
    system%stress_loc_lr_scr_energy = E_lr_scr
    system%stress_loc = -strs_sum
    if(fh_grad_dump >= 0) close(fh_grad_dump)
    if(yn_out_stress_numerics == 'y') then
      call calc_local_sr_gspace_diagnostics(system, pp, fg, info, ppg, cutoff_g2)
    else
      system%stress_loc_sr_gspace_dvg_maxdiff = 0d0
      system%stress_loc_sr_gspace_dvg_meandiff = 0d0
      system%stress_loc_sr_gspace_dvg_g_at_max = 0d0
    end if

  contains

    subroutine open_loc_sr_grad_sampled_dump()
      implicit none
      character(16) :: rank_label

      write(rank_label,'(I6.6)') info%id_r
      file_grad_dump = trim(base_directory)//trim(sysname)//'_loc_sr_grad_sampled_rank'//trim(rank_label)//'.data'
      fh_grad_dump = open_filehandle(trim(file_grad_dump), status='replace')
      write(fh_grad_dump,'(a)') '# Local-SR sampled current-vs-legacy dV_sr/dG^2 dump'
      write(fh_grad_dump,'(a)') '# rank means icomm_r rank; only id_k=0 and id_o=0 write files'
      write(fh_grad_dump,'(a)') '# rank ix iy iz gx gy gz g_abs rho_e_re rho_e_im dVsr_dG2_current_re dVsr_dG2_current_im'
      write(fh_grad_dump,'(a)') '# dVsr_dG2_legacy_re dVsr_dG2_legacy_im'
      write(fh_grad_dump,'(a)') '# pressure_current_gpa pressure_legacy_gpa delta_pressure_gpa'
    end subroutine open_loc_sr_grad_sampled_dump

    subroutine dump_loc_sr_grad_sampled_points(ix, iy, iz, g, g_abs, rho_e, dVsr_dG2_current, dVsr_dG2_legacy, &
         pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa)
      implicit none
      integer, intent(in) :: ix, iy, iz
      real(8), intent(in) :: g(3), g_abs, pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa
      complex(8), intent(in) :: rho_e, dVsr_dG2_current, dVsr_dG2_legacy

      write(fh_grad_dump,*) info%id_r, ix, iy, iz, g(1), g(2), g(3), g_abs, &
           dble(rho_e), aimag(rho_e), dble(dVsr_dG2_current), aimag(dVsr_dG2_current), &
           dble(dVsr_dG2_legacy), aimag(dVsr_dG2_legacy), &
           pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa
    end subroutine dump_loc_sr_grad_sampled_points
  end subroutine calc_stress_loc

  subroutine calc_local_sr_gspace_diagnostics(system, pp, fg, info, ppg, cutoff_g2)
    use communication, only: comm_get_max_array1d_double, comm_summation
    use salmon_global, only: kion
    use structures
    implicit none
    type(s_dft_system),      intent(inout) :: system
    type(s_pp_info),         intent(in)    :: pp
    type(s_reciprocal_grid), intent(in)    :: fg
    type(s_parallel_info),   intent(in)    :: info
    type(s_pp_grid),         intent(in)    :: ppg
    real(8),                 intent(in)    :: cutoff_g2
    integer :: ix, iy, iz, ik, ia
    real(8) :: gvec(3), gmag, h, deriv_num, diff_abs
    real(8) :: diff_sum_local, diff_sum_global, count_local, count_global
    real(8) :: maxdiff_in(1), maxdiff_out(1), g_at_max_local, g_at_max_global, g_pick_local
    logical :: species_used

    diff_sum_local = 0d0
    count_local = 0d0
    maxdiff_in(1) = 0d0
    g_at_max_local = 0d0

    do ik = 1, size(pp%zps)
      species_used = .false.
      do ia = 1, system%nion
        if(kion(ia) == ik) then
          species_used = .true.
          exit
        end if
      end do
      if(.not. species_used) cycle

      do iz = lbound(fg%if_Gzero,3), ubound(fg%if_Gzero,3)
      do iy = lbound(fg%if_Gzero,2), ubound(fg%if_Gzero,2)
      do ix = lbound(fg%if_Gzero,1), ubound(fg%if_Gzero,1)
        if(fg%if_Gzero(ix,iy,iz)) cycle
        gvec(1) = fg%vec_G(1,ix,iy,iz)
        gvec(2) = fg%vec_G(2,ix,iy,iz)
        gvec(3) = fg%vec_G(3,ix,iy,iz)
        gmag = sqrt(gvec(1)**2 + gvec(2)**2 + gvec(3)**2)
        if(gmag**2 > cutoff_g2) cycle

        h = max(1d-4, 1d-2 * gmag)
        deriv_num = numerical_dvg_dg2_from_vg(pp, ik, gmag, h)
        diff_abs = abs(ppg%dVG_ion_dG2(ix,iy,iz,ik) - deriv_num)
        diff_sum_local = diff_sum_local + diff_abs
        count_local = count_local + 1d0
        if(diff_abs > maxdiff_in(1)) then
          maxdiff_in(1) = diff_abs
          g_at_max_local = gmag
        end if
      end do
      end do
      end do
    end do

    call comm_get_max_array1d_double(maxdiff_in, maxdiff_out, 1, info%icomm_r)
    call comm_summation(diff_sum_local, diff_sum_global, info%icomm_r)
    call comm_summation(count_local, count_global, info%icomm_r)
    g_pick_local = 0d0
    if(maxdiff_in(1) > 0d0 .and. abs(maxdiff_in(1) - maxdiff_out(1)) <= max(1d-12, 1d-10 * maxdiff_out(1))) then
      g_pick_local = g_at_max_local
    end if
    call comm_summation(g_pick_local, g_at_max_global, info%icomm_r)

    system%stress_loc_sr_gspace_dvg_maxdiff = maxdiff_out(1)
    if(count_global > 0d0) then
      system%stress_loc_sr_gspace_dvg_meandiff = diff_sum_global / count_global
    else
      system%stress_loc_sr_gspace_dvg_meandiff = 0d0
    end if
    system%stress_loc_sr_gspace_dvg_g_at_max = g_at_max_global
  end subroutine calc_local_sr_gspace_diagnostics

  real(8) function numerical_dvg_dg2_from_vg(pp, ik, gmag, h) result(dvg_num)
    use structures
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer,         intent(in) :: ik
    real(8),         intent(in) :: gmag, h
    real(8) :: g_lo, g_hi, v_lo, v_hi, denom

    g_lo = max(0d0, gmag - h)
    g_hi = gmag + h
    v_lo = eval_local_sr_vg_from_table(pp, ik, g_lo)
    v_hi = eval_local_sr_vg_from_table(pp, ik, g_hi)
    denom = g_hi**2 - g_lo**2
    if(denom <= 0d0) then
      dvg_num = 0d0
    else
      dvg_num = (v_hi - v_lo) / denom
    end if
  end function numerical_dvg_dg2_from_vg

  real(8) function eval_local_sr_vg_from_table(pp, ik, gmag) result(vg)
    use structures
    use salmon_global, only: npt_loc_sr_aux_2pi
    use prep_pp_sub, only: integrate_local_sr_vg_shell_stress
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer,         intent(in) :: ik
    real(8),         intent(in) :: gmag

    vg = integrate_local_sr_vg_shell_stress(pp, ik, gmag, npt_loc_sr_aux_2pi)
  end function eval_local_sr_vg_from_table

  subroutine calc_stress_nl(system, pp, info, mg, stencil, ppg, tpsi, energy, field_state)
    use structures
    use math_constants,     only: zi
    use stencil_sub,        only: calc_gradient_psi
    use nonlocal_potential, only: calc_uVpsi, calc_uVpsi_rdivided
    use communication,      only: comm_summation
    use salmon_global,      only: yn_out_stress_details, stress_l_decomp, kion
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
    integer :: ispec, nspecies
    real(8) :: rtmp, kAc(3), strs(3,3), strs_sum(3,3), V, nl_energy_part
    real(8) :: contrib
    complex(8) :: psi_r, w(3), r_uVpsi_b(3), uVpsi_ilma
    complex(8), allocatable :: gtpsi(:,:,:,:), uVpsibox(:,:,:,:,:), uVpsibox2(:,:,:,:,:)
    integer, allocatable :: ll_of_ilma(:), species_of_ilma(:)
    real(8), allocatable :: strs_species_l(:,:,:,:), strs_species_l_sum(:,:,:,:), e_nl_species_l(:,:), e_nl_species_l_sum(:,:)
    logical :: want_l_species

    V = system%det_a
    im = 1
    strs = 0d0
    want_l_species = (trim(stress_l_decomp) == 'species')

    if(want_l_species) then
      call build_nl_l_channel_map(pp, kion, ppg%nlma, ll_of_ilma, lmax_nl, species_of_ilma)
      nspecies = size(pp%zps)
      allocate(strs_species_l(1:nspecies,0:lmax_nl,3,3), strs_species_l_sum(1:nspecies,0:lmax_nl,3,3))
      allocate(e_nl_species_l(1:nspecies,0:lmax_nl), e_nl_species_l_sum(1:nspecies,0:lmax_nl))
      strs_species_l = 0d0
      e_nl_species_l = 0d0
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

      call accumulate_stress_nl_orbital(system, info, mg, ppg, tpsi, field_state, ik, io, ispin, im, gtpsi, &
           uVpsibox2, rtmp, want_l_species, ll_of_ilma, species_of_ilma, strs, strs_species_l, e_nl_species_l)
    end do
    end do
    end do

    deallocate(gtpsi)
    if(allocated(uVpsibox)) deallocate(uVpsibox)
    if(allocated(uVpsibox2)) deallocate(uVpsibox2)

    call comm_summation(strs, strs_sum, 9, info%icomm_rko)
    if(want_l_species) then
      call comm_summation(strs_species_l, strs_species_l_sum, size(strs_species_l), info%icomm_rko)
      call comm_summation(e_nl_species_l, e_nl_species_l_sum, size(e_nl_species_l), info%icomm_rko)
    end if
    strs_sum = strs_sum / V
    do a = 1, 3
      strs_sum(a,a) = strs_sum(a,a) + energy%E_ion_nloc / V
    end do
    system%stress_nl = -strs_sum

    if(want_l_species) then
      allocate(system%stress_nl_species_l(1:nspecies,0:lmax_nl,3,3))
      system%stress_nl_species_l = -strs_species_l_sum / V
      do ispec = 1, nspecies
        do ll = 0, lmax_nl
          do a = 1, 3
            system%stress_nl_species_l(ispec,ll,a,a) = system%stress_nl_species_l(ispec,ll,a,a) - e_nl_species_l_sum(ispec,ll) / V
          end do
        end do
      end do

      deallocate(strs_species_l, strs_species_l_sum, e_nl_species_l, e_nl_species_l_sum, ll_of_ilma, species_of_ilma)
    end if
  end subroutine calc_stress_nl

  subroutine accumulate_stress_nl_orbital(system, info, mg, ppg, tpsi, field_state, ik, io, ispin, im, gtpsi, &
                                          uVpsibox2, rtmp, want_l_species, ll_of_ilma, species_of_ilma, &
                                          strs, strs_species_l, e_nl_species_l)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dft_system),         intent(in)    :: system
    type(s_parallel_info),      intent(in)    :: info
    type(s_rgrid),              intent(in)    :: mg
    type(s_pp_grid),            intent(in)    :: ppg
    type(s_orbital),            intent(inout) :: tpsi
    type(s_stress_field_state), intent(in)    :: field_state
    integer,                    intent(in)    :: ik, io, ispin, im
    complex(8),                 intent(in)    :: gtpsi(3, mg%is_array(1):mg%ie_array(1), mg%is_array(2):mg%ie_array(2), mg%is_array(3):mg%ie_array(3))
    complex(8),                 intent(in)    :: uVpsibox2(system%nspin, info%io_s:info%io_e, info%ik_s:info%ik_e, info%im_s:info%im_e, ppg%nlma)
    real(8),                    intent(in)    :: rtmp
    logical,                    intent(in)    :: want_l_species
    integer,                    intent(in), optional :: ll_of_ilma(:), species_of_ilma(:)
    real(8),                    intent(inout) :: strs(3,3)
    real(8),                    intent(inout), optional :: strs_species_l(:,:,:,:), e_nl_species_l(:,:)
    integer :: ilma, ia, j, a, b, ix, iy, iz, ll, ispec
    integer :: nspecies_local, lmax_local
    real(8) :: contrib, nl_energy_part
    real(8) :: kAc(3), kAc_uniform(3)
    real(8) :: strs_thread(3,3)
    real(8), allocatable :: strs_species_thread(:,:,:,:), e_nl_species_thread(:,:)
    complex(8) :: psi_r, w_a, r_uVpsi_b(3), uVpsi_ilma, projector

    kAc_uniform = system%vec_k(:,ik) + field_state%Ac_uniform(:)
    if (want_l_species) then
      if (.not. present(ll_of_ilma) .or. .not. present(species_of_ilma) .or. .not. present(strs_species_l) .or. .not. present(e_nl_species_l)) then
        call fail_stress("species-l nonlocal stress requires l-channel maps")
      end if
    end if

!$omp parallel default(none) &
!$omp shared(system,mg,ppg,tpsi,field_state,ik,io,ispin,im,gtpsi,uVpsibox2,rtmp,want_l_species,ll_of_ilma,species_of_ilma) &
!$omp shared(strs,strs_species_l,e_nl_species_l,kAc_uniform) &
!$omp private(ilma,ia,j,a,b,ix,iy,iz,ll,ispec,contrib,nl_energy_part,kAc,psi_r,w_a,r_uVpsi_b,uVpsi_ilma,projector) &
!$omp private(strs_thread,strs_species_thread,e_nl_species_thread,nspecies_local,lmax_local)
    strs_thread = 0d0
    if (want_l_species) then
      nspecies_local = size(strs_species_l, 1)
      lmax_local = ubound(strs_species_l, 2)
      allocate(strs_species_thread(1:nspecies_local, 0:lmax_local, 3, 3))
      allocate(e_nl_species_thread(1:nspecies_local, 0:lmax_local))
      strs_species_thread = 0d0
      e_nl_species_thread = 0d0
    end if
!$omp do schedule(static)
    do ilma = 1, ppg%nlma
      ia = ppg%ia_tbl(ilma)
      uVpsi_ilma = uVpsibox2(ispin, io, ik, im, ilma)
      if (want_l_species) then
        ll = ll_of_ilma(ilma)
        ispec = species_of_ilma(ilma)
        nl_energy_part = rtmp * dble(conjg(uVpsi_ilma) * uVpsi_ilma) / ppg%rinv_uvu(ilma)
        e_nl_species_thread(ispec,ll) = e_nl_species_thread(ispec,ll) + nl_energy_part
      end if
      do a = 1, 3
        r_uVpsi_b = (0d0, 0d0)
        do j = 1, ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          if (field_state%use_micro_ac) then
            kAc(1) = system%vec_k(1,ik) + system%Ac_micro%v(1,ix,iy,iz)
            kAc(2) = system%vec_k(2,ik) + system%Ac_micro%v(2,ix,iy,iz)
            kAc(3) = system%vec_k(3,ik) + system%Ac_micro%v(3,ix,iy,iz)
          else
            kAc = kAc_uniform
          end if
          psi_r = tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          w_a = gtpsi(a,ix,iy,iz) + zi * kAc(a) * psi_r
          projector = conjg(ppg%zekr_uV(j,ilma,ik))
          do b = 1, 3
            r_uVpsi_b(b) = r_uVpsi_b(b) + ppg%rxyz(b,j,ia) * projector * w_a
          end do
        end do
        do b = 1, 3
          contrib = 2d0 * rtmp * dble(conjg(uVpsi_ilma) * r_uVpsi_b(b))
          strs_thread(a,b) = strs_thread(a,b) + contrib
          if (want_l_species) strs_species_thread(ispec,ll,a,b) = strs_species_thread(ispec,ll,a,b) + contrib
        end do
      end do
    end do
!$omp end do
!$omp critical(stress_nl_orbital_accum)
    strs = strs + strs_thread
    if (want_l_species) then
      strs_species_l = strs_species_l + strs_species_thread
      e_nl_species_l = e_nl_species_l + e_nl_species_thread
    end if
!$omp end critical(stress_nl_orbital_accum)
    if (want_l_species) then
      deallocate(strs_species_thread, e_nl_species_thread)
    end if
!$omp end parallel
  end subroutine accumulate_stress_nl_orbital

  subroutine build_nl_l_channel_map(pp, kion, nlma, ll_of_ilma, lmax_nl, species_of_ilma)
    use structures
    implicit none
    type(s_pp_info),      intent(in)  :: pp
    integer,              intent(in)  :: kion(:)
    integer,              intent(in)  :: nlma
    integer, allocatable, intent(out) :: ll_of_ilma(:)
    integer,              intent(out) :: lmax_nl
    integer, allocatable, intent(out), optional :: species_of_ilma(:)
    integer :: ia, ik, ll, l, l0, m, ilma

    lmax_nl = 0
    do ia = 1, size(kion)
      ik = kion(ia)
      lmax_nl = max(lmax_nl, pp%mlps(ik))
    end do

    allocate(ll_of_ilma(nlma))
    ll_of_ilma = -1
    if(present(species_of_ilma)) then
      allocate(species_of_ilma(nlma))
      species_of_ilma = 0
    end if

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
            if(present(species_of_ilma)) species_of_ilma(ilma) = ik
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
    real(8) :: cutoff_r2, exp_ewald_pref
    complex(8) :: SG

    V = system%det_a
    strs_G = 0d0
    strs_R = 0d0
    strs_G_grad = 0d0
    strs_G_diag = 0d0
    strs_G_self = 0d0
    E_ewa_G_loc = 0d0
    E_ewa_R_loc = 0d0
    cutoff_r2 = cutoff_r**2
    exp_ewald_pref = 2d0 * sqrt(aEwald/pi)
    Qtot = 0d0
    zps2 = 0d0
    do ia = 1, system%nion
      Qtot = Qtot + pp%zps(kion(ia))
      zps2 = zps2 + pp%zps(kion(ia))**2
    end do

    call get_g_bounds(fg, ig_s, ig_e)
    !$omp parallel do collapse(2) default(none) &
    !$omp   private(ix,iy,iz,ia,a,b,g,G2,SG,fact) &
    !$omp   shared(fg,ig_s,ig_e,system,pp,kion,V,aEwald) &
    !$omp   reduction(+:strs_G_grad,strs_G_diag,E_ewa_G_loc)
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
    !$omp end parallel do

    strs_G = strs_G_grad + strs_G_diag
    call sum_stress_tensor(info%icomm_r, strs_G, strs_G_sum)
    call sum_stress_tensor(info%icomm_r, strs_G_grad, strs_G_grad_sum)
    call sum_stress_tensor(info%icomm_r, strs_G_diag, strs_G_diag_sum)
    do a = 1, 3
      strs_G_sum(a,a) = strs_G_sum(a,a) + pi * Qtot**2 / (2d0*aEwald*V**2)
      strs_G_self(a,a) = pi * Qtot**2 / (2d0*aEwald*V**2)
    end do

    !$omp parallel do default(none) &
    !$omp   private(iia,ia,ipair,ib,a,b,r,rab,rr,r_abs,fact) &
    !$omp   shared(info,ewald,system,pp,kion,cutoff_r2,V,aEwald,exp_ewald_pref) &
    !$omp   reduction(+:strs_R,E_ewa_R_loc)
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
        if(rr > cutoff_r2) cycle
        r_abs = sqrt(rr)
        fact = -0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib)) / (V * r_abs**3) &
             * (erfc_salmon(sqrt(aEwald)*r_abs) + exp_ewald_pref*r_abs*exp(-aEwald*rr))
        E_ewa_R_loc = E_ewa_R_loc + 0.5d0 * pp%zps(kion(ia)) * pp%zps(kion(ib)) &
                   * erfc_salmon(sqrt(aEwald)*r_abs) / r_abs
        do b = 1, 3
        do a = 1, 3
          strs_R(a,b) = strs_R(a,b) + fact * rab(a) * rab(b)
        end do
        end do
      end do
    end do
    !$omp end parallel do

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
  subroutine calc_stress_loc_sr_rs(system, pp, info, mg, ppg, rho_s, srg_scalar)
    use structures
    use communication,  only: comm_summation
    use filesystem,     only: open_filehandle
    use salmon_global,  only: base_directory, kion, sysname, &
                              yn_out_loc_sr_rs_sampled_dump, yn_out_loc_sr_rs_subdiv_probe, &
                              num_loc_sr_rs_subdiv_probe, mode_loc_sr_rs_subdiv_probe_rho
    use inputoutput,    only: au_pressure_gpa
    use prep_pp_sub, only: eval_local_sr_shared_u_stress
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_grid),       intent(in)    :: ppg
    type(s_scalar),        intent(in)    :: rho_s(:)
    type(s_sendrecv_grid), intent(inout), optional :: srg_scalar
    integer  :: ia, ib, ik, j, a, b, ix, iy, iz, i1, i2, i3, intr_legacy, intr_stress, ispin, nspin, bin_idx
    integer  :: probe_subdiv_n
    integer  :: fh_sampled_dump
    real(8)  :: s(3), r_abs, r_inv, rho_val, dvsr_dr_current, dvsr_dr_legacy
    real(8)  :: u_r, du_r, r1, r2, r3
    real(8)  :: V, hvol
    real(8)  :: contrib, contrib_legacy, pressure_probe_gpa
    real(8)  :: probe_tensor(3,3)
    real(8)  :: strs(3,3), strs_sum(3,3), strs_probe(3,3), strs_probe_sum(3,3)
    real(8)  :: strs_bins(3,3,4), strs_bins_sum(3,3,4)
    real(8), pointer :: probe_rho_box(:,:,:)
    logical  :: dump_sampled, legacy_ok, probe_rs_subdiv, probe_rho_trilinear
    character(256) :: file_sampled_dump
    character(16) :: probe_rho_mode

    V = system%det_a
    hvol = system%Hvol
    nspin = system%nspin
    strs = 0d0
    strs_probe = 0d0
    strs_bins = 0d0
    nullify(probe_rho_box)
    probe_rs_subdiv = (yn_out_loc_sr_rs_subdiv_probe == 'y')
    probe_subdiv_n = 1
    if(probe_rs_subdiv) probe_subdiv_n = num_loc_sr_rs_subdiv_probe
    probe_rho_mode = trim(mode_loc_sr_rs_subdiv_probe_rho)
    probe_rho_trilinear = probe_rs_subdiv .and. probe_rho_mode == 'trilinear'
    if(probe_rho_trilinear) then
      call ensure_r2scan_stress_workspace(mg)
      probe_rho_box => stress_rho_box_work
      probe_rho_box = 0d0
      !$omp parallel do collapse(2) default(none) private(ix,iy,iz,ispin) shared(probe_rho_box,rho_s,mg,nspin)
      do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        do ispin = 1, nspin
          probe_rho_box(ix,iy,iz) = probe_rho_box(ix,iy,iz) + rho_s(ispin)%f(ix,iy,iz)
        end do
      end do
      end do
      end do
      !$omp end parallel do
      if(info%if_divide_rspace) then
        if(.not. present(srg_scalar)) call fail_stress("trilinear rho probe requires srg_scalar")
        call update_overlap_real8(srg_scalar, mg, probe_rho_box)
      end if
    end if
    dump_sampled = (yn_out_loc_sr_rs_sampled_dump == 'y' .and. info%id_k == 0 .and. info%id_o == 0)
    fh_sampled_dump = -1
    if(dump_sampled) call open_loc_sr_rs_sampled_dump()

    do ia = 1, system%nion
      ik = kion(ia)
      i1 = 0
      i2 = 0
      i3 = 0
      do ib = 1, pp%nrloc(ik)
        if(pp%rad(ib,ik) <= 0d0) cycle
        if(i1 == 0) then
          i1 = ib
        else if(i2 == 0) then
          i2 = ib
        else
          i3 = ib
          exit
        end if
      end do
      if(i1 == 0 .or. i2 == 0 .or. i3 == 0) stop 'calc_stress_loc_sr_rs: need three positive radial points'
      r1 = pp%rad(i1,ik)
      r2 = pp%rad(i2,ik)
      r3 = pp%rad(i3,ik)

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

        call find_radial_index(r_abs, pp%rad_u_sr_stress(:,ik), pp%nr_u_sr_stress(ik), intr_stress)
        if(intr_stress < 1 .or. intr_stress >= pp%nr_u_sr_stress(ik)) cycle

        call eval_local_sr_shared_u_stress(pp,ik,r_abs,u_r,du_r,intr_stress)
        ! V_sr(r) = u(r) / r, so V'_sr(r) = (u'(r) * r - u(r)) / r^2
        dvsr_dr_current = (du_r * r_abs - u_r) * r_inv**2
        contrib = -rho_val * dvsr_dr_current * r_inv * hvol
        if(probe_rs_subdiv) then
          call calc_loc_sr_rs_subdiv_probe_point(ik, ix, iy, iz, s, rho_val, probe_subdiv_n, probe_tensor, pressure_probe_gpa)
        else
          probe_tensor = 0d0
          do b = 1, 3
          do a = 1, 3
            probe_tensor(a,b) = contrib * s(a) * s(b)
          end do
          end do
          pressure_probe_gpa = contrib * r_abs * r_abs * au_pressure_gpa / (3d0 * V)
        end if

        bin_idx = 4
        if(r_abs < r1) then
          bin_idx = 1
        else if(r_abs < r2) then
          bin_idx = 2
        else if(r_abs < r3) then
          bin_idx = 3
        end if

        if(dump_sampled) then
          call find_radial_index(r_abs, pp%rad(:,ik), pp%nrloc(ik), intr_legacy)
          call legacy_loc_sr_dvsr_from_table_currentdump(pp, ik, r_abs, intr_legacy, &
               dvsr_dr_legacy, legacy_ok)
          if(legacy_ok) then
            contrib_legacy = -rho_val * dvsr_dr_legacy * r_inv * hvol
          else
            dvsr_dr_legacy = 0d0
            contrib_legacy = 0d0
          end if
          call dump_loc_sr_rs_sampled_points(ia, ik, j, intr_legacy, ix, iy, iz, s, bin_idx, r_abs, rho_val, &
               dvsr_dr_current, dvsr_dr_legacy, contrib, contrib_legacy, pressure_probe_gpa, legacy_ok)
        end if

        ! Accumulate: -rho * V'_sr * s_a * s_b / |s| * hvol
        do b = 1, 3
        do a = 1, 3
          strs(a,b) = strs(a,b) + contrib * s(a) * s(b)
          strs_probe(a,b) = strs_probe(a,b) + probe_tensor(a,b)
          strs_bins(a,b,bin_idx) = strs_bins(a,b,bin_idx) + contrib * s(a) * s(b)
        end do
        end do
      end do
    end do

    ! Use icomm_r (not icomm_rko): this is a density-only sum with no k/orbital loops.
    call comm_summation(strs, strs_sum, 9, info%icomm_r)
    call comm_summation(strs_probe, strs_probe_sum, 9, info%icomm_r)
    call comm_summation(strs_bins, strs_bins_sum, size(strs_bins), info%icomm_r)
    system%stress_loc_sr_rs = -strs_sum / V
    system%stress_loc_sr_rs_subdiv_probe = -strs_probe_sum / V
    do ib = 1, 4
      system%stress_loc_sr_rs_bins(:,:,ib) = -strs_bins_sum(:,:,ib) / V
    end do
    if(fh_sampled_dump >= 0) close(fh_sampled_dump)
    call calc_stress_loc_sr_rs_legacy_compare(system, pp, info, mg, ppg, rho_s)

  contains

    subroutine calc_loc_sr_rs_subdiv_probe_point(ik, ix_center, iy_center, iz_center, s_center, rho_center, probe_subdiv_n, tensor_probe, pressure_probe_gpa)
      implicit none
      integer, intent(in) :: ik, ix_center, iy_center, iz_center, probe_subdiv_n
      real(8), intent(in) :: s_center(3), rho_center
      real(8), intent(out) :: tensor_probe(3,3), pressure_probe_gpa
      integer :: isx, isy, isz, a, b
      real(8) :: fx, fy, fz
      real(8) :: du_shift, dv_shift, dw_shift
      real(8) :: s_probe(3), r_probe, r_probe_inv, r_probe2
      real(8) :: u_probe, du_probe, dvsr_dr_probe, contrib_probe_sub, rho_probe_sub
      real(8) :: subdiv_count_inv, subcell_weight

      tensor_probe = 0d0
      pressure_probe_gpa = 0d0
      subdiv_count_inv = 1d0 / dble(probe_subdiv_n)
      subcell_weight = hvol * subdiv_count_inv**3

      do isx = 1, probe_subdiv_n
        fx = (dble(isx) - 0.5d0) * subdiv_count_inv
        du_shift = (fx - 0.5d0) * system%hgs(1)
        do isy = 1, probe_subdiv_n
          fy = (dble(isy) - 0.5d0) * subdiv_count_inv
          dv_shift = (fy - 0.5d0) * system%hgs(2)
          do isz = 1, probe_subdiv_n
            fz = (dble(isz) - 0.5d0) * subdiv_count_inv
            dw_shift = (fz - 0.5d0) * system%hgs(3)
            s_probe(1) = s_center(1) + system%rmatrix_a(1,1) * du_shift + system%rmatrix_a(1,2) * dv_shift + &
                         system%rmatrix_a(1,3) * dw_shift
            s_probe(2) = s_center(2) + system%rmatrix_a(2,1) * du_shift + system%rmatrix_a(2,2) * dv_shift + &
                         system%rmatrix_a(2,3) * dw_shift
            s_probe(3) = s_center(3) + system%rmatrix_a(3,1) * du_shift + system%rmatrix_a(3,2) * dv_shift + &
                         system%rmatrix_a(3,3) * dw_shift
            r_probe2 = s_probe(1)**2 + s_probe(2)**2 + s_probe(3)**2
            if(r_probe2 < 1d-24) cycle
            r_probe = sqrt(r_probe2)
            if(r_probe >= pp%rps(ik)) cycle
            r_probe_inv = 1d0 / r_probe
            call eval_local_sr_shared_u_stress(pp, ik, r_probe, u_probe, du_probe)
            dvsr_dr_probe = (du_probe * r_probe - u_probe) * r_probe_inv**2
            rho_probe_sub = rho_center
            if(probe_rho_trilinear) rho_probe_sub = sample_probe_rho_trilinear(ix_center, iy_center, iz_center, fx, fy, fz)
            contrib_probe_sub = -rho_probe_sub * dvsr_dr_probe * r_probe_inv * subcell_weight
            do b = 1, 3
            do a = 1, 3
              tensor_probe(a,b) = tensor_probe(a,b) + contrib_probe_sub * s_probe(a) * s_probe(b)
            end do
            end do
            pressure_probe_gpa = pressure_probe_gpa + contrib_probe_sub * r_probe2 * au_pressure_gpa / (3d0 * V)
          end do
        end do
      end do
    end subroutine calc_loc_sr_rs_subdiv_probe_point

    real(8) function sample_probe_rho_trilinear(ix_center, iy_center, iz_center, fx, fy, fz)
      implicit none
      integer, intent(in) :: ix_center, iy_center, iz_center
      real(8), intent(in) :: fx, fy, fz
      integer :: ix0, ix1, iy0, iy1, iz0, iz1
      real(8) :: tx, ty, tz

      call select_probe_rho_axis(ix_center, fx, 1, ix0, ix1, tx)
      call select_probe_rho_axis(iy_center, fy, 2, iy0, iy1, ty)
      call select_probe_rho_axis(iz_center, fz, 3, iz0, iz1, tz)

      sample_probe_rho_trilinear = &
           (1d0 - tx) * ((1d0 - ty) * ((1d0 - tz) * probe_rho_box(ix0,iy0,iz0) + tz * probe_rho_box(ix0,iy0,iz1)) + &
                         ty         * ((1d0 - tz) * probe_rho_box(ix0,iy1,iz0) + tz * probe_rho_box(ix0,iy1,iz1))) + &
           tx         * ((1d0 - ty) * ((1d0 - tz) * probe_rho_box(ix1,iy0,iz0) + tz * probe_rho_box(ix1,iy0,iz1)) + &
                         ty         * ((1d0 - tz) * probe_rho_box(ix1,iy1,iz0) + tz * probe_rho_box(ix1,iy1,iz1)))
    end function sample_probe_rho_trilinear

    subroutine select_probe_rho_axis(ic, frac_cell, idir, i0, i1, t)
      implicit none
      integer, intent(in) :: ic, idir
      real(8), intent(in) :: frac_cell
      integer, intent(out) :: i0, i1
      real(8), intent(out) :: t
      integer :: ngrid_axis

      if(frac_cell < 0.5d0) then
        i0 = ic - 1
        i1 = ic
        t = frac_cell + 0.5d0
      else
        i0 = ic
        i1 = ic + 1
        t = frac_cell - 0.5d0
      end if

      if(.not. info%if_divide_rspace) then
        ngrid_axis = mg%num(idir)
        if(i0 < mg%is(idir)) i0 = i0 + ngrid_axis
        if(i0 > mg%ie(idir)) i0 = i0 - ngrid_axis
        if(i1 < mg%is(idir)) i1 = i1 + ngrid_axis
        if(i1 > mg%ie(idir)) i1 = i1 - ngrid_axis
      end if
    end subroutine select_probe_rho_axis

    subroutine open_loc_sr_rs_sampled_dump()
      implicit none
      character(16) :: rank_label
      write(rank_label,'(I6.6)') info%id_r
      file_sampled_dump = trim(base_directory)//trim(sysname)//'_loc_sr_rs_sampled_rank' &
                        // trim(rank_label)//'.data'
      fh_sampled_dump = open_filehandle(trim(file_sampled_dump), status='replace')
      write(fh_sampled_dump,'(a)') '# Local-SR sampled current-vs-legacy dV_sr/dr dump'
      write(fh_sampled_dump,'(a)') '# rank means icomm_r rank; only id_k=0 and id_o=0 write files'
      write(fh_sampled_dump,'(a)') '# rank ia ik j intr ix iy iz bin legacy_ok sx sy sz r_abs rho dvsr_dr_current dvsr_dr_legacy'
      write(fh_sampled_dump,'(a)') '# delta_dvsr_dr pressure_current_gpa pressure_legacy_gpa delta_pressure_gpa'
      if(probe_rs_subdiv) then
        write(fh_sampled_dump,'(a,i0)') '# probe_subdiv_n = ', probe_subdiv_n
        write(fh_sampled_dump,'(a,a)') '# probe_rho_mode = ', trim(probe_rho_mode)
        write(fh_sampled_dump,'(a)') '# pressure_probe_gpa delta_pressure_probe_current_gpa'
      end if
    end subroutine open_loc_sr_rs_sampled_dump

    subroutine dump_loc_sr_rs_sampled_points(ia, ik, j, intr, ix, iy, iz, s, bin_idx, r_abs, rho_val, &
         dvsr_dr_current, dvsr_dr_legacy, contrib_current, contrib_legacy, pressure_probe_gpa, legacy_ok)
      implicit none
      integer, intent(in) :: ia, ik, j, intr, ix, iy, iz, bin_idx
      real(8), intent(in) :: s(3), r_abs, rho_val, dvsr_dr_current, dvsr_dr_legacy
      real(8), intent(in) :: contrib_current, contrib_legacy, pressure_probe_gpa
      logical, intent(in) :: legacy_ok
      integer :: legacy_ok_int
      real(8) :: delta_dvsr_dr, pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa
      real(8) :: delta_pressure_probe_current_gpa

      delta_dvsr_dr = dvsr_dr_current - dvsr_dr_legacy
      pressure_current_gpa = contrib_current * r_abs * r_abs * au_pressure_gpa / (3d0 * V)
      pressure_legacy_gpa = contrib_legacy * r_abs * r_abs * au_pressure_gpa / (3d0 * V)
      delta_pressure_gpa = pressure_current_gpa - pressure_legacy_gpa
      delta_pressure_probe_current_gpa = pressure_probe_gpa - pressure_current_gpa
      legacy_ok_int = 0
      if(legacy_ok) legacy_ok_int = 1

      if(probe_rs_subdiv) then
        write(fh_sampled_dump,*) &
             info%id_r, ia, ik, j, intr, ix, iy, iz, bin_idx, &
             legacy_ok_int, s(1), s(2), s(3), r_abs, rho_val, dvsr_dr_current, dvsr_dr_legacy, delta_dvsr_dr, &
             pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa, pressure_probe_gpa, delta_pressure_probe_current_gpa
      else
        write(fh_sampled_dump,*) &
             info%id_r, ia, ik, j, intr, ix, iy, iz, bin_idx, &
             legacy_ok_int, s(1), s(2), s(3), r_abs, rho_val, dvsr_dr_current, dvsr_dr_legacy, delta_dvsr_dr, &
             pressure_current_gpa, pressure_legacy_gpa, delta_pressure_gpa
      end if
    end subroutine dump_loc_sr_rs_sampled_points

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

    subroutine legacy_loc_sr_dvsr_from_table_currentdump(pp, ik, r, intr, dvsr_dr_legacy, ok)
      implicit none
      type(s_pp_info), intent(in) :: pp
      integer, intent(in) :: ik, intr
      real(8), intent(in) :: r
      real(8), intent(out) :: dvsr_dr_legacy
      logical, intent(out) :: ok
      real(8) :: r_lo, r_hi, dv_lo, dv_hi, ratio, dvloc_r

      ok = .false.
      dvsr_dr_legacy = 0d0
      if(intr < 2 .or. intr >= pp%nrloc(ik)) return

      r_lo = pp%rad(intr, ik)
      r_hi = pp%rad(intr+1, ik)
      if(r_hi <= r_lo) return

      dv_lo = pp%dvloctbl(intr, ik)
      dv_hi = pp%dvloctbl(intr+1, ik)
      ratio = (r - r_lo) / (r_hi - r_lo)
      dvloc_r = dv_lo + ratio * (dv_hi - dv_lo)

      dvsr_dr_legacy = dvloc_r - dble(pp%zps(ik)) / (r * r)
      ok = .true.
    end subroutine legacy_loc_sr_dvsr_from_table_currentdump

  end subroutine calc_stress_loc_sr_rs

  subroutine calc_stress_loc_sr_rs_legacy_compare(system, pp, info, mg, ppg, rho_s)
    use structures
    use communication,  only: comm_summation
    use salmon_global,  only: kion
    implicit none
    type(s_dft_system),    intent(inout) :: system
    type(s_pp_info),       intent(in)    :: pp
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg
    type(s_pp_grid),       intent(in)    :: ppg
    type(s_scalar),        intent(in)    :: rho_s(:)
    integer  :: ia, ib, ik, j, a, b, ix, iy, iz, i1, i2, i3, intr, ispin, nspin, bin_idx
    real(8)  :: s(3), r_abs, r_inv, rho_val, dvsr_dr_legacy, r1, r2, r3
    real(8)  :: V, hvol, contrib_legacy
    real(8)  :: strs_legacy(3,3), strs_legacy_sum(3,3), strs_legacy_bins(3,3,4), strs_legacy_bins_sum(3,3,4)
    logical  :: legacy_ok

    V = system%det_a
    hvol = system%Hvol
    nspin = system%nspin
    strs_legacy = 0d0
    strs_legacy_bins = 0d0

    do ia = 1, system%nion
      ik = kion(ia)
      i1 = 0
      i2 = 0
      i3 = 0
      do ib = 1, pp%nrloc(ik)
        if(pp%rad(ib,ik) <= 0d0) cycle
        if(i1 == 0) then
          i1 = ib
        else if(i2 == 0) then
          i2 = ib
        else
          i3 = ib
          exit
        end if
      end do
      if(i1 == 0 .or. i2 == 0 .or. i3 == 0) stop 'calc_stress_loc_sr_rs_legacy_compare: need three positive radial points'
      r1 = pp%rad(i1,ik)
      r2 = pp%rad(i2,ik)
      r3 = pp%rad(i3,ik)

      do j = 1, ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        s(1) = ppg%rxyz(1,j,ia)
        s(2) = ppg%rxyz(2,j,ia)
        s(3) = ppg%rxyz(3,j,ia)
        r_abs = sqrt(s(1)**2 + s(2)**2 + s(3)**2)
        if(r_abs < 1d-12) cycle
        r_inv = 1d0 / r_abs

        rho_val = 0d0
        do ispin = 1, nspin
          rho_val = rho_val + rho_s(ispin)%f(ix,iy,iz)
        end do

        call find_radial_index(r_abs, pp%rad(:,ik), pp%nrloc(ik), intr)
        call legacy_loc_sr_dvsr_from_table(pp, ik, r_abs, intr, dvsr_dr_legacy, legacy_ok)
        if(.not. legacy_ok) cycle
        contrib_legacy = -rho_val * dvsr_dr_legacy * r_inv * hvol

        bin_idx = 4
        if(r_abs < r1) then
          bin_idx = 1
        else if(r_abs < r2) then
          bin_idx = 2
        else if(r_abs < r3) then
          bin_idx = 3
        end if

        do b = 1, 3
        do a = 1, 3
          strs_legacy(a,b) = strs_legacy(a,b) + contrib_legacy * s(a) * s(b)
          strs_legacy_bins(a,b,bin_idx) = strs_legacy_bins(a,b,bin_idx) + contrib_legacy * s(a) * s(b)
        end do
        end do
      end do
    end do

    call comm_summation(strs_legacy, strs_legacy_sum, 9, info%icomm_r)
    call comm_summation(strs_legacy_bins, strs_legacy_bins_sum, size(strs_legacy_bins), info%icomm_r)
    system%stress_loc_sr_rs_legacy = -strs_legacy_sum / V
    do ib = 1, 4
      system%stress_loc_sr_rs_legacy_bins(:,:,ib) = -strs_legacy_bins_sum(:,:,ib) / V
    end do

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

    subroutine legacy_loc_sr_dvsr_from_table(pp, ik, r, intr, dvsr_dr_legacy, ok)
      implicit none
      type(s_pp_info), intent(in) :: pp
      integer, intent(in) :: ik, intr
      real(8), intent(in) :: r
      real(8), intent(out) :: dvsr_dr_legacy
      logical, intent(out) :: ok
      real(8) :: r_lo, r_hi, dv_lo, dv_hi, ratio, dvloc_r

      ok = .false.
      dvsr_dr_legacy = 0d0
      if(intr < 2 .or. intr >= pp%nrloc(ik)) return

      r_lo = pp%rad(intr, ik)
      r_hi = pp%rad(intr+1, ik)
      if(r_hi <= r_lo) return

      dv_lo = pp%dvloctbl(intr, ik)
      dv_hi = pp%dvloctbl(intr+1, ik)
      ratio = (r - r_lo) / (r_hi - r_lo)
      dvloc_r = dv_lo + ratio * (dv_hi - dv_lo)

      dvsr_dr_legacy = dvloc_r - dble(pp%zps(ik)) / (r * r)
      ok = .true.
    end subroutine legacy_loc_sr_dvsr_from_table

  end subroutine calc_stress_loc_sr_rs_legacy_compare

end module stress_sub
