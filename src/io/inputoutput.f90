!
!  Copyright 2017-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

module inputoutput
  use phys_constants
  use salmon_global
  implicit none
!Physical constant
  real(8),parameter :: au_time_fs = 0.02418884326505d0
  real(8),parameter :: au_energy_ev = 27.21138505d0
  real(8),parameter :: au_length_aa = 0.52917721067d0


  integer :: fh_variables_log
  integer :: fh_namelist
!  integer, parameter :: fh_atomic_spiecies = 902
  integer :: fh_atomic_coor
  integer :: fh_reentrance
  integer :: fh_atomic_red_coor
  logical :: if_nml_coor, if_nml_red_coor

  integer :: inml_calculation
  integer :: inml_control
  integer :: inml_units
  integer :: inml_parallel
  integer :: inml_system
  integer :: inml_pseudo
  integer :: inml_functional
  integer :: inml_rgrid
  integer :: inml_kgrid
  integer :: inml_tgrid
  integer :: inml_propagation
  integer :: inml_scf
  integer :: inml_emfield
  integer :: inml_multiscale
  integer :: inml_maxwell
  integer :: inml_analysis
  integer :: inml_poisson
  integer :: inml_ewald
  integer :: inml_opt
  integer :: inml_md
  integer :: inml_group_fundamental
  integer :: inml_group_hartree
  integer :: inml_group_others
  integer :: inml_code

!Input/Output units
  integer :: iflag_unit_time
  integer,parameter :: ntype_unit_time_au = 0
  integer,parameter :: ntype_unit_time_fs = 1
  real(8) :: utime_to_au, utime_from_au

  integer :: iflag_unit_length
  integer,parameter :: ntype_unit_length_au = 0
  integer,parameter :: ntype_unit_length_aa = 1
  real(8) :: ulength_to_au, ulength_from_au

  integer :: iflag_unit_energy
  integer,parameter :: ntype_unit_energy_au = 0
  integer,parameter :: ntype_unit_energy_ev = 1
  real(8) :: uenergy_to_au, uenergy_from_au

  integer :: iflag_unit_charge
  integer,parameter :: ntype_unit_charge_au = 0
  real(8) :: ucharge_to_au, ucharge_from_au



  type unit_t
     character(32) :: name
     real(8)       :: conv
  end type unit_t

  type(unit_t) :: t_unit_length
  type(unit_t) :: t_unit_length_inv
  type(unit_t) :: t_unit_energy
  type(unit_t) :: t_unit_energy_inv
  type(unit_t) :: t_unit_time
  type(unit_t) :: t_unit_time_inv
  type(unit_t) :: t_unit_spectrum_dipole
  type(unit_t) :: t_unit_spectrum_dipole_square
  type(unit_t) :: t_unit_current
  type(unit_t) :: t_unit_spectrum_current
  type(unit_t) :: t_unit_spectrum_current_square
  type(unit_t) :: t_unit_ac
  type(unit_t) :: t_unit_elec
  type(unit_t) :: t_unit_spectrum_elec
  type(unit_t) :: t_unit_spectrum_elec_square
  type(unit_t) :: t_unit_polarizability
  type(unit_t) :: t_unit_conductivity

contains
  subroutine read_input
    implicit none

    call read_stdin
    call read_input_common ! Should be renamed properly later
    call read_atomic_coordinates
    call dump_input_common ! Should be renamed properly later
    call check_bad_input

  end subroutine read_input


  subroutine read_stdin
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: get_filehandle
    implicit none
    integer :: cur
    integer :: ret = 0
    character(256) :: buff, text

    if (comm_is_root(nproc_id_global)) then

      fh_namelist = get_filehandle()
      open(fh_namelist, file='.namelist.tmp', status='replace')
     !open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')

      if_nml_coor =.false.
      fh_atomic_coor = get_filehandle()
      open(fh_atomic_coor, file='.atomic_coor.tmp', status='replace')

      if_nml_red_coor = .false.
      fh_atomic_red_coor = get_filehandle()
      open(fh_atomic_red_coor, file='.atomic_red_coor.tmp', status='replace')

      fh_reentrance = get_filehandle()
      open(fh_reentrance, file='.reenetrance.tmp', status='replace')

      cur = fh_namelist
      do while (.true.)
        read(*, '(a)', iostat=ret) buff
        if (ret < 0) then
          exit
        else
          text = trim(adjustl(buff))
          ! Comment lines
          if (text(1:1) == '!') cycle
         !! Beginning of 'atomic_species' part
         !if (text == '&atomic_spiecies') then
         !  cur = fh_atomic_spiecies
         !  cycle
         !end if
          ! Beginning of 'atomic_positions' part
          if (text == '&atomic_coor') then
            cur = fh_atomic_coor
            if_nml_coor =.true.
            cycle
          end if
          if (text == '&atomic_red_coor') then
            cur = fh_atomic_red_coor
            if_nml_red_coor = .true.
            cycle
          end if
          ! Beginning of 'atomic_species' part
          if (text == '&reentrance') then
            cur = fh_reentrance
            cycle
          end if
          ! End of 'atomic_(spiecies|positions)' part
          if ((text == '/') .and. (cur /= fh_namelist)) then
            cur = fh_namelist
            cycle
          end if

          write(cur, '(a)') trim(text)
        end if
      end do

      close(fh_namelist)
      close(fh_atomic_coor)
      close(fh_atomic_red_coor)
     !close(fh_atomic_spiecies)
      close(fh_reentrance)
    end if

   !call comm_sync_all()

    return
  end subroutine read_stdin

  subroutine read_input_common
    use parallelization
    use communication
    use filesystem, only: get_filehandle
    use misc_routines, only: string_lowercase
    implicit none
    integer :: ii
    real(8) :: norm

    namelist/calculation/ &
      & theory, &
      & yn_md,  &
      & yn_opt

    namelist/control/ &
      & sysname, &
      & base_directory, &
      & output_buffer_interval, &
      & yn_restart, &
      & directory_read_data, &
      & yn_self_checkpoint,  &
      & checkpoint_interval, &
      & yn_reset_step_restart, &
      & read_gs_restart_data,  &
      & write_gs_restart_data, &
      & time_shutdown,         &
      & method_wf_distributor, &
      & nblock_wf_distribute,  &
      & yn_gbp,            &  !should change name or temporary
      & yn_gbp_fourier0,   &  !should change name or temporary
      & read_gs_dns_cube,  &  !remove later (but this is used currently)
      & write_gs_wfn_k,    &  !remove later (but this is used currently)
      & write_rt_wfn_k        !remove later (but this is used currently)

    namelist/units/ &
      & unit_system

    namelist/parallel/ &
      & nproc_k, &
      & nproc_ob, &
      & nproc_rgrid, &
      & yn_ffte, &
      & yn_scalapack, &
      & yn_scalapack_red_mem, &
      & yn_eigenexa, &
      & process_allocation

    namelist/system/ &
      & yn_periodic, &
      & ispin, &
      & al, &
      & al_vec1,al_vec2,al_vec3, &
      & isym, &
      & crystal_structure, &
      & nstate, &
      & nstate_spin, &
      & nelec, &
      & nelec_spin, &
      & temperature, &
      & temperature_k, &
      & nelem, &
      & natom, &
      & file_atom_coor, &
      & file_atom_red_coor

    namelist/pseudo/ &
      & file_pseudo, &
      & lmax_ps, &
      & lloc_ps, &
      & izatom, &
      & yn_psmask, &
      & alpha_mask, &
      & gamma_mask, &
      & eta_mask

    namelist/functional/ &
      & xc, &
      & cname, &
      & xname, &
#ifdef USE_LIBXC
      & alibx, &
      & alibc, &
      & alibxc, &
#endif
      & cval

    namelist/rgrid/ &
      & dl, &
      & num_rgrid

    namelist/kgrid/ &
      & num_kgrid, &
      & file_kw

    namelist/tgrid/ &
      & nt, &
      & dt, &
      & gram_schmidt_interval

    namelist/propagation/ &
      & n_hamil, &
      & propagator, &
      & yn_fix_func

    namelist/scf/ &
      & method_init_wf, &
      & iseed_number_change, &
      & method_min, &
      & ncg, &
      & ncg_init, &
      & method_mixing, &
      & mixrate, &
      & nmemory_mb, &
      & alpha_mb, &
      & nmemory_p, &
      & beta_p, &
      & yn_auto_mixing, &
      & update_mixing_ratio, &
      & nscf, &
      & yn_subspace_diagonalization, &
      & convergence, &
      & threshold, &
      & iditer_notemperature, &
      & step_initial_mix_zero, &
      & conv_gap_mix_zero

    namelist/emfield/ &
      & trans_longi, &
      & ae_shape1, &
      & file_input1, &
      & e_impulse, &
      & E_amplitude1, &
      & I_wcm2_1, &
      & tw1, &
      & omega1, &
      & epdir_re1, &
      & epdir_im1, &
      & phi_cep1, &
      & ae_shape2, &
      & E_amplitude2, &
      & I_wcm2_2, &
      & tw2, &
      & omega2, &
      & epdir_re2, &
      & epdir_im2, &
      & phi_cep2, &
      & t1_t2, &
      & t1_start, &
      & num_dipole_source , &
      & vec_dipole_source , &
      & cood_dipole_source , &
      & rad_dipole_source, &
      & cutoff_G2_emfield

    namelist/multiscale/ &
      & fdtddim, &
      & twod_shape, &
      & nx_m, &
      & ny_m, &
      & nz_m, &
      & hx_m, &
      & hy_m, &
      & hz_m, &
      & nksplit, &
      & nxysplit, &
      & nxvacl_m, &
      & nxvacr_m, &
      & nx_origin_m, &
      & ny_origin_m, &
      & nz_origin_m, &
      & file_macropoint, &
      & num_macropoint,  &
      & set_ini_coor_vel,&
      & nmacro_write_group

    namelist/maxwell/    &
      & al_em,           &
      & dl_em,           &
      & dt_em,           &
      & nt_em,           &
      & boundary_em,     &
      & shape_file,      &
      & media_num,       &
      & media_type,      &
      & epsilon_em,      &
      & mu_em,           &
      & sigma_em,        &
      & pole_num_ld,     &
      & omega_p_ld,      &
      & f_ld,            &
      & gamma_ld,        &
      & omega_ld,        &
      & wave_input,      &
      & ek_dir1,         &
      & source_loc1,     &
      & ek_dir2,         &
      & source_loc2,     &
      & obs_num_em,      &
      & obs_samp_em,     &
      & obs_loc_em,      &
      & yn_obs_plane_em, &
      & yn_wf_em

    namelist/analysis/ &
      & projection_option, &
      & out_projection_step, &
      & nenergy, &
      & de, &
      & yn_out_psi, &
      & yn_out_dos, &
      & yn_out_dos_set_fe_origin, &
      & out_dos_start, &
      & out_dos_end, &
      & out_dos_nenergy, &
      & out_dos_width, &
      & out_dos_function, &
      & yn_out_pdos, &
      & yn_out_dns, &
      & yn_out_dns_rt, &
      & yn_out_dns_ac_je, &
      & out_dns_rt_step, &
      & out_dns_ac_je_step, &
      & out_old_dns, &
      & yn_out_dns_trans, &
      & out_dns_trans_energy, &
      & yn_out_elf, &
      & yn_out_elf_rt, &
      & out_elf_rt_step, &
      & yn_out_estatic_rt, &
      & out_estatic_rt_step, &
      & yn_out_rvf_rt, &
      & out_rvf_rt_step, &
      & yn_out_tm, &
      & out_ms_step, &
      & format_voxel_data, &
      & nsplit_voxel_data, &
      & timer_process

    namelist/poisson/ &
      & layout_multipole, &
      & num_multipole_xyz, &
      & threshold_cg

    namelist/ewald/ &
      & newald, &
      & aewald, &
      & cutoff_r, &
      & cutoff_r_buff, &
      & cutoff_g

    namelist/opt/ &
      & nopt, &
      & max_step_len_adjust, &
      & convrg_opt_fmax

    namelist/md/ &
      & ensemble, &
      & thermostat, &
      & step_velocity_scaling, &
      & step_update_ps, &
      & temperature0_ion_k, &
      & yn_set_ini_velocity, &
      & file_ini_velocity, &
      & thermostat_tau, &
      & yn_stop_system_momt

    namelist/group_fundamental/ &
      & iditer_nosubspace_diag, &
      & ntmg, &
      & idisnum, &
      & iwrite_projection, &
      & itwproj, &
      & iwrite_projnum, &
      & itcalc_ene

    namelist/group_hartree/ &
      & lmax_lmp

    namelist/group_others/ &
      & iswitch_orbital_mesh, &
      & iflag_psicube, &
      & num_projection, &
      & iwrite_projection_ob, &
      & iwrite_projection_k, &
      & filename_pot, &
      & iwrite_external, &
      & iflag_intelectron, &
      & num_dip2, &
      & dip2boundary, &
      & dip2center, &
      & itotntime2, &
      & iwdenoption, &
      & iwdenstep, &
      & iflag_estatic

    namelist/code/ &
      & yn_want_stencil_hand_vectorization, &
      & yn_want_communication_overlapping, &
      & stencil_openmp_mode, &
      & current_openmp_mode, &
      & force_openmp_mode

!! == default for &unit ==
    unit_system='au'
!! =======================

    if (comm_is_root(nproc_id_global)) then
      open(fh_namelist, file='.namelist.tmp', status='old')
      read(fh_namelist, nml=units, iostat=inml_units)
      rewind(fh_namelist)
      close(fh_namelist)
    end if

    call comm_bcast(unit_system,nproc_group_global)

    select case(unit_system)
    case('au','a.u.','A_eV_fs')
      continue
    case default
      stop 'invalid unit_system'
    end select

    select case(unit_system)
    case('au','a.u.')
      unit_time='au'
      unit_length='au'
      unit_energy='au'
      unit_charge='au'
    case('A_eV_fs')
      unit_time='fs'
      unit_length='AA'
      unit_energy='eV'
      unit_charge='au'
    end select

    call initialize_inputoutput_units

!! == default for &calculation
    theory              = 'tddft'
    yn_md               = 'n'
    yn_opt              = 'n'
!! == default for &control
    sysname               = 'default'
    base_directory        = './'
    output_buffer_interval= -1
    yn_restart            = 'n'
    directory_read_data   = 'restart/'
    yn_self_checkpoint    = 'n'
    checkpoint_interval   = 0
    yn_reset_step_restart = 'n'
    read_gs_restart_data  = 'all'
    write_gs_restart_data = 'all'
    time_shutdown         = -1d0
    method_wf_distributor = 'single'
    nblock_wf_distribute = 16
    yn_gbp        = 'n'
    yn_gbp_fourier0 = 'n' ! temporary
    !remove later
    read_gs_dns_cube = 'n'
    write_gs_wfn_k   = 'n'
    write_rt_wfn_k   = 'n'
!! == default for &parallel
    nproc_k              = 0
    nproc_ob             = 0
    nproc_rgrid          = 0
    yn_ffte              = 'n'
    yn_scalapack         = 'n'
    yn_scalapack_red_mem = 'n'
    yn_eigenexa          = 'n'
    process_allocation   = 'grid_sequential'
!! == default for &system
    yn_periodic        = 'n'
    ispin              = 0
    al                 = 0d0
    al_vec1            = 0d0
    al_vec2            = 0d0
    al_vec3            = 0d0
    isym               = 1
    crystal_structure  = 'none'
    nstate             = 0
    nstate_spin(:)     = 0
    nelec              = 0
    nelec_spin (:)     = 0
    temperature        = -1d0
    temperature_k      = -1d0
    nelem              = 0
    natom              = 0
    file_atom_coor     = 'none'
    file_atom_red_coor = 'none'
!! == default for &pseudo
    file_pseudo = 'none'
    lmax_ps     = -1
    lloc_ps     = -1
    izatom      = -1
    yn_psmask   = 'n'
    alpha_mask  = 0.8d0
    gamma_mask  = 1.8d0
    eta_mask    = 15d0
!! == default for &functional
    xc    = 'none'
    ! xcname = 'PZ'
    xname = 'none'
    cname = 'none'
    alibx = 'none'
    alibc = 'none'
    alibxc= 'none'
    cval  = -1d0
!! == default for &rgrid
    dl        = 0d0
    num_rgrid = 0
!! == default for &kgrid
    num_kgrid = 1
    file_kw   = 'none'
!! == default for &tgrid
    nt = 0
    dt = 0
    gram_schmidt_interval = 0
!! == default for &propagation
    n_hamil     = 4
    propagator  = 'middlepoint'
    yn_fix_func = 'n'
!! == default for &scf
    method_init_wf = 'gauss'
    iseed_number_change  =  0
    method_min    = 'cg'
    ncg           = 4
    ncg_init      = 4
    method_mixing = 'broyden'
    mixrate       = 0.5d0
    nmemory_mb    = 8
    alpha_mb      = 0.75d0
    nmemory_p     = 4
    beta_p        = 0.75d0
    yn_auto_mixing = 'n'
    update_mixing_ratio = 3.d0
    nscf          = 0
    yn_subspace_diagonalization = 'y'
    convergence   = 'rho_dne'
    threshold     = -1d0  !a.u. (default value for 'rho_dne'is given later)
    iditer_notemperature = 10
    step_initial_mix_zero= -1
    conv_gap_mix_zero    = 99999d0*uenergy_from_au

!! == default for &emfield
    trans_longi    = 'tr'
    ae_shape1      = 'none'
    file_input1    = ''
    e_impulse      = 1d-2*uenergy_from_au/ulength_from_au*utime_from_au ! a.u.
    E_amplitude1   = 0d0
    I_wcm2_1       = -1d0
    tw1            = 0d0
    omega1         = 0d0
    epdir_re1      = (/1d0,0d0,0d0/)
    epdir_im1      = 0d0
    phi_cep1       = 0d0
    ae_shape2      = 'none'
    E_amplitude2   = 0d0
    I_wcm2_2       = -1d0
    tw2            = 0d0
    omega2         = 0d0
    epdir_re2      = (/1d0,0d0,0d0/)
    epdir_im2      = 0d0
    phi_cep2       = 0d0
    t1_t2          = 0d0
    t1_start       = 0d0
    num_dipole_source  = 0
    vec_dipole_source  = 0d0
    cood_dipole_source = 0d0
    rad_dipole_source  = 2d0 ! a.u.
    cutoff_G2_emfield  = -1d0

!! == default for &multiscale
    fdtddim    = '1d'
    twod_shape = 'periodic'
    nx_m       = 1
    ny_m       = 1
    nz_m       = 1
    hx_m       = 0d0
    hy_m       = 0d0
    hz_m       = 0d0
    nksplit    = 0
    nxysplit   = 0
    nxvacl_m   = 1
    nxvacr_m   = 1
    nx_origin_m = 1
    ny_origin_m = 1
    nz_origin_m = 1
    file_macropoint = ''
    set_ini_coor_vel= 'n'
    nmacro_write_group= -1

!! == default for &maxwell
    al_em(:)           = 0d0
    dl_em(:)           = 0d0
    dt_em              = 0d0
    nt_em              = 0
    boundary_em(:,:)   = 'default'
    shape_file         = 'none'
    media_num          = 0
    media_type(:)      = 'vacuum'
    epsilon_em(:)      = 1d0
    mu_em(:)           = 1d0
    sigma_em(:)        = 0d0
    pole_num_ld(:)     = 1
    omega_p_ld(:)      = 0d0
    f_ld(:,:)          = 0d0
    gamma_ld(:,:)      = 0d0
    omega_ld(:,:)      = 0d0
    wave_input         = 'none'
    ek_dir1(:)         = 0d0
    source_loc1(:)     = 0d0
    ek_dir2(:)         = 0d0
    source_loc2(:)     = 0d0
    obs_num_em         = 0
    obs_samp_em        = 1
    obs_loc_em(:,:)    = 0d0
    yn_obs_plane_em(:) = 'n'
    yn_wf_em           = 'y'

!! == default for &analysis
    projection_option   = 'no'
    out_projection_step = 100
    nenergy             = 1000
    de                  = (0.01d0/au_energy_ev)*uenergy_from_au  ! eV
    yn_out_psi          = 'n'
    yn_out_dos          = 'n'
    yn_out_dos_set_fe_origin = 'n'
    out_dos_start       = -1.d10 / au_energy_ev * uenergy_from_au
    out_dos_end         = +1.d10 / au_energy_ev * uenergy_from_au
    out_dos_nenergy     = 601
    out_dos_width       = 0.1d0 / au_energy_ev * uenergy_from_au
    out_dos_function    = 'gaussian'
    yn_out_pdos         = 'n'
    yn_out_dns          = 'n'
    yn_out_dns_rt       = 'n'
    yn_out_dns_ac_je    = 'n'
    out_dns_rt_step     = 50
    out_dns_ac_je_step  = 50
    out_old_dns         = 'n'
    yn_out_dns_trans    = 'n'
    out_dns_trans_energy= 1.55d0 / au_energy_ev * uenergy_from_au  ! eV

    yn_out_elf          = 'n'
    yn_out_elf_rt       = 'n'
    out_elf_rt_step     = 50
    yn_out_estatic_rt   = 'n'
    out_estatic_rt_step = 50
    yn_out_rvf_rt       = 'n'
    out_rvf_rt_step     = 10
    yn_out_tm           = 'n'
    out_ms_step         = 100
    format_voxel_data   = 'cube'
    nsplit_voxel_data   = 1
    timer_process       = 'n'

!! == default for &poisson
    layout_multipole  = 3
    num_multipole_xyz = 0
    threshold_cg      = 1.d-15*uenergy_from_au**2*ulength_from_au**3 ! a.u., 1.d-15 a.u. = ! 1.10d-13 eV**2*AA**3
!! == default for &ewald
    newald = 4
    aewald = 0.5d0
    cutoff_r = -1d0
    cutoff_r_buff =  2.0d0
    cutoff_g = -1d0
!! == default for &opt
    nopt                = 100
    max_step_len_adjust =  -1d0 ![au] (no adjust if negative number)
    convrg_opt_fmax     =  1d-3
!! == default for &md
    ensemble              = 'nve'
    thermostat            = 'nose-hoover'
    step_velocity_scaling = -1
    step_update_ps        = 10
    temperature0_ion_k    = 298.15d0
    yn_set_ini_velocity   = 'n'
    file_ini_velocity     = 'none'
    thermostat_tau        =  41.34d0/utime_to_au  !=1[fs]: test value
    yn_stop_system_momt   = 'n'
!! == default for &group_fundamental
    iditer_nosubspace_diag = 10
    ntmg                   = 1
    idisnum                = (/1,2/)
    iwrite_projection      = 0
    itwproj                = -1
    iwrite_projnum         = 0
    itcalc_ene             = 10
!! == default for &group_hartree
    lmax_lmp = 4
!! == default for &group_others
    iswitch_orbital_mesh = 0
    iflag_psicube        = 0
    num_projection       = 1
    do ii=1,200
      iwrite_projection_ob(ii) = ii
    end do
    iwrite_projection_k(1:200) = 1
    filename_pot               = 'pot'
    iwrite_external            = 0
    iflag_intelectron          = 0
    num_dip2                   = 1
    dip2boundary(1:100)        = 0.d0*ulength_from_au ! a.u.
    dip2center(1:100)          = 0.d0*ulength_from_au ! a.u.
    itotntime2                 = 0
    iwdenoption                = 0
    iwdenstep                  = 0
    iflag_estatic              = 0
!! == default for code
    yn_want_stencil_hand_vectorization = 'y'
    yn_want_communication_overlapping  = 'n'
    stencil_openmp_mode = 'auto'
    current_openmp_mode = 'auto'
    force_openmp_mode   = 'auto'

    if (comm_is_root(nproc_id_global)) then
      fh_namelist = get_filehandle()
      open(fh_namelist, file='.namelist.tmp', status='old')

      read(fh_namelist, nml=calculation, iostat=inml_calculation)
      rewind(fh_namelist)

      read(fh_namelist, nml=control, iostat=inml_control)
      rewind(fh_namelist)

      read(fh_namelist, nml=parallel, iostat=inml_parallel)
      rewind(fh_namelist)

      read(fh_namelist, nml=system, iostat=inml_system)
      rewind(fh_namelist)

      read(fh_namelist, nml=pseudo, iostat=inml_pseudo)
      rewind(fh_namelist)

      read(fh_namelist, nml=functional, iostat=inml_functional)
      rewind(fh_namelist)

      read(fh_namelist, nml=rgrid, iostat=inml_rgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=kgrid, iostat=inml_kgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=tgrid, iostat=inml_tgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=propagation, iostat=inml_propagation)
      rewind(fh_namelist)

      read(fh_namelist, nml=scf, iostat=inml_scf)
      rewind(fh_namelist)

      read(fh_namelist, nml=emfield, iostat=inml_emfield)
      rewind(fh_namelist)

      read(fh_namelist, nml=multiscale, iostat=inml_multiscale)
      rewind(fh_namelist)

      read(fh_namelist, nml=maxwell, iostat=inml_maxwell)
      rewind(fh_namelist)

      read(fh_namelist, nml=analysis, iostat=inml_analysis)
      rewind(fh_namelist)

      read(fh_namelist, nml=poisson, iostat=inml_poisson)
      rewind(fh_namelist)

      read(fh_namelist, nml=ewald, iostat=inml_ewald)
      rewind(fh_namelist)

      read(fh_namelist, nml=opt, iostat=inml_opt)
      rewind(fh_namelist)

      read(fh_namelist, nml=md, iostat=inml_md)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_fundamental, iostat=inml_group_fundamental)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_hartree, iostat=inml_group_hartree)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_others, iostat=inml_group_others)
      rewind(fh_namelist)

      read(fh_namelist, nml=code, iostat=inml_code)
      rewind(fh_namelist)

      close(fh_namelist)
    end if


!! convert lowercase
    call string_lowercase(theory)


! Broad cast
!! == bcast for &calculation
    call comm_bcast(theory             ,nproc_group_global)
    call comm_bcast(yn_md              ,nproc_group_global)
    call comm_bcast(yn_opt             ,nproc_group_global)

!! == bcast for &control
    call comm_bcast(sysname         ,nproc_group_global)
    call comm_bcast(base_directory  ,nproc_group_global)
    ii = len_trim(base_directory)
    if(base_directory(ii:ii).ne.'/') &
       base_directory = trim(base_directory)//'/'
    call comm_bcast(output_buffer_interval,nproc_group_global)
    call comm_bcast(yn_restart            ,nproc_group_global)
    call comm_bcast(directory_read_data   ,nproc_group_global)
    ii = len_trim(directory_read_data)
    if(directory_read_data(ii:ii).ne.'/') &
       directory_read_data = trim(directory_read_data)//'/'
    call comm_bcast(yn_self_checkpoint    ,nproc_group_global)
    call comm_bcast(checkpoint_interval   ,nproc_group_global)
    call comm_bcast(yn_reset_step_restart ,nproc_group_global)
    call comm_bcast(read_gs_restart_data  ,nproc_group_global)
    call comm_bcast(write_gs_restart_data ,nproc_group_global)
    call comm_bcast(time_shutdown         ,nproc_group_global)
    call comm_bcast(method_wf_distributor ,nproc_group_global)
    call comm_bcast(nblock_wf_distribute  ,nproc_group_global)
    call comm_bcast(yn_gbp                ,nproc_group_global)
    call comm_bcast(yn_gbp_fourier0       ,nproc_group_global) ! temporary
    !remove later
    call comm_bcast(read_gs_dns_cube,nproc_group_global)
    call comm_bcast(write_gs_wfn_k  ,nproc_group_global)
    call comm_bcast(write_rt_wfn_k  ,nproc_group_global)

!! == bcast for &parallel
    call comm_bcast(nproc_k             ,nproc_group_global)
    call comm_bcast(nproc_ob            ,nproc_group_global)
    call comm_bcast(nproc_rgrid         ,nproc_group_global)
    call comm_bcast(yn_ffte             ,nproc_group_global)
    call comm_bcast(yn_scalapack        ,nproc_group_global)
    call comm_bcast(yn_scalapack_red_mem ,nproc_group_global)
    call comm_bcast(yn_eigenexa         ,nproc_group_global)
    call comm_bcast(process_allocation  ,nproc_group_global)
!! == bcast for &system
    call comm_bcast(yn_periodic,nproc_group_global)
    if(yn_periodic=='y') iperiodic=3
    if(yn_periodic=='n') iperiodic=0
    call comm_bcast(ispin      ,nproc_group_global)
    call comm_bcast(al         ,nproc_group_global)
    al = al * ulength_to_au
    call comm_bcast(al_vec1            ,nproc_group_global)
    call comm_bcast(al_vec2            ,nproc_group_global)
    call comm_bcast(al_vec3            ,nproc_group_global)
    al_vec1 = al_vec1 * ulength_to_au
    al_vec2 = al_vec2 * ulength_to_au
    al_vec3 = al_vec3 * ulength_to_au
    call comm_bcast(isym               ,nproc_group_global)
    call comm_bcast(crystal_structure  ,nproc_group_global)
    call comm_bcast(nstate             ,nproc_group_global)
    call comm_bcast(nstate_spin        ,nproc_group_global)
    call comm_bcast(nelec              ,nproc_group_global)
    call comm_bcast(nelec_spin         ,nproc_group_global)
    call comm_bcast(temperature        ,nproc_group_global)
    call comm_bcast(temperature_k      ,nproc_group_global)
    if(temperature_k>=0d0) temperature = temperature_k * kB_au ! Kelvin --> atomic units
    call comm_bcast(nelem              ,nproc_group_global)
    call comm_bcast(natom              ,nproc_group_global)
    call comm_bcast(file_atom_coor     ,nproc_group_global)
    call comm_bcast(file_atom_red_coor ,nproc_group_global)
!! == bcast for &pseudo
    call comm_bcast(file_pseudo  ,nproc_group_global)
    call comm_bcast(lmax_ps      ,nproc_group_global)
    call comm_bcast(lloc_ps      ,nproc_group_global)
    call comm_bcast(izatom       ,nproc_group_global)
    call comm_bcast(yn_psmask,nproc_group_global)
    call comm_bcast(alpha_mask   ,nproc_group_global)
    call comm_bcast(gamma_mask   ,nproc_group_global)
    call comm_bcast(eta_mask     ,nproc_group_global)
!! == bcast for &functional

#ifdef USE_LIBXC
    if (alibxc .ne. 'none') xc =  trim(alibxc)
    if (alibx .ne. 'none') xname =  trim(alibx)
    if (alibc .ne. 'none') cname =  trim(alibc)
#endif
    call comm_bcast(xc           ,nproc_group_global)
    call comm_bcast(cname        ,nproc_group_global)
    call comm_bcast(xname        ,nproc_group_global)
#ifdef USE_LIBXC
    call comm_bcast(alibxc       ,nproc_group_global)
    call comm_bcast(alibx        ,nproc_group_global)
    call comm_bcast(alibc        ,nproc_group_global)
#endif
    call comm_bcast(cval         ,nproc_group_global)
!! == bcast for &rgrid
    call comm_bcast(dl,nproc_group_global)
    dl = dl * ulength_to_au
    call comm_bcast(num_rgrid,nproc_group_global)
!! == bcast for &kgrid
    call comm_bcast(num_kgrid,nproc_group_global)
    call comm_bcast(file_kw  ,nproc_group_global)
!! == bcast for &tgrid
    call comm_bcast(nt,nproc_group_global)
    call comm_bcast(dt,nproc_group_global)
    dt = dt * utime_to_au
    call comm_bcast(gram_schmidt_interval,nproc_group_global)
!! == bcast for &propagation
    call comm_bcast(n_hamil    ,nproc_group_global)
    call comm_bcast(propagator ,nproc_group_global)
    call comm_bcast(yn_fix_func,nproc_group_global)
!! == bcast for &scf
    call comm_bcast(method_init_wf          ,nproc_group_global)
    call comm_bcast(iseed_number_change     ,nproc_group_global)
    call comm_bcast(method_min              ,nproc_group_global)
    call comm_bcast(ncg                     ,nproc_group_global)
    call comm_bcast(ncg_init                ,nproc_group_global)
    call comm_bcast(method_mixing           ,nproc_group_global)
    call comm_bcast(mixrate                 ,nproc_group_global)
    call comm_bcast(nmemory_mb              ,nproc_group_global)
    call comm_bcast(alpha_mb                ,nproc_group_global)
    call comm_bcast(nmemory_p               ,nproc_group_global)
    call comm_bcast(beta_p                  ,nproc_group_global)
    call comm_bcast(yn_auto_mixing          ,nproc_group_global)
    call comm_bcast(update_mixing_ratio     ,nproc_group_global)
    call comm_bcast(nscf                    ,nproc_group_global)
    call comm_bcast(yn_subspace_diagonalization,nproc_group_global)
    call comm_bcast(convergence             ,nproc_group_global)
    call comm_bcast(threshold               ,nproc_group_global)
    if(threshold.lt.-0.5d0) then !<-- when not input value
      select case(convergence)
      case('rho_dne')
         threshold = 1d-17  !default (a.u.)
      case('norm_rho','norm_rho_dng')
         threshold = -1d0   !default (a.u.)
      case('norm_pot','norm_pot_dng')
         threshold = -1d0   !default (a.u.)
      end select
    else
      select case(convergence)
      case('rho_dne')
         continue
      case('norm_rho','norm_rho_dng')
         threshold = threshold / (ulength_to_au)**6
      case('norm_pot','norm_pot_dng')
         threshold = threshold * (uenergy_to_au)**2 / (ulength_to_au)**6
      end select
    endif
    call comm_bcast(iditer_notemperature    ,nproc_group_global)
    call comm_bcast(step_initial_mix_zero   ,nproc_group_global)
    call comm_bcast(conv_gap_mix_zero       ,nproc_group_global)
    conv_gap_mix_zero = conv_gap_mix_zero * uenergy_to_au

!! == bcast for &emfield
    call comm_bcast(trans_longi,nproc_group_global)
    call comm_bcast(ae_shape1  ,nproc_group_global)
    call comm_bcast(file_input1,nproc_group_global)
    call comm_bcast(e_impulse,nproc_group_global)
    e_impulse = e_impulse *uenergy_to_au/ulength_to_au*utime_to_au
    call comm_bcast(E_amplitude1 ,nproc_group_global)
    E_amplitude1 = E_amplitude1*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(I_wcm2_1,nproc_group_global)
    call comm_bcast(tw1  ,nproc_group_global)
    tw1 = tw1 * utime_to_au
    call comm_bcast(omega1,nproc_group_global)
    omega1 = omega1 * uenergy_to_au
    call comm_bcast(epdir_re1 ,nproc_group_global)
    call comm_bcast(epdir_im1 ,nproc_group_global)
    norm = sqrt(sum(epdir_re1(:)**2+epdir_im1(:)**2))
    if ( norm > 0.0d0 ) then
    epdir_re1 = epdir_re1 / norm
    epdir_im1 = epdir_im1 / norm
    end if
    call comm_bcast(phi_cep1  ,nproc_group_global)
    call comm_bcast(ae_shape2 ,nproc_group_global)
    call comm_bcast(E_amplitude2,nproc_group_global)
    E_amplitude2 = E_amplitude2*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(I_wcm2_2,nproc_group_global)
    call comm_bcast(tw2  ,nproc_group_global)
    tw2 = tw2 * utime_to_au
    call comm_bcast(omega2,nproc_group_global)
    omega2 = omega2 * uenergy_to_au
    call comm_bcast(epdir_re2,nproc_group_global)
    call comm_bcast(epdir_im2,nproc_group_global)
    norm = sqrt(sum(epdir_re2(:)**2+epdir_im2(:)**2))
    epdir_re2 = epdir_re2 / norm
    epdir_im2 = epdir_im2 / norm
    call comm_bcast(phi_cep2 ,nproc_group_global)
    call comm_bcast(t1_t2    ,nproc_group_global)
    t1_t2 = t1_t2 * utime_to_au
    call comm_bcast(t1_start ,nproc_group_global)
    t1_start = t1_start * utime_to_au
    call comm_bcast(num_dipole_source,nproc_group_global)
    call comm_bcast(vec_dipole_source,nproc_group_global)
    vec_dipole_source = vec_dipole_source * ulength_to_au
    call comm_bcast(cood_dipole_source,nproc_group_global)
    cood_dipole_source = cood_dipole_source * ulength_to_au
    call comm_bcast(rad_dipole_source,nproc_group_global)
    rad_dipole_source = rad_dipole_source * ulength_to_au
    call comm_bcast(cutoff_G2_emfield,nproc_group_global)
    cutoff_G2_emfield = cutoff_G2_emfield * uenergy_to_au

!! == bcast for &multiscale
    call comm_bcast(fdtddim   ,nproc_group_global)
    call comm_bcast(twod_shape,nproc_group_global)
    call comm_bcast(nx_m      ,nproc_group_global)
    call comm_bcast(ny_m      ,nproc_group_global)
    call comm_bcast(nz_m      ,nproc_group_global)
    call comm_bcast(hx_m      ,nproc_group_global)
    hx_m = hx_m * ulength_to_au
    call comm_bcast(hy_m      ,nproc_group_global)
    hy_m = hy_m * ulength_to_au
    call comm_bcast(hz_m      ,nproc_group_global)
    hz_m = hz_m * ulength_to_au
    call comm_bcast(nksplit   ,nproc_group_global)
    call comm_bcast(nxysplit  ,nproc_group_global)
    call comm_bcast(nxvacl_m  ,nproc_group_global)
    call comm_bcast(nxvacr_m  ,nproc_group_global)
    call comm_bcast(nx_origin_m,nproc_group_global)
    call comm_bcast(ny_origin_m,nproc_group_global)
    call comm_bcast(nz_origin_m,nproc_group_global)
    call comm_bcast(file_macropoint, nproc_group_global)
    call comm_bcast(num_macropoint,  nproc_group_global)
    call comm_bcast(set_ini_coor_vel,nproc_group_global)
    call comm_bcast(nmacro_write_group,nproc_group_global)

!! == bcast for &maxwell
    call comm_bcast(al_em           ,nproc_group_global)
    al_em = al_em * ulength_to_au
    call comm_bcast(dl_em           ,nproc_group_global)
    dl_em = dl_em * ulength_to_au
    call comm_bcast(dt_em           ,nproc_group_global)
    dt_em = dt_em * utime_to_au
    call comm_bcast(nt_em           ,nproc_group_global)
    call comm_bcast(boundary_em     ,nproc_group_global)
    call comm_bcast(shape_file      ,nproc_group_global)
    call comm_bcast(media_num       ,nproc_group_global)
    call comm_bcast(media_type      ,nproc_group_global)
    call comm_bcast(epsilon_em      ,nproc_group_global)
    call comm_bcast(mu_em           ,nproc_group_global)
    call comm_bcast(sigma_em        ,nproc_group_global)
    call comm_bcast(pole_num_ld     ,nproc_group_global)
    call comm_bcast(omega_p_ld      ,nproc_group_global)
    omega_p_ld = omega_p_ld * uenergy_to_au
    call comm_bcast(f_ld            ,nproc_group_global)
    call comm_bcast(gamma_ld        ,nproc_group_global)
    gamma_ld = gamma_ld * uenergy_to_au
    call comm_bcast(omega_ld        ,nproc_group_global)
    omega_ld = omega_ld * uenergy_to_au
    call comm_bcast(yn_wf_em        ,nproc_group_global)
    call comm_bcast(wave_input,nproc_group_global)
    call comm_bcast(ek_dir1         ,nproc_group_global)
    call comm_bcast(source_loc1     ,nproc_group_global)
    source_loc1 = source_loc1 * ulength_to_au
    call comm_bcast(ek_dir2         ,nproc_group_global)
    call comm_bcast(source_loc2     ,nproc_group_global)
    source_loc2 = source_loc2 * ulength_to_au
    call comm_bcast(obs_num_em      ,nproc_group_global)
    call comm_bcast(obs_samp_em     ,nproc_group_global)
    call comm_bcast(obs_loc_em      ,nproc_group_global)
    obs_loc_em = obs_loc_em * ulength_to_au
    call comm_bcast(yn_obs_plane_em ,nproc_group_global)

!! == bcast for &analysis
    call comm_bcast(projection_option   ,nproc_group_global)
    call comm_bcast(out_projection_step ,nproc_group_global)
    call comm_bcast(nenergy             ,nproc_group_global)
    call comm_bcast(de                  ,nproc_group_global)
    de = de * uenergy_to_au
    call comm_bcast(yn_out_psi          ,nproc_group_global)
    call comm_bcast(yn_out_dos          ,nproc_group_global)
    call comm_bcast(yn_out_dos_set_fe_origin ,nproc_group_global)
    call comm_bcast(out_dos_start       ,nproc_group_global)
    out_dos_start = out_dos_start * uenergy_to_au
    call comm_bcast(out_dos_end         ,nproc_group_global)
    out_dos_end = out_dos_end * uenergy_to_au
    call comm_bcast(out_dos_nenergy     ,nproc_group_global)
    call comm_bcast(out_dos_width       ,nproc_group_global)
    out_dos_width = out_dos_width * uenergy_to_au
    call comm_bcast(out_dos_function    ,nproc_group_global)
    call comm_bcast(yn_out_pdos         ,nproc_group_global)
    call comm_bcast(yn_out_dns          ,nproc_group_global)
    call comm_bcast(yn_out_dns_rt       ,nproc_group_global)
    call comm_bcast(yn_out_dns_ac_je    ,nproc_group_global)
    call comm_bcast(out_dns_rt_step     ,nproc_group_global)
    call comm_bcast(out_dns_ac_je_step  ,nproc_group_global)
    call comm_bcast(out_old_dns         ,nproc_group_global)
    call comm_bcast(yn_out_dns_trans    ,nproc_group_global)
    call comm_bcast(out_dns_trans_energy,nproc_group_global)
    out_dns_trans_energy = out_dns_trans_energy * uenergy_to_au
    call comm_bcast(yn_out_elf          ,nproc_group_global)
    call comm_bcast(yn_out_elf_rt       ,nproc_group_global)
    call comm_bcast(out_elf_rt_step     ,nproc_group_global)
    call comm_bcast(yn_out_estatic_rt   ,nproc_group_global)
    call comm_bcast(out_estatic_rt_step ,nproc_group_global)
    call comm_bcast(yn_out_rvf_rt       ,nproc_group_global)
    call comm_bcast(out_rvf_rt_step     ,nproc_group_global)
    call comm_bcast(yn_out_tm           ,nproc_group_global)
    call comm_bcast(out_ms_step         ,nproc_group_global)
    call comm_bcast(format_voxel_data   ,nproc_group_global)
    call comm_bcast(nsplit_voxel_data   ,nproc_group_global)
    call comm_bcast(timer_process       ,nproc_group_global)

!! == bcast for &poisson
    call comm_bcast(layout_multipole  ,nproc_group_global)
    call comm_bcast(num_multipole_xyz ,nproc_group_global)
    call comm_bcast(threshold_cg      ,nproc_group_global)
    threshold_cg = threshold_cg * (uenergy_to_au)**2 * (ulength_to_au)**3
!! == bcast for &ewald
    call comm_bcast(newald        ,nproc_group_global)
    call comm_bcast(aewald        ,nproc_group_global)
    call comm_bcast(cutoff_r      ,nproc_group_global)
    call comm_bcast(cutoff_r_buff ,nproc_group_global)
    call comm_bcast(cutoff_g      ,nproc_group_global)
    cutoff_r      = cutoff_r * ulength_to_au
    cutoff_r_buff = cutoff_r_buff * ulength_to_au
    cutoff_g      = cutoff_g / ulength_to_au
!! == bcast for &opt
    call comm_bcast(nopt                ,nproc_group_global)
    call comm_bcast(max_step_len_adjust ,nproc_group_global)
    call comm_bcast(convrg_opt_fmax     ,nproc_group_global)
!! == bcast for &md
    call comm_bcast(ensemble               ,nproc_group_global)
    call comm_bcast(thermostat             ,nproc_group_global)
    call comm_bcast(step_velocity_scaling  ,nproc_group_global)
    call comm_bcast(step_update_ps         ,nproc_group_global)
    call comm_bcast(temperature0_ion_k     ,nproc_group_global)
    call comm_bcast(yn_set_ini_velocity    ,nproc_group_global)
    call comm_bcast(file_ini_velocity      ,nproc_group_global)
    call comm_bcast(thermostat_tau         ,nproc_group_global)
    thermostat_tau = thermostat_tau * utime_to_au
    call comm_bcast(yn_stop_system_momt    ,nproc_group_global)
!! == bcast for &group_fundamental
    call comm_bcast(iditer_nosubspace_diag,nproc_group_global)
    call comm_bcast(ntmg                  ,nproc_group_global)
    call comm_bcast(idisnum               ,nproc_group_global)
    call comm_bcast(iwrite_projection     ,nproc_group_global)
    call comm_bcast(itwproj               ,nproc_group_global)
    call comm_bcast(iwrite_projnum        ,nproc_group_global)
    call comm_bcast(itcalc_ene            ,nproc_group_global)
!! == bcast for &group_hartree
    call comm_bcast(lmax_lmp,nproc_group_global)
!! == bcast for &group_others
    call comm_bcast(iswitch_orbital_mesh,nproc_group_global)
    call comm_bcast(iflag_psicube       ,nproc_group_global)
    call comm_bcast(num_projection      ,nproc_group_global)
    call comm_bcast(iwrite_projection_ob,nproc_group_global)
    call comm_bcast(iwrite_projection_k ,nproc_group_global)
    call comm_bcast(filename_pot        ,nproc_group_global)
    call comm_bcast(iwrite_external     ,nproc_group_global)
    call comm_bcast(iflag_intelectron   ,nproc_group_global)
    call comm_bcast(num_dip2            ,nproc_group_global)
    call comm_bcast(dip2boundary        ,nproc_group_global)
    dip2boundary  = dip2boundary  * ulength_to_au
    call comm_bcast(dip2center          ,nproc_group_global)
    dip2center    = dip2center    * ulength_to_au
    call comm_bcast(itotntime2          ,nproc_group_global)
    call comm_bcast(iwdenoption         ,nproc_group_global)
    call comm_bcast(iwdenstep           ,nproc_group_global)
    call comm_bcast(iflag_estatic       ,nproc_group_global)
!! == bcast for code
    call comm_bcast(yn_want_stencil_hand_vectorization     ,nproc_group_global)
    call comm_bcast(yn_want_communication_overlapping      ,nproc_group_global)
    call comm_bcast(stencil_openmp_mode                    ,nproc_group_global)
    call comm_bcast(current_openmp_mode                    ,nproc_group_global)
    call comm_bcast(force_openmp_mode                      ,nproc_group_global)
  end subroutine read_input_common

  subroutine read_atomic_coordinates
    use parallelization
    use communication
    use filesystem, only: get_filehandle
    use salmon_global, only: directory_read_data,yn_restart
    use checkpoint_restart_sub, only: generate_restart_directory_name
    character(256) :: filename_tmp,char_atom, gdir,wdir
    integer :: icount,i
    logical :: if_error, if_cartesian

    if (comm_is_root(nproc_id_global)) then

      if_error     = .false.
      if_cartesian = .true.
      iflag_atom_coor = ntype_atom_coor_none
      icount = 0
      if(file_atom_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = trim(file_atom_coor)
        iflag_atom_coor = ntype_atom_coor_cartesian
      end if

      if(file_atom_red_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = trim(file_atom_coor)
        iflag_atom_coor = ntype_atom_coor_reduced
      end if

      if(if_nml_coor)then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = '.atomic_coor.tmp'
        iflag_atom_coor = ntype_atom_coor_cartesian
      end if

      if(if_nml_red_coor)then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = '.atomic_red_coor.tmp'
        iflag_atom_coor = ntype_atom_coor_reduced
      end if

      if(icount==0 .and. (yn_restart == 'y' .or. &
         index(theory,'tddft_')/=0 .or. index(theory,'_tddft')/=0 ) ) then

        if (comm_is_root(nproc_id_global))then
           write(*,"(A)") '  Atomic coordinate is read from restart directory'
        end if

        call generate_restart_directory_name(directory_read_data,gdir,wdir)

        icount = icount + 1
        if_cartesian = .true.
        if(yn_self_checkpoint == 'y') then
           filename_tmp = trim(gdir)//"rank_000000/atomic_coor.txt"
        else
           filename_tmp = trim(gdir)//"atomic_coor.txt"
        endif
        iflag_atom_coor = ntype_atom_coor_cartesian
      end if

    end if

    call comm_bcast(icount,nproc_group_global)
    call comm_bcast(if_cartesian,nproc_group_global)
    call comm_bcast(iflag_atom_coor,nproc_group_global)

    if(0 < natom .and. icount/=1)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(A)")'Error in input: The following inputs are incompatible.'
         write(*,"(A)")'file_atom_coor, file_atom_red_coor, &atomic_coor, and &atomic_red_coor.'
       end if
       call end_parallel
       stop
    end if

    if( (.not.if_cartesian) .and. iperiodic == 0)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(A)")'Error in input: Reduced coordinate is invalid for isolated systems.'
       end if
       call end_parallel
       stop
    end if

    allocate(atom_name(natom))
    allocate(Rion(3,natom), Rion_red(3,natom),kion(natom), flag_opt_atom(natom))
    Rion     = 0d0
    Rion_red = 0d0
    kion = 0
    flag_opt_atom = 'n'

    if (0 < natom) then

      if (comm_is_root(nproc_id_global))then
        fh_atomic_coor = get_filehandle()
        open(fh_atomic_coor, file=filename_tmp, status='old')
        select case(iflag_atom_coor)
        case(ntype_atom_coor_cartesian)
           do i=1,natom
              if(yn_opt == 'y')then
                 read(fh_atomic_coor, *) char_atom, Rion(:,i), kion(i), flag_opt_atom(i)
              else
                 read(fh_atomic_coor, *) char_atom, Rion(:,i), kion(i)
              end if
              atom_name(i) = char_atom
           end do
           Rion = Rion*ulength_to_au

        case(ntype_atom_coor_reduced)
           do i=1,natom
              if(yn_opt == 'y')then
                 read(fh_atomic_coor, *) char_atom, Rion_red(:,i), kion(i), flag_opt_atom(i)
              else
                 read(fh_atomic_coor, *) char_atom, Rion_red(:,i), kion(i)
              end if
              atom_name(i) = char_atom
           end do

        end select
        close(fh_atomic_coor)

      end if

      call comm_bcast(Rion,nproc_group_global)
      call comm_bcast(Rion_red,nproc_group_global)
      call comm_bcast(kion,nproc_group_global)
      call comm_bcast(flag_opt_atom,nproc_group_global)
      call comm_bcast(atom_name,nproc_group_global)
    end if ! if 0 < natom


  end subroutine read_atomic_coordinates

  subroutine initialize_inputoutput_units
    implicit none


! Unit for time
    select case(unit_time)
    case('au','a.u.')
      utime_to_au   = 1d0
      utime_from_au = 1d0
      iflag_unit_time = ntype_unit_time_au
    case('fs','femtosecond')
      utime_to_au   = 1d0/au_time_fs
      utime_from_au = au_time_fs
      iflag_unit_time = ntype_unit_time_fs
    case default
      stop "Invalid unit for time."
    end select

! Unit for length
    select case(unit_length)
    case('au','a.u.')
      ulength_to_au   = 1d0
      ulength_from_au = 1d0
      iflag_unit_length = ntype_unit_length_au
    case('AA','angstrom','Angstrom')
      ulength_to_au   = 1d0/au_length_aa
      ulength_from_au = au_length_aa
      iflag_unit_length = ntype_unit_length_aa
    case default
      stop "Invalid unit for length."
    end select

! Unit for energy
    select case(unit_energy)
    case('au','a.u.')
      uenergy_to_au   = 1d0
      uenergy_from_au = 1d0
      iflag_unit_energy = ntype_unit_energy_au
    case('ev','eV')
      uenergy_to_au   = 1d0/au_energy_ev
      uenergy_from_au = au_energy_ev
      iflag_unit_energy = ntype_unit_energy_ev
    case default
      stop "Invalid unit for energy."
    end select

! Unit for charge
    select case(unit_charge)
    case('au','a.u.')
      ucharge_to_au   = 1d0
      ucharge_from_au = 1d0
      iflag_unit_charge = ntype_unit_charge_au
    case default
      stop "Invalid unit for charge."
    end select

!! prepare type(unit_t) :: t_unit_length,t_unit_length_inv
    t_unit_length%conv = ulength_from_au
    t_unit_length_inv%conv = 1d0/ulength_from_au
    if(iflag_unit_length == ntype_unit_length_aa)then
      t_unit_length%name     = 'Angstrom'
      t_unit_length_inv%name = '1/Angstrom'
    else
      t_unit_length%name     = 'a.u.'
      t_unit_length_inv%name = 'a.u.'
      t_unit_length%conv = 1d0
      t_unit_length_inv%conv = 1d0
    end if

!! prepare type(unit_t) :: t_unit_energy,t_unit_energy_inv
    t_unit_energy%conv = uenergy_from_au
    t_unit_energy_inv%conv = 1d0/uenergy_from_au
    if(iflag_unit_energy == ntype_unit_energy_ev)then
      t_unit_energy%name     = 'eV'
      t_unit_energy_inv%name = '1/eV'
    else
      t_unit_energy%name     = 'a.u.'
      t_unit_energy_inv%name = 'a.u.'
      t_unit_energy%conv = 1d0
      t_unit_energy_inv%conv = 1d0
    end if

!! prepare type(unit_t) :: t_unit_time,t_unit_time_inv
    t_unit_time%conv = utime_from_au
    t_unit_time_inv%conv = 1d0/utime_from_au
    if(iflag_unit_time == ntype_unit_time_fs)then
      t_unit_time%name     = 'fs'
      t_unit_time_inv%name = '1/fs'
    else
      t_unit_time%name     = 'a.u.'
      t_unit_time_inv%name = 'a.u.'
      t_unit_time%conv = 1d0
      t_unit_time_inv%conv = 1d0
    end if

!! prepare type(unit_t) :: t_unit_spectrum_dipole
    t_unit_spectrum_dipole%conv = utime_from_au*ulength_from_au
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_spectrum_dipole%name  = 'fs*Angstrom'
    else
      t_unit_spectrum_dipole%name  = 'a.u.'
      t_unit_spectrum_dipole%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_spectrum_dipole_square
    t_unit_spectrum_dipole_square%conv = utime_from_au**2*ulength_from_au**2
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_spectrum_dipole_square%name  = 'fs^2*Angstrom^2'
    else
      t_unit_spectrum_dipole_square%name  = 'a.u.'
      t_unit_spectrum_dipole_square%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_current
    t_unit_current%conv = (ulength_from_au/utime_from_au)/ulength_from_au**3
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_current%name  = '1/fs*Angstrom^2'
    else
      t_unit_current%name  = 'a.u.'
      t_unit_current%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_spectrum_current
    t_unit_spectrum_current%conv = 1d0/ulength_from_au**2
    if(iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_spectrum_current%name  = '1/Angstrom^2'
    else
      t_unit_spectrum_current%name  = 'a.u.'
      t_unit_spectrum_current%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_spectrum_current_square
    t_unit_spectrum_current_square%conv = 1d0/ulength_from_au**4
    if(iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_spectrum_current_square%name  = '1/Angstrom^4'
    else
      t_unit_spectrum_current_square%name  = 'a.u.'
      t_unit_spectrum_current_square%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_ac
    t_unit_ac%conv = utime_from_au*uenergy_from_au/ulength_from_au/ucharge_from_au
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_ac%name     = 'fs*V/Angstrom'
    else
      t_unit_ac%name     = 'a.u.'
      t_unit_ac%conv     = 1d0
    end if

    !! prepare type(unit_t) :: t_unit_elec
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_elec%name     = 'V/Angstrom'
      t_unit_elec%conv     = 51.42206707d0
    else
      t_unit_elec%name     = 'a.u.'
      t_unit_elec%conv     = 1d0
    end if

    !! prepare type(unit_t) :: t_unit_spectrum_elec
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_spectrum_elec%name     = 'fs*V/Angstrom'
      t_unit_spectrum_elec%conv     = utime_from_au*51.42206707d0
    else
      t_unit_spectrum_elec%name     = 'a.u.'
      t_unit_spectrum_elec%conv     = 1d0
    end if

    !! prepare type(unit_t) :: t_unit_spectrum_elec_square
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_spectrum_elec_square%name     = 'fs^2*V^2/Angstrom^2'
      t_unit_spectrum_elec_square%conv     = utime_from_au**2*51.42206707d0**2
    else
      t_unit_spectrum_elec_square%name     = 'a.u.'
      t_unit_spectrum_elec_square%conv     = 1d0
    end if

!! prepare type(unit_t) :: t_unit_polarizability
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_polarizability%name  = 'Augstrom^2/V'
      t_unit_polarizability%conv = ulength_from_au**2/51.42206707d0
    else
      t_unit_polarizability%name  = 'a.u.'
      t_unit_polarizability%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_conductivity
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_conductivity%name  = '1/fs*V*Angstrom'
      t_unit_conductivity%conv = 1.d0/utime_from_au/51.42206707d0/ulength_from_au
    else
      t_unit_conductivity%name  = 'a.u.'
      t_unit_conductivity%conv  = 1d0
    end if

  end subroutine initialize_inputoutput_units

  subroutine dump_input_common
    use parallelization
    use communication
    use filesystem, only: get_filehandle
    use filesystem, only: atomic_create_directory
    implicit none
    integer :: i,j,ierr_nml
    ierr_nml = 0

    if (comm_is_root(nproc_id_global)) then

      fh_variables_log = get_filehandle()
      open(fh_variables_log,file='variables.log')

      if(inml_calculation >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'calculation', inml_calculation
      write(fh_variables_log, '("#",4X,A,"=",A)') 'theory', theory
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_md', yn_md
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_opt', yn_opt

      if(inml_control >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'control', inml_control
      write(fh_variables_log, '("#",4X,A,"=",A)') 'sysname', trim(sysname)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'base_directory', trim(base_directory)
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'output_buffer_interval', output_buffer_interval
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_restart', yn_restart
      write(fh_variables_log, '("#",4X,A,"=",A)') 'directory_read_data', trim(directory_read_data)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_self_checkpoint', yn_self_checkpoint
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'checkpoint_interval', checkpoint_interval
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_reset_step_restart', yn_reset_step_restart
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_gs_restart_data', trim(read_gs_restart_data)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_gs_restart_data', trim(write_gs_restart_data)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'time_shutdown', time_shutdown
      write(fh_variables_log, '("#",4X,A,"=",A)') 'method_wf_distributor', method_wf_distributor
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nblock_wf_distribute', nblock_wf_distribute
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_gbp', yn_gbp
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_gbp_fourier0', yn_gbp_fourier0 ! temporary
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_gs_dns_cube', trim(read_gs_dns_cube)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_gs_wfn_k', trim(write_gs_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_rt_wfn_k', trim(write_rt_wfn_k)


      if(inml_units >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'units', inml_units
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_system', unit_system

      if(inml_parallel >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'parallel', inml_parallel
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_k', nproc_k
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_ob', nproc_ob
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_rgrid(1)', nproc_rgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_rgrid(2)', nproc_rgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_rgrid(3)', nproc_rgrid(3)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_ffte', yn_ffte
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_scalapack', yn_scalapack
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_scalapack_red_mem', yn_scalapack_red_mem
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_eigenexa', yn_eigenexa
      write(fh_variables_log, '("#",4X,A,"=",A)') 'process_allocation', process_allocation

      if(inml_system >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'system', inml_system
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_periodic', yn_periodic
      write(fh_variables_log, '("#",4X,A,"=",I1)') 'ispin', ispin
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(1)', al(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(2)', al(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(3)', al(3)
      write(fh_variables_log, '("#",4X,A,"=",3ES12.5)') 'al_vec1(1:3)', al_vec1(1:3)
      write(fh_variables_log, '("#",4X,A,"=",3ES12.5)') 'al_vec2(1:3)', al_vec2(1:3)
      write(fh_variables_log, '("#",4X,A,"=",3ES12.5)') 'al_vec3(1:3)', al_vec3(1:3)
      write(fh_variables_log, '("#",4X,A,"=",I1)') 'isym', isym
      write(fh_variables_log, '("#",4X,A,"=",A)') 'crystal_structure', crystal_structure
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nstate', nstate
      write(fh_variables_log, '("#",4X,A,"=",I4,2x,I4)') 'nstate_spin(1:2)', nstate_spin
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nelec', nelec
      write(fh_variables_log, '("#",4X,A,"=",I4,2x,I4)') 'nelec_spin(1:2)', nelec_spin
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature', temperature
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature_k', temperature_k
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nelem', nelem
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'natom', natom
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_coor', trim(file_atom_coor)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_red_coor', trim(file_atom_red_coor)

      if(inml_pseudo >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'pseudo', inml_pseudo

      do i = 1,nelem
        write(fh_variables_log, '("#",4X,A,I2,A,"=",A)') 'file_pseudo(',i,')', trim(file_pseudo(i))
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'lmax_ps(',i,')', lmax_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'lloc_ps(',i,')', lloc_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'izatom(',i,')', izatom(i)
      end do
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_psmask', yn_psmask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'alpha_mask', alpha_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'gamma_mask', gamma_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'eta_mask', eta_mask

      if(inml_functional >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'functional', inml_functional
      write(fh_variables_log, '("#",4X,A,"=",A)') 'xc', trim(xc)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'xname', trim(xname)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'cname', trim(cname)
#ifdef USE_LIBXC
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibxc', trim(alibxc)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibx', trim(alibx)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibc', trim(alibc)
#endif
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cval', cval

      if(inml_rgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'rgrid', inml_rgrid
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(1)', dl(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(2)', dl(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(3)', dl(3)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(1)', num_rgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(2)', num_rgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(3)', num_rgrid(3)

      if(inml_kgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'kgrid', inml_kgrid
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(1)', num_kgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(2)', num_kgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(3)', num_kgrid(3)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_kw', trim(file_kw)

      if(inml_tgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'tgrid', inml_tgrid
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'nt', nt
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dt', dt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'gram_schmidt_interval', gram_schmidt_interval

      if(inml_propagation >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'propagation', inml_propagation
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'n_hamil', n_hamil
      write(fh_variables_log, '("#",4X,A,"=",A)') 'propagator', trim(propagator)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_fix_func', yn_fix_func

      if(inml_scf >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'scf', inml_scf
      write(fh_variables_log, '("#",4X,A,"=",A)') 'method_init_wf', method_init_wf
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iseed_number_change', iseed_number_change
      write(fh_variables_log, '("#",4X,A,"=",A)') 'method_min', method_min
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'ncg', ncg
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'ncg_init', ncg_init
      write(fh_variables_log, '("#",4X,A,"=",A)') 'method_mixing', method_mixing
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'mixrate', mixrate
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nmemory_mb', nmemory_mb
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'alpha_mb', alpha_mb
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nmemory_p', nmemory_p
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'beta_p', beta_p
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_auto_mixing', yn_auto_mixing
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'update_mixing_ratio', update_mixing_ratio
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nscf', nscf
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_subspace_diagonalization', yn_subspace_diagonalization
      write(fh_variables_log, '("#",4X,A,"=",A)') 'convergence', convergence
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold', threshold
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'iditer_notemperature', iditer_notemperature
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'step_initial_mix_zero', step_initial_mix_zero
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'conv_gap_mix_zero', conv_gap_mix_zero

      if(inml_emfield >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'emfield', inml_emfield
      write(fh_variables_log, '("#",4X,A,"=",A)') 'trans_longi', trans_longi
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ae_shape1', ae_shape1
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_input1', file_input1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'e_impulse', e_impulse
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'E_amplitude1', E_amplitude1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'I_wcm2_1', I_wcm2_1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'tw1', tw1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'omega1', omega1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(1)', epdir_re1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(2)', epdir_re1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(3)', epdir_re1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(1)', epdir_im1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(2)', epdir_im1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(3)', epdir_im1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'phi_cep1', phi_cep1
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ae_shape2', ae_shape2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'E_amplitude2', E_amplitude2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'I_wcm2_2', I_wcm2_2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'tw2', tw2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'omega2', omega2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(1)', epdir_re2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(2)', epdir_re2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(3)', epdir_re2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(1)', epdir_im2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(2)', epdir_im2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(3)', epdir_im2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'phi_cep2', phi_cep2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 't1_t2', t1_t2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 't1_start', t1_start
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_dipole_source', num_dipole_source
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(1,1)', vec_dipole_source(1,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(2,1)', vec_dipole_source(2,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(3,1)', vec_dipole_source(3,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(1,2)', vec_dipole_source(1,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(2,2)', vec_dipole_source(2,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vec_dipole_source(3,2)', vec_dipole_source(3,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(1,1)', cood_dipole_source(1,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(2,1)', cood_dipole_source(2,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(3,1)', cood_dipole_source(3,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(1,2)', cood_dipole_source(1,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(2,2)', cood_dipole_source(2,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cood_dipole_source(3,2)', cood_dipole_source(3,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rad_dipole_source', rad_dipole_source
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cutoff_G2_emfield', cutoff_G2_emfield

      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'multiscale', inml_multiscale
      write(fh_variables_log, '("#",4X,A,"=",A)') 'fdtddim', fdtddim
      write(fh_variables_log, '("#",4X,A,"=",A)') 'twod_shape', twod_shape
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nx_m', nx_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ny_m', ny_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nz_m', nz_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hx_m', hx_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hy_m', hy_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hz_m', hz_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nksplit', nksplit
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nxysplit', nxysplit
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nxvacl_m', nxvacl_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nxvacr_m', nxvacr_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nx_origin_m', nx_origin_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'ny_origin_m', ny_origin_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nz_origin_m', nz_origin_m
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_macropoint', trim(file_macropoint)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'num_macropoint', num_macropoint
      write(fh_variables_log, '("#",4X,A,"=",A)') 'set_ini_coor_vel', set_ini_coor_vel
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nmacro_write_group', nmacro_write_group

      if(inml_maxwell >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'maxwell', inml_maxwell
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(1)', al_em(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(2)', al_em(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(3)', al_em(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(1)', dl_em(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(2)', dl_em(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(3)', dl_em(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dt_em', dt_em
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'nt_em', nt_em
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(1,1)', boundary_em(1,1)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(1,2)', boundary_em(1,2)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(2,1)', boundary_em(2,1)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(2,2)', boundary_em(2,2)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(3,1)', boundary_em(3,1)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'boundary_em(3,2)', boundary_em(3,2)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'shape_file', trim(shape_file)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'media_num', media_num
      do i = 0,media_num
        write(fh_variables_log, '("#",4X,A,I3,A,"=",A)')      'media_type(',i,')', media_type(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'epsilon_em(',i,')', epsilon_em(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'mu_em(',i,')', mu_em(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'sigma_em(',i,')', sigma_em(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",I6)')     'pole_num_ld(',i,')', pole_num_ld(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'omega_p_ld(',i,')', omega_p_ld(i)
        do j = 1,pole_num_ld(i)
          write(fh_variables_log, '("#",4X,A,I3,A,I3,A,"=",ES12.5)') 'f_ld(',i,',',j,')', f_ld(i,j)
          write(fh_variables_log, '("#",4X,A,I3,A,I3,A,"=",ES12.5)') 'gamma_ld(',i,',',j,')', gamma_ld(i,j)
          write(fh_variables_log, '("#",4X,A,I3,A,I3,A,"=",ES12.5)') 'omega_ld(',i,',',j,')', omega_ld(i,j)
        end do
      end do
      write(fh_variables_log, '("#",4X,A,"=",A)')      'wave_input', wave_input
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(1)', ek_dir1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(2)', ek_dir1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(3)', ek_dir1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(1)', source_loc1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(2)', source_loc1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(3)', source_loc1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(1)', ek_dir2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(2)', ek_dir2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(3)', ek_dir2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(1)', source_loc2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(2)', source_loc2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(3)', source_loc2(3)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'obs_num_em', obs_num_em
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'obs_samp_em', obs_samp_em
      if(obs_num_em==0) then
        write(fh_variables_log, '("#",4X,A,"=",3ES14.5)') 'obs_loc_em', obs_loc_em(1,1),obs_loc_em(1,2),obs_loc_em(1,3)
        write(fh_variables_log, '("#",4X,A,"=",A)')       'yn_obs_plane_em', yn_obs_plane_em(1)
      else
        do i = 1,obs_num_em
          write(fh_variables_log, '("#",4X,A,I3,A,"=",3ES14.5)') 'obs_loc_em(',i,',:)', obs_loc_em(i,:)
          write(fh_variables_log, '("#",4X,A,I3,A,"=",A)')       'yn_obs_plane_em(',i,')', yn_obs_plane_em(i)
        end do
      end if
      write(fh_variables_log, '("#",4X,A,"=",A)')      'yn_wf_em', yn_wf_em

      if(inml_analysis >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'analysis', inml_analysis
      write(fh_variables_log, '("#",4X,A,"=",A)') 'projection_option', projection_option
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_projection_step', out_projection_step
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'nenergy', nenergy
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'de', de
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_psi', yn_out_psi
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dos', yn_out_dos
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dos_set_fe_origin', yn_out_dos_set_fe_origin
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_start', out_dos_start
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_end', out_dos_end
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_dos_nenergy', out_dos_nenergy
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_width', out_dos_width
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dos_function', out_dos_function
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_pdos', yn_out_pdos
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dns', yn_out_dns
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dns_rt', yn_out_dns_rt
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dns_ac_je', yn_out_dns_ac_je
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_dns_rt_step', out_dns_rt_step
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_dns_ac_je_step', out_dns_ac_je_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_old_dns', out_old_dns
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_dns_trans', yn_out_dns_trans
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dns_trans_energy', out_dns_trans_energy
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_elf', yn_out_elf
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_elf_rt', yn_out_elf_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_elf_rt_step', out_elf_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_estatic_rt', yn_out_estatic_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_estatic_rt_step', out_estatic_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_rvf_rt', yn_out_rvf_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_rvf_rt_step', out_rvf_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_out_tm', yn_out_tm
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_ms_step', out_ms_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'format_voxel_data', format_voxel_data
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'nsplit_voxel_data', nsplit_voxel_data
      write(fh_variables_log, '("#",4X,A,"=",A)') 'timer_process', timer_process

      if(inml_poisson >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'poisson', inml_poisson
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'layout_multipole', layout_multipole
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_multipole_xyz(1)', num_multipole_xyz(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_multipole_xyz(2)', num_multipole_xyz(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_multipole_xyz(3)', num_multipole_xyz(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold_cg', threshold_cg

      if(inml_ewald >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'ewald', inml_ewald
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'newald', newald
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'aewald', aewald
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cutoff_r', cutoff_r
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cutoff_r_buff', cutoff_r_buff
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cutoff_g', cutoff_g

      if(inml_opt >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'opt', inml_opt
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nopt', nopt
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'max_step_len_adjust', max_step_len_adjust
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'convrg_opt_fmax',convrg_opt_fmax
      if(inml_md >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'md', inml_md
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ensemble', ensemble
      write(fh_variables_log, '("#",4X,A,"=",A)') 'thermostat', thermostat
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'step_velocity_scaling', step_velocity_scaling
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'step_update_ps', step_update_ps
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature0_ion_k', temperature0_ion_k
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_set_ini_velocity', yn_set_ini_velocity
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_ini_velocity', trim(file_ini_velocity)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'thermostat_tau', thermostat_tau
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_stop_system_momt', yn_stop_system_momt

      if(inml_group_fundamental >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_fundamental', inml_group_fundamental
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iditer_nosubspace_diag', iditer_nosubspace_diag
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ntmg', ntmg
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'idisnum(1)', idisnum(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'idisnum(2)', idisnum(2)
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwrite_projection', iwrite_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itwproj', itwproj
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projnum', iwrite_projnum
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itcalc_ene', itcalc_ene

      if(inml_group_hartree >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_hartree', inml_group_hartree
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'lmax_lmp', lmax_lmp

      if(inml_group_others >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_others', inml_group_others
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iswitch_orbital_mesh', iswitch_orbital_mesh
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_psicube', iflag_psicube
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_projection', num_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_projection', num_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_ob(1)', iwrite_projection_ob(1)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_ob(2)', iwrite_projection_ob(2)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_k(1)', iwrite_projection_k(1)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_k(2)', iwrite_projection_k(2)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'filename_pot', filename_pot
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwrite_external', iwrite_external
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_intelectron', iflag_intelectron
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_dip2', num_dip2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2boundary(1)', dip2boundary(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2boundary(2)', dip2boundary(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2center(1)', dip2center(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2center(2)', dip2center(2)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itotntime2', itotntime2
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwdenoption', iwdenoption
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwdenstep', iwdenstep
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_estatic', iflag_estatic


      select case(iflag_atom_coor)
      case(ntype_atom_coor_cartesian)
        write(fh_variables_log, '("#namelist: ",A)') 'atom_coor'
        do i = 1,natom
          write(fh_variables_log, '("#",4X,A,I4,A,"=",3ES14.5)') 'Rion(',i,')', Rion(1:3,i)
        end do
      case(ntype_atom_coor_reduced)
        write(fh_variables_log, '("#namelist: ",A)') 'atom_red_coor'
        do i = 1,natom
          write(fh_variables_log, '("#",4X,A,I4,A,"=",3ES14.5)') 'Rion_red(',i,')', Rion_red(1:3,i)
        end do
      case default
      end select

      if(inml_code >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'code', inml_code
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_want_stencil_hand_vectorization', yn_want_stencil_hand_vectorization
      write(fh_variables_log, '("#",4X,A,"=",A)') 'yn_want_communication_overlapping', yn_want_communication_overlapping
      write(fh_variables_log, '("#",4X,A,"=",A)') 'stencil_openmp_mode', stencil_openmp_mode
      write(fh_variables_log, '("#",4X,A,"=",A)') 'current_openmp_mode', current_openmp_mode
      write(fh_variables_log, '("#",4X,A,"=",A)') 'force_openmp_mode',   force_openmp_mode

      close(fh_variables_log)

    end if

    call comm_bcast(ierr_nml,nproc_group_global)
    if(ierr_nml > 0)then
      if (comm_is_root(nproc_id_global)) write(*,"(I4,2x,A)")ierr_nml,'error(s) in input.'
      call end_parallel
      stop
    end if

    !(create output directory)
    if(base_directory(1:3).ne."./ ") then
      call atomic_create_directory(base_directory &
                                  ,nproc_group_global,nproc_id_global)
    endif

  end subroutine dump_input_common

  subroutine check_bad_input
    use parallelization
    use communication
    implicit none
    integer :: round_phi
    real(8) :: udp_phi  ! udp: under dicimal point

    !! Add wrong input keyword or wrong/unavailable input combinations here
    !! (now only a few)

    ! to correct 'Y' and 'N' to be 'y' and 'n', also error check too
    call yn_argument_check(yn_md)
    call yn_argument_check(yn_opt)
    call yn_argument_check(yn_restart)
    call yn_argument_check(yn_self_checkpoint)
    call yn_argument_check(yn_reset_step_restart)
    call yn_argument_check(yn_ffte)
    call yn_argument_check(yn_scalapack)
    call yn_argument_check(yn_scalapack_red_mem)
    call yn_argument_check(yn_eigenexa)
    call yn_argument_check(yn_periodic)
    call yn_argument_check(yn_psmask)
    call yn_argument_check(yn_fix_func)
    call yn_argument_check(yn_auto_mixing)
    call yn_argument_check(yn_subspace_diagonalization)
    call yn_argument_check(yn_out_psi)
    call yn_argument_check(yn_out_dos)
    call yn_argument_check(yn_out_dos_set_fe_origin)
    call yn_argument_check(yn_out_pdos)
    call yn_argument_check(yn_out_dns)
    call yn_argument_check(yn_out_dns_rt)
    call yn_argument_check(yn_out_dns_ac_je)
    call yn_argument_check(yn_out_dns_trans)
    call yn_argument_check(yn_out_elf)
    call yn_argument_check(yn_out_elf_rt)
    call yn_argument_check(yn_out_estatic_rt)
    call yn_argument_check(yn_out_rvf_rt)
    call yn_argument_check(yn_out_tm)
    call yn_argument_check(yn_set_ini_velocity)
    call yn_argument_check(yn_stop_system_momt)
    call yn_argument_check(yn_want_stencil_hand_vectorization)
    call yn_argument_check(yn_want_communication_overlapping)
    call yn_argument_check(yn_gbp)
    call yn_argument_check(yn_gbp_fourier0) ! temporary

    select case(method_wf_distributor)
    case ('single','slice') ; continue
    case default            ; stop 'method_wf_distributor must be single or slice'
    end select

    select case(method_init_wf)
    case ('gauss','random') ; continue
    case ('gauss2','gauss3','gauss4','gauss5','gauss10') ; continue
    case default            ; stop 'method_init_wf must be gauss or random'
    end select

    select case(convergence)
    case('rho_dne')
      continue
    case('norm_rho','norm_rho_dng')
      if(threshold<-1.d-12)then
        if (comm_is_root(nproc_id_global)) then
          write(*,*) 'set threshold when convergence is norm_rho or norm_rho_dng.'
        endif
        call end_parallel
      end if
    case('norm_pot','norm_pot_dng')
      if(threshold<-1.d-12)then
        if (comm_is_root(nproc_id_global)) then
          write(*,*) 'set threshold when convergence is norm_pot or norm_rho_pot.'
        endif
        call end_parallel
      end if
    case default
      if (comm_is_root(nproc_id_global)) then
        write(*,*) 'check a keyword of convergence.'
      endif
      call end_parallel
    end select

    if(yn_out_dos=='y'.or.yn_out_pdos=='y')then
      select case(out_dos_function)
      case("gaussian","lorentzian")
        continue
      case default
        stop 'set out_dos_meshotd to "gaussian" or "lorentzian"'
      end select
    end if

    if (yn_eigenexa == 'y') then
#ifdef USE_EIGENEXA
#else
      stop 'EigenExa does not supported, please reconfiguration and rebuild it.'
#endif
    end if

    if (yn_scalapack == 'y') then
#ifdef USE_SCALAPACK
#else
      stop 'ScaLAPACK does not supported, please reconfiguration and rebuild it.'
#endif
    end if

    if (yn_eigenexa == 'y' .and. yn_scalapack == 'y') then
      stop "both yn_scalapack and yn_eigenexa is specified 'y'"
    end if

  ! for main_tddft
    select case(theory)
    case('tddft_response','tddft_pulse','single_scale_maxwell_tddft')

      round_phi=int((phi_cep1-0.25d0)*2.d0)
      udp_phi=(phi_cep1-0.25d0)*2.d0-round_phi
      if(ae_shape1=="Ecos2".and.abs(udp_phi)>=1.d-12)then
        stop "phi_cep1 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape1."
      end if

      round_phi=int((phi_cep2-0.25d0)*2.d0)
      udp_phi=(phi_cep2-0.25d0)*2.d0-round_phi
      if(ae_shape2=="Ecos2".and.abs(udp_phi)>=1.d-12)then
        stop "phi_cep2 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape2."
      end if

      select case(ae_shape1)
      case("impulse","Ecos2","Acos2","input")
        continue
      case default
        stop 'set ae_shape1 to "impulse", "Ecos2", or "Acos2"'
      end select

      select case(ae_shape2)
      case("none","impulse","Ecos2","Acos2")
        continue
      case default
        stop 'set ae_shape2 to "none", "impulse", "Ecos2", or "Acos2"'
      end select
    end select

  end subroutine check_bad_input

  subroutine stop_by_bad_input2(inp1,inp2,inp3)
    use parallelization
    use communication
    implicit none
    character(*) :: inp1
    character(*) :: inp2
    character(*),optional :: inp3
    if (comm_is_root(nproc_id_global)) then
      write(*,*) ' Bad input combination: '
      if(present(inp3))then
        write(*,*) ' check keywords of ',trim(inp1),' and ',trim(inp2),' and ',trim(inp3)
      else
        write(*,*) ' check keywords of ',trim(inp1),' and ',trim(inp2)
      end if
    endif
    call end_parallel
    stop
  end subroutine stop_by_bad_input2

  subroutine yn_argument_check(str)
    use misc_routines, only: string_lowercase
    use parallelization, only: end_parallel
    implicit none
    character(*), intent(inout) :: str

    call string_lowercase(str)

    if (str /= 'y' .and. str /= 'n') then
      write (*,*) "Bad input: yn_* option only accepts 'y' or 'n'."
      call end_parallel
      stop
    end if
  end subroutine yn_argument_check

end module inputoutput
