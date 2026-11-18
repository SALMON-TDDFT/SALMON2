! SALMONpy model-owned initialization values for the standard scalar TDDFT profile.
! Input values are captured after parsing/unit conversion; numerical/MPI owners are excluded.
module salmonpy_model_inputs
  use salmon_global, only: maxmki, iflag_atom_coor, &
    theory, &
    yn_md, &
    yn_opt, &
    yn_dc, &
    yn_conventional_from_dcdft, &
    sysname, &
    base_directory, &
    yn_restart, &
    directory_read_data, &
    yn_self_checkpoint, &
    checkpoint_interval, &
    yn_reset_step_restart, &
    read_gs_restart_data, &
    write_gs_restart_data, &
    time_shutdown, &
    method_wf_distributor, &
    nblock_wf_distribute, &
    write_gs_wfn_k, &
    write_rt_wfn_k, &
    unit_system, &
    unit_time, &
    unit_length, &
    unit_energy, &
    unit_charge, &
    nproc_k, &
    nproc_ob, &
    nproc_rgrid, &
    yn_ffte, &
    yn_fftw, &
    yn_scalapack, &
    yn_gramschmidt_blas, &
    yn_eigenexa, &
    yn_diagonalization_red_mem, &
    yn_nccl_reduction, &
    process_allocation, &
    iperiodic, &
    yn_periodic, &
    spin, &
    al, &
    al_vec1, &
    al_vec2, &
    al_vec3, &
    nstate, &
    nelec, &
    nelec_spin, &
    temperature, &
    temperature_k, &
    nelem, &
    natom, &
    file_atom_coor, &
    file_atom_red_coor, &
    yn_spinorbit, &
    yn_symmetry, &
    absorbing_boundary, &
    imagnary_potential_w0, &
    imagnary_potential_dr, &
    file_pseudo, &
    lmax_ps, &
    lloc_ps, &
    izatom, &
    yn_psmask, &
    yn_pseudo_atomic_orbital, &
    alpha_mask, &
    gamma_mask, &
    eta_mask, &
    xc, &
    xname, &
    hse_omega, &
    cname, &
    alibx, &
    alibc, &
    alibxc, &
    cval, &
    dl, &
    num_rgrid, &
    num_kgrid, &
    file_kw, &
    dk_shift, &
    yn_gamma_centered, &
    nt, &
    dt, &
    gram_schmidt_interval, &
    n_hamil, &
    propagator, &
    yn_fix_func, &
    yn_predictor_corrector, &
    method_init_wf, &
    iseed_number_change, &
    method_min, &
    ncg, &
    ncg_init, &
    method_mixing, &
    mixrate, &
    nmemory_mb, &
    alpha_mb, &
    nmemory_p, &
    beta_p, &
    yn_auto_mixing, &
    update_mixing_ratio, &
    nscf, &
    yn_subspace_diagonalization, &
    convergence, &
    threshold, &
    nscf_init_redistribution, &
    nscf_init_no_diagonal, &
    nscf_init_mix_zero, &
    conv_gap_mix_zero, &
    method_init_density, &
    magdir_atom, &
    yn_preconditioning, &
    alpha_pre, &
    trans_longi, &
    ae_shape1, &
    file_input1, &
    e_impulse, &
    E_amplitude1, &
    I_wcm2_1, &
    tw1, &
    omega1, &
    epdir_re1, &
    epdir_im1, &
    phi_cep1, &
    ae_shape2, &
    E_amplitude2, &
    I_wcm2_2, &
    tw2, &
    omega2, &
    epdir_re2, &
    epdir_im2, &
    phi_cep2, &
    t1_t2, &
    t1_start, &
    num_dipole_source, &
    vec_dipole_source, &
    cood_dipole_source, &
    rad_dipole_source, &
    projection_option, &
    threshold_projection, &
    nenergy, &
    de, &
    out_rt_energy_step, &
    yn_out_psi, &
    yn_out_dos, &
    yn_out_dos_set_fe_origin, &
    out_dos_start, &
    out_dos_end, &
    out_dos_nenergy, &
    out_dos_width, &
    out_dos_function, &
    yn_out_pdos, &
    yn_out_dns, &
    yn_out_dns_rt, &
    yn_out_dns_ac_je, &
    yn_out_micro_je, &
    out_dns_rt_step, &
    out_dns_ac_je_step, &
    out_micro_je_step, &
    out_old_dns, &
    yn_out_dns_trans, &
    out_dns_trans_energy, &
    yn_out_elf, &
    yn_out_elf_rt, &
    out_elf_rt_step, &
    yn_out_estatic_rt, &
    out_estatic_rt_step, &
    yn_out_rvf_rt, &
    out_rvf_rt_step, &
    yn_out_tm, &
    yn_out_tm_bin, &
    yn_out_gs_sgm_eps, &
    out_gs_sgm_eps_mu_nu, &
    out_gs_sgm_eps_width, &
    out_projection_step, &
    out_ms_step, &
    format_voxel_data, &
    nsplit_voxel_data, &
    yn_lr_w0_correction, &
    yn_out_intraband_current, &
    yn_out_current_decomposed, &
    out_current_decomposed_step, &
    out_rt_spin_step, &
    yn_out_mag_decomposed_rt, &
    yn_out_mag_micro_rt, &
    yn_out_spin_current_decomposed, &
    yn_out_spin_current_micro, &
    yn_out_rt_energy_components, &
    yn_out_perflog, &
    format_perflog, &
    layout_multipole, &
    num_multipole_xyz, &
    lmax_multipole, &
    threshold_cg, &
    method_poisson, &
    newald, &
    aewald, &
    cutoff_r, &
    cutoff_r_buff, &
    cutoff_g, &
    yn_jm, &
    yn_charge_neutral_jm, &
    yn_output_dns_jm, &
    shape_file_jm, &
    num_jm, &
    rs_bohr_jm, &
    sphere_nion_jm, &
    sphere_loc_jm, &
    kion, &
    Rion, &
    Rion_red, &
    flag_opt_atom, &
    atom_name, &
    yn_want_stencil_hand_vectorization, &
    yn_want_communication_overlapping, &
    stencil_openmp_mode, &
    current_openmp_mode, &
    force_openmp_mode, &
    dm_unfold_option, &
    num_lkgrid, &
    num_skgrid, &
    no_pr, &
    out_dm_unfold_step, &
    yn_out_mom_distr_gs, &
    yn_out_mom_distr_rt, &
    out_mom_distr_rt_step, &
    nq_mom, &
    dq_mom
  implicit none
  private
  type, public :: s_model_inputs
    private
    logical :: ready=.false.
    integer :: iflag_atom_coor
    character(32) :: theory
    character(1) :: yn_md
    character(1) :: yn_opt
    character(1) :: yn_dc
    character(1) :: yn_conventional_from_dcdft
    character(256) :: sysname
    character(256) :: base_directory
    character(1) :: yn_restart
    character(256) :: directory_read_data
    character(1) :: yn_self_checkpoint
    integer :: checkpoint_interval
    character(1) :: yn_reset_step_restart
    character(256) :: read_gs_restart_data
    character(256) :: write_gs_restart_data
    real(8) :: time_shutdown
    character(20) :: method_wf_distributor
    integer :: nblock_wf_distribute
    character(1) :: write_gs_wfn_k
    character(1) :: write_rt_wfn_k
    character(16) :: unit_system
    character(16) :: unit_time
    character(16) :: unit_length
    character(16) :: unit_energy
    character(16) :: unit_charge
    integer :: nproc_k
    integer :: nproc_ob
    integer :: nproc_rgrid(3)
    character(1) :: yn_ffte
    character(1) :: yn_fftw
    character(1) :: yn_scalapack
    character(1) :: yn_gramschmidt_blas
    character(1) :: yn_eigenexa
    character(1) :: yn_diagonalization_red_mem
    character(1) :: yn_nccl_reduction
    character(32) :: process_allocation
    integer :: iperiodic
    character(1) :: yn_periodic
    character(16) :: spin
    real(8) :: al(3)
    real(8) :: al_vec1(3)
    real(8) :: al_vec2(3)
    real(8) :: al_vec3(3)
    integer :: nstate
    integer :: nelec
    integer :: nelec_spin(2)
    real(8) :: temperature
    real(8) :: temperature_k
    integer :: nelem
    integer :: natom
    character(256) :: file_atom_coor
    character(256) :: file_atom_red_coor
    character(1) :: yn_spinorbit
    character(3) :: yn_symmetry
    character(16) :: absorbing_boundary
    real(8) :: imagnary_potential_w0
    real(8) :: imagnary_potential_dr
    character(256) :: file_pseudo(maxmki)
    integer :: lmax_ps(maxmki)
    integer :: lloc_ps(maxmki)
    integer :: izatom(maxmki)
    character(1) :: yn_psmask
    character(1) :: yn_pseudo_atomic_orbital
    real(8) :: alpha_mask
    real(8) :: gamma_mask
    real(8) :: eta_mask
    character(64) :: xc
    character(64) :: xname
    real(8) :: hse_omega
    character(64) :: cname
    character(64) :: alibx
    character(64) :: alibc
    character(64) :: alibxc
    real(8) :: cval
    real(8) :: dl(3)
    integer :: num_rgrid(3)
    integer :: num_kgrid(3)
    character(256) :: file_kw
    real(8) :: dk_shift(3)
    character(1) :: yn_gamma_centered
    integer :: nt
    real(8) :: dt
    integer :: gram_schmidt_interval
    integer :: n_hamil
    character(16) :: propagator
    character(1) :: yn_fix_func
    character(1) :: yn_predictor_corrector
    character(8) :: method_init_wf
    integer :: iseed_number_change
    character(8) :: method_min
    integer :: ncg
    integer :: ncg_init
    character(16) :: method_mixing
    real(8) :: mixrate
    integer :: nmemory_mb
    real(8) :: alpha_mb
    integer :: nmemory_p
    real(8) :: beta_p
    character(1) :: yn_auto_mixing
    real(8) :: update_mixing_ratio
    integer :: nscf
    character(1) :: yn_subspace_diagonalization
    character(16) :: convergence
    real(8) :: threshold
    integer :: nscf_init_redistribution
    integer :: nscf_init_no_diagonal
    integer :: nscf_init_mix_zero
    real(8) :: conv_gap_mix_zero
    character(16) :: method_init_density
    real(8) :: magdir_atom(100)
    character(1) :: yn_preconditioning
    real(8) :: alpha_pre
    character(2) :: trans_longi
    character(16) :: ae_shape1
    character(256) :: file_input1
    real(8) :: e_impulse
    real(8) :: E_amplitude1
    real(8) :: I_wcm2_1
    real(8) :: tw1
    real(8) :: omega1
    real(8) :: epdir_re1(3)
    real(8) :: epdir_im1(3)
    real(8) :: phi_cep1
    character(16) :: ae_shape2
    real(8) :: E_amplitude2
    real(8) :: I_wcm2_2
    real(8) :: tw2
    real(8) :: omega2
    real(8) :: epdir_re2(3)
    real(8) :: epdir_im2(3)
    real(8) :: phi_cep2
    real(8) :: t1_t2
    real(8) :: t1_start
    integer :: num_dipole_source
    real(8) :: vec_dipole_source(3,2)
    real(8) :: cood_dipole_source(3,2)
    real(8) :: rad_dipole_source
    character(2) :: projection_option
    real(8) :: threshold_projection
    integer :: nenergy
    real(8) :: de
    integer :: out_rt_energy_step
    character(1) :: yn_out_psi
    character(1) :: yn_out_dos
    character(1) :: yn_out_dos_set_fe_origin
    real(8) :: out_dos_start
    real(8) :: out_dos_end
    integer :: out_dos_nenergy
    real(8) :: out_dos_width
    character(16) :: out_dos_function
    character(1) :: yn_out_pdos
    character(1) :: yn_out_dns
    character(1) :: yn_out_dns_rt
    character(1) :: yn_out_dns_ac_je
    character(1) :: yn_out_micro_je
    integer :: out_dns_rt_step
    integer :: out_dns_ac_je_step
    integer :: out_micro_je_step
    character(1) :: out_old_dns
    character(1) :: yn_out_dns_trans
    real(8) :: out_dns_trans_energy
    character(1) :: yn_out_elf
    character(1) :: yn_out_elf_rt
    integer :: out_elf_rt_step
    character(1) :: yn_out_estatic_rt
    integer :: out_estatic_rt_step
    character(1) :: yn_out_rvf_rt
    integer :: out_rvf_rt_step
    character(1) :: yn_out_tm
    character(1) :: yn_out_tm_bin
    character(1) :: yn_out_gs_sgm_eps
    integer :: out_gs_sgm_eps_mu_nu(2)
    real(8) :: out_gs_sgm_eps_width
    integer :: out_projection_step
    integer :: out_ms_step
    character(16) :: format_voxel_data
    integer :: nsplit_voxel_data
    character(1) :: yn_lr_w0_correction
    character(1) :: yn_out_intraband_current
    character(1) :: yn_out_current_decomposed
    integer :: out_current_decomposed_step
    integer :: out_rt_spin_step
    character(1) :: yn_out_mag_decomposed_rt
    character(1) :: yn_out_mag_micro_rt
    character(1) :: yn_out_spin_current_decomposed
    character(1) :: yn_out_spin_current_micro
    character(1) :: yn_out_rt_energy_components
    character(1) :: yn_out_perflog
    character(6) :: format_perflog
    integer :: layout_multipole
    integer :: num_multipole_xyz(3)
    integer :: lmax_multipole
    real(8) :: threshold_cg
    character(9) :: method_poisson
    integer :: newald
    real(8) :: aewald
    real(8) :: cutoff_r
    real(8) :: cutoff_r_buff
    real(8) :: cutoff_g
    character(1) :: yn_jm
    character(1) :: yn_charge_neutral_jm
    character(1) :: yn_output_dns_jm
    character(256) :: shape_file_jm
    integer :: num_jm
    real(8) :: rs_bohr_jm(200)
    integer :: sphere_nion_jm(200)
    real(8) :: sphere_loc_jm(200,3)
    integer,allocatable :: kion(:)
    real(8),allocatable :: Rion(:,:)
    real(8),allocatable :: Rion_red(:,:)
    character(1),allocatable :: flag_opt_atom(:)
    character(256),allocatable :: atom_name(:)
    character(1) :: yn_want_stencil_hand_vectorization
    character(1) :: yn_want_communication_overlapping
    character(10) :: stencil_openmp_mode
    character(10) :: current_openmp_mode
    character(10) :: force_openmp_mode
    character(256) :: dm_unfold_option
    integer :: num_lkgrid(3)
    integer :: num_skgrid(3)
    integer :: no_pr
    integer :: out_dm_unfold_step
    character(1) :: yn_out_mom_distr_gs
    character(1) :: yn_out_mom_distr_rt
    integer :: out_mom_distr_rt_step
    integer :: nq_mom
    real(8) :: dq_mom
  contains
    procedure :: capture => capture_model_inputs
    procedure :: activate => activate_model_inputs
  end type
contains
  subroutine capture_model_inputs(self)
    class(s_model_inputs),intent(out) :: self
    self%theory=theory
    self%yn_md=yn_md
    self%yn_opt=yn_opt
    self%yn_dc=yn_dc
    self%yn_conventional_from_dcdft=yn_conventional_from_dcdft
    self%sysname=sysname
    self%base_directory=base_directory
    self%yn_restart=yn_restart
    self%directory_read_data=directory_read_data
    self%yn_self_checkpoint=yn_self_checkpoint
    self%checkpoint_interval=checkpoint_interval
    self%yn_reset_step_restart=yn_reset_step_restart
    self%read_gs_restart_data=read_gs_restart_data
    self%write_gs_restart_data=write_gs_restart_data
    self%time_shutdown=time_shutdown
    self%method_wf_distributor=method_wf_distributor
    self%nblock_wf_distribute=nblock_wf_distribute
    self%write_gs_wfn_k=write_gs_wfn_k
    self%write_rt_wfn_k=write_rt_wfn_k
    self%unit_system=unit_system
    self%unit_time=unit_time
    self%unit_length=unit_length
    self%unit_energy=unit_energy
    self%unit_charge=unit_charge
    self%nproc_k=nproc_k
    self%nproc_ob=nproc_ob
    self%nproc_rgrid=nproc_rgrid
    self%yn_ffte=yn_ffte
    self%yn_fftw=yn_fftw
    self%yn_scalapack=yn_scalapack
    self%yn_gramschmidt_blas=yn_gramschmidt_blas
    self%yn_eigenexa=yn_eigenexa
    self%yn_diagonalization_red_mem=yn_diagonalization_red_mem
    self%yn_nccl_reduction=yn_nccl_reduction
    self%process_allocation=process_allocation
    self%iperiodic=iperiodic
    self%yn_periodic=yn_periodic
    self%spin=spin
    self%al=al
    self%al_vec1=al_vec1
    self%al_vec2=al_vec2
    self%al_vec3=al_vec3
    self%nstate=nstate
    self%nelec=nelec
    self%nelec_spin=nelec_spin
    self%temperature=temperature
    self%temperature_k=temperature_k
    self%nelem=nelem
    self%natom=natom
    self%file_atom_coor=file_atom_coor
    self%file_atom_red_coor=file_atom_red_coor
    self%yn_spinorbit=yn_spinorbit
    self%yn_symmetry=yn_symmetry
    self%absorbing_boundary=absorbing_boundary
    self%imagnary_potential_w0=imagnary_potential_w0
    self%imagnary_potential_dr=imagnary_potential_dr
    self%file_pseudo=file_pseudo
    self%lmax_ps=lmax_ps
    self%lloc_ps=lloc_ps
    self%izatom=izatom
    self%yn_psmask=yn_psmask
    self%yn_pseudo_atomic_orbital=yn_pseudo_atomic_orbital
    self%alpha_mask=alpha_mask
    self%gamma_mask=gamma_mask
    self%eta_mask=eta_mask
    self%xc=xc
    self%xname=xname
    self%hse_omega=hse_omega
    self%cname=cname
    self%alibx=alibx
    self%alibc=alibc
    self%alibxc=alibxc
    self%cval=cval
    self%dl=dl
    self%num_rgrid=num_rgrid
    self%num_kgrid=num_kgrid
    self%file_kw=file_kw
    self%dk_shift=dk_shift
    self%yn_gamma_centered=yn_gamma_centered
    self%nt=nt
    self%dt=dt
    self%gram_schmidt_interval=gram_schmidt_interval
    self%n_hamil=n_hamil
    self%propagator=propagator
    self%yn_fix_func=yn_fix_func
    self%yn_predictor_corrector=yn_predictor_corrector
    self%method_init_wf=method_init_wf
    self%iseed_number_change=iseed_number_change
    self%method_min=method_min
    self%ncg=ncg
    self%ncg_init=ncg_init
    self%method_mixing=method_mixing
    self%mixrate=mixrate
    self%nmemory_mb=nmemory_mb
    self%alpha_mb=alpha_mb
    self%nmemory_p=nmemory_p
    self%beta_p=beta_p
    self%yn_auto_mixing=yn_auto_mixing
    self%update_mixing_ratio=update_mixing_ratio
    self%nscf=nscf
    self%yn_subspace_diagonalization=yn_subspace_diagonalization
    self%convergence=convergence
    self%threshold=threshold
    self%nscf_init_redistribution=nscf_init_redistribution
    self%nscf_init_no_diagonal=nscf_init_no_diagonal
    self%nscf_init_mix_zero=nscf_init_mix_zero
    self%conv_gap_mix_zero=conv_gap_mix_zero
    self%method_init_density=method_init_density
    self%magdir_atom=magdir_atom
    self%yn_preconditioning=yn_preconditioning
    self%alpha_pre=alpha_pre
    self%trans_longi=trans_longi
    self%ae_shape1=ae_shape1
    self%file_input1=file_input1
    self%e_impulse=e_impulse
    self%E_amplitude1=E_amplitude1
    self%I_wcm2_1=I_wcm2_1
    self%tw1=tw1
    self%omega1=omega1
    self%epdir_re1=epdir_re1
    self%epdir_im1=epdir_im1
    self%phi_cep1=phi_cep1
    self%ae_shape2=ae_shape2
    self%E_amplitude2=E_amplitude2
    self%I_wcm2_2=I_wcm2_2
    self%tw2=tw2
    self%omega2=omega2
    self%epdir_re2=epdir_re2
    self%epdir_im2=epdir_im2
    self%phi_cep2=phi_cep2
    self%t1_t2=t1_t2
    self%t1_start=t1_start
    self%num_dipole_source=num_dipole_source
    self%vec_dipole_source=vec_dipole_source
    self%cood_dipole_source=cood_dipole_source
    self%rad_dipole_source=rad_dipole_source
    self%projection_option=projection_option
    self%threshold_projection=threshold_projection
    self%nenergy=nenergy
    self%de=de
    self%out_rt_energy_step=out_rt_energy_step
    self%yn_out_psi=yn_out_psi
    self%yn_out_dos=yn_out_dos
    self%yn_out_dos_set_fe_origin=yn_out_dos_set_fe_origin
    self%out_dos_start=out_dos_start
    self%out_dos_end=out_dos_end
    self%out_dos_nenergy=out_dos_nenergy
    self%out_dos_width=out_dos_width
    self%out_dos_function=out_dos_function
    self%yn_out_pdos=yn_out_pdos
    self%yn_out_dns=yn_out_dns
    self%yn_out_dns_rt=yn_out_dns_rt
    self%yn_out_dns_ac_je=yn_out_dns_ac_je
    self%yn_out_micro_je=yn_out_micro_je
    self%out_dns_rt_step=out_dns_rt_step
    self%out_dns_ac_je_step=out_dns_ac_je_step
    self%out_micro_je_step=out_micro_je_step
    self%out_old_dns=out_old_dns
    self%yn_out_dns_trans=yn_out_dns_trans
    self%out_dns_trans_energy=out_dns_trans_energy
    self%yn_out_elf=yn_out_elf
    self%yn_out_elf_rt=yn_out_elf_rt
    self%out_elf_rt_step=out_elf_rt_step
    self%yn_out_estatic_rt=yn_out_estatic_rt
    self%out_estatic_rt_step=out_estatic_rt_step
    self%yn_out_rvf_rt=yn_out_rvf_rt
    self%out_rvf_rt_step=out_rvf_rt_step
    self%yn_out_tm=yn_out_tm
    self%yn_out_tm_bin=yn_out_tm_bin
    self%yn_out_gs_sgm_eps=yn_out_gs_sgm_eps
    self%out_gs_sgm_eps_mu_nu=out_gs_sgm_eps_mu_nu
    self%out_gs_sgm_eps_width=out_gs_sgm_eps_width
    self%out_projection_step=out_projection_step
    self%out_ms_step=out_ms_step
    self%format_voxel_data=format_voxel_data
    self%nsplit_voxel_data=nsplit_voxel_data
    self%yn_lr_w0_correction=yn_lr_w0_correction
    self%yn_out_intraband_current=yn_out_intraband_current
    self%yn_out_current_decomposed=yn_out_current_decomposed
    self%out_current_decomposed_step=out_current_decomposed_step
    self%out_rt_spin_step=out_rt_spin_step
    self%yn_out_mag_decomposed_rt=yn_out_mag_decomposed_rt
    self%yn_out_mag_micro_rt=yn_out_mag_micro_rt
    self%yn_out_spin_current_decomposed=yn_out_spin_current_decomposed
    self%yn_out_spin_current_micro=yn_out_spin_current_micro
    self%yn_out_rt_energy_components=yn_out_rt_energy_components
    self%yn_out_perflog=yn_out_perflog
    self%format_perflog=format_perflog
    self%layout_multipole=layout_multipole
    self%num_multipole_xyz=num_multipole_xyz
    self%lmax_multipole=lmax_multipole
    self%threshold_cg=threshold_cg
    self%method_poisson=method_poisson
    self%newald=newald
    self%aewald=aewald
    self%cutoff_r=cutoff_r
    self%cutoff_r_buff=cutoff_r_buff
    self%cutoff_g=cutoff_g
    self%yn_jm=yn_jm
    self%yn_charge_neutral_jm=yn_charge_neutral_jm
    self%yn_output_dns_jm=yn_output_dns_jm
    self%shape_file_jm=shape_file_jm
    self%num_jm=num_jm
    self%rs_bohr_jm=rs_bohr_jm
    self%sphere_nion_jm=sphere_nion_jm
    self%sphere_loc_jm=sphere_loc_jm
    if(allocated(kion))self%kion=kion
    if(allocated(Rion))self%Rion=Rion
    if(allocated(Rion_red))self%Rion_red=Rion_red
    if(allocated(flag_opt_atom))self%flag_opt_atom=flag_opt_atom
    if(allocated(atom_name))self%atom_name=atom_name
    self%yn_want_stencil_hand_vectorization=yn_want_stencil_hand_vectorization
    self%yn_want_communication_overlapping=yn_want_communication_overlapping
    self%stencil_openmp_mode=stencil_openmp_mode
    self%current_openmp_mode=current_openmp_mode
    self%force_openmp_mode=force_openmp_mode
    self%dm_unfold_option=dm_unfold_option
    self%num_lkgrid=num_lkgrid
    self%num_skgrid=num_skgrid
    self%no_pr=no_pr
    self%out_dm_unfold_step=out_dm_unfold_step
    self%yn_out_mom_distr_gs=yn_out_mom_distr_gs
    self%yn_out_mom_distr_rt=yn_out_mom_distr_rt
    self%out_mom_distr_rt_step=out_mom_distr_rt_step
    self%nq_mom=nq_mom
    self%dq_mom=dq_mom
    self%iflag_atom_coor=iflag_atom_coor
    self%ready=.true.
  end subroutine
  subroutine activate_model_inputs(self,status)
    class(s_model_inputs),intent(in) :: self
    integer,intent(out) :: status
    status=1
    if(.not.self%ready)return
    theory=self%theory
    yn_md=self%yn_md
    yn_opt=self%yn_opt
    yn_dc=self%yn_dc
    yn_conventional_from_dcdft=self%yn_conventional_from_dcdft
    sysname=self%sysname
    base_directory=self%base_directory
    yn_restart=self%yn_restart
    directory_read_data=self%directory_read_data
    yn_self_checkpoint=self%yn_self_checkpoint
    checkpoint_interval=self%checkpoint_interval
    yn_reset_step_restart=self%yn_reset_step_restart
    read_gs_restart_data=self%read_gs_restart_data
    write_gs_restart_data=self%write_gs_restart_data
    time_shutdown=self%time_shutdown
    method_wf_distributor=self%method_wf_distributor
    nblock_wf_distribute=self%nblock_wf_distribute
    write_gs_wfn_k=self%write_gs_wfn_k
    write_rt_wfn_k=self%write_rt_wfn_k
    unit_system=self%unit_system
    unit_time=self%unit_time
    unit_length=self%unit_length
    unit_energy=self%unit_energy
    unit_charge=self%unit_charge
    nproc_k=self%nproc_k
    nproc_ob=self%nproc_ob
    nproc_rgrid=self%nproc_rgrid
    yn_ffte=self%yn_ffte
    yn_fftw=self%yn_fftw
    yn_scalapack=self%yn_scalapack
    yn_gramschmidt_blas=self%yn_gramschmidt_blas
    yn_eigenexa=self%yn_eigenexa
    yn_diagonalization_red_mem=self%yn_diagonalization_red_mem
    yn_nccl_reduction=self%yn_nccl_reduction
    process_allocation=self%process_allocation
    iperiodic=self%iperiodic
    yn_periodic=self%yn_periodic
    spin=self%spin
    al=self%al
    al_vec1=self%al_vec1
    al_vec2=self%al_vec2
    al_vec3=self%al_vec3
    nstate=self%nstate
    nelec=self%nelec
    nelec_spin=self%nelec_spin
    temperature=self%temperature
    temperature_k=self%temperature_k
    nelem=self%nelem
    natom=self%natom
    file_atom_coor=self%file_atom_coor
    file_atom_red_coor=self%file_atom_red_coor
    yn_spinorbit=self%yn_spinorbit
    yn_symmetry=self%yn_symmetry
    absorbing_boundary=self%absorbing_boundary
    imagnary_potential_w0=self%imagnary_potential_w0
    imagnary_potential_dr=self%imagnary_potential_dr
    file_pseudo=self%file_pseudo
    lmax_ps=self%lmax_ps
    lloc_ps=self%lloc_ps
    izatom=self%izatom
    yn_psmask=self%yn_psmask
    yn_pseudo_atomic_orbital=self%yn_pseudo_atomic_orbital
    alpha_mask=self%alpha_mask
    gamma_mask=self%gamma_mask
    eta_mask=self%eta_mask
    xc=self%xc
    xname=self%xname
    hse_omega=self%hse_omega
    cname=self%cname
    alibx=self%alibx
    alibc=self%alibc
    alibxc=self%alibxc
    cval=self%cval
    dl=self%dl
    num_rgrid=self%num_rgrid
    num_kgrid=self%num_kgrid
    file_kw=self%file_kw
    dk_shift=self%dk_shift
    yn_gamma_centered=self%yn_gamma_centered
    nt=self%nt
    dt=self%dt
    gram_schmidt_interval=self%gram_schmidt_interval
    n_hamil=self%n_hamil
    propagator=self%propagator
    yn_fix_func=self%yn_fix_func
    yn_predictor_corrector=self%yn_predictor_corrector
    method_init_wf=self%method_init_wf
    iseed_number_change=self%iseed_number_change
    method_min=self%method_min
    ncg=self%ncg
    ncg_init=self%ncg_init
    method_mixing=self%method_mixing
    mixrate=self%mixrate
    nmemory_mb=self%nmemory_mb
    alpha_mb=self%alpha_mb
    nmemory_p=self%nmemory_p
    beta_p=self%beta_p
    yn_auto_mixing=self%yn_auto_mixing
    update_mixing_ratio=self%update_mixing_ratio
    nscf=self%nscf
    yn_subspace_diagonalization=self%yn_subspace_diagonalization
    convergence=self%convergence
    threshold=self%threshold
    nscf_init_redistribution=self%nscf_init_redistribution
    nscf_init_no_diagonal=self%nscf_init_no_diagonal
    nscf_init_mix_zero=self%nscf_init_mix_zero
    conv_gap_mix_zero=self%conv_gap_mix_zero
    method_init_density=self%method_init_density
    magdir_atom=self%magdir_atom
    yn_preconditioning=self%yn_preconditioning
    alpha_pre=self%alpha_pre
    trans_longi=self%trans_longi
    ae_shape1=self%ae_shape1
    file_input1=self%file_input1
    e_impulse=self%e_impulse
    E_amplitude1=self%E_amplitude1
    I_wcm2_1=self%I_wcm2_1
    tw1=self%tw1
    omega1=self%omega1
    epdir_re1=self%epdir_re1
    epdir_im1=self%epdir_im1
    phi_cep1=self%phi_cep1
    ae_shape2=self%ae_shape2
    E_amplitude2=self%E_amplitude2
    I_wcm2_2=self%I_wcm2_2
    tw2=self%tw2
    omega2=self%omega2
    epdir_re2=self%epdir_re2
    epdir_im2=self%epdir_im2
    phi_cep2=self%phi_cep2
    t1_t2=self%t1_t2
    t1_start=self%t1_start
    num_dipole_source=self%num_dipole_source
    vec_dipole_source=self%vec_dipole_source
    cood_dipole_source=self%cood_dipole_source
    rad_dipole_source=self%rad_dipole_source
    projection_option=self%projection_option
    threshold_projection=self%threshold_projection
    nenergy=self%nenergy
    de=self%de
    out_rt_energy_step=self%out_rt_energy_step
    yn_out_psi=self%yn_out_psi
    yn_out_dos=self%yn_out_dos
    yn_out_dos_set_fe_origin=self%yn_out_dos_set_fe_origin
    out_dos_start=self%out_dos_start
    out_dos_end=self%out_dos_end
    out_dos_nenergy=self%out_dos_nenergy
    out_dos_width=self%out_dos_width
    out_dos_function=self%out_dos_function
    yn_out_pdos=self%yn_out_pdos
    yn_out_dns=self%yn_out_dns
    yn_out_dns_rt=self%yn_out_dns_rt
    yn_out_dns_ac_je=self%yn_out_dns_ac_je
    yn_out_micro_je=self%yn_out_micro_je
    out_dns_rt_step=self%out_dns_rt_step
    out_dns_ac_je_step=self%out_dns_ac_je_step
    out_micro_je_step=self%out_micro_je_step
    out_old_dns=self%out_old_dns
    yn_out_dns_trans=self%yn_out_dns_trans
    out_dns_trans_energy=self%out_dns_trans_energy
    yn_out_elf=self%yn_out_elf
    yn_out_elf_rt=self%yn_out_elf_rt
    out_elf_rt_step=self%out_elf_rt_step
    yn_out_estatic_rt=self%yn_out_estatic_rt
    out_estatic_rt_step=self%out_estatic_rt_step
    yn_out_rvf_rt=self%yn_out_rvf_rt
    out_rvf_rt_step=self%out_rvf_rt_step
    yn_out_tm=self%yn_out_tm
    yn_out_tm_bin=self%yn_out_tm_bin
    yn_out_gs_sgm_eps=self%yn_out_gs_sgm_eps
    out_gs_sgm_eps_mu_nu=self%out_gs_sgm_eps_mu_nu
    out_gs_sgm_eps_width=self%out_gs_sgm_eps_width
    out_projection_step=self%out_projection_step
    out_ms_step=self%out_ms_step
    format_voxel_data=self%format_voxel_data
    nsplit_voxel_data=self%nsplit_voxel_data
    yn_lr_w0_correction=self%yn_lr_w0_correction
    yn_out_intraband_current=self%yn_out_intraband_current
    yn_out_current_decomposed=self%yn_out_current_decomposed
    out_current_decomposed_step=self%out_current_decomposed_step
    out_rt_spin_step=self%out_rt_spin_step
    yn_out_mag_decomposed_rt=self%yn_out_mag_decomposed_rt
    yn_out_mag_micro_rt=self%yn_out_mag_micro_rt
    yn_out_spin_current_decomposed=self%yn_out_spin_current_decomposed
    yn_out_spin_current_micro=self%yn_out_spin_current_micro
    yn_out_rt_energy_components=self%yn_out_rt_energy_components
    yn_out_perflog=self%yn_out_perflog
    format_perflog=self%format_perflog
    layout_multipole=self%layout_multipole
    num_multipole_xyz=self%num_multipole_xyz
    lmax_multipole=self%lmax_multipole
    threshold_cg=self%threshold_cg
    method_poisson=self%method_poisson
    newald=self%newald
    aewald=self%aewald
    cutoff_r=self%cutoff_r
    cutoff_r_buff=self%cutoff_r_buff
    cutoff_g=self%cutoff_g
    yn_jm=self%yn_jm
    yn_charge_neutral_jm=self%yn_charge_neutral_jm
    yn_output_dns_jm=self%yn_output_dns_jm
    shape_file_jm=self%shape_file_jm
    num_jm=self%num_jm
    rs_bohr_jm=self%rs_bohr_jm
    sphere_nion_jm=self%sphere_nion_jm
    sphere_loc_jm=self%sphere_loc_jm
    if(allocated(self%kion))then
      kion=self%kion
    else if(allocated(kion))then
      deallocate(kion)
    endif
    if(allocated(self%Rion))then
      Rion=self%Rion
    else if(allocated(Rion))then
      deallocate(Rion)
    endif
    if(allocated(self%Rion_red))then
      Rion_red=self%Rion_red
    else if(allocated(Rion_red))then
      deallocate(Rion_red)
    endif
    if(allocated(self%flag_opt_atom))then
      flag_opt_atom=self%flag_opt_atom
    else if(allocated(flag_opt_atom))then
      deallocate(flag_opt_atom)
    endif
    if(allocated(self%atom_name))then
      atom_name=self%atom_name
    else if(allocated(atom_name))then
      deallocate(atom_name)
    endif
    yn_want_stencil_hand_vectorization=self%yn_want_stencil_hand_vectorization
    yn_want_communication_overlapping=self%yn_want_communication_overlapping
    stencil_openmp_mode=self%stencil_openmp_mode
    current_openmp_mode=self%current_openmp_mode
    force_openmp_mode=self%force_openmp_mode
    dm_unfold_option=self%dm_unfold_option
    num_lkgrid=self%num_lkgrid
    num_skgrid=self%num_skgrid
    no_pr=self%no_pr
    out_dm_unfold_step=self%out_dm_unfold_step
    yn_out_mom_distr_gs=self%yn_out_mom_distr_gs
    yn_out_mom_distr_rt=self%yn_out_mom_distr_rt
    out_mom_distr_rt_step=self%out_mom_distr_rt_step
    nq_mom=self%nq_mom
    dq_mom=self%dq_mom
    iflag_atom_coor=self%iflag_atom_coor
    status=0
  end subroutine
end module
