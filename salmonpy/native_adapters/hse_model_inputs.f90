! Value-only fast-PR model settings. Never store native owners here.
module salmonpy_hse_model_inputs
use salmon_global
use inputoutput, only: unit_t, &
  iflag_unit_time, &
  utime_to_au, &
  utime_from_au, &
  iflag_unit_length, &
  ulength_to_au, &
  ulength_from_au, &
  iflag_unit_energy, &
  uenergy_to_au, &
  uenergy_from_au, &
  iflag_unit_charge, &
  ucharge_to_au, &
  ucharge_from_au, &
  t_unit_length, &
  t_unit_length_inv, &
  t_unit_energy, &
  t_unit_energy_inv, &
  t_unit_time, &
  t_unit_time_inv, &
  t_unit_spectrum_dipole, &
  t_unit_spectrum_dipole_square, &
  t_unit_current, &
  t_unit_spectrum_current, &
  t_unit_spectrum_current_square, &
  t_unit_ac, &
  t_unit_elec, &
  t_unit_spectrum_elec, &
  t_unit_spectrum_elec_square, &
  t_unit_polarizability, &
  t_unit_conductivity
implicit none
private
public :: s_hse_model_inputs
type :: s_hse_model_inputs
  integer :: ipsfileform(maxmki)
  character(16) :: ps_format(maxmki)
  integer :: iflag_atom_coor
  logical :: quiet
  character(16) :: calc_mode
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
  real(8) :: pbeh_coulomb_radius
  real(8) :: rvv10_b
  real(8) :: rvv10_c
  character(8) :: rvv10_fft
  integer :: rvv10_nq
  real(8) :: hse_omega
  character(1) :: yn_hse_wannier
  character(1) :: yn_exx_dc_mlwf
  real(8) :: exx_pre_scf_threshold
  integer :: exx_pre_scf_steps
  logical :: exx_pre_scf_active
  character(8) :: exx_ace_support
  character(8) :: exx_pair_screening
  character(8) :: exx_factor_mode
  character(16) :: exx_k_algorithm
  character(120) :: exx_k_algorithm_reason
  character(1) :: yn_exx_spinor_pair_fft
  character(1) :: yn_exx_fft_orbits
  character(1) :: yn_exx_k_pair_symmetry
  character(1) :: yn_hse_sr_pair_symmetry
  real(8) :: hse_sr_tolerance
  real(8) :: exx_pair_tolerance
  character(8) :: exx_kpoint_backend
  character(8) :: exx_local_backend
  integer :: exx_gpu_batch_size
  character(8) :: exx_local_fft
  integer :: exx_factor_exact_interval
  integer :: exx_factor_warmup_steps
  integer :: exx_factor_history_frames
  integer :: exx_factor_n_teach
  integer :: exx_mlwf_interval
  integer :: exx_mlwf_maxiter
  real(8) :: exx_mlwf_radius
  real(8) :: exx_mlwf_norm_fraction
  real(8) :: exx_mlwf_tolerance
  integer :: hse_q_reduction(3)
  integer :: exx_worker_tile_divisor
  integer :: exx_worker_rprime_groups
  integer :: hse_block_rows
  character(64) :: hse_fft_layout
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
  character(32) :: method_singlescale
  real(8) :: cutoff_G2_emfield
  character(1) :: yn_symmetrized_stencil
  character(1) :: yn_put_wall_z_boundary
  real(8) :: wall_height
  real(8) :: wall_width
  character(16) :: fdtddim
  character(16) :: twod_shape
  integer :: nx_m
  integer :: ny_m
  integer :: nz_m
  real(8) :: hx_m
  real(8) :: hy_m
  real(8) :: hz_m
  integer :: nksplit
  integer :: nxysplit
  integer :: nxvacl_m
  integer :: nxvacr_m
  integer :: nxvac_m(1:2)
  integer :: nyvac_m(1:2)
  integer :: nzvac_m(1:2)
  integer :: nx_origin_m
  integer :: ny_origin_m
  integer :: nz_origin_m
  integer :: out_ms_region_ix_m(1:2)
  integer :: out_ms_region_iy_m(1:2)
  integer :: out_ms_region_iz_m(1:2)
  character(256) :: file_macropoint
  character(1) :: set_ini_coor_vel
  integer :: nmacro_write_group
  integer :: nmacro_chunk
  real(8) :: rmat_ms(3, 3)
  real(8) :: al_em(3)
  real(8) :: dl_em(3)
  integer :: num_rgrid_em(3)
  real(8) :: at_em
  real(8) :: dt_em
  integer :: nt_em
  character(8) :: boundary_em(3,2)
  character(256) :: shape_file
  integer :: media_num
  character(16) :: media_type(0:200)
  real(8) :: epsilon_em(0:200)
  real(8) :: mu_em(0:200)
  real(8) :: sigma_em(0:200)
  integer :: pole_num_ld(0:200)
  real(8) :: omega_p_ld(0:200)
  real(8) :: f_ld(0:200,1:100)
  real(8) :: gamma_ld(0:200,1:100)
  real(8) :: omega_ld(0:200,1:100)
  character(16) :: wave_input
  real(8) :: ek_dir1(3)
  real(8) :: source_loc1(3)
  real(8) :: gbeam_sigma_plane1(3)
  real(8) :: gbeam_sigma_line1(3)
  real(8) :: ek_dir2(3)
  real(8) :: source_loc2(3)
  real(8) :: gbeam_sigma_plane2(3)
  real(8) :: gbeam_sigma_line2(3)
  integer :: obs_num_em
  integer :: obs_samp_em
  real(8) :: obs_loc_em(200,3)
  real(8) :: obs_plane_ene_em(200,100)
  character(1) :: yn_obs_plane_em(200)
  character(1) :: yn_obs_plane_integral_em(200)
  character(1) :: yn_wf_em
  real(8) :: film_thickness
  integer :: media_id_pml(3,2)
  integer :: media_id_source1
  integer :: media_id_source2
  real(8) :: bloch_k_em(3)
  character(4) :: bloch_real_imag_em(3)
  integer :: ase_num_em
  real(8) :: ase_ene_min_em
  real(8) :: ase_ene_max_em
  real(8) :: ase_wav_min_em
  real(8) :: ase_wav_max_em
  integer :: ase_smedia_id_em
  real(8) :: ase_box_cent_em(3)
  real(8) :: ase_box_size_em(3)
  integer :: art_num_em
  real(8) :: art_ene_min_em
  real(8) :: art_ene_max_em
  real(8) :: art_wav_min_em
  real(8) :: art_wav_max_em
  integer :: art_smedia_id_em
  real(8) :: art_plane_bot_em(3)
  real(8) :: art_plane_top_em(3)
  character(1) :: yn_make_shape
  character(1) :: yn_output_shape
  character(1) :: yn_copy_x
  character(1) :: yn_copy_y
  character(1) :: yn_copy_z
  character(6) :: rot_type
  integer :: n_s
  character(32) :: typ_s(1000)
  integer :: id_s(1000)
  real(8) :: inf_s(1000,10)
  real(8) :: ori_s(1000,3)
  real(8) :: rot_s(1000,3)
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
  integer :: nopt
  real(8) :: max_step_len_adjust
  real(8) :: convrg_opt_fmax
  character(5) :: method_opt
  real(8) :: step_steep
  real(8) :: step_fire
  character(10) :: ensemble
  character(20) :: thermostat
  integer :: step_velocity_scaling
  integer :: step_update_ps
  real(8) :: temperature0_ion_k
  character(1) :: yn_set_ini_velocity
  character(256) :: file_ini_velocity
  real(8) :: thermostat_tau
  character(1) :: yn_stop_system_momt
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
  character(3) :: lattice
  integer :: nref_band
  real(8) :: tol_esp_diff
  integer :: num_of_segments
  integer :: ndiv_segment(max_num_of_segments)
  real(8) :: kpt(3,max_num_of_segments+1)
  character(1) :: kpt_label(max_num_of_segments+1)
  character(1) :: yn_vnl_correction
  integer :: num_sbe
  character(256) :: sysname_sbe(1:200)
  integer :: nk_sbe(1:200)
  integer :: nstate_sbe(1:200)
  integer :: nelec_sbe(1:200)
  real(8) :: al_sbe(3,200)
  real(8) :: al_vec1_sbe(3,200)
  real(8) :: al_vec2_sbe(3,200)
  real(8) :: al_vec3_sbe(3,200)
  integer :: norder_correction
  character(16) :: gauge_sbe
  real(8) :: t_2
  integer :: am_s
  integer :: num_fragment(3)
  integer :: num_rgrid_buffer(3)
  integer :: nproc_rgrid_tot(3)
  character(256) :: file_atom_coor_frag
  character(1) :: yn_out_dc_fragment_coor
  real(8) :: xi_dc
  character(1) :: yn_dc_lcfo
  character(1) :: yn_dc_lcfo_diag
  character(16) :: lcfo_eigensolver
  integer :: lcfo_diag_chefsi_filter_degree
  integer :: lcfo_diag_chefsi_filter_chunk_size
  integer :: lcfo_diag_chefsi_max_cycle
  real(8) :: lcfo_diag_chefsi_residual_tolerance
  integer :: nstate_frag
  real(8) :: energy_cut
  real(8) :: lambda_cut
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
  integer :: iflag_unit_time
  real(8) :: utime_to_au
  real(8) :: utime_from_au
  integer :: iflag_unit_length
  real(8) :: ulength_to_au
  real(8) :: ulength_from_au
  integer :: iflag_unit_energy
  real(8) :: uenergy_to_au
  real(8) :: uenergy_from_au
  integer :: iflag_unit_charge
  real(8) :: ucharge_to_au
  real(8) :: ucharge_from_au
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
  procedure :: capture => capture_inputs
  procedure :: activate => activate_inputs
end type
contains
subroutine capture_inputs(self)
  class(s_hse_model_inputs),intent(inout) :: self
  self%ipsfileform=ipsfileform
  self%ps_format=ps_format
  self%iflag_atom_coor=iflag_atom_coor
  self%quiet=quiet
  self%calc_mode=calc_mode
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
  self%pbeh_coulomb_radius=pbeh_coulomb_radius
  self%rvv10_b=rvv10_b
  self%rvv10_c=rvv10_c
  self%rvv10_fft=rvv10_fft
  self%rvv10_nq=rvv10_nq
  self%hse_omega=hse_omega
  self%yn_hse_wannier=yn_hse_wannier
  self%yn_exx_dc_mlwf=yn_exx_dc_mlwf
  self%exx_pre_scf_threshold=exx_pre_scf_threshold
  self%exx_pre_scf_steps=exx_pre_scf_steps
  self%exx_pre_scf_active=exx_pre_scf_active
  self%exx_ace_support=exx_ace_support
  self%exx_pair_screening=exx_pair_screening
  self%exx_factor_mode=exx_factor_mode
  self%exx_k_algorithm=exx_k_algorithm
  self%exx_k_algorithm_reason=exx_k_algorithm_reason
  self%yn_exx_spinor_pair_fft=yn_exx_spinor_pair_fft
  self%yn_exx_fft_orbits=yn_exx_fft_orbits
  self%yn_exx_k_pair_symmetry=yn_exx_k_pair_symmetry
  self%yn_hse_sr_pair_symmetry=yn_hse_sr_pair_symmetry
  self%hse_sr_tolerance=hse_sr_tolerance
  self%exx_pair_tolerance=exx_pair_tolerance
  self%exx_kpoint_backend=exx_kpoint_backend
  self%exx_local_backend=exx_local_backend
  self%exx_gpu_batch_size=exx_gpu_batch_size
  self%exx_local_fft=exx_local_fft
  self%exx_factor_exact_interval=exx_factor_exact_interval
  self%exx_factor_warmup_steps=exx_factor_warmup_steps
  self%exx_factor_history_frames=exx_factor_history_frames
  self%exx_factor_n_teach=exx_factor_n_teach
  self%exx_mlwf_interval=exx_mlwf_interval
  self%exx_mlwf_maxiter=exx_mlwf_maxiter
  self%exx_mlwf_radius=exx_mlwf_radius
  self%exx_mlwf_norm_fraction=exx_mlwf_norm_fraction
  self%exx_mlwf_tolerance=exx_mlwf_tolerance
  self%hse_q_reduction=hse_q_reduction
  self%exx_worker_tile_divisor=exx_worker_tile_divisor
  self%exx_worker_rprime_groups=exx_worker_rprime_groups
  self%hse_block_rows=hse_block_rows
  self%hse_fft_layout=hse_fft_layout
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
  self%method_singlescale=method_singlescale
  self%cutoff_G2_emfield=cutoff_G2_emfield
  self%yn_symmetrized_stencil=yn_symmetrized_stencil
  self%yn_put_wall_z_boundary=yn_put_wall_z_boundary
  self%wall_height=wall_height
  self%wall_width=wall_width
  self%fdtddim=fdtddim
  self%twod_shape=twod_shape
  self%nx_m=nx_m
  self%ny_m=ny_m
  self%nz_m=nz_m
  self%hx_m=hx_m
  self%hy_m=hy_m
  self%hz_m=hz_m
  self%nksplit=nksplit
  self%nxysplit=nxysplit
  self%nxvacl_m=nxvacl_m
  self%nxvacr_m=nxvacr_m
  self%nxvac_m=nxvac_m
  self%nyvac_m=nyvac_m
  self%nzvac_m=nzvac_m
  self%nx_origin_m=nx_origin_m
  self%ny_origin_m=ny_origin_m
  self%nz_origin_m=nz_origin_m
  self%out_ms_region_ix_m=out_ms_region_ix_m
  self%out_ms_region_iy_m=out_ms_region_iy_m
  self%out_ms_region_iz_m=out_ms_region_iz_m
  self%file_macropoint=file_macropoint
  self%set_ini_coor_vel=set_ini_coor_vel
  self%nmacro_write_group=nmacro_write_group
  self%nmacro_chunk=nmacro_chunk
  self%rmat_ms=rmat_ms
  self%al_em=al_em
  self%dl_em=dl_em
  self%num_rgrid_em=num_rgrid_em
  self%at_em=at_em
  self%dt_em=dt_em
  self%nt_em=nt_em
  self%boundary_em=boundary_em
  self%shape_file=shape_file
  self%media_num=media_num
  self%media_type=media_type
  self%epsilon_em=epsilon_em
  self%mu_em=mu_em
  self%sigma_em=sigma_em
  self%pole_num_ld=pole_num_ld
  self%omega_p_ld=omega_p_ld
  self%f_ld=f_ld
  self%gamma_ld=gamma_ld
  self%omega_ld=omega_ld
  self%wave_input=wave_input
  self%ek_dir1=ek_dir1
  self%source_loc1=source_loc1
  self%gbeam_sigma_plane1=gbeam_sigma_plane1
  self%gbeam_sigma_line1=gbeam_sigma_line1
  self%ek_dir2=ek_dir2
  self%source_loc2=source_loc2
  self%gbeam_sigma_plane2=gbeam_sigma_plane2
  self%gbeam_sigma_line2=gbeam_sigma_line2
  self%obs_num_em=obs_num_em
  self%obs_samp_em=obs_samp_em
  self%obs_loc_em=obs_loc_em
  self%obs_plane_ene_em=obs_plane_ene_em
  self%yn_obs_plane_em=yn_obs_plane_em
  self%yn_obs_plane_integral_em=yn_obs_plane_integral_em
  self%yn_wf_em=yn_wf_em
  self%film_thickness=film_thickness
  self%media_id_pml=media_id_pml
  self%media_id_source1=media_id_source1
  self%media_id_source2=media_id_source2
  self%bloch_k_em=bloch_k_em
  self%bloch_real_imag_em=bloch_real_imag_em
  self%ase_num_em=ase_num_em
  self%ase_ene_min_em=ase_ene_min_em
  self%ase_ene_max_em=ase_ene_max_em
  self%ase_wav_min_em=ase_wav_min_em
  self%ase_wav_max_em=ase_wav_max_em
  self%ase_smedia_id_em=ase_smedia_id_em
  self%ase_box_cent_em=ase_box_cent_em
  self%ase_box_size_em=ase_box_size_em
  self%art_num_em=art_num_em
  self%art_ene_min_em=art_ene_min_em
  self%art_ene_max_em=art_ene_max_em
  self%art_wav_min_em=art_wav_min_em
  self%art_wav_max_em=art_wav_max_em
  self%art_smedia_id_em=art_smedia_id_em
  self%art_plane_bot_em=art_plane_bot_em
  self%art_plane_top_em=art_plane_top_em
  self%yn_make_shape=yn_make_shape
  self%yn_output_shape=yn_output_shape
  self%yn_copy_x=yn_copy_x
  self%yn_copy_y=yn_copy_y
  self%yn_copy_z=yn_copy_z
  self%rot_type=rot_type
  self%n_s=n_s
  self%typ_s=typ_s
  self%id_s=id_s
  self%inf_s=inf_s
  self%ori_s=ori_s
  self%rot_s=rot_s
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
  self%nopt=nopt
  self%max_step_len_adjust=max_step_len_adjust
  self%convrg_opt_fmax=convrg_opt_fmax
  self%method_opt=method_opt
  self%step_steep=step_steep
  self%step_fire=step_fire
  self%ensemble=ensemble
  self%thermostat=thermostat
  self%step_velocity_scaling=step_velocity_scaling
  self%step_update_ps=step_update_ps
  self%temperature0_ion_k=temperature0_ion_k
  self%yn_set_ini_velocity=yn_set_ini_velocity
  self%file_ini_velocity=file_ini_velocity
  self%thermostat_tau=thermostat_tau
  self%yn_stop_system_momt=yn_stop_system_momt
  self%yn_jm=yn_jm
  self%yn_charge_neutral_jm=yn_charge_neutral_jm
  self%yn_output_dns_jm=yn_output_dns_jm
  self%shape_file_jm=shape_file_jm
  self%num_jm=num_jm
  self%rs_bohr_jm=rs_bohr_jm
  self%sphere_nion_jm=sphere_nion_jm
  self%sphere_loc_jm=sphere_loc_jm
  if(allocated(kion))then
    self%kion=kion
  else if(allocated(self%kion))then
    deallocate(self%kion)
  endif
  if(allocated(Rion))then
    self%Rion=Rion
  else if(allocated(self%Rion))then
    deallocate(self%Rion)
  endif
  if(allocated(Rion_red))then
    self%Rion_red=Rion_red
  else if(allocated(self%Rion_red))then
    deallocate(self%Rion_red)
  endif
  if(allocated(flag_opt_atom))then
    self%flag_opt_atom=flag_opt_atom
  else if(allocated(self%flag_opt_atom))then
    deallocate(self%flag_opt_atom)
  endif
  if(allocated(atom_name))then
    self%atom_name=atom_name
  else if(allocated(self%atom_name))then
    deallocate(self%atom_name)
  endif
  self%yn_want_stencil_hand_vectorization=yn_want_stencil_hand_vectorization
  self%yn_want_communication_overlapping=yn_want_communication_overlapping
  self%stencil_openmp_mode=stencil_openmp_mode
  self%current_openmp_mode=current_openmp_mode
  self%force_openmp_mode=force_openmp_mode
  self%lattice=lattice
  self%nref_band=nref_band
  self%tol_esp_diff=tol_esp_diff
  self%num_of_segments=num_of_segments
  self%ndiv_segment=ndiv_segment
  self%kpt=kpt
  self%kpt_label=kpt_label
  self%yn_vnl_correction=yn_vnl_correction
  self%num_sbe=num_sbe
  self%sysname_sbe=sysname_sbe
  self%nk_sbe=nk_sbe
  self%nstate_sbe=nstate_sbe
  self%nelec_sbe=nelec_sbe
  self%al_sbe=al_sbe
  self%al_vec1_sbe=al_vec1_sbe
  self%al_vec2_sbe=al_vec2_sbe
  self%al_vec3_sbe=al_vec3_sbe
  self%norder_correction=norder_correction
  self%gauge_sbe=gauge_sbe
  self%t_2=t_2
  self%am_s=am_s
  self%num_fragment=num_fragment
  self%num_rgrid_buffer=num_rgrid_buffer
  self%nproc_rgrid_tot=nproc_rgrid_tot
  self%file_atom_coor_frag=file_atom_coor_frag
  self%yn_out_dc_fragment_coor=yn_out_dc_fragment_coor
  self%xi_dc=xi_dc
  self%yn_dc_lcfo=yn_dc_lcfo
  self%yn_dc_lcfo_diag=yn_dc_lcfo_diag
  self%lcfo_eigensolver=lcfo_eigensolver
  self%lcfo_diag_chefsi_filter_degree=lcfo_diag_chefsi_filter_degree
  self%lcfo_diag_chefsi_filter_chunk_size=lcfo_diag_chefsi_filter_chunk_size
  self%lcfo_diag_chefsi_max_cycle=lcfo_diag_chefsi_max_cycle
  self%lcfo_diag_chefsi_residual_tolerance=lcfo_diag_chefsi_residual_tolerance
  self%nstate_frag=nstate_frag
  self%energy_cut=energy_cut
  self%lambda_cut=lambda_cut
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
  self%iflag_unit_time=iflag_unit_time
  self%utime_to_au=utime_to_au
  self%utime_from_au=utime_from_au
  self%iflag_unit_length=iflag_unit_length
  self%ulength_to_au=ulength_to_au
  self%ulength_from_au=ulength_from_au
  self%iflag_unit_energy=iflag_unit_energy
  self%uenergy_to_au=uenergy_to_au
  self%uenergy_from_au=uenergy_from_au
  self%iflag_unit_charge=iflag_unit_charge
  self%ucharge_to_au=ucharge_to_au
  self%ucharge_from_au=ucharge_from_au
  self%t_unit_length=t_unit_length
  self%t_unit_length_inv=t_unit_length_inv
  self%t_unit_energy=t_unit_energy
  self%t_unit_energy_inv=t_unit_energy_inv
  self%t_unit_time=t_unit_time
  self%t_unit_time_inv=t_unit_time_inv
  self%t_unit_spectrum_dipole=t_unit_spectrum_dipole
  self%t_unit_spectrum_dipole_square=t_unit_spectrum_dipole_square
  self%t_unit_current=t_unit_current
  self%t_unit_spectrum_current=t_unit_spectrum_current
  self%t_unit_spectrum_current_square=t_unit_spectrum_current_square
  self%t_unit_ac=t_unit_ac
  self%t_unit_elec=t_unit_elec
  self%t_unit_spectrum_elec=t_unit_spectrum_elec
  self%t_unit_spectrum_elec_square=t_unit_spectrum_elec_square
  self%t_unit_polarizability=t_unit_polarizability
  self%t_unit_conductivity=t_unit_conductivity
end subroutine
subroutine activate_inputs(self)
  class(s_hse_model_inputs),intent(inout) :: self
  ipsfileform=self%ipsfileform
  ps_format=self%ps_format
  iflag_atom_coor=self%iflag_atom_coor
  quiet=self%quiet
  calc_mode=self%calc_mode
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
  pbeh_coulomb_radius=self%pbeh_coulomb_radius
  rvv10_b=self%rvv10_b
  rvv10_c=self%rvv10_c
  rvv10_fft=self%rvv10_fft
  rvv10_nq=self%rvv10_nq
  hse_omega=self%hse_omega
  yn_hse_wannier=self%yn_hse_wannier
  yn_exx_dc_mlwf=self%yn_exx_dc_mlwf
  exx_pre_scf_threshold=self%exx_pre_scf_threshold
  exx_pre_scf_steps=self%exx_pre_scf_steps
  exx_pre_scf_active=self%exx_pre_scf_active
  exx_ace_support=self%exx_ace_support
  exx_pair_screening=self%exx_pair_screening
  exx_factor_mode=self%exx_factor_mode
  exx_k_algorithm=self%exx_k_algorithm
  exx_k_algorithm_reason=self%exx_k_algorithm_reason
  yn_exx_spinor_pair_fft=self%yn_exx_spinor_pair_fft
  yn_exx_fft_orbits=self%yn_exx_fft_orbits
  yn_exx_k_pair_symmetry=self%yn_exx_k_pair_symmetry
  yn_hse_sr_pair_symmetry=self%yn_hse_sr_pair_symmetry
  hse_sr_tolerance=self%hse_sr_tolerance
  exx_pair_tolerance=self%exx_pair_tolerance
  exx_kpoint_backend=self%exx_kpoint_backend
  exx_local_backend=self%exx_local_backend
  exx_gpu_batch_size=self%exx_gpu_batch_size
  exx_local_fft=self%exx_local_fft
  exx_factor_exact_interval=self%exx_factor_exact_interval
  exx_factor_warmup_steps=self%exx_factor_warmup_steps
  exx_factor_history_frames=self%exx_factor_history_frames
  exx_factor_n_teach=self%exx_factor_n_teach
  exx_mlwf_interval=self%exx_mlwf_interval
  exx_mlwf_maxiter=self%exx_mlwf_maxiter
  exx_mlwf_radius=self%exx_mlwf_radius
  exx_mlwf_norm_fraction=self%exx_mlwf_norm_fraction
  exx_mlwf_tolerance=self%exx_mlwf_tolerance
  hse_q_reduction=self%hse_q_reduction
  exx_worker_tile_divisor=self%exx_worker_tile_divisor
  exx_worker_rprime_groups=self%exx_worker_rprime_groups
  hse_block_rows=self%hse_block_rows
  hse_fft_layout=self%hse_fft_layout
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
  method_singlescale=self%method_singlescale
  cutoff_G2_emfield=self%cutoff_G2_emfield
  yn_symmetrized_stencil=self%yn_symmetrized_stencil
  yn_put_wall_z_boundary=self%yn_put_wall_z_boundary
  wall_height=self%wall_height
  wall_width=self%wall_width
  fdtddim=self%fdtddim
  twod_shape=self%twod_shape
  nx_m=self%nx_m
  ny_m=self%ny_m
  nz_m=self%nz_m
  hx_m=self%hx_m
  hy_m=self%hy_m
  hz_m=self%hz_m
  nksplit=self%nksplit
  nxysplit=self%nxysplit
  nxvacl_m=self%nxvacl_m
  nxvacr_m=self%nxvacr_m
  nxvac_m=self%nxvac_m
  nyvac_m=self%nyvac_m
  nzvac_m=self%nzvac_m
  nx_origin_m=self%nx_origin_m
  ny_origin_m=self%ny_origin_m
  nz_origin_m=self%nz_origin_m
  out_ms_region_ix_m=self%out_ms_region_ix_m
  out_ms_region_iy_m=self%out_ms_region_iy_m
  out_ms_region_iz_m=self%out_ms_region_iz_m
  file_macropoint=self%file_macropoint
  set_ini_coor_vel=self%set_ini_coor_vel
  nmacro_write_group=self%nmacro_write_group
  nmacro_chunk=self%nmacro_chunk
  rmat_ms=self%rmat_ms
  al_em=self%al_em
  dl_em=self%dl_em
  num_rgrid_em=self%num_rgrid_em
  at_em=self%at_em
  dt_em=self%dt_em
  nt_em=self%nt_em
  boundary_em=self%boundary_em
  shape_file=self%shape_file
  media_num=self%media_num
  media_type=self%media_type
  epsilon_em=self%epsilon_em
  mu_em=self%mu_em
  sigma_em=self%sigma_em
  pole_num_ld=self%pole_num_ld
  omega_p_ld=self%omega_p_ld
  f_ld=self%f_ld
  gamma_ld=self%gamma_ld
  omega_ld=self%omega_ld
  wave_input=self%wave_input
  ek_dir1=self%ek_dir1
  source_loc1=self%source_loc1
  gbeam_sigma_plane1=self%gbeam_sigma_plane1
  gbeam_sigma_line1=self%gbeam_sigma_line1
  ek_dir2=self%ek_dir2
  source_loc2=self%source_loc2
  gbeam_sigma_plane2=self%gbeam_sigma_plane2
  gbeam_sigma_line2=self%gbeam_sigma_line2
  obs_num_em=self%obs_num_em
  obs_samp_em=self%obs_samp_em
  obs_loc_em=self%obs_loc_em
  obs_plane_ene_em=self%obs_plane_ene_em
  yn_obs_plane_em=self%yn_obs_plane_em
  yn_obs_plane_integral_em=self%yn_obs_plane_integral_em
  yn_wf_em=self%yn_wf_em
  film_thickness=self%film_thickness
  media_id_pml=self%media_id_pml
  media_id_source1=self%media_id_source1
  media_id_source2=self%media_id_source2
  bloch_k_em=self%bloch_k_em
  bloch_real_imag_em=self%bloch_real_imag_em
  ase_num_em=self%ase_num_em
  ase_ene_min_em=self%ase_ene_min_em
  ase_ene_max_em=self%ase_ene_max_em
  ase_wav_min_em=self%ase_wav_min_em
  ase_wav_max_em=self%ase_wav_max_em
  ase_smedia_id_em=self%ase_smedia_id_em
  ase_box_cent_em=self%ase_box_cent_em
  ase_box_size_em=self%ase_box_size_em
  art_num_em=self%art_num_em
  art_ene_min_em=self%art_ene_min_em
  art_ene_max_em=self%art_ene_max_em
  art_wav_min_em=self%art_wav_min_em
  art_wav_max_em=self%art_wav_max_em
  art_smedia_id_em=self%art_smedia_id_em
  art_plane_bot_em=self%art_plane_bot_em
  art_plane_top_em=self%art_plane_top_em
  yn_make_shape=self%yn_make_shape
  yn_output_shape=self%yn_output_shape
  yn_copy_x=self%yn_copy_x
  yn_copy_y=self%yn_copy_y
  yn_copy_z=self%yn_copy_z
  rot_type=self%rot_type
  n_s=self%n_s
  typ_s=self%typ_s
  id_s=self%id_s
  inf_s=self%inf_s
  ori_s=self%ori_s
  rot_s=self%rot_s
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
  nopt=self%nopt
  max_step_len_adjust=self%max_step_len_adjust
  convrg_opt_fmax=self%convrg_opt_fmax
  method_opt=self%method_opt
  step_steep=self%step_steep
  step_fire=self%step_fire
  ensemble=self%ensemble
  thermostat=self%thermostat
  step_velocity_scaling=self%step_velocity_scaling
  step_update_ps=self%step_update_ps
  temperature0_ion_k=self%temperature0_ion_k
  yn_set_ini_velocity=self%yn_set_ini_velocity
  file_ini_velocity=self%file_ini_velocity
  thermostat_tau=self%thermostat_tau
  yn_stop_system_momt=self%yn_stop_system_momt
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
  lattice=self%lattice
  nref_band=self%nref_band
  tol_esp_diff=self%tol_esp_diff
  num_of_segments=self%num_of_segments
  ndiv_segment=self%ndiv_segment
  kpt=self%kpt
  kpt_label=self%kpt_label
  yn_vnl_correction=self%yn_vnl_correction
  num_sbe=self%num_sbe
  sysname_sbe=self%sysname_sbe
  nk_sbe=self%nk_sbe
  nstate_sbe=self%nstate_sbe
  nelec_sbe=self%nelec_sbe
  al_sbe=self%al_sbe
  al_vec1_sbe=self%al_vec1_sbe
  al_vec2_sbe=self%al_vec2_sbe
  al_vec3_sbe=self%al_vec3_sbe
  norder_correction=self%norder_correction
  gauge_sbe=self%gauge_sbe
  t_2=self%t_2
  am_s=self%am_s
  num_fragment=self%num_fragment
  num_rgrid_buffer=self%num_rgrid_buffer
  nproc_rgrid_tot=self%nproc_rgrid_tot
  file_atom_coor_frag=self%file_atom_coor_frag
  yn_out_dc_fragment_coor=self%yn_out_dc_fragment_coor
  xi_dc=self%xi_dc
  yn_dc_lcfo=self%yn_dc_lcfo
  yn_dc_lcfo_diag=self%yn_dc_lcfo_diag
  lcfo_eigensolver=self%lcfo_eigensolver
  lcfo_diag_chefsi_filter_degree=self%lcfo_diag_chefsi_filter_degree
  lcfo_diag_chefsi_filter_chunk_size=self%lcfo_diag_chefsi_filter_chunk_size
  lcfo_diag_chefsi_max_cycle=self%lcfo_diag_chefsi_max_cycle
  lcfo_diag_chefsi_residual_tolerance=self%lcfo_diag_chefsi_residual_tolerance
  nstate_frag=self%nstate_frag
  energy_cut=self%energy_cut
  lambda_cut=self%lambda_cut
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
  iflag_unit_time=self%iflag_unit_time
  utime_to_au=self%utime_to_au
  utime_from_au=self%utime_from_au
  iflag_unit_length=self%iflag_unit_length
  ulength_to_au=self%ulength_to_au
  ulength_from_au=self%ulength_from_au
  iflag_unit_energy=self%iflag_unit_energy
  uenergy_to_au=self%uenergy_to_au
  uenergy_from_au=self%uenergy_from_au
  iflag_unit_charge=self%iflag_unit_charge
  ucharge_to_au=self%ucharge_to_au
  ucharge_from_au=self%ucharge_from_au
  t_unit_length=self%t_unit_length
  t_unit_length_inv=self%t_unit_length_inv
  t_unit_energy=self%t_unit_energy
  t_unit_energy_inv=self%t_unit_energy_inv
  t_unit_time=self%t_unit_time
  t_unit_time_inv=self%t_unit_time_inv
  t_unit_spectrum_dipole=self%t_unit_spectrum_dipole
  t_unit_spectrum_dipole_square=self%t_unit_spectrum_dipole_square
  t_unit_current=self%t_unit_current
  t_unit_spectrum_current=self%t_unit_spectrum_current
  t_unit_spectrum_current_square=self%t_unit_spectrum_current_square
  t_unit_ac=self%t_unit_ac
  t_unit_elec=self%t_unit_elec
  t_unit_spectrum_elec=self%t_unit_spectrum_elec
  t_unit_spectrum_elec_square=self%t_unit_spectrum_elec_square
  t_unit_polarizability=self%t_unit_polarizability
  t_unit_conductivity=self%t_unit_conductivity
end subroutine
end module
