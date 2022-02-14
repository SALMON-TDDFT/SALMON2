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
module salmon_global
  implicit none

!Parameters for pseudo-potential
  integer, parameter :: maxmki=10
   !shinohara
  integer :: ipsfileform(maxmki)   ! file format for pseudo potential
  character(16)  :: ps_format(maxmki)
! List of pseudopotential file formats
  integer,parameter :: n_Yabana_Bertsch_psformat = 1 !.rps
  integer,parameter :: n_ABINIT_psformat = 2 ! .pspnc
  integer,parameter :: n_FHI_psformat = 3 ! .cpi
  integer,parameter :: n_ABINITFHI_psformat = 4 ! .fhi

! Flag for atomic coordinate type
  integer :: iflag_atom_coor
  integer,parameter :: ntype_atom_coor_none      = 0
  integer,parameter :: ntype_atom_coor_cartesian = 1
  integer,parameter :: ntype_atom_coor_reduced   = 2

! Flag for suppress standard outputs 
  logical :: quiet

! For band
  integer,parameter :: max_num_of_segments = 10

  character(16)  :: calc_mode      !old input variable, but used as a flag; move later

!Input variables
!! &calculation
  character(32)  :: theory
  character(1)   :: yn_md
  character(1)   :: yn_opt
!! &control
  character(256) :: sysname
  character(256) :: base_directory
  character(1)   :: yn_restart
  character(256) :: directory_read_data
  character(1)   :: yn_self_checkpoint
  integer        :: checkpoint_interval
  character(1)   :: yn_reset_step_restart
  character(256) :: read_gs_restart_data
  character(256) :: write_gs_restart_data
  real(8)        :: time_shutdown
  character(20)  :: method_wf_distributor
  integer        :: nblock_wf_distribute
  !remove later
  character(1)   :: read_gs_dns_cube
  character(1)   :: write_gs_wfn_k
  character(1)   :: write_rt_wfn_k

!! &units
  character(16)  :: unit_system
  character(16)  :: unit_time
  character(16)  :: unit_length
  character(16)  :: unit_energy
  character(16)  :: unit_charge

!! &parallel
  integer        :: nproc_k
  integer        :: nproc_ob
  integer        :: nproc_rgrid(3)
  character(1)   :: yn_ffte
#ifdef USE_FFTW
  character(1)   :: yn_fftw
#endif
  character(1)   :: yn_scalapack
  character(1)   :: yn_gramschmidt_blas
  character(1)   :: yn_eigenexa
  character(1)   :: yn_diagonalization_red_mem
  character(32)  :: process_allocation

!! &system
  integer        :: iperiodic  !this is old keyword but still defined here
  character(1)   :: yn_periodic
  character(16)  :: spin
  real(8)        :: al(3)
  real(8)        :: al_vec1(3),al_vec2(3),al_vec3(3)
  integer        :: nstate
  integer        :: nelec
  integer        :: nelec_spin(2)
  real(8)        :: temperature
  real(8)        :: temperature_k
  integer        :: nelem
  integer        :: natom
  character(256) :: file_atom_coor
  character(256) :: file_atom_red_coor
  character(1)   :: yn_spinorbit
  character(3)   :: yn_symmetry
  character(16)  :: absorbing_boundary
  real(8)        :: imagnary_potential_w0
  real(8)        :: imagnary_potential_dr

!! &pseudo
  character(256) :: file_pseudo(maxmki)
  integer        :: lmax_ps(maxmki)
  integer        :: lloc_ps(maxmki)
  integer        :: izatom(maxmki)
  character(1)   :: yn_psmask
  real(8)        :: alpha_mask
  real(8)        :: gamma_mask
  real(8)        :: eta_mask

!! &functional
  character(64)  :: xc !, xcname
  character(64)  :: xname
  character(64)  :: cname
  character(64)  :: alibx
  character(64)  :: alibc
  character(64)  :: alibxc
  real(8)        :: cval

!! &rgrid
  real(8)        :: dl(3)
  integer        :: num_rgrid(3)

!! &kgrid
  integer        :: num_kgrid(3)
  character(256) :: file_kw

!! &tgrid
  integer        :: nt
  real(8)        :: dt
  integer        :: gram_schmidt_interval

!! &propagation
  integer        :: n_hamil
  character(16)  :: propagator
  character(1)   :: yn_fix_func
  character(1)   :: yn_predictor_corrector

!! &scf
  character(8)   :: method_init_wf
  integer        :: iseed_number_change
  character(8)   :: method_min
  integer        :: ncg,ncg_init
  character(16)  :: method_mixing
  real(8)        :: mixrate
  integer        :: nmemory_mb
  real(8)        :: alpha_mb
  integer        :: nmemory_p
  real(8)        :: beta_p
  character(1)   :: yn_auto_mixing
  real(8)        :: update_mixing_ratio
  integer        :: nscf
  character(1)   :: yn_subspace_diagonalization
  character(16)  :: convergence
  real(8)        :: threshold
  integer        :: nscf_init_redistribution
  integer        :: nscf_init_no_diagonal
  integer        :: nscf_init_mix_zero
  real(8)        :: conv_gap_mix_zero
  character(16)  :: method_init_density
  real(8)        :: magdir_atom(100)

!! &emfield
  character(2)   :: trans_longi
  character(16)  :: ae_shape1
  character(256) :: file_input1
  real(8)        :: e_impulse
  real(8)        :: E_amplitude1
  real(8)        :: I_wcm2_1
  real(8)        :: tw1
  real(8)        :: omega1
  real(8)        :: epdir_re1(3)
  real(8)        :: epdir_im1(3)
  real(8)        :: phi_cep1
  character(16)  :: ae_shape2
  real(8)        :: E_amplitude2
  real(8)        :: I_wcm2_2
  real(8)        :: tw2
  real(8)        :: omega2
  real(8)        :: epdir_re2(3)
  real(8)        :: epdir_im2(3)
  real(8)        :: phi_cep2
  real(8)        :: t1_t2
  real(8)        :: t1_start
  integer        :: num_dipole_source
  real(8)        :: vec_dipole_source(3,2)
  real(8)        :: cood_dipole_source(3,2)
  real(8)        :: rad_dipole_source

!! &singlescale
  character(32)  :: method_singlescale
  real(8)        :: cutoff_G2_emfield
  character(1)   :: yn_symmetrized_stencil
  character(1)   :: yn_put_wall_z_boundary
  real(8)        :: wall_height
  real(8)        :: wall_width

!! &multiscale
  character(16)  :: fdtddim
  character(16)  :: twod_shape
  integer        :: nx_m
  integer        :: ny_m
  integer        :: nz_m
  real(8)        :: hx_m
  real(8)        :: hy_m
  real(8)        :: hz_m
  integer        :: nksplit !! TODO: remove this variable
  integer        :: nxysplit !! TODO: remove this variable
  ! The input variables nxvac(l|r)_m do not recommend to use,
  ! However I tempolary remain them for the reason of the compatibility.
  ! Please use  n(x|y|z)_origin_m to provide the same functionality.
  integer        :: nxvacl_m
  integer        :: nxvacr_m
  integer        :: nx_origin_m
  integer        :: ny_origin_m
  integer        :: nz_origin_m
  character(100) :: file_macropoint
  character(1)   :: set_ini_coor_vel
  integer        :: nmacro_write_group
  integer        :: nmacro_chunk

!! &maxwell
  real(8)        :: al_em(3)
  real(8)        :: dl_em(3)
  integer        :: num_rgrid_em(3)
  real(8)        :: dt_em
  integer        :: nt_em
  character(8)   :: boundary_em(3,2)
  character(256) :: shape_file
  integer        :: media_num
  character(16)  :: media_type(0:200)
  real(8)        :: epsilon_em(0:200)
  real(8)        :: mu_em(0:200)
  real(8)        :: sigma_em(0:200)
  integer        :: pole_num_ld(0:200)
  real(8)        :: omega_p_ld(0:200)
  real(8)        :: f_ld(0:200,1:100)
  real(8)        :: gamma_ld(0:200,1:100)
  real(8)        :: omega_ld(0:200,1:100)
  character(16)  :: wave_input
  real(8)        :: ek_dir1(3)
  real(8)        :: source_loc1(3)
  real(8)        :: ek_dir2(3)
  real(8)        :: source_loc2(3)
  integer        :: obs_num_em
  integer        :: obs_samp_em
  real(8)        :: obs_loc_em(200,3)
  real(8)        :: obs_plane_ene_em(200,100)
  character(1)   :: yn_obs_plane_em(200)
  character(1)   :: yn_obs_plane_integral_em(200)
  character(1)   :: yn_wf_em
  real(8)        :: film_thickness
  integer        :: media_id_pml(3,2)
  integer        :: media_id_source1
  integer        :: media_id_source2
  character(1)   :: yn_make_shape
  character(1)   :: yn_output_shape
  character(1)   :: yn_copy_x
  character(1)   :: yn_copy_y
  character(1)   :: yn_copy_z
  character(6)   :: rot_type
  integer        :: n_s
  character(32)  :: typ_s(1000)
  integer        :: id_s(1000)
  real(8)        :: inf_s(1000,10)
  real(8)        :: ori_s(1000,3)
  real(8)        :: rot_s(1000,3)
  
!! &analysis
  character(2)   :: projection_option
  integer        :: nenergy
  real(8)        :: de
  integer        :: out_rt_energy_step
  character(1)   :: yn_out_psi
  character(1)   :: yn_out_dos
  character(1)   :: yn_out_dos_set_fe_origin
  real(8)        :: out_dos_start
  real(8)        :: out_dos_end
  integer        :: out_dos_nenergy
  real(8)        :: out_dos_width
  character(16)  :: out_dos_function
  character(1)   :: yn_out_pdos
  character(1)   :: yn_out_dns
  character(1)   :: yn_out_dns_rt
  character(1)   :: yn_out_dns_ac_je
  integer        :: out_dns_rt_step
  integer        :: out_dns_ac_je_step
  character(1)   :: out_old_dns
  character(1)   :: yn_out_dns_trans
  real(8)        :: out_dns_trans_energy
  character(1)   :: yn_out_elf
  character(1)   :: yn_out_elf_rt
  integer        :: out_elf_rt_step
  character(1)   :: yn_out_estatic_rt
  integer        :: out_estatic_rt_step
  character(1)   :: yn_out_rvf_rt
  integer        :: out_rvf_rt_step
  character(1)   :: yn_out_tm
  character(1)   :: yn_out_gs_sgm_eps
  integer        :: out_gs_sgm_eps_mu_nu(2)
  real(8)        :: out_gs_sgm_eps_width
  integer        :: out_projection_step
  integer        :: out_ms_step
  character(16)  :: format_voxel_data
  integer        :: nsplit_voxel_data
  character(1)   :: yn_lr_w0_correction
  integer        :: out_magnetization_step
  character(1)   :: yn_out_perflog
  character(6)   :: format_perflog ! 'stdout','text','csv'
  
!! &poisson
  integer        :: layout_multipole
  integer        :: num_multipole_xyz(3)
  integer        :: lmax_multipole
  real(8)        :: threshold_cg
  character(2)   :: method_poisson ! 'cg','ft'

!! &ewald
  integer        :: newald
  real(8)        :: aewald
  real(8)        :: cutoff_r
  real(8)        :: cutoff_r_buff
  real(8)        :: cutoff_g

!! &opt
  integer        :: nopt
  real(8)        :: max_step_len_adjust
  real(8)        :: convrg_opt_fmax

!! &md
  character(10)  :: ensemble
  character(20)  :: thermostat
  integer        :: step_velocity_scaling
  integer        :: step_update_ps
  real(8)        :: temperature0_ion_k
  character(1)   :: yn_set_ini_velocity
  character(256) :: file_ini_velocity
  real(8)        :: thermostat_tau
  character(1)   :: yn_stop_system_momt

!! &jellium
  character(1)   :: yn_jm
  character(1)   :: yn_charge_neutral_jm
  character(1)   :: yn_output_dns_jm
  character(256) :: shape_file_jm
  integer        :: num_jm
  real(8)        :: rs_bohr_jm(200)
  integer        :: sphere_nion_jm(200)
  real(8)        :: sphere_loc_jm(200,3)

!! &atomic_coor
!! &atomic_red_coor
integer,allocatable :: kion(:)
real(8),allocatable :: Rion(:,:)
real(8),allocatable :: Rion_red(:,:)
character(1),allocatable :: flag_opt_atom(:)
character(256),allocatable :: atom_name(:)

!! &code
  character(1) :: yn_want_stencil_hand_vectorization
  character(1) :: yn_want_communication_overlapping
  character(10) :: stencil_openmp_mode  ! 'auto', 'orbital', 'rgrid'
  character(10) :: current_openmp_mode  ! 'auto', 'orbital', 'rgrid'
  character(10) :: force_openmp_mode    ! 'auto', 'orbital', 'rgrid'

!! &band
  character(3) :: lattice
  integer :: nref_band
  real(8) :: tol_esp_diff
  integer :: num_of_segments
  integer :: ndiv_segment(max_num_of_segments)
  real(8) :: kpt(3,max_num_of_segments+1)
  character(1) :: kpt_label(max_num_of_segments+1)

end module salmon_global
