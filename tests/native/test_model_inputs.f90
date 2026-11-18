program test_model_inputs
  use salmon_global
  use salmonpy_model_inputs, only: s_model_inputs
  implicit none
  type(s_model_inputs) :: first,second,clone,empty
  integer :: status
  call empty%activate(status)
  if(status==0)error stop "uncaptured input accepted"
  call fill(1,2)
  call first%capture()
  clone=first
  call fill(2,3)
  call second%capture()
  call first%activate(status)
  if(status/=0)error stop "activate first failed"
  if(natom/=2.or.nelec/=1.or.xc/="a".or.dt/=1d0.or.iflag_atom_coor/=1)error stop "first scalar ownership"
  if(any(shape(Rion)/=[3,2]).or.any(Rion/=1d0).or.any(kion/=1))error stop "first array ownership"
  if(any(num_kgrid/=1).or.any(dl/=1d0).or.any(file_pseudo/="a"))error stop "first model input"
  Rion=9d0;kion=9;xc="changed";dt=9d0
  call second%activate(status)
  if(status/=0)error stop "activate second failed"
  if(natom/=3.or.nelec/=2.or.xc/="b".or.dt/=2d0)error stop "second scalar ownership"
  if(any(shape(Rion)/=[3,3]).or.any(Rion/=2d0).or.any(kion/=2))error stop "second array ownership"
  call clone%activate(status)
  if(status/=0.or.any(Rion/=1d0).or.natom/=2)error stop "deep assignment alias"
  deallocate(Rion,Rion_red,kion,flag_opt_atom,atom_name)
  call empty%capture()
  call fill(2,3)
  call empty%activate(status)
  if(status/=0.or.allocated(Rion).or.allocated(kion))error stop "unallocated input ownership"
  call first%activate(status)
  if(status/=0.or.any(Rion/=1d0))error stop "input lifetime independence"
  print *,"model input ownership PASS"
contains
  subroutine fill(value,atoms)
    integer,intent(in) :: value,atoms
    theory=merge('a','b',value==1)
    yn_md=merge('a','b',value==1)
    yn_opt=merge('a','b',value==1)
    yn_dc=merge('a','b',value==1)
    yn_conventional_from_dcdft=merge('a','b',value==1)
    sysname=merge('a','b',value==1)
    base_directory=merge('a','b',value==1)
    yn_restart=merge('a','b',value==1)
    directory_read_data=merge('a','b',value==1)
    yn_self_checkpoint=merge('a','b',value==1)
    checkpoint_interval=value
    yn_reset_step_restart=merge('a','b',value==1)
    read_gs_restart_data=merge('a','b',value==1)
    write_gs_restart_data=merge('a','b',value==1)
    time_shutdown=real(value,8)
    method_wf_distributor=merge('a','b',value==1)
    nblock_wf_distribute=value
    write_gs_wfn_k=merge('a','b',value==1)
    write_rt_wfn_k=merge('a','b',value==1)
    unit_system=merge('a','b',value==1)
    unit_time=merge('a','b',value==1)
    unit_length=merge('a','b',value==1)
    unit_energy=merge('a','b',value==1)
    unit_charge=merge('a','b',value==1)
    nproc_k=value
    nproc_ob=value
    nproc_rgrid=value
    yn_ffte=merge('a','b',value==1)
    yn_fftw=merge('a','b',value==1)
    yn_scalapack=merge('a','b',value==1)
    yn_gramschmidt_blas=merge('a','b',value==1)
    yn_eigenexa=merge('a','b',value==1)
    yn_diagonalization_red_mem=merge('a','b',value==1)
    yn_nccl_reduction=merge('a','b',value==1)
    process_allocation=merge('a','b',value==1)
    iperiodic=value
    yn_periodic=merge('a','b',value==1)
    spin=merge('a','b',value==1)
    al=real(value,8)
    al_vec1=real(value,8)
    al_vec2=real(value,8)
    al_vec3=real(value,8)
    nstate=value
    nelec=value
    nelec_spin=value
    temperature=real(value,8)
    temperature_k=real(value,8)
    nelem=value
    natom=value
    file_atom_coor=merge('a','b',value==1)
    file_atom_red_coor=merge('a','b',value==1)
    yn_spinorbit=merge('a','b',value==1)
    yn_symmetry=merge('a','b',value==1)
    absorbing_boundary=merge('a','b',value==1)
    imagnary_potential_w0=real(value,8)
    imagnary_potential_dr=real(value,8)
    file_pseudo=merge('a','b',value==1)
    lmax_ps=value
    lloc_ps=value
    izatom=value
    yn_psmask=merge('a','b',value==1)
    yn_pseudo_atomic_orbital=merge('a','b',value==1)
    alpha_mask=real(value,8)
    gamma_mask=real(value,8)
    eta_mask=real(value,8)
    xc=merge('a','b',value==1)
    xname=merge('a','b',value==1)
    hse_omega=real(value,8)
    cname=merge('a','b',value==1)
    alibx=merge('a','b',value==1)
    alibc=merge('a','b',value==1)
    alibxc=merge('a','b',value==1)
    cval=real(value,8)
    dl=real(value,8)
    num_rgrid=value
    num_kgrid=value
    file_kw=merge('a','b',value==1)
    dk_shift=real(value,8)
    yn_gamma_centered=merge('a','b',value==1)
    nt=value
    dt=real(value,8)
    gram_schmidt_interval=value
    n_hamil=value
    propagator=merge('a','b',value==1)
    yn_fix_func=merge('a','b',value==1)
    yn_predictor_corrector=merge('a','b',value==1)
    method_init_wf=merge('a','b',value==1)
    iseed_number_change=value
    method_min=merge('a','b',value==1)
    ncg=value
    ncg_init=value
    method_mixing=merge('a','b',value==1)
    mixrate=real(value,8)
    nmemory_mb=value
    alpha_mb=real(value,8)
    nmemory_p=value
    beta_p=real(value,8)
    yn_auto_mixing=merge('a','b',value==1)
    update_mixing_ratio=real(value,8)
    nscf=value
    yn_subspace_diagonalization=merge('a','b',value==1)
    convergence=merge('a','b',value==1)
    threshold=real(value,8)
    nscf_init_redistribution=value
    nscf_init_no_diagonal=value
    nscf_init_mix_zero=value
    conv_gap_mix_zero=real(value,8)
    method_init_density=merge('a','b',value==1)
    magdir_atom=real(value,8)
    yn_preconditioning=merge('a','b',value==1)
    alpha_pre=real(value,8)
    trans_longi=merge('a','b',value==1)
    ae_shape1=merge('a','b',value==1)
    file_input1=merge('a','b',value==1)
    e_impulse=real(value,8)
    E_amplitude1=real(value,8)
    I_wcm2_1=real(value,8)
    tw1=real(value,8)
    omega1=real(value,8)
    epdir_re1=real(value,8)
    epdir_im1=real(value,8)
    phi_cep1=real(value,8)
    ae_shape2=merge('a','b',value==1)
    E_amplitude2=real(value,8)
    I_wcm2_2=real(value,8)
    tw2=real(value,8)
    omega2=real(value,8)
    epdir_re2=real(value,8)
    epdir_im2=real(value,8)
    phi_cep2=real(value,8)
    t1_t2=real(value,8)
    t1_start=real(value,8)
    num_dipole_source=value
    vec_dipole_source=real(value,8)
    cood_dipole_source=real(value,8)
    rad_dipole_source=real(value,8)
    projection_option=merge('a','b',value==1)
    threshold_projection=real(value,8)
    nenergy=value
    de=real(value,8)
    out_rt_energy_step=value
    yn_out_psi=merge('a','b',value==1)
    yn_out_dos=merge('a','b',value==1)
    yn_out_dos_set_fe_origin=merge('a','b',value==1)
    out_dos_start=real(value,8)
    out_dos_end=real(value,8)
    out_dos_nenergy=value
    out_dos_width=real(value,8)
    out_dos_function=merge('a','b',value==1)
    yn_out_pdos=merge('a','b',value==1)
    yn_out_dns=merge('a','b',value==1)
    yn_out_dns_rt=merge('a','b',value==1)
    yn_out_dns_ac_je=merge('a','b',value==1)
    yn_out_micro_je=merge('a','b',value==1)
    out_dns_rt_step=value
    out_dns_ac_je_step=value
    out_micro_je_step=value
    out_old_dns=merge('a','b',value==1)
    yn_out_dns_trans=merge('a','b',value==1)
    out_dns_trans_energy=real(value,8)
    yn_out_elf=merge('a','b',value==1)
    yn_out_elf_rt=merge('a','b',value==1)
    out_elf_rt_step=value
    yn_out_estatic_rt=merge('a','b',value==1)
    out_estatic_rt_step=value
    yn_out_rvf_rt=merge('a','b',value==1)
    out_rvf_rt_step=value
    yn_out_tm=merge('a','b',value==1)
    yn_out_tm_bin=merge('a','b',value==1)
    yn_out_gs_sgm_eps=merge('a','b',value==1)
    out_gs_sgm_eps_mu_nu=value
    out_gs_sgm_eps_width=real(value,8)
    out_projection_step=value
    out_ms_step=value
    format_voxel_data=merge('a','b',value==1)
    nsplit_voxel_data=value
    yn_lr_w0_correction=merge('a','b',value==1)
    yn_out_intraband_current=merge('a','b',value==1)
    yn_out_current_decomposed=merge('a','b',value==1)
    out_current_decomposed_step=value
    out_rt_spin_step=value
    yn_out_mag_decomposed_rt=merge('a','b',value==1)
    yn_out_mag_micro_rt=merge('a','b',value==1)
    yn_out_spin_current_decomposed=merge('a','b',value==1)
    yn_out_spin_current_micro=merge('a','b',value==1)
    yn_out_rt_energy_components=merge('a','b',value==1)
    yn_out_perflog=merge('a','b',value==1)
    format_perflog=merge('a','b',value==1)
    layout_multipole=value
    num_multipole_xyz=value
    lmax_multipole=value
    threshold_cg=real(value,8)
    method_poisson=merge('a','b',value==1)
    newald=value
    aewald=real(value,8)
    cutoff_r=real(value,8)
    cutoff_r_buff=real(value,8)
    cutoff_g=real(value,8)
    yn_jm=merge('a','b',value==1)
    yn_charge_neutral_jm=merge('a','b',value==1)
    yn_output_dns_jm=merge('a','b',value==1)
    shape_file_jm=merge('a','b',value==1)
    num_jm=value
    rs_bohr_jm=real(value,8)
    sphere_nion_jm=value
    sphere_loc_jm=real(value,8)
    yn_want_stencil_hand_vectorization=merge('a','b',value==1)
    yn_want_communication_overlapping=merge('a','b',value==1)
    stencil_openmp_mode=merge('a','b',value==1)
    current_openmp_mode=merge('a','b',value==1)
    force_openmp_mode=merge('a','b',value==1)
    dm_unfold_option=merge('a','b',value==1)
    num_lkgrid=value
    num_skgrid=value
    no_pr=value
    out_dm_unfold_step=value
    yn_out_mom_distr_gs=merge('a','b',value==1)
    yn_out_mom_distr_rt=merge('a','b',value==1)
    out_mom_distr_rt_step=value
    nq_mom=value
    dq_mom=real(value,8)
    iflag_atom_coor=value
    natom=atoms
    if(allocated(Rion))deallocate(Rion,Rion_red,kion,flag_opt_atom,atom_name)
    allocate(Rion(3,atoms),Rion_red(3,atoms),kion(atoms),flag_opt_atom(atoms),atom_name(atoms))
    Rion=real(value,8);Rion_red=real(value,8);kion=value;flag_opt_atom="n";atom_name="Si"
  end subroutine
end program
