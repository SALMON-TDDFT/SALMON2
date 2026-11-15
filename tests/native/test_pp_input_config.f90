! Loader ownership probe; no assertion of full TDDFT/pseudopotential-grid independence.
program probe
 use structures,only:s_pp_input_config,s_pp_info,s_dft_system
 use salmon_pp,only:capture_pp_input,read_pslfile
 use salmon_global
 use parallelization,only:nproc_id_global,nproc_group_global
 implicit none
 type(s_pp_input_config)::a,b
 type(s_pp_info)::legacy,explicit,alternate,repeated
 type(s_dft_system)::sys
 nproc_id_global=0;nproc_group_global=0
 iperiodic=3;nelem=1;izatom(1)=14;lloc_ps(1)=2;lmax_ps(1)=-1
 file_pseudo(1)='Si_rps.dat';yn_psmask='n';yn_spinorbit='n';yn_pseudo_atomic_orbital='n'
 alpha_mask=0.8d0;eta_mask=15d0;gamma_mask=1.8d0;method_init_density='wf';quiet=.true.;base_directory='./'
 allocate(sys%mass(1));sys%hgs=0.5d0
 call capture_pp_input(a)
 call read_pslfile(sys,legacy)
 ! Poison legacy global input after capture. Explicit path must never consult it.
 nelem=2;iperiodic=0;izatom=-1;lloc_ps=-1;lmax_ps=8;file_pseudo='missing_file'
 ps_format='POISON';yn_psmask='?';yn_spinorbit='y';yn_pseudo_atomic_orbital='y'
 alpha_mask=-1d0;eta_mask=-1d0;gamma_mask=-1d0;method_init_density='pp';base_directory='missing/'
 call read_pslfile(sys,explicit,a)
 call same(legacy,explicit)
 if(sys%mass(1)/=28d0.or.explicit%nrmax/=3000.or.any(ps_format/='POISON'))stop 1
 b=a;b%lloc_ps=1;b%izatom=6
 call read_pslfile(sys,alternate,b)
 if(alternate%lref(1)/=1.or.alternate%rmass(1)/=12d0)stop 2
 if(all(alternate%vloctbl==explicit%vloctbl))stop 3
 b%file_pseudo='caller_mutated';b%lloc_ps=0
 if(alternate%input_config%file_pseudo(1)/='Si_rps.dat'.or.alternate%input_config%lloc_ps(1)/=1)stop 4
 call read_pslfile(sys,repeated,a)
 call same(explicit,repeated)
 if(explicit%input_config%lloc_ps(1)/=2.or.explicit%rmass(1)/=28d0)stop 5
 print *,'Pseudopotential input ownership PASS'
contains
 subroutine same(x,y)
 type(s_pp_info),intent(in)::x,y
 if(any(x%mr/=y%mr).or.any(x%mlps/=y%mlps).or.any(x%nproj/=y%nproj))stop 10
 if(any(x%rmass/=y%rmass).or.any(x%lref/=y%lref).or.any(x%ps_format/=y%ps_format))stop 11
 if(any(x%rad/=y%rad).or.any(x%radnl/=y%radnl).or.any(x%vloctbl/=y%vloctbl))stop 12
 if(any(x%dvloctbl/=y%dvloctbl).or.any(x%udvtbl/=y%udvtbl).or.any(x%dudvtbl/=y%dudvtbl))stop 13
 if(any(x%anorm/=y%anorm).or.any(x%inorm/=y%inorm).or.any(x%rps/=y%rps))stop 14
 if(any(x%rho_nlcc_tbl/=y%rho_nlcc_tbl).or.any(x%tau_nlcc_tbl/=y%tau_nlcc_tbl))stop 15
 if(any(x%udvtbl_so/=y%udvtbl_so).or.any(x%anorm_so/=y%anorm_so))stop 16
 if(any(x%upptbl_ao/=y%upptbl_ao).or.any(x%upptbl_ao_l/=y%upptbl_ao_l))stop 17
 end subroutine
end program
