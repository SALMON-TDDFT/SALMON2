! Scalar grid/atom ownership and geometry-update regression; synthetic small system.
program probe
 use structures
 use prep_pp_sub,only:cache_jxyz,init_ps,dealloc_init_ps,pp_grid_config_from_input
 use salmon_pp,only:read_pslfile,calc_nlcc
 use salmon_global
 use parallelization,only:setup_parallel,end_parallel
 use communication,only:COMM_GROUP_SELF
 use timer,only:timer_initialize
 use hamiltonian,only:update_kvector_nonlocalpt
 use plusU_global,only:PLUS_U_ON
 implicit none
 type(s_pp_grid)::a,b,c,cache_a,cache_b
 type(s_pp_grid_config)::cfg
 type(s_pp_input_config)::radial
 type(s_pp_info)::pp
 type(s_dft_system)::sa,sb,sc
 type(s_rgrid)::lg,mg
 type(s_parallel_info)::info
 type(s_reciprocal_grid)::fg
 type(s_poisson)::poisson
 type(s_scalar)::va,vb,vc
 type(s_pp_nlcc)::na,nb,nc
 real(8)::ra(3,2),rb(3,1)
 real(8),allocatable::saved_v(:,:,:),saved_rho(:,:,:)
 integer::i,j
 character(24)::mode
 call get_command_argument(1,mode)
 call setup_parallel
 call timer_initialize
 ra=0d0;rb=1d0;natom=9
 call cache_jxyz(cache_a,ra)
 if(size(cache_a%mps_old)/=2.or.size(cache_a%jxyz_old,3)/=2)stop 1
 natom=17
 call cache_jxyz(cache_b,rb)
 if(size(cache_b%mps_old)/=1.or.size(cache_b%rxyz_old,3)/=1)stop 2
 if(any(cache_a%rion_old/=ra))stop 3
 ! Two radial species, with different local channels, from the licensed Si fixture.
 radial%iperiodic=0;radial%nelem=2
 radial%izatom=[14,6];radial%lloc_ps=[2,1];radial%lmax_ps=[-1,-1]
 allocate(radial%file_pseudo(2));radial%file_pseudo='Si_rps.dat'
 radial%quiet=.true.;radial%base_directory='./';radial%method_init_density='wf'
 allocate(sa%mass(2));sa%hgs=.5d0
 call read_pslfile(sa,pp,radial)
 sa%nion=1;sa%kion=[1];allocate(sa%rion(3,1));sa%rion=0d0
 sa%primitive_a=0d0;sa%rmatrix_a=0d0;sa%hvol=product(sa%hgs)
 do i=1,3
  sa%primitive_a(i,i)=4d0;sa%rmatrix_a(i,i)=1d0
 enddo
 sb=sa;sb%nion=2;sb%kion=[2,1];deallocate(sb%rion);allocate(sb%rion(3,2))
 sb%rion(:,1)=[.25d0,0d0,0d0];sb%rion(:,2)=[-.25d0,.2d0,0d0]
 lg%is=-3;lg%ie=4;lg%num=8
 allocate(lg%coordinate(-3:4,3))
 do i=-3,4
  lg%coordinate(i,:)=(real(i,8)-.5d0)*sa%hgs
 enddo
 mg=lg
 allocate(va%f(-3:4,-3:4,-3:4),vb%f(-3:4,-3:4,-3:4),vc%f(-3:4,-3:4,-3:4))
 info%ia_s=1;info%ia_e=1;info%icomm_r=COMM_GROUP_SELF;info%icomm_ko=COMM_GROUP_SELF
 a%nps=0;b%nps=0;c%nps=0
 iperiodic=0;natom=1;kion=[1];nelem=2;quiet=.true.;yn_spinorbit='n'
 method_poisson='cg';yn_ffte='n';yn_fftw='n'
 call pp_grid_config_from_input(cfg)
 call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
 if(trim(mode)=='bad_grid')then
  sa%hgs=sa%hgs*2d0
  call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
  error stop 'Unexpected accepted grid'
 endif
 if(trim(mode)=='bad_layout')then
  pp%mlps=[2,2,2]
  call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
  error stop 'Unexpected accepted layout'
 endif
 if(trim(mode)=='bad_mapping')then
  sa%kion=[2]
  call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
  error stop 'Unexpected accepted mapping'
 endif
 if(trim(mode)=='bad_source')then
  pp%input_config%file_pseudo='different-file'
  call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
  error stop 'Unexpected accepted source'
 endif
 saved_v=va%f
 ! Explicit B ignores poisoned globals; A owns its original arrays/settings.
 iperiodic=3;natom=19;nelem=17;kion=9;yn_spinorbit='y';method_poisson='ft';yn_ffte='y';yn_fftw='y'
 info%ia_e=2
 call init_ps(lg,mg,sb,info,fg,poisson,pp,b,vb,cfg)
 if(size(b%mps)/=2.or.size(a%mps)/=1.or.any(va%f/=saved_v))stop 4
 if(any(b%atom_species/=[2,1]).or.maxval(abs(vb%f-va%f))<1d-6)stop 5
 ! Periodic k-phase dispatch must also ignore global SOC/+U flags.
 allocate(sa%vec_k(3,1),sb%vec_k(3,1));sa%vec_k=.2d0;sb%vec_k=.3d0
 sa%vec_Ac=.1d0;sb%vec_Ac=.05d0;PLUS_U_ON=.true.
 call update_kvector_nonlocalpt(1,1,sa,a)
 call update_kvector_nonlocalpt(1,1,sb,b)
 do j=1,a%nlma
  do i=1,a%mps(a%ia_tbl(j))
   if(abs(a%zekr_uv(i,j,1)-exp(cmplx(0d0,-dot_product(sa%vec_k(:,1)+sa%vec_Ac, &
     a%rxyz(:,i,a%ia_tbl(j))),8))*a%uv(i,j))>1d-13)stop 11
  enddo
 enddo
 cfg%method_poisson='ft'
 if(b%input_config%method_poisson/='cg')stop 6
 cfg=a%input_config
 ! Existing geometry-update sequence must equal a new state at the moved position.
 sa%rion(1,1)=.1d0;sc=sa;info%ia_e=1
 call dealloc_init_ps(a)
 call init_ps(lg,mg,sa,info,fg,poisson,pp,a,va)
 call init_ps(lg,mg,sc,info,fg,poisson,pp,c,vc,cfg)
 if(any(va%f/=vc%f).or.any(a%mps/=c%mps).or.any(a%jxyz/=c%jxyz))stop 7
 if(any(a%uv/=c%uv).or.any(a%zekr_uv/=c%zekr_uv))stop 8
 ! Synthetic nonzero NLCC tables isolate species/boundary ownership.
 pp%flag_nlcc=.true.
 do j=1,2
  do i=1,pp%nrmax
   pp%rho_nlcc_tbl(i,j)=real(j,8)*exp(-pp%rad(i,j)**2)
   pp%tau_nlcc_tbl(i,j)=.5d0*pp%rho_nlcc_tbl(i,j)
  enddo
 enddo
 call calc_nlcc(pp,sa,mg,na);saved_rho=na%rho_nlcc
 call calc_nlcc(pp,sb,mg,nb)
 call calc_nlcc(pp,sa,mg,nc)
 if(any(nc%rho_nlcc/=saved_rho).or.any(na%tau_nlcc/=nc%tau_nlcc))stop 9
 if(maxval(abs(nb%rho_nlcc-saved_rho))<1d-6.or.maxval(saved_rho)<=0d0)stop 10
 print *,'Scalar PP grid/config/mapping/update and NLCC ownership PASS'
 call end_parallel
end program
