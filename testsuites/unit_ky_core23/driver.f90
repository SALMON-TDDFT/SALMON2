! Read through the production MPI-broadcast path and construct frozen core grids.
program test_reader_grid
#ifdef TEST_MPI
  use mpi,only:MPI_COMM_SELF
#endif
  use structures
  use salmon_pp, only: init_pp,calc_nlcc
  use input_pp_sub, only: input_pp
  use nc_mgga23, only: nc_mgga23_refresh,nc_mgga23_add_action
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use salmon_global, only: nelem,izatom,file_pseudo,lloc_ps,lmax_ps,yn_psmask,yn_spinorbit, &
    method_init_density,quiet,base_directory,kion,iperiodic
  use parallelization, only: setup_parallel,end_parallel,nproc_id_global
  use math_constants, only: set_math_constants
  implicit none
  type(s_pp_info) :: pp
  type(s_pp_nlcc) :: ppn
  type(s_dft_system) :: sys
  type(s_rgrid) :: rg
  type(s_parallel_info) :: info
  type(s_orbital) :: psi,hpsi
  type(s_stencil) :: stencil
  type(s_sendrecv_grid) :: scalar_halo,orb_halo
  real(8) :: exc,action_norm,band
  character(256) :: filename,output,mode
  integer :: u,i,l,x,y,z
  call setup_parallel
  call set_math_constants
  call get_command_argument(1,filename)
  call get_command_argument(2,mode)
  nelem=1;izatom=34;file_pseudo(1)=filename;lloc_ps=2;lmax_ps=-1
  yn_psmask='n';yn_spinorbit='y';method_init_density='wf';quiet=.true.;base_directory='./'
  call init_pp(pp,3000,8,.false.)
  call input_pp(pp,0.4d0,0.4d0,0.4d0)
  sys%nion=1;allocate(kion(1));kion(1)=1;iperiodic=0;sys%hgs=0.4d0
  allocate(sys%Rion(3,1));sys%Rion=0d0
  sys%primitive_a=0d0;sys%rmatrix_a=0d0
  do i=1,3
    sys%primitive_a(i,i)=20d0;sys%rmatrix_a(i,i)=1d0
  enddo
  if(trim(mode)=='--periodic') then
    iperiodic=3;sys%primitive_a=0d0
    do i=1,3
      sys%primitive_a(i,i)=1d0
    enddo
  endif
  rg%is=1;rg%ie=[3,2,2];rg%num=rg%ie-rg%is+1
  call calc_nlcc(pp,sys,rg,ppn)
  exc=0d0;action_norm=0d0;band=0d0
  if(ppn%core23_complete) then
    ! Bounded synthetic reader-to-gKS coupling check, not a Ni material state.
    rg%is_array=rg%is;rg%ie_array=rg%ie
    allocate(rg%idx(-3:rg%ie(1)+4),rg%idy(-3:rg%ie(2)+4),rg%idz(-3:rg%ie(3)+4))
    do i=-3,rg%ie(1)+4
      rg%idx(i)=modulo(i-1,rg%ie(1))+1
    enddo
    do i=-3,rg%ie(2)+4
      rg%idy(i)=modulo(i-1,rg%ie(2))+1
    enddo
    do i=-3,rg%ie(3)+4
      rg%idz(i)=modulo(i-1,rg%ie(3))+1
    enddo
    info%im_s=1;info%im_e=1;info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=2
    info%if_divide_rspace=.false.;info%icomm_r=0;info%icomm_ko=0
#ifdef TEST_MPI
    info%icomm_r=MPI_COMM_SELF;info%icomm_ko=MPI_COMM_SELF
#endif
    sys%nspin=2;sys%hvol=product(sys%hgs);sys%vec_Ac=0d0;sys%rmatrix_B=0d0
    do i=1,3
      sys%rmatrix_B(i,i)=1d0
    enddo
    allocate(sys%vec_k(3,1),sys%wtk(1),sys%rocc(2,1,2))
    sys%vec_k(:,1)=[0.07d0,-0.03d0,0.01d0];sys%wtk=1d0;sys%rocc=0.7d0
    allocate(psi%zwf(3,2,2,2,2,1,1))
    do z=1,2
    do y=1,2
    do x=1,3
      psi%zwf(x,y,z,1,1,1,1)=cmplx(0.31d0+0.01d0*cos(real(x,8)),0.01d0*sin(real(y+z,8)),8)
      psi%zwf(x,y,z,2,1,1,1)=cmplx(0.17d0,0.01d0*cos(real(x+y,8)),8)
      psi%zwf(x,y,z,1,2,1,1)=cmplx(0.23d0,0.01d0*sin(real(y,8)),8)
      psi%zwf(x,y,z,2,2,1,1)=cmplx(-0.21d0,0.01d0*cos(real(z,8)),8)
    enddo
    enddo
    enddo
    hpsi=psi;hpsi%zwf=0d0
    stencil%if_orthogonal=.true.;stencil%coef_nab=0d0;stencil%coef_nab(1,:)=0.5d0/sys%hgs
    stencil%coef_lap=0d0;stencil%coef_lap(1,:)=1d0/sys%hgs**2;stencil%coef_lap0=sum(stencil%coef_lap)
    call nc_mgga23_refresh(sys,rg,info,psi,stencil,scalar_halo,orb_halo,ppn,exc)
    call nc_mgga23_add_action(sys,rg,info,psi,hpsi,stencil,orb_halo)
    action_norm=sum(abs(hpsi%zwf));band=sys%nc_mgga23_payload%band_expectation
    if(.not.ieee_is_finite(exc).or..not.ieee_is_finite(action_norm).or..not.ieee_is_finite(band)) &
      error stop 'loaded core produced nonfinite gKS operator'
    if(action_norm<=0d0) error stop 'loaded core gKS action unexpectedly zero'
  endif
  write(output,'(a,i0,a)')'reader.rank',nproc_id_global,'.data'
  open(newunit=u,file=trim(output),status='replace')
  write(u,*)pp%core23_present(1),ppn%core23_complete,pp%mr(1),pp%nrps(1),pp%mlps(1)
  do i=1,pp%mr(1)+1
    write(u,'(100es27.16e3)')pp%rad(i,1),pp%rho_nlcc_tbl(i,1),pp%core23_tau_tbl(i,1), &
      pp%core23_current_tbl(i,1),(pp%udvtbl(i,l,1),pp%udvtbl_so(i,l,1),l=0,pp%mlps(1))
  enddo
  do z=rg%is(3),rg%ie(3)
  do y=rg%is(2),rg%ie(2)
  do x=rg%is(1),rg%ie(1)
    write(u,'(100es27.16e3)')ppn%rho_nlcc(x,y,z),ppn%core23_tau(x,y,z),ppn%core23_current(:,:,x,y,z)
  enddo
  enddo
  enddo
  write(u,'(3es27.16e3)')exc,action_norm,band
  close(u)
  call end_parallel
end program
