! Synthetic occupied spinors: no ground-state calculation or collected GS data.
program test_core_current
#ifdef TEST_MPI
  use mpi
#endif
  use structures
  use nc_mgga23
  use hamiltonian, only: hpsi
  use taylor_sub, only: taylor
  use salmon_global, only: yn_jm,yn_periodic,yn_want_communication_overlapping,absorbing_boundary
  use code_optimization, only: stencil_is_parallelized_by_omp
  use timer, only: timer_initialize
  use inputoutput, only: n_hamil,t_unit_time,t_unit_current
  use salmon_xc, only: init_xc,exchange_correlation
  use noncollinear_module, only: calc_dm_noncollinear
  use write_sub, only: write_rt_spin
  use parallelization, only: nproc_id_global
  use salmon_global, only: yn_spinorbit,base_directory,sysname,dt, &
    yn_out_mag_decomposed_rt,yn_out_mag_micro_rt,yn_out_spin_current_decomposed,yn_out_spin_current_micro
  use density_matrix, only: calc_current
  implicit none
  integer,parameter :: ng=8
  type(s_dft_system) :: sys
  type(s_rgrid) :: mg
  type(s_parallel_info) :: info
  type(s_stencil) :: st
  type(s_orbital) :: psi,hp,base,delta,phi,hphi,hfull,hbare,kinetic,coarse,fine
  type(s_pp_grid) :: ppg
  type(s_pp_info) :: pp
  type(s_xc_functional) :: func
  type(s_ofile) :: ofl
  type(s_scalar) :: vlocal(2),rho(2)
  type(s_rt) :: rt
  type(s_pp_nlcc) :: core
  type(s_sendrecv_grid) :: scalar_halo,orb_halo
  real(8) :: theta,e0,ep,em,eps,fd,analytic,band,err,erot,rt_error,spin_error,norm_error
  complex(8) :: lhs,rhs,u(2,2),q(2)
  integer :: x,y,z,i,io,rank,nrank,mpierr
  real(8) :: before(3,2,1),after(3,2,1),expected(3),reduced
  character(32) :: mode
  rank=0;nrank=1
#ifdef TEST_MPI
  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,mpierr)
  if(nrank/=2) stop 'MPI test requires exactly two ranks'
#endif
  mg%is=1;mg%ie=ng;mg%num=ng;mg%is_array=1;mg%ie_array=ng
  allocate(mg%idx(-3:ng+4),mg%idy(-3:ng+4),mg%idz(-3:ng+4))
  do i=-3,ng+4
    mg%idx(i)=modulo(i-1,ng)+1;mg%idy(i)=mg%idx(i);mg%idz(i)=mg%idx(i)
  end do
  info%im_s=1;info%im_e=1;info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=2
  info%if_divide_rspace=.false.;info%icomm_r=0;info%icomm_ko=0
#ifdef TEST_MPI
  info%io_s=rank+1;info%io_e=rank+1
  info%icomm_r=MPI_COMM_SELF;info%icomm_ko=MPI_COMM_WORLD
#endif
  info%icomm_rko=info%icomm_ko;nproc_id_global=rank
  sys%no=2;sys%nk=1;sys%ngrid=ng**3
  sys%nspin=2;sys%hvol=0.31d0
  allocate(sys%vec_k(3,1),sys%wtk(1),sys%rocc(2,1,2))
  sys%vec_k(:,1)=[0.17d0,-0.13d0,0.09d0];sys%wtk=1d0
  sys%rocc(:,1,1)=[0.8d0,0.57d0];sys%rocc(:,1,2)=sys%rocc(:,1,1)
  sys%rmatrix_B=reshape([1d0,0.12d0,0.04d0,0d0,1.1d0,-0.08d0,0.03d0,0d0,0.95d0],[3,3])
  st%if_orthogonal=.true.;st%coef_nab=0d0;st%coef_nab(1,:)=0.5d0
  st%coef_lap=0d0;st%coef_lap(1,:)=1d0;st%coef_lap0=3d0
  allocate(psi%zwf(ng,ng,ng,2,info%io_s:info%io_e,1,1), &
    delta%zwf(ng,ng,ng,2,info%io_s:info%io_e,1,1))
  allocate(core%rho_nlcc(ng,ng,ng))
  do z=1,ng
  do y=1,ng
  do x=1,ng
    theta=2d0*acos(-1d0)/ng
    core%rho_nlcc(x,y,z)=0d0
    do io=info%io_s,info%io_e
      psi%zwf(x,y,z,1,io,1,1)=cmplx(0.34d0+0.03d0*cos(theta*x*io),0.02d0*sin(theta*y),8)
      psi%zwf(x,y,z,2,io,1,1)=cmplx((-1d0)**io*0.19d0+0.02d0*sin(theta*y),0.03d0*cos(theta*z*io),8)
      delta%zwf(x,y,z,1,io,1,1)=cmplx(0.02d0*cos(theta*(x+z)),0.03d0*sin(theta*y*io),8)
      delta%zwf(x,y,z,2,io,1,1)=cmplx(0.03d0*sin(theta*(y+z)),0.02d0*cos(theta*x*io),8)
    end do
  end do
  end do
  end do
  call timer_initialize
  yn_jm='y';yn_spinorbit='y'
  sys%vec_Ac=0d0
  sys%nc_mgga23_payload%active=.false.
  call calc_current(sys,mg,st,info,orb_halo,psi,ppg,before)
  sys%nc_mgga23_payload%active=.true.
  sys%nc_mgga23_payload%core_deda=[0.17d0,-0.23d0,0.31d0]
  expected=sys%nc_mgga23_payload%core_deda/(sys%ngrid*sys%hvol)
  call calc_current(sys,mg,st,info,orb_halo,psi,ppg,after)
  if(maxval(abs(sum(after(:,:,1)-before(:,:,1),dim=2)-expected))>1d-13) &
    error stop 'FAIL missing/duplicate/sign core current correction'
  if(maxval(abs(after(:,2,1)-before(:,2,1)))>1d-13) error stop 'FAIL double spin core correction'
  print *, 'PASS core current added once with native positive dE/dAc convention'
end program
