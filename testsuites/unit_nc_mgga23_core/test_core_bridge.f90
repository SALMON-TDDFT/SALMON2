! Synthetic occupied spinors: no ground-state calculation or collected GS data.
program test_core_bridge
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
  real(8) :: deda0(3),ac0(3),eoff
  integer :: irun,c
  real(8) :: reduced
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
  call get_command_argument(1,mode)
  if(trim(mode)=='--nonzero-core') core%rho_nlcc=0.08d0
  if(trim(mode)=='--nonzero-core-tau') then
    allocate(core%tau_nlcc(ng,ng,ng));core%tau_nlcc=0.01d0
  end if
  info%numo=info%io_e-info%io_s+1;info%numk=1;info%numm=1
  yn_jm='y';yn_periodic='y';yn_want_communication_overlapping='n';absorbing_boundary='n'
  stencil_is_parallelized_by_omp=.false.;sys%vec_Ac=0d0
  call timer_initialize
  do i=1,2
    allocate(vlocal(i)%f(ng,ng,ng),rho(i)%f(ng,ng,ng));vlocal(i)%f=0d0;rho(i)%f=0d0
  end do
  base=psi
  if(trim(mode)=='--invalid-metadata') then
    core%rho_nlcc=.08d0
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,e0)
    stop 'unexpected accepted incomplete metadata'
  endif
  if(trim(mode)=='--missing-tau') then
    core%core23_complete=.true.;core%rho_nlcc=.08d0
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,e0)
    stop 'unexpected accepted missing V2 tau'
  endif
  if(trim(mode)=='--bad-layout') then
    core%core23_complete=.true.;core%rho_nlcc=.08d0
    allocate(core%core23_tau(ng-1,ng,ng),core%core23_current(4,3,ng,ng,ng))
    core%core23_tau=.3d0;core%core23_current=0d0
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,e0)
    stop 'unexpected accepted bad V2 layout'
  endif
  allocate(core%core23_tau(ng,ng,ng),core%core23_current(4,3,ng,ng,ng))
  core%core23_complete=.true.
  sys%vec_Ac=[.07d0,-.03d0,.02d0];ac0=sys%vec_Ac
  do irun=1,2
    psi=base
    core%core23_current=0d0
    do z=1,ng
    do y=1,ng
    do x=1,ng
      core%rho_nlcc(x,y,z)=.08d0*(1d0+.2d0*cos(theta*x))
      core%core23_tau(x,y,z)=.4d0*(1d0+.1d0*cos(theta*y))
      if(irun==2) then
        core%core23_current(2,2,x,y,z)=.002d0*sin(theta*z)
        core%core23_current(3,1,x,y,z)=-.002d0*sin(theta*z)
      endif
    enddo
    enddo
    enddo
    if(trim(mode)=='--bad-core') core%core23_tau(1,1,1)=-1d0
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,e0)
    deda0=sys%nc_mgga23_payload%core_deda
    hp=psi;hp%zwf=0d0
    call nc_mgga23_add_action(sys,mg,info,psi,hp,st,orb_halo)
    analytic=0d0;band=0d0
    do io=info%io_s,info%io_e
      analytic=analytic+2d0*real(sum(conjg(delta%zwf(:,:,:,:,io,1,1))*hp%zwf(:,:,:,:,io,1,1)),8) &
        *sys%rocc(io,1,1)*sys%hvol
      band=band+real(sum(conjg(psi%zwf(:,:,:,:,io,1,1))*hp%zwf(:,:,:,:,io,1,1)),8) &
        *sys%rocc(io,1,1)*sys%hvol
    enddo
#ifdef TEST_MPI
    call MPI_Allreduce(analytic,reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr);analytic=reduced
    call MPI_Allreduce(band,reduced,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr);band=reduced
#endif
    if(abs(band-sys%nc_mgga23_payload%band_expectation)>1d-9) stop 'frozen core included in band expectation'
    phi=delta;hphi=phi;hphi%zwf=0d0
    call nc_mgga23_add_action(sys,mg,info,phi,hphi,st,orb_halo)
    lhs=sum(conjg(phi%zwf)*hp%zwf);rhs=sum(conjg(hphi%zwf)*psi%zwf)
    if(abs(lhs-rhs)>1d-10) stop 'core bridge not Hermitian'
    eps=1d-5
    psi%zwf=base%zwf+eps*delta%zwf
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,ep)
    psi%zwf=base%zwf-eps*delta%zwf
    call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,em)
    fd=(ep-em)/(2*eps)
    if(abs(fd-analytic)>2d-7*max(1d0,abs(fd),abs(analytic))) stop 'core orbital derivative mismatch'
    psi=base
    eps=1d-6
    do c=1,3
      sys%vec_Ac=ac0;sys%vec_Ac(c)=ac0(c)+eps
      call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,ep)
      sys%vec_Ac(c)=ac0(c)-eps
      call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,em)
      if(abs((ep-em)/(2*eps)-deda0(c))>2d-7) stop 'core dE/dAc mismatch'
    enddo
    sys%vec_Ac=ac0
    if(irun==1) then
      ! Scalar reference core has no spin torque; uniform phase cannot change this.
      call evolve(base,.001d0,5,fine)
      spin_error=maxval(abs(magnetization(fine)-magnetization(base)))
      if(spin_error>1d-7) stop 'scalar-core spin conservation failed'
    endif
    print *, 'PASS full core bridge orbital FD, band, Hermiticity, Ac FD; model=',irun
  enddo
  core%rho_nlcc=0d0;core%core23_tau=0d0;core%core23_current=0d0
  psi=base;core%core23_complete=.false.
  call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,eoff)
  core%core23_complete=.true.
  call nc_mgga23_refresh(sys,mg,info,psi,st,scalar_halo,orb_halo,core,e0)
  if(abs(eoff-e0)>1d-13) stop 'zero-core regression failed'
  print *, 'PASS zero-core regression'
#ifdef TEST_MPI
  call MPI_Finalize(mpierr)
#endif
contains
  function magnetization(p) result(m)
    type(s_orbital),intent(in) :: p
    real(8) :: m(3),w,mglobal(3)
    complex(8) :: cross
    integer :: io
    m=0d0
    do io=info%io_s,info%io_e
      w=sys%rocc(io,1,1)*sys%hvol
      cross=sum(conjg(p%zwf(:,:,:,1,io,1,1))*p%zwf(:,:,:,2,io,1,1))
      m(1)=m(1)+2d0*w*real(cross,8)
      m(2)=m(2)+2d0*w*aimag(cross)
      m(3)=m(3)+w*sum(abs(p%zwf(:,:,:,1,io,1,1))**2-abs(p%zwf(:,:,:,2,io,1,1))**2)
    end do
#ifdef TEST_MPI
    call MPI_Allreduce(m,mglobal,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr);m=mglobal
#endif
  end function
  subroutine propagate(p,dt,out)
    type(s_orbital),intent(in) :: p
    real(8),intent(in) :: dt
    type(s_orbital),intent(out) :: out
    type(s_orbital) :: scratch,work
    complex(8) :: coeff
    integer :: j
    n_hamil=6
    if(.not.allocated(rt%zc)) allocate(rt%zc(n_hamil))
    coeff=(1d0,0d0)
    do j=1,n_hamil
      coeff=coeff*cmplx(0d0,-dt,8)/j;rt%zc(j)=coeff
    end do
    scratch=p;work=p;out=p
    call taylor(mg,sys,info,st,orb_halo,scratch,out,work,ppg,vlocal,rt)
  end subroutine
  subroutine evolve(p,dt,steps,out)
    type(s_orbital),intent(in) :: p
    real(8),intent(in) :: dt
    integer,intent(in) :: steps
    type(s_orbital),intent(out) :: out
    type(s_orbital) :: state,mid,next
    real(8) :: e
    integer :: it
    state=p
    do it=1,steps
      call nc_mgga23_refresh(sys,mg,info,state,st,scalar_halo,orb_halo,core,e)
      call propagate(state,0.5d0*dt,mid)
      call nc_mgga23_refresh(sys,mg,info,mid,st,scalar_halo,orb_halo,core,e)
      call propagate(state,dt,next)
      state=next
    end do
    out=state
  end subroutine
end program
