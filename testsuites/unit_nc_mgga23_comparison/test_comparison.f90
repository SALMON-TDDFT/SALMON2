! Synthetic, deliberately non-ground-state trajectories; no material or GS data.
program test_comparison
  use structures
  use taylor_sub, only: taylor
  use salmon_xc, only: init_xc,exchange_correlation
  use noncollinear_module, only: calc_dm_noncollinear,rot_dm_noncollinear, &
    snapshot_xc_noncollinear,average_xc_noncollinear,calc_magnetization
  use write_sub, only: write_rt_spin
  use salmon_global, only: yn_jm,yn_periodic,yn_want_communication_overlapping,absorbing_boundary, &
    yn_spinorbit,base_directory,sysname,dt,xc,yn_md,yn_opt,yn_dc, &
    yn_out_mag_decomposed_rt,yn_out_mag_micro_rt,yn_out_spin_current_decomposed,yn_out_spin_current_micro
  use inputoutput, only: n_hamil,t_unit_time,t_unit_current
  use code_optimization, only: stencil_is_parallelized_by_omp
  use timer, only: timer_initialize
  use parallelization, only: nproc_id_global
  implicit none
  integer,parameter::ng=8,steps=20
  type(s_dft_system)::sys
  type(s_rgrid)::mg
  type(s_parallel_info)::info
  type(s_stencil)::st
  type(s_orbital)::psi,delta,base,pred,next,state
  type(s_pp_grid)::ppg
  type(s_pp_info)::pp
  type(s_xc_functional)::func
  type(s_ofile)::ofl
  type(s_scalar)::vlocal(2),rho(2)
  type(s_rt)::rt
  type(s_pp_nlcc)::core
  type(s_sendrecv_grid)::scalar_halo,orb_halo
  type(s_nc_mgga23_payload)::initial_payload
  complex(8),allocatable::initial_matrix(:,:,:,:,:)
  real(8)::theta,e,spin0(3),spin(3),max_spin_error,max_norm_error,norm0,norm_current,time
  integer::x,y,z,i,io,rank,run,it,unit
  rank=0
  mg%is=1;mg%ie=ng;mg%num=ng;mg%is_array=1;mg%ie_array=ng
  allocate(mg%idx(-3:ng+4),mg%idy(-3:ng+4),mg%idz(-3:ng+4))
  do i=-3,ng+4
    mg%idx(i)=modulo(i-1,ng)+1;mg%idy(i)=mg%idx(i);mg%idz(i)=mg%idx(i)
  end do
  info%im_s=1;info%im_e=1;info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=2
  info%if_divide_rspace=.false.;info%icomm_r=0;info%icomm_ko=0
  info%icomm_rko=info%icomm_ko;nproc_id_global=rank
  sys%no=2;sys%nk=1;sys%ngrid=ng**3
  sys%nspin=2;sys%hvol=0.31d0
  allocate(sys%vec_k(3,1),sys%wtk(1),sys%rocc(2,1,2))
  sys%vec_k(:,1)=0d0;sys%wtk=1d0
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
  info%numo=info%io_e-info%io_s+1;info%numk=1;info%numm=1
  yn_jm='n';yn_periodic='y';yn_want_communication_overlapping='n';absorbing_boundary='n'
  stencil_is_parallelized_by_omp=.false.;sys%vec_Ac=0d0
  call timer_initialize
  do i=1,2
    allocate(vlocal(i)%f(ng,ng,ng),rho(i)%f(ng,ng,ng));vlocal(i)%f=0d0;rho(i)%f=0d0
  end do

  yn_spinorbit='y';yn_md='n';yn_opt='n';yn_dc='n';pp%flag_nlcc=.false.
  base_directory='./';dt=.002d0
  yn_out_mag_decomposed_rt='n';yn_out_mag_micro_rt='n'
  yn_out_spin_current_decomposed='n';yn_out_spin_current_micro='n'
  t_unit_time%conv=1d0;t_unit_time%name='au'
  t_unit_current%conv=1d0;t_unit_current%name='au'
  allocate(ppg%ia_tbl_so(0));base=psi
  open(newunit=unit,file='validation.txt',status='replace')
  write(unit,'(a)') 'SYNTHETIC non-ground-state spinors; no Ni/material, no GS, no ionic SOC.'
  write(unit,'(a)') 'Periodic 8^3 points, 2 occupied spinors, same initial state for both functionals.'
  write(unit,'(a)') 'Thus this is an explicit functional-quench/operator comparison, not own-GS response.'
  write(unit,'(a)') 'dt=0.002 au; 20 steps; Taylor order=6; midpoint predictor/corrector XC.'
  write(unit,'(a)') 'Scalar local field 0.02*sin(5*t)*cos(2*pi*x/8); no spin torque.'
  do run=1,2
    if(run==1) then
      xc='pz';sysname='LSDA'
    else
      xc='nc_mgga23';sysname='MGGA23'
    endif
    call init_xc(func,'noncollinear',0d0,xcname=trim(xc))
    state=base
    call refresh(state,0d0)
    call calc_magnetization(sys,mg,info,spin0)
    norm0=weighted_norm(state)
    max_spin_error=0d0;max_norm_error=0d0
    call write_rt_spin(-1,ofl,sys,mg,mg,info,st,ppg,state)
    call write_rt_spin(0,ofl,sys,mg,mg,info,st,ppg,state)
    do it=1,steps
      time=(it-.5d0)*dt
      call refresh(state,time)
      if(sys%nc_mgga23_payload%active) then
        initial_payload=sys%nc_mgga23_payload
      else
        call snapshot_xc_noncollinear(initial_matrix)
      endif
      call propagate(state,pred)
      call refresh(pred,time)
      if(sys%nc_mgga23_payload%active) then
        sys%nc_mgga23_payload%local=.5d0*(initial_payload%local+sys%nc_mgga23_payload%local)
        sys%nc_mgga23_payload%tau=.5d0*(initial_payload%tau+sys%nc_mgga23_payload%tau)
        sys%nc_mgga23_payload%current=.5d0*(initial_payload%current+sys%nc_mgga23_payload%current)
      else
        call average_xc_noncollinear(initial_matrix)
      endif
      call propagate(state,next)
      state=next
      call refresh(state,it*dt)
      call calc_magnetization(sys,mg,info,spin)
      norm_current=weighted_norm(state)
      max_spin_error=max(max_spin_error,maxval(abs(spin-spin0)))
      max_norm_error=max(max_norm_error,abs(norm_current-norm0))
      call write_rt_spin(it,ofl,sys,mg,mg,info,st,ppg,state)
    enddo
    close(ofl%fh_rt_spin)
    write(unit,'(a,2(1x,es24.16))') trim(sysname)//' max spin,norm drift:',max_spin_error,max_norm_error
    if(max_spin_error>1d-6.or.max_norm_error>1d-7) stop 'spin/norm conservation failure'
  enddo
  close(unit)
  print *, 'PASS synthetic LSDA and MGGA23 trajectories, 21 initial-inclusive samples each'
contains
  subroutine refresh(p,time)
    type(s_orbital),intent(inout)::p
    real(8),intent(in)::time
    integer::x,y,z,sp
    call calc_dm_noncollinear(p,sys,info,mg)
    call rot_dm_noncollinear(rho,sys,mg)
    call exchange_correlation(sys,func,mg,scalar_halo,orb_halo,rho,pp,core,info,p,st,vlocal,e)
    do sp=1,2
    do z=1,ng
    do y=1,ng
    do x=1,ng
      vlocal(sp)%f(x,y,z)=vlocal(sp)%f(x,y,z)+.02d0*sin(5d0*time)*cos(theta*x)
    enddo
    enddo
    enddo
    enddo
  end subroutine
  real(8) function weighted_norm(p) result(value)
    type(s_orbital),intent(in)::p
    integer::io
    value=0d0
    do io=1,2
      value=value+sys%rocc(io,1,1)*sys%hvol*sum(abs(p%zwf(:,:,:,:,io,1,1))**2)
    enddo
  end function
  subroutine propagate(p,out)
    type(s_orbital),intent(in)::p
    type(s_orbital),intent(out)::out
    type(s_orbital)::scratch,work
    complex(8)::coeff
    integer::j
    n_hamil=6
    if(.not.allocated(rt%zc)) allocate(rt%zc(n_hamil))
    coeff=(1d0,0d0)
    do j=1,n_hamil
      coeff=coeff*cmplx(0d0,-dt,8)/j;rt%zc(j)=coeff
    enddo
    scratch=p;work=p;out=p
    call taylor(mg,sys,info,st,orb_halo,scratch,out,work,ppg,vlocal,rt)
  end subroutine
end program
