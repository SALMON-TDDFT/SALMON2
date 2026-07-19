program test_dg_wpw_lcfo_potential_map_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  use dg_wpw_lcfo_potential_map,only:s_dg_wpw_potential_map_state,s_dg_wpw_operator_candidate,&
    initialize_dg_wpw_potential_map_state,&
    run_dg_wpw_lcfo_potential_map,build_dg_wpw_raw_density,redistribute_dg_wpw_core_density
  use dg_wpw_lcfo_potential_map,only:build_dg_wpw_core_density
  use dg_wpw_lcfo_potential_map,only:s_dg_wpw_grid_route,prepare_dg_wpw_core_density_route,&
    return_dg_wpw_core_potential
  use dg_wpw_lcfo_potential_map,only:evaluate_and_run_dg_wpw_lcfo_potential_map
  implicit none
  type(s_dg_wpw_potential_map_state)::state
  integer::ierr,rank,info,stage
  real(8)::rho0(2)=[0.4d0,0.6d0],v0(2)=[1d0,2d0],raw(2)=[0.5d0,0.7d0],vc(2)=[1.5d0,2.5d0]
  real(8)::saved_rho(2),saved_v(2),saved_mix(2),saved_res(2),saved_energy
  integer(int64)::saved_fp
  complex(8),allocatable::saved_ww(:),saved_wp(:),saved_pp(:)
  integer::saved_epoch,saved_iteration
  integer::grid_ids(2),bad_ids(2)
  complex(8)::windows(2,1),coefficients(1,1)
  integer::owned_ids(2)
  integer::block_owned_ids(2),destinations(2)
  real(8)::occupations(1),weights(2),rho_local(2),rho_owned(2),charge_local,charge_global
  real(8)::potential_owned(2),potential_core(2)
  type(s_dg_wpw_grid_route)::grid_route
  complex(8)::w_points(1,2),p_points(1,2),cw_support(1,1),cp_support(1,1)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  windows(:,1)=[cmplx(1d0,0d0,8),cmplx(0d0,1d0,8)]
  coefficients(1,1)=cmplx(1d0,0d0,8);occupations=1d0;weights=0.5d0
  call build_dg_wpw_raw_density(windows,coefficients,occupations,weights,rho_local,charge_local,info)
  if(info/=0.or.maxval(abs(rho_local-1d0))>1d-14.or.abs(charge_local-1d0)>1d-14)&
    call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  w_points=reshape([(1d0,0d0),(2d0,0d0)],[1,2]);p_points=reshape([(0.5d0,0d0),(1d0,0d0)],[1,2])
  cw_support=(1d0,0d0);cp_support=(2d0,0d0)
  call build_dg_wpw_core_density(w_points,p_points,cw_support,cp_support,occupations,weights,&
    rho_local,charge_local,info)
  if(info/=0.or.maxval(abs(rho_local-[4d0,16d0]))>1d-14.or.abs(charge_local-10d0)>1d-14)&
    call MPI_Abort(MPI_COMM_WORLD,28,ierr)
  rho_local=1d0;charge_local=1d0
  grid_ids=[2*rank+1,2*rank+2]
  owned_ids=[rank+1,rank+3]
  call redistribute_dg_wpw_core_density(MPI_COMM_WORLD,grid_ids,rho_local,weights,4,owned_ids,rho_owned,&
    charge_global,info)
  if(info/=0.or.maxval(abs(rho_owned-1d0))>1d-14.or.abs(charge_global-2d0)>1d-14)&
    call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  call prepare_dg_wpw_core_density_route(MPI_COMM_WORLD,grid_ids,rho_local,weights,4,owned_ids,&
    grid_route,rho_owned,charge_global,info)
  potential_owned=10d0*dble(owned_ids)
  call return_dg_wpw_core_potential(grid_route,owned_ids,potential_owned,potential_core,info)
  if(info/=0.or.maxval(abs(potential_core-10d0*dble(grid_ids)))>1d-14)call MPI_Abort(MPI_COMM_WORLD,29,ierr)
  block_owned_ids=[2*rank+1,2*rank+2];destinations=[(grid_ids(1)-1)/2,(grid_ids(2)-1)/2]
  call prepare_dg_wpw_core_density_route(MPI_COMM_WORLD,grid_ids,rho_local,weights,4,block_owned_ids,&
    grid_route,rho_owned,charge_global,info,destinations=destinations)
  potential_owned=10d0*dble(block_owned_ids)
  call return_dg_wpw_core_potential(grid_route,block_owned_ids,potential_owned,potential_core,info)
  if(info/=0.or.maxval(abs(potential_core-10d0*dble(grid_ids)))>1d-14)call MPI_Abort(MPI_COMM_WORLD,30,ierr)
  bad_ids=grid_ids;if(rank==1)bad_ids(1)=2
  call redistribute_dg_wpw_core_density(MPI_COMM_WORLD,bad_ids,rho_local,weights,4,owned_ids,rho_owned,&
    charge_global,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  call initialize_dg_wpw_potential_map_state(state,MPI_COMM_WORLD,7,rho0,v0,3d0,991_int64,[0.1d0,0.2d0],&
    [0.3d0,0.4d0],4,info)
  if(info/=0)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  do stage=1,3
    saved_rho=state%rho;saved_v=state%potential;saved_energy=state%energy
    saved_fp=state%operator_fingerprint;saved_mix=state%mixer_history;saved_res=state%residual_history
    saved_epoch=state%publication_epoch;saved_iteration=state%iteration
    call run_dg_wpw_lcfo_potential_map(state,7,raw,0.5d0,vc,4d0,1234_int64,merge(stage,0,rank==0),info)
    if(info==0.or.any(state%rho/=saved_rho).or.any(state%potential/=saved_v).or.state%energy/=saved_energy.or.&
      state%operator_fingerprint/=saved_fp.or.any(state%mixer_history/=saved_mix).or.&
      any(state%residual_history/=saved_res).or.state%publication_epoch/=saved_epoch.or.&
      state%iteration/=saved_iteration)call MPI_Abort(MPI_COMM_WORLD,10+stage,ierr)
  enddo
  call run_dg_wpw_lcfo_potential_map(state,6,raw,0.5d0,vc,4d0,1234_int64,0,info)
  if(info==0.or.state%publication_epoch/=7)call MPI_Abort(MPI_COMM_WORLD,20,ierr)
  call run_dg_wpw_lcfo_potential_map(state,7,raw,0.5d0,vc,4d0,1234_int64,0,info,&
    potential_epoch=6,energy_epoch=7,operator_epoch=7)
  if(info==0.or.state%publication_epoch/=7.or.state%iteration/=4)call MPI_Abort(MPI_COMM_WORLD,23,ierr)
  call run_dg_wpw_lcfo_potential_map(state,7,raw,0.5d0,vc,4d0,1234_int64,0,info)
  if(info/=0.or.state%publication_epoch/=8.or.state%iteration/=5)call MPI_Abort(MPI_COMM_WORLD,21,ierr)
  if(maxval(abs(state%rho-[0.45d0,0.65d0]))>1d-14.or.any(state%potential/=vc).or.&
    state%energy/=4d0.or.state%operator_fingerprint/=1234_int64)call MPI_Abort(MPI_COMM_WORLD,22,ierr)
  saved_rho=state%rho;saved_v=state%potential;saved_energy=state%energy;saved_fp=state%operator_fingerprint
  saved_epoch=state%publication_epoch;saved_iteration=state%iteration
  call evaluate_and_run_dg_wpw_lcfo_potential_map(state,8,[0.55d0,0.75d0],0.5d0,&
    evaluate_candidate,rebuild_candidate,info)
  if(info/=0.or.state%publication_epoch/=9.or.state%iteration/=6.or.&
    maxval(abs(state%rho-[0.5d0,0.7d0]))>1d-14.or.maxval(abs(state%potential-[1.5d0,2.1d0]))>1d-14.or.&
    abs(state%energy-2.22d0)>1d-14.or.state%operator_fingerprint/=3600_int64.or.&
    any(state%operator_ww/=[(1d0,8d0),(2d0,8d0)]).or.any(state%operator_wp/=[(3d0,8d0)]).or.&
    any(state%operator_pp/=[(4d0,8d0)]))call MPI_Abort(MPI_COMM_WORLD,24,ierr)
  saved_rho=state%rho;saved_v=state%potential;saved_energy=state%energy;saved_fp=state%operator_fingerprint
  saved_ww=state%operator_ww;saved_wp=state%operator_wp;saved_pp=state%operator_pp
  saved_epoch=state%publication_epoch;saved_iteration=state%iteration
  call evaluate_and_run_dg_wpw_lcfo_potential_map(state,9,[0.6d0,0.8d0],0.5d0,&
    evaluate_failure,rebuild_candidate,info)
  if(info==0.or.any(state%rho/=saved_rho).or.any(state%potential/=saved_v).or.&
    state%energy/=saved_energy.or.state%operator_fingerprint/=saved_fp.or.&
    any(state%operator_ww/=saved_ww).or.any(state%operator_wp/=saved_wp).or.any(state%operator_pp/=saved_pp).or.&
    state%publication_epoch/=saved_epoch.or.state%iteration/=saved_iteration)call MPI_Abort(MPI_COMM_WORLD,25,ierr)
  call evaluate_and_run_dg_wpw_lcfo_potential_map(state,9,[0.6d0,0.8d0],0.5d0,&
    evaluate_candidate,rebuild_failure,info)
  if(info==0.or.any(state%rho/=saved_rho).or.any(state%potential/=saved_v).or.&
    state%energy/=saved_energy.or.state%operator_fingerprint/=saved_fp.or.&
    any(state%operator_ww/=saved_ww).or.any(state%operator_wp/=saved_wp).or.any(state%operator_pp/=saved_pp).or.&
    state%publication_epoch/=saved_epoch.or.state%iteration/=saved_iteration)call MPI_Abort(MPI_COMM_WORLD,26,ierr)
  call evaluate_and_run_dg_wpw_lcfo_potential_map(state,9,[0.6d0,0.8d0],0.5d0,&
    evaluate_candidate,rebuild_missing_blocks,info)
  if(info==0.or.any(state%operator_ww/=saved_ww).or.any(state%operator_wp/=saved_wp).or.&
    any(state%operator_pp/=saved_pp).or.state%publication_epoch/=saved_epoch)call MPI_Abort(MPI_COMM_WORLD,27,ierr)
  if(rank==0)print '(a)','PASS transactional LCFO potential-map publication and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine evaluate_candidate(rho,epoch,potential,energy,info)
    real(8),intent(in)::rho(:);integer,intent(in)::epoch
    real(8),intent(out)::potential(:),energy;integer,intent(out)::info
    potential=3d0*rho;energy=sum(potential*rho);info=0
  end subroutine
  subroutine evaluate_failure(rho,epoch,potential,energy,info)
    real(8),intent(in)::rho(:);integer,intent(in)::epoch
    real(8),intent(out)::potential(:),energy;integer,intent(out)::info
    potential=0d0;energy=0d0;info=merge(1,0,rank==0)
  end subroutine
  subroutine rebuild_candidate(potential,epoch,candidate,info)
    real(8),intent(in)::potential(:);integer,intent(in)::epoch
    type(s_dg_wpw_operator_candidate),intent(out)::candidate;integer,intent(out)::info
    candidate%fingerprint=int(nint(1000d0*sum(potential)),int64)
    candidate%ww=[cmplx(1d0,dble(epoch),8),cmplx(2d0,dble(epoch),8)]
    candidate%wp=[cmplx(3d0,dble(epoch),8)];candidate%pp=[cmplx(4d0,dble(epoch),8)]
    candidate%valid=.true.;info=merge(0,1,epoch==8)
  end subroutine
  subroutine rebuild_failure(potential,epoch,candidate,info)
    real(8),intent(in)::potential(:);integer,intent(in)::epoch
    type(s_dg_wpw_operator_candidate),intent(out)::candidate;integer,intent(out)::info
    candidate%fingerprint=9999_int64;candidate%ww=[cmplx(99d0,0d0,8)]
    candidate%wp=[cmplx(99d0,0d0,8)];candidate%pp=[cmplx(99d0,0d0,8)]
    candidate%valid=.true.;info=merge(1,0,rank==1)
  end subroutine
  subroutine rebuild_missing_blocks(potential,epoch,candidate,info)
    real(8),intent(in)::potential(:);integer,intent(in)::epoch
    type(s_dg_wpw_operator_candidate),intent(out)::candidate;integer,intent(out)::info
    candidate%fingerprint=7777_int64;candidate%valid=.true.;info=0
  end subroutine
end program
