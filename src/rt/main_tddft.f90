!
!  Copyright 2019-2020 SALMON developers
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

!=======================================================================

#include "config.h"

subroutine main_tddft
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
use math_constants, only: pi
#ifdef USE_MPI
use mpi, only: MPI_Comm_rank,MPI_Bcast,MPI_INTEGER,MPI_Allreduce,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,&
  MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_SUCCESS
use mpi, only: MPI_MAX
#endif
use salmon_global
use structures
use parallelization, only: adjust_elapse_time, nproc_group_global
use communication, only: comm_is_root, comm_sync_all, comm_bcast, comm_summation
use salmon_xc, only: finalize_xc
use timer
use write_sub, only: write_response_0d,write_response_3d,write_pulse_0d,write_pulse_3d, &
  write_dg_polarization_data, write_dg_polarization_response_3d
use initialization_rt_sub
use checkpoint_restart_sub
use jellium, only: check_condition_jm
use rt_angular_momentum, only: write_local_angular_momentum_xy, flush_local_angular_momentum_xy
use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
use dg_overlapping_wannier_checkpoint, only: s_dg_overlapping_wannier_checkpoint, &
  read_dg_overlapping_wannier_checkpoint
use rt_dg_overlapping_wannier, only: s_dg_overlapping_wannier_rt_state, &
  initialize_dg_overlapping_wannier_rt,advance_dg_overlapping_wannier_rt,&
  write_dg_overlapping_wannier_rt_restart,read_dg_overlapping_wannier_rt_restart,&
  evaluate_dg_overlapping_wannier_observables,&
  write_dg_overlapping_wannier_rt_observable_sample
use em_field, only: calc_Ac_ext_t
use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state,initialize_rt_dg_hybrid_from_checkpoint
use rt_dg_hybrid_density_update,only:update_rt_dg_hybrid_density,reconstruct_rt_dg_hybrid_density
use rt_dg_hybrid_stationarity,only:s_rt_dg_hybrid_stationarity_reference,&
  s_rt_dg_hybrid_stationarity_receipt,initialize_rt_dg_hybrid_stationarity,&
  evaluate_rt_dg_hybrid_stationarity
use dg_hybrid_total_energy,only:evaluate_dg_hybrid_fixed_energy
use rt_dg_hybrid_length_gauge,only:propagate_rt_dg_hybrid_length_gauge
use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange
use dg_overlapping_wannier_construction,only:redistribute_dg_row_owned_real_field_to_requests
use hartree_sub,only:hartree
use salmon_xc,only:exchange_correlation_density
use hamiltonian,only:update_vlocal
use Total_Energy,only:calc_Total_Energy_periodic
use plusU_global,only:PLUS_U_ON
use nvtx
use parallelization, only: nproc_id_global
implicit none

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_parallel_info) :: info
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_ewald_ion_ion) :: ewald
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_scalar) :: Vpsl
type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox
type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_scalar
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_singlescale) :: singlescale

integer :: Mit, itt
logical :: is_checkpoint_iter, is_shutdown_time, is_checkpoint
type(s_rt_dg_hybrid_state) :: hybrid_state

if(yn_dg_overlapping_wannier_rt=='y')then
  call run_dg_overlapping_wannier_coefficient_rt()
  return
endif

!check condition for using jellium model
if(yn_jm=='y') call check_condition_jm

call timer_begin(LOG_TOTAL)

if(yn_rt_dg_hybrid_continuation=='y')then
  call initialization_rt_dg_hybrid( Mit, system, energy, ewald, rt, md, &
                          singlescale, stencil, fg, poisson, lg, mg, info, xc_func, ofl, &
                          srg, srg_scalar, rho, rho_jm, rho_s, &
                          V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl, pp, ppg, ppn )
  call run_dg_hybrid_continuation_rt()
  return
endif

call initialization_rt( Mit, system, energy, ewald, rt, md, &
                        singlescale,  &
                        stencil, fg, poisson,  &
                        lg, mg,   &
                        info,  &
                        xc_func, ofl,  &
                        srg, srg_scalar,  &
                        spsi_in, spsi_out, tpsi, rho, rho_jm, rho_s,  &
                        V_local, Vbox, Vh, Vh_stock1, Vh_stock2, Vxc, Vpsl,&
                        pp, ppg, ppn )

#ifdef __FUJITSU
call fapp_start('time_evol',1,0) ! performance profiling
#endif

call print_header()

if (yn_out_lcm_rt == 'y') then
  call write_local_chern_marker_xy(Mit, mg, system, info, spsi_in)
end if
if (yn_out_lz_rt == 'y') then
  if (.not. singlescale%flag_use) stop 'yn_out_lz_rt=y requires theory=single_scale_maxwell_tddft'
  call write_local_angular_momentum_xy(Mit, lg, mg, system, info, singlescale, spsi_in)
end if

if (iperiodic == 3) then
  call write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, 'full-initial-density')
end if

#ifdef USE_OPENACC
!$acc enter data copyin(rt, rt%zc)
!$acc enter data copyin(mg, mg%is, mg%ie)
!$acc enter data copyin(poisson)
!$acc enter data copyin(fg)
!$acc enter data copyin(lg)
!$acc enter data copyin(Vh, Vxc, Vpsl)
!$acc enter data copyin(ewald, pp, ppg)
#endif

call comm_sync_all
call timer_enable_sub
call timer_begin(LOG_RT_ITERATION)

! === Standard real-space RT time evolution ===
TE : do itt=Mit+1,nt
    call nvtxStartRange('main loop', itt)

  if(mod(itt,2)==1)then
    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
     & ,srg,srg_scalar,pp,ppg,ppn,spsi_in,spsi_out,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
	  else
	    call time_evolution_step(Mit,nt,itt,lg,mg,system,rt,info,stencil,xc_func &
	     & ,srg,srg_scalar,pp,ppg,ppn,spsi_out,spsi_in,tpsi,rho,rho_jm,rho_s,V_local,Vbox,Vh,Vh_stock1,Vh_stock2,Vxc &
	     & ,Vpsl,fg,energy,ewald,md,ofl,poisson,singlescale)
	  end if

      if (yn_out_lcm_rt == 'y') then
        if (mod(itt, out_lcm_rt_step) == 0) then
          if (mod(itt,2) == 1) then
            call write_local_chern_marker_xy(itt, mg, system, info, spsi_out)
          else
            call write_local_chern_marker_xy(itt, mg, system, info, spsi_in)
          end if
        end if
      end if
      if (yn_out_lz_rt == 'y') then
        if (mod(itt, out_lz_rt_step) == 0) then
          if (mod(itt,2) == 1) then
            call write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale, spsi_out)
          else
            call write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale, spsi_in)
          end if
        end if
      end if

	  is_checkpoint_iter = (checkpoint_interval >= 1) .and. (mod(itt,checkpoint_interval) == 0)
  is_shutdown_time   = (time_shutdown > 0d0) .and. (adjust_elapse_time(timer_now(LOG_TOTAL)) > time_shutdown)

  is_checkpoint = is_checkpoint_iter .or. is_shutdown_time
  call nvtxStartRange('comm_bcast', __LINE__)
  call comm_bcast(is_checkpoint,nproc_group_global)
  call nvtxEndRange

  call nvtxEndRange
  if(is_checkpoint) then
    if (is_shutdown_time .and. comm_is_root(info%id_rko)) then
      print *, 'shutdown the calculation, iter =', itt
    end if

    call timer_begin(LOG_CHECKPOINT_SYNC)
    call timer_begin(LOG_CHECKPOINT_SELF)
    if (mod(itt,2)==1) then
      call checkpoint_rt(lg,mg,system,info,spsi_out,itt,rt,Vh_stock1,Vh_stock2,singlescale)
    else
      call checkpoint_rt(lg,mg,system,info,spsi_in, itt,rt,Vh_stock1,Vh_stock2,singlescale)
    endif
    call timer_end(LOG_CHECKPOINT_SELF)
    call comm_sync_all
    call timer_end(LOG_CHECKPOINT_SYNC)

    if (is_shutdown_time) then
      exit TE
    end if
  endif

end do TE

if (yn_out_lz_rt == 'y') then
  call flush_local_angular_momentum_xy(system)
end if

call timer_end(LOG_RT_ITERATION)
call timer_disable_sub

#ifdef __FUJITSU
call fapp_stop('time_evol',1,0) ! performance profiling
#endif

close(030) ! laser


!--------------------------------- end of time-evolution

!------------ Writing part -----------


call timer_begin(LOG_WRITE_RT_RESULTS)

!
select case(iperiodic)
case(0)
  if(theory=="tddft_response")then
    call write_response_0d(ofl,rt)
  else
    call write_pulse_0d(ofl,rt)
  end if
case(3)
  if(theory=="tddft_response")then
    call write_response_3d(ofl,rt)
    if (yn_dg_length_gauge == 'y') call write_dg_polarization_response_3d(ofl)
  else
    call write_pulse_3d(ofl,rt)
  end if
end select

if(comm_is_root(nproc_id_global))then
  close(ofl%fh_rt)  ! Close _rt.data file
end if

call timer_end(LOG_WRITE_RT_RESULTS)
call timer_end(LOG_TOTAL)

if(write_rt_wfn_k=='y')then
  call checkpoint_rt(lg,mg,system,info,spsi_out,Mit,rt,Vh_stock1,Vh_stock2,singlescale,ofl%dir_out_restart)
end if

call finalize_xc(xc_func)

contains

subroutine run_dg_hybrid_continuation_rt()
  type(s_rt_dg_sparse_exchange)::metric_exchange,operator_exchange
  type(s_rt_dg_hybrid_stationarity_reference)::stationarity_reference
  type(s_rt_dg_hybrid_stationarity_receipt)::stationarity_receipt
  complex(8),allocatable::next(:),initial_hamiltonian(:)
  real(8),allocatable::vector_potential_samples(:,:)
  real(8)::electric_field(3),previous_polarization(3),periods(3),metric_norm,orbital_energy,polarization(3),&
    local_defect,global_defect,local_scale,global_scale,current_total_energy,current_electron_count,&
    current_hamiltonian_residual,stationarity_tolerances(5)
  integer(8)::workspace,fingerprint
  integer::step,orbital,iterations,ierr,update_count,local_bad,global_bad
  logical::ok,stationarity_enabled
  character(256)::message
  call initialize_rt_dg_hybrid_from_checkpoint(nproc_group_global,'./hybrid_dg_ground_state.chk',theory,&
    iperiodic==3,system%nspin,yn_spinorbit=='y',PLUS_U_ON,yn_hse=='y',yn_fix_func=='y',yn_jm=='y',&
    xc_func%xctype,hybrid_state,ok,message)
  if(.not.ok)then;write(0,'(a)')trim(message);error stop 'hybrid DG RT initialization failed';endif
  allocate(initial_hamiltonian,source=hybrid_state%operators%hamiltonian_values)
  update_count=0
  call update_rt_dg_hybrid_density(nproc_group_global,hybrid_state,hybrid_state%density,project_salmon_local_rows,ok,message)
  update_count=update_count+1
  local_bad=merge(0,1,ok);local_defect=huge(1d0);local_scale=1d0
  global_bad=1;global_defect=huge(1d0);global_scale=1d0
  if(ok)then
    local_defect=maxval(abs(hybrid_state%operators%hamiltonian_values-initial_hamiltonian))
    local_scale=max(1d0,maxval(abs(initial_hamiltonian)))
  endif
  call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,nproc_group_global,ierr)
  if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
    nproc_group_global,ierr)
  if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
    nproc_group_global,ierr)
  if(ierr/=MPI_SUCCESS.or.global_bad/=0.or..not.ieee_is_finite(global_defect).or.&
      .not.ieee_is_finite(global_scale).or..not.(global_defect<=1d-10*global_scale))&
    error stop 'hybrid DG RT initial Hamiltonian reconstruction failed'
  call evaluate_hybrid_rt_physical_invariants(current_total_energy,current_electron_count,&
    current_hamiltonian_residual,ok,message)
  if(.not.ok)then;write(0,'(a)')trim(message);error stop 'hybrid DG RT initial physical invariants failed';endif
  if(.not.allocated(hybrid_state%energy_receipt))error stop 'hybrid DG RT physical energy receipt is absent'
  if(nproc_id_global==0)write(*,'(a,i0,2(a,es16.8))')'[HYBRID-RT-HANDOFF] payload_fingerprint=',&
    hybrid_state%payload_fingerprint,' operator_symmetry=',hybrid_state%startup_operator_covariance,&
    ' projector_symmetry=',hybrid_state%startup_projector_covariance
  stationarity_enabled=size(hybrid_state%energy_receipt)==7.and.any(hybrid_state%energy_receipt/=0d0)
  if(stationarity_enabled)then
    if(abs(current_total_energy-hybrid_state%energy_receipt(1))>&
        1d-10*max(1d0,abs(hybrid_state%energy_receipt(1))))&
      error stop 'hybrid DG RT initial physical energy does not match the checkpoint'
  endif
  call initialize_rt_dg_hybrid_stationarity(nproc_group_global,hybrid_state%owned_row_ids,&
    hybrid_state%density,current_total_energy,hybrid_state%coefficients,&
    apply_hybrid_metric_to_coefficients(),hybrid_state%occupations,current_electron_count,&
    current_hamiltonian_residual,stationarity_reference,ok,message)
  if(.not.ok)error stop 'hybrid DG RT stationarity reference failed'
  stationarity_tolerances=[dg_dc_gs_final_density_tolerance,dg_dc_gs_final_orbital_tolerance,&
    dg_dc_gs_final_orbital_tolerance,dg_dc_gs_electron_count_tolerance,dg_dc_gs_final_orbital_tolerance]
  call build_rt_dg_sparse_exchange(nproc_group_global,hybrid_state%global_count,hybrid_state%metric%fingerprint,&
    hybrid_state%metric%owned_row_ids,hybrid_state%metric%column_ids,metric_exchange,ok,message)
  if(.not.ok)error stop 'hybrid DG RT metric exchange setup failed'
  call build_rt_dg_sparse_exchange(nproc_group_global,hybrid_state%global_count,hybrid_state%operator_structure_fingerprint,&
    hybrid_state%operators%owned_row_ids,hybrid_state%operators%column_ids,operator_exchange,ok,message)
  if(.not.ok)error stop 'hybrid DG RT operator exchange setup failed'
  allocate(vector_potential_samples(3,0:nt+1));call calc_Ac_ext_t(0d0,dt,0,nt+1,vector_potential_samples)
  previous_polarization=0d0
  periods=[max(1d0,sqrt(sum(system%primitive_a(:,1)**2))),max(1d0,sqrt(sum(system%primitive_a(:,2)**2))),&
    max(1d0,sqrt(sum(system%primitive_a(:,3)**2)))]
  do step=1,nt
    call reconstruct_rt_dg_hybrid_density(nproc_group_global,hybrid_state,ok,message)
    if(.not.ok)error stop 'hybrid DG RT density reconstruction failed'
    call update_rt_dg_hybrid_density(nproc_group_global,hybrid_state,hybrid_state%density,project_salmon_local_rows,ok,message)
    update_count=update_count+1
    if(.not.ok)error stop 'hybrid DG RT Hartree/XC update failed'
    electric_field=-(vector_potential_samples(:,step)-vector_potential_samples(:,step-1))/dt
    if(stationarity_enabled.and.maxval(abs(electric_field))<=10d0*epsilon(1d0))then
      call evaluate_hybrid_rt_physical_invariants(current_total_energy,current_electron_count,&
        current_hamiltonian_residual,ok,message)
      if(.not.ok)error stop 'hybrid DG RT stationarity invariant evaluation failed'
      call evaluate_rt_dg_hybrid_stationarity(nproc_group_global,stationarity_reference,&
        hybrid_state%density,current_total_energy,hybrid_state%coefficients,&
        apply_hybrid_metric_to_coefficients(),current_electron_count,current_hamiltonian_residual,&
        stationarity_tolerances,stationarity_receipt,ok,message)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'hybrid DG RT zero-field stationarity failed';endif
      if(nproc_id_global==0)write(*,'(a,i0,5(a,es16.8))')'[HYBRID-RT-STATIONARITY] step=',step,&
        ' density=',stationarity_receipt%density_drift,' energy=',stationarity_receipt%energy_drift,&
        ' projector=',stationarity_receipt%projector_drift,' electron=',stationarity_receipt%electron_drift,&
        ' h_residual=',stationarity_receipt%hamiltonian_residual
    endif
    do orbital=1,hybrid_state%noccupied
      call propagate_rt_dg_hybrid_length_gauge(nproc_group_global,hybrid_state%metric,hybrid_state%operators,&
        hybrid_state%coefficients(:,orbital),electric_field,dt,1d-12,24,previous_polarization,periods,next,&
        metric_norm,orbital_energy,polarization,iterations,workspace,fingerprint,ok,message,&
        metric_exchange,operator_exchange)
      if(.not.ok)then;write(0,'(a)')trim(message);error stop 'hybrid DG RT propagation failed';endif
      hybrid_state%coefficients(:,orbital)=next
    enddo
    previous_polarization=polarization
  enddo
  if(update_count/=nt+1)error stop 'hybrid DG RT density update schedule violated'
end subroutine run_dg_hybrid_continuation_rt

function apply_hybrid_metric_to_coefficients() result(s_coefficients)
    complex(8),allocatable::s_coefficients(:,:),global_coefficients(:,:)
    integer::i,j,edge,ierr_local
    allocate(global_coefficients(hybrid_state%global_count,hybrid_state%noccupied),&
      s_coefficients(size(hybrid_state%owned_row_ids),hybrid_state%noccupied))
    global_coefficients=(0d0,0d0);s_coefficients=(0d0,0d0)
    do i=1,size(hybrid_state%owned_row_ids)
      global_coefficients(int(hybrid_state%owned_row_ids(i)),:)=hybrid_state%coefficients(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,size(global_coefficients),MPI_DOUBLE_COMPLEX,MPI_SUM,&
      nproc_group_global,ierr_local)
    if(ierr_local/=MPI_SUCCESS)error stop 'hybrid DG RT metric coefficient redistribution failed'
    do i=1,size(hybrid_state%owned_row_ids)
      do edge=hybrid_state%metric%row_offsets(i),hybrid_state%metric%row_offsets(i+1)-1
        j=hybrid_state%metric%column_ids(edge)
        s_coefficients(i,:)=s_coefficients(i,:)+hybrid_state%metric%values(edge)*global_coefficients(j,:)
      enddo
    enddo
end function apply_hybrid_metric_to_coefficients

subroutine evaluate_hybrid_rt_physical_invariants(total,electron_count,h_residual,evaluate_ok,evaluate_message)
    real(8),intent(out)::total,electron_count,h_residual
    logical,intent(out)::evaluate_ok
    character(*),intent(out)::evaluate_message
    complex(8),allocatable::global_coefficients(:,:),s_coefficients(:,:),h_coefficients(:,:)
    real(8)::kinetic_energy,nonlocal_energy,local_norms(2),global_norms(2)
    integer::i,j,edge,state_index,ierr_local
    call evaluate_dg_hybrid_fixed_energy(nproc_group_global,hybrid_state%owned_row_ids,&
      hybrid_state%coefficients,hybrid_state%occupations,hybrid_state%kinetic_rows,&
      hybrid_state%sipg_rows,hybrid_state%nonlocal_rows,kinetic_energy,nonlocal_energy,evaluate_ok,evaluate_message)
    if(.not.evaluate_ok)return
    energy%E_kin=kinetic_energy;energy%E_ion_nloc=nonlocal_energy
    call calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,.false.,energy)
    total=energy%E_tot;electron_count=sum(hybrid_state%occupations)
    allocate(global_coefficients(hybrid_state%global_count,hybrid_state%noccupied),&
      s_coefficients(size(hybrid_state%owned_row_ids),hybrid_state%noccupied),&
      h_coefficients(size(hybrid_state%owned_row_ids),hybrid_state%noccupied))
    global_coefficients=(0d0,0d0);s_coefficients=(0d0,0d0);h_coefficients=(0d0,0d0)
    do i=1,size(hybrid_state%owned_row_ids)
      global_coefficients(int(hybrid_state%owned_row_ids(i)),:)=hybrid_state%coefficients(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,size(global_coefficients),MPI_DOUBLE_COMPLEX,MPI_SUM,&
      nproc_group_global,ierr_local)
    if(ierr_local/=MPI_SUCCESS)then;evaluate_ok=.false.;evaluate_message='hybrid coefficient redistribution failed';return;endif
    do i=1,size(hybrid_state%owned_row_ids)
      do edge=hybrid_state%metric%row_offsets(i),hybrid_state%metric%row_offsets(i+1)-1
        j=hybrid_state%metric%column_ids(edge)
        s_coefficients(i,:)=s_coefficients(i,:)+hybrid_state%metric%values(edge)*global_coefficients(j,:)
      enddo
      do edge=hybrid_state%operators%row_offsets(i),hybrid_state%operators%row_offsets(i+1)-1
        j=hybrid_state%operators%column_ids(edge)
        h_coefficients(i,:)=h_coefficients(i,:)+hybrid_state%operators%hamiltonian_values(edge)*global_coefficients(j,:)
      enddo
    enddo
    local_norms=0d0
    do state_index=1,hybrid_state%noccupied
      local_norms(1)=local_norms(1)+sum(abs(h_coefficients(:,state_index)-&
        hybrid_state%eigenvalues(state_index)*s_coefficients(:,state_index))**2)
      local_norms(2)=local_norms(2)+sum(abs(h_coefficients(:,state_index))**2)
    enddo
    call MPI_Allreduce(local_norms,global_norms,2,MPI_DOUBLE_PRECISION,MPI_SUM,nproc_group_global,ierr_local)
    h_residual=sqrt(global_norms(1))/max(1d0,sqrt(global_norms(2)))
    evaluate_ok=ierr_local==MPI_SUCCESS.and.ieee_is_finite(total).and.ieee_is_finite(h_residual)
    if(evaluate_ok)then;evaluate_message='';else;evaluate_message='nonfinite hybrid RT physical invariants';endif
end subroutine evaluate_hybrid_rt_physical_invariants

subroutine project_salmon_local_rows(row_ids,grid_ids,density,local_rows,callback_ok,callback_message)
    integer(8),intent(in)::row_ids(:),grid_ids(:)
    real(8),intent(in)::density(:)
    complex(8),intent(out)::local_rows(:,:)
    logical,intent(out)::callback_ok
    character(*),intent(out)::callback_message
    real(8),allocatable::density_on_grid(:),potential_on_basis_grid(:),potential_source(:)
    complex(8),allocatable::row_contribution(:),reduced_row(:)
    integer(8),allocatable::local_grid_ids(:)
    integer(8)::workspace_peak
    integer::p,j,ix,iy,iz,local_grid_count,ierr,row,owner,row_position,rank,nproc,row_failed,global_row_failed
    logical::redistribution_ok
    character(256)::redistribution_message
    call MPI_Comm_rank(nproc_group_global,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;callback_ok=.false.;callback_message='physical callback rank query failed';return;endif
    call MPI_Comm_size(nproc_group_global,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;callback_ok=.false.;callback_message='physical callback size query failed';return;endif
    local_grid_count=product(mg%ie-mg%is+1)
    allocate(local_grid_ids(local_grid_count),potential_source(local_grid_count));p=0
    do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
      p=p+1;local_grid_ids(p)=int(ix,8)+int(lg%num(1),8)*(int(iy-1,8)+int(lg%num(2),8)*int(iz-1,8))
    enddo;enddo;enddo
    call redistribute_dg_row_owned_real_field_to_requests(nproc_group_global,int(product(lg%num),8),grid_ids,density,&
      local_grid_ids,density_on_grid,workspace_peak,redistribution_ok,redistribution_message)
    if(.not.redistribution_ok)then
      callback_ok=.false.;callback_message='physical density redistribution failed: '//trim(redistribution_message);return
    endif
    p=0
    do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
      p=p+1;rho_s(1)%f(ix,iy,iz)=density_on_grid(p)
    enddo;enddo;enddo
    rho%f=rho_s(1)%f
    call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,Vh)
    call exchange_correlation_density(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,stencil,Vxc,energy%E_xc)
    call update_vlocal(mg,system%nspin,Vh,Vpsl,Vxc,V_local)
    p=0
    do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
      p=p+1;potential_source(p)=V_local(1)%f(ix,iy,iz)
    enddo;enddo;enddo
    call redistribute_dg_row_owned_real_field_to_requests(nproc_group_global,int(product(lg%num),8),local_grid_ids,&
      potential_source,grid_ids,potential_on_basis_grid,workspace_peak,redistribution_ok,redistribution_message)
    if(.not.redistribution_ok)then
      callback_ok=.false.;callback_message='physical potential redistribution failed: '//trim(redistribution_message);return
    endif
    allocate(row_contribution(hybrid_state%global_count),reduced_row(hybrid_state%global_count));local_rows=(0d0,0d0)
    global_row_failed=0
    do row=1,hybrid_state%global_count
      row_contribution=(0d0,0d0)
      do p=1,size(grid_ids);do j=1,hybrid_state%global_count
        row_contribution(j)=row_contribution(j)+hybrid_state%grid_weights(p)*&
          conjg(hybrid_state%basis_values(row,p))*hybrid_state%basis_values(j,p)*potential_on_basis_grid(p)
      enddo;enddo
      owner=mod(hybrid_state%global_count-row,nproc);reduced_row=(0d0,0d0)
      call MPI_Reduce(row_contribution,reduced_row,hybrid_state%global_count,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,&
        nproc_group_global,ierr)
      row_failed=merge(0,1,ierr==MPI_SUCCESS)
      call MPI_Allreduce(row_failed,global_row_failed,1,MPI_INTEGER,MPI_MAX,nproc_group_global,ierr)
      if(ierr/=MPI_SUCCESS.or.global_row_failed/=0)exit
      if(rank==owner)then
        row_position=findloc(row_ids,int(row,8),dim=1)
        if(row_position>0)local_rows(row_position,:)=reduced_row
      endif
    enddo
    callback_ok=ierr==MPI_SUCCESS.and.global_row_failed==0
    if(callback_ok)then;callback_message='';else;callback_message='local-potential projection failed';endif
end subroutine project_salmon_local_rows

subroutine run_dg_overlapping_wannier_coefficient_rt()
  type(s_dg_overlapping_wannier_checkpoint)::checkpoint
  type(s_dg_overlapping_wannier_rt_state)::state
  complex(8),allocatable::coefficients(:,:)
  real(8),allocatable::vector_potential_samples(:,:)
  real(8)::electric_field(3),vector_potential(3),acceptance_gates(6),&
    polarization(3),current(3),cell_volume
  integer,allocatable::row_ids(:)
  integer::step,rank,ierr
  logical::ok,reusable
  character(256)::message
#ifdef USE_MPI
  call MPI_Comm_rank(nproc_group_global,rank,ierr)
#else
  rank=0
#endif
  acceptance_gates=[dg_dc_gs_final_density_tolerance,dg_dc_gs_final_orbital_tolerance,&
    10d0*dg_dc_gs_final_orbital_tolerance,dg_dc_gs_electron_count_tolerance,&
    1d0/dg_dc_metric_rank_tolerance,dg_ow_symmetry_tolerance]
  call read_dg_overlapping_wannier_checkpoint(nproc_group_global,'./overlapping_wannier_gs',&
    0,0,0_8,0_8,acceptance_gates,checkpoint,reusable,ok,message)
  if(.not.ok.or..not.reusable)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'accepted V3 overlapping-Wannier checkpoint is required'
  endif
  if(trim(checkpoint%field_coupling_convention)/='cell_wrapped_length_velocity')&
    error stop 'unsupported overlapping-Wannier field convention'
  cell_volume=product(al)
  if(cell_volume<=0d0)error stop 'invalid overlapping-Wannier RT cell volume'
  ! The V3 reader certifies that every retained tail covers each physical
  ! periodic-grid id at least once; overlapping buffers may repeat IDs.
  ! With basis updates forbidden, every
  ! coefficient combination remains in that closed periodic support, so
  ! a nonzero representational tail escape is structurally impossible.
  allocate(row_ids(size(checkpoint%overlap_row_ids)));row_ids=int(checkpoint%overlap_row_ids)
  allocate(coefficients,source=checkpoint%coefficients)
  call initialize_dg_overlapping_wannier_rt(nproc_group_global,row_ids,checkpoint%overlap,&
    checkpoint%hamiltonian0,checkpoint%position,checkpoint%velocity,&
    checkpoint%basis_generation,checkpoint%geometry_generation,checkpoint%basis_fingerprint,&
    checkpoint%operator_fingerprint,checkpoint%hamiltonian_fingerprint,&
    checkpoint%observable_fingerprint,checkpoint%field_coupling_convention,&
    checkpoint%basis_generation,checkpoint%geometry_generation,&
    checkpoint%basis_fingerprint,checkpoint%operator_fingerprint,coefficients,state,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT initialization failed'
  endif
  if(yn_dg_overlapping_wannier_rt_restart=='y')then
    call read_dg_overlapping_wannier_rt_restart(nproc_group_global,&
      './overlapping_wannier_rt.restart',coefficients,state,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a)')trim(message)
      error stop 'overlapping-Wannier coefficient RT restart failed'
    endif
  endif
  electric_field=0d0
  call evaluate_dg_overlapping_wannier_observables(nproc_group_global,coefficients,&
    checkpoint%occupations,cell_volume,state,polarization,current,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT observable evaluation failed'
  endif
  call write_dg_overlapping_wannier_rt_observable_sample(nproc_group_global,&
    './overlapping_wannier_rt_observables.dat',electric_field,polarization,current,&
    cell_volume,state,yn_dg_overlapping_wannier_rt_restart=='y',ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'overlapping-Wannier coefficient RT observable publication failed'
  endif
  allocate(vector_potential_samples(3,0:nt+1))
  call calc_Ac_ext_t(0d0,dt,0,nt+1,vector_potential_samples)
  do step=state%step+1,nt
    if(step==1.and.state%step==0.and.trim(ae_shape1)=='impulse')then
      ! The SALMON impulse is a vector-potential jump at t=0.  The value just
      ! before the first coefficient interval is zero, not the already-jumped
      ! sample stored at index zero.
      electric_field=-vector_potential_samples(:,step)/dt
    else
      electric_field=-(vector_potential_samples(:,step)-vector_potential_samples(:,step-1))/dt
    endif
    vector_potential=0d0
    call advance_dg_overlapping_wannier_rt(nproc_group_global,dt,electric_field,&
      vector_potential,coefficients,state,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'coefficient RT failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT propagation failed'
    endif
    call evaluate_dg_overlapping_wannier_observables(nproc_group_global,coefficients,&
      checkpoint%occupations,cell_volume,state,polarization,current,ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'observable evaluation failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT observable evaluation failed'
    endif
    call write_dg_overlapping_wannier_rt_observable_sample(nproc_group_global,&
      './overlapping_wannier_rt_observables.dat',electric_field,polarization,current,&
      cell_volume,state,yn_dg_overlapping_wannier_rt_restart=='y',ok,message)
    if(.not.ok)then
      if(rank==0)write(0,'(a,i0,2a)')'observable publication failed at step ',step,': ',trim(message)
      error stop 'overlapping-Wannier coefficient RT observable publication failed'
    endif
  enddo
  call write_dg_overlapping_wannier_rt_restart(nproc_group_global,&
    './overlapping_wannier_rt.restart',coefficients,state,ok,message)
  if(.not.ok)then
    if(rank==0)write(0,'(a)')trim(message)
    error stop 'cannot publish overlapping-Wannier coefficient RT restart'
  endif
end subroutine

subroutine write_initial_density_probe(system, info, mg, rho, rho_s, Vh, Vxc, Vpsl, label)
  use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_scalar
  use parallelization, only: nproc_id_global
  use communication, only: comm_summation, comm_get_max, comm_is_root
  implicit none
  type(s_dft_system),    intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_rgrid),         intent(in) :: mg
  type(s_scalar),        intent(in) :: rho, Vh, Vpsl
  type(s_scalar),        intent(in) :: rho_s(system%nspin), Vxc(system%nspin)
  character(*),          intent(in) :: label
  integer :: ix, iy, iz, ispin
  real(8) :: rho_vh_local, rho_vh_sum
  real(8) :: rho_vxc_local, rho_vxc_sum
  real(8) :: rho_vpsl_local, rho_vpsl_sum
  real(8) :: rho2_local, rho2_sum
  real(8) :: rho_max_local(1), rho_max_sum(1)

  rho_vh_local = 0.0d0
  rho_vxc_local = 0.0d0
  rho_vpsl_local = 0.0d0
  rho2_local = 0.0d0
  rho_max_local(1) = 0.0d0

  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        rho_max_local(1) = max(rho_max_local(1), rho%f(ix, iy, iz))
        rho2_local = rho2_local + rho%f(ix, iy, iz) * rho%f(ix, iy, iz)
        rho_vh_local = rho_vh_local + rho%f(ix, iy, iz) * Vh%f(ix, iy, iz)
        rho_vpsl_local = rho_vpsl_local + rho%f(ix, iy, iz) * Vpsl%f(ix, iy, iz)
        do ispin = 1, system%nspin
          rho_vxc_local = rho_vxc_local + rho_s(ispin)%f(ix, iy, iz) * Vxc(ispin)%f(ix, iy, iz)
        end do
      end do
    end do
  end do

  rho_vh_local = rho_vh_local * system%Hvol
  rho_vxc_local = rho_vxc_local * system%Hvol
  rho_vpsl_local = rho_vpsl_local * system%Hvol
  rho2_local = rho2_local * system%Hvol

  call comm_summation(rho_vh_local, rho_vh_sum, info%icomm_r)
  call comm_summation(rho_vxc_local, rho_vxc_sum, info%icomm_r)
  call comm_summation(rho_vpsl_local, rho_vpsl_sum, info%icomm_r)
  call comm_summation(rho2_local, rho2_sum, info%icomm_r)
  call comm_get_max(rho_max_local, rho_max_sum, 1, info%icomm_r)

  if (comm_is_root(nproc_id_global)) then
    write(*,'(1x,a,a,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
      '        ', trim(label), ': rhoVh=', rho_vh_sum, ' rhoVxc=', rho_vxc_sum, &
      ' rhoVpsl=', rho_vpsl_sum, ' rho2=', rho2_sum, ' rhomax=', rho_max_sum(1)
    flush(6)
  end if
end subroutine write_initial_density_probe

subroutine write_local_chern_marker_xy(itt, mg, system, info, psi_fin)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
  use communication, only: comm_is_root, comm_summation
  use rt_local_chern_marker, only: compute_local_chern_marker_from_orbital
  use rt_local_chern_marker_soi, only: compute_local_chern_marker_from_orbital_soi => compute_local_chern_marker_from_orbital
  use filesystem, only: create_directory
  use inputoutput, only: t_unit_length
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, sysname, yn_spinorbit
  implicit none
  integer, intent(in) :: itt
  type(s_rgrid), intent(in) :: mg
  type(s_dft_system), intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_orbital), intent(in) :: psi_fin
  real(8), allocatable :: marker(:,:,:)
  real(8), allocatable :: marker_xy(:,:), marker_xy_full_local(:,:), marker_xy_full(:,:)
  character(256) :: filename, filenum, map_directory
  integer :: ix, iy, iz, iunit, nx, ny

  allocate(marker(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
  allocate(marker_xy(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2)))
  if (yn_spinorbit == 'y') then
    call compute_local_chern_marker_from_orbital_soi(mg, system, info, psi_fin, marker)
  else
    call compute_local_chern_marker_from_orbital(mg, system, info, psi_fin, marker)
  end if

  marker_xy(:,:) = 0.0d0
  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      do ix = mg%is(1), mg%ie(1)
        marker_xy(ix,iy) = marker_xy(ix,iy) + marker(ix,iy,iz) * system%hgs(3)
      end do
    end do
  end do

  nx = maxval(mg%ie_all(1,:))
  ny = maxval(mg%ie_all(2,:))
  allocate(marker_xy_full_local(nx, ny), marker_xy_full(nx, ny))
  marker_xy_full_local(:,:) = 0.0d0
  marker_xy_full(:,:) = 0.0d0
  do iy = mg%is(2), mg%ie(2)
    do ix = mg%is(1), mg%ie(1)
      marker_xy_full_local(ix,iy) = marker_xy(ix,iy)
    end do
  end do
  call comm_summation(marker_xy_full_local, marker_xy_full, nx*ny, info%icomm_r)

  if (comm_is_root(nproc_id_global)) then
    write(filenum, '(i6.6)') itt
    map_directory = trim(base_directory)//trim(sysname)//'_lcm_xy/'
    call create_directory(trim(map_directory))
    filename = trim(map_directory)//trim(sysname)//'_lcm_xy_'//trim(adjustl(filenum))//'.data'
    open(newunit=iunit, file=trim(filename), status='replace', action='write')
    write(iunit,'(a)') '# Local Chern marker integrated along z'
    write(iunit,'(a)') '# x: x coordinate'
    write(iunit,'(a)') '# y: y coordinate'
    write(iunit,'(a)') '# local_chern_marker_zint: local Chern marker integrated over z'
    write(iunit,'(a,a,a,a,a)') '# 1:x[', trim(t_unit_length%name), '] 2:y[', &
      trim(t_unit_length%name), '] 3:local_chern_marker_zint[none]'
    do iy = 1, ny
      do ix = 1, nx
        write(iunit,'(3(1x,es24.16))') dble(ix-1) * system%hgs(1) * t_unit_length%conv, &
                                       dble(iy-1) * system%hgs(2) * t_unit_length%conv, marker_xy_full(ix,iy)
      end do
      write(iunit,*)
    end do
    close(iunit)
  end if
  deallocate(marker_xy_full, marker_xy_full_local, marker_xy, marker)
end subroutine write_local_chern_marker_xy

subroutine print_header()
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use salmon_global, only: iperiodic, yn_jm
  implicit none
  !(header of standard output)
  if(comm_is_root(nproc_id_global))then
    write(*,*)
    select case(iperiodic)
    case(0)
      write(*,'(1x,a10,a11,a48,a15,a18,a10)') &
                  "time-step ", "time[fs]",   &
                  "Dipole moment(xyz)[A]"     &
                ,"electrons", "Total energy[eV]", "iterVh"
    case(3)
      if(yn_jm=='n')then
        write(*,'(1x,a10,a11,a48,a15,a18)')   &
                    "time-step", "time[fs] ", &
                    "Current(xyz)[a.u.]",     &
                    "electrons", "Total energy[eV] "
      else
        write(*,'(1x,a10,a11,a48,a15,a18)')   &
                    "time-step", "time[fs] ", &
                    "Current(xyz)[a.u.]",     &
                    "electrons"
      end if
    end select
    write(*,'("#",7("----------"))')
  endif
end subroutine print_header

end subroutine main_tddft
