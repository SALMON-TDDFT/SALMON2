! Test-only legacy reference; excluded from the SALMON production build.
#include "config.h"
module dg_hybrid_scf
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine dg_hybrid_update_potential(input_density,ok)
      import real64
      real(real64),intent(in)::input_density(:);logical,intent(out)::ok
    end subroutine dg_hybrid_update_potential
    subroutine dg_hybrid_assemble_hamiltonian(iteration,ok)
      integer,intent(in)::iteration;logical,intent(out)::ok
    end subroutine dg_hybrid_assemble_hamiltonian
    subroutine dg_hybrid_solve_occupied_states(iteration,band_energy_indicator,eigensystem_residual,electron_count_defect,&
        symmetry_defect,ok)
      import real64
      integer,intent(in)::iteration
      real(real64),intent(out)::band_energy_indicator,eigensystem_residual,electron_count_defect,symmetry_defect
      logical,intent(out)::ok
    end subroutine dg_hybrid_solve_occupied_states
    subroutine dg_hybrid_reconstruct_density(output_density,ok)
      import real64
      real(real64),intent(out)::output_density(:);logical,intent(out)::ok
    end subroutine dg_hybrid_reconstruct_density
    subroutine dg_hybrid_mix_density(iteration,input_density,output_density,reset_history,reduce_rate,mixed_density,ok)
      import real64
      integer,intent(in)::iteration
      real(real64),intent(in)::input_density(:),output_density(:)
      logical,intent(in)::reset_history,reduce_rate
      real(real64),intent(out)::mixed_density(:);logical,intent(out)::ok
    end subroutine dg_hybrid_mix_density
  end interface
  public::run_dg_hybrid_self_consistent_ground_state
contains
  subroutine run_dg_hybrid_self_consistent_ground_state(comm,global_point_count,point_ids,initial_density,&
      hybrid_basis_fingerprint,metric_fingerprint,update_potential,assemble_hamiltonian,solve_occupied_states,&
      reconstruct_density,mix_density,maximum_iterations,density_tolerance,energy_tolerance,eigensystem_tolerance,&
      physical_tolerance,converged_density,iterations,density_residual,energy_residual,eigensystem_residual,&
      electron_count_defect,symmetry_defect,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,maximum_iterations
    integer(int64),intent(in)::point_ids(:),hybrid_basis_fingerprint,metric_fingerprint
    real(real64),intent(in)::initial_density(:),density_tolerance,energy_tolerance,eigensystem_tolerance,physical_tolerance
    procedure(dg_hybrid_update_potential)::update_potential
    procedure(dg_hybrid_assemble_hamiltonian)::assemble_hamiltonian
    procedure(dg_hybrid_solve_occupied_states)::solve_occupied_states
    procedure(dg_hybrid_reconstruct_density)::reconstruct_density
    procedure(dg_hybrid_mix_density)::mix_density
    real(real64),allocatable,intent(out)::converged_density(:)
    integer,intent(out)::iterations
    real(real64),intent(out)::density_residual,energy_residual,eigensystem_residual,electron_count_defect,symmetry_defect
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,nlocal,local_bad,global_bad,minimum_integer,maximum_integer,allocation_status,i
    integer,allocatable::ownership_count(:)
    integer(int64)::minimum_bits,maximum_bits,bits,local_hash,global_hash,entry_hash,quantized
    real(real64),allocatable::current_density(:),output_density(:),mixed_density(:),delta(:),previous_delta(:)
    real(real64)::band_energy_indicator,previous_band_energy,previous_density_residual,local_value,global_value,local_dot,global_dot,&
      previous_norm,current_norm,quantization_scale,quantization_limit
    logical::callback_ok,reset_history,reduce_rate,converged
    ok=.false.;message='';iterations=0;fingerprint=0_int64;density_residual=huge(1d0);energy_residual=huge(1d0)
    eigensystem_residual=huge(1d0);electron_count_defect=huge(1d0);symmetry_defect=huge(1d0)
    nlocal=size(point_ids);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='hybrid SCF communicator failed';return;endif
    call agree_integer(global_point_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing hybrid SCF extent';return;endif
    call agree_integer(maximum_iterations,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing hybrid SCF iteration cap';return;endif
    call agree_receipt(hybrid_basis_fingerprint,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing hybrid basis fingerprint';return;endif
    call agree_receipt(metric_fingerprint,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing hybrid metric fingerprint';return;endif
    call agree_real(density_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing density tolerance';return;endif
    call agree_real(energy_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing energy tolerance';return;endif
    call agree_real(eigensystem_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing eigensystem tolerance';return;endif
    call agree_real(physical_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing physical tolerance';return;endif
    if(global_point_count<1.or.maximum_iterations<1.or.size(initial_density)/=nlocal)local_bad=1
    if(any(point_ids<1_int64).or.any(point_ids>int(max(0,global_point_count),int64)))local_bad=1
    if(.not.all(ieee_is_finite(initial_density)).or.any(initial_density<0d0))local_bad=1
    if(.not.valid_tolerance(density_tolerance).or..not.valid_tolerance(energy_tolerance).or.&
      .not.valid_tolerance(eigensystem_tolerance).or..not.valid_tolerance(physical_tolerance))local_bad=1
    if(hybrid_basis_fingerprint==0_int64.or.metric_fingerprint==0_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid SCF contract';return;endif
    allocate(ownership_count(global_point_count),current_density(nlocal),output_density(nlocal),mixed_density(nlocal),&
      delta(nlocal),previous_delta(nlocal),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid SCF workspace';return;endif
    ownership_count=0;do i=1,nlocal;ownership_count(int(point_ids(i)))=ownership_count(int(point_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_point_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;call cleanup();message='hybrid SCF points are not owned exactly once';return;endif
    current_density=initial_density;previous_delta=0d0;previous_band_energy=huge(1d0)
    previous_density_residual=huge(1d0);converged=.false.
    do iterations=1,maximum_iterations
      call update_potential(current_density,callback_ok);call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid SCF potential update failed';return;endif
      call assemble_hamiltonian(iterations,callback_ok);call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid SCF Hamiltonian assembly failed';return;endif
      call solve_occupied_states(iterations,band_energy_indicator,eigensystem_residual,electron_count_defect,symmetry_defect,callback_ok)
      local_bad=merge(0,1,callback_ok.and.ieee_is_finite(band_energy_indicator).and.ieee_is_finite(eigensystem_residual).and.&
        ieee_is_finite(electron_count_defect).and.ieee_is_finite(symmetry_defect))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid occupied solve failed';return;endif
      call reconstruct_density(output_density,callback_ok)
      local_bad=merge(0,1,callback_ok.and.all(ieee_is_finite(output_density)).and.all(output_density>=0d0))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid density reconstruction failed';return;endif
      delta=output_density-current_density;local_value=sum(delta**2)
      call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      density_residual=sqrt(global_value/real(global_point_count,real64))
      if(iterations==1)then
        energy_residual=huge(1d0)
      else
        ! This is an occupied band-energy stability indicator, not a total-DFT-energy difference.
        energy_residual=abs(band_energy_indicator-previous_band_energy)
      endif
      converged=iterations>1.and.density_residual<=density_tolerance.and.energy_residual<=energy_tolerance.and.&
        eigensystem_residual<=eigensystem_tolerance.and.electron_count_defect<=physical_tolerance.and.&
        symmetry_defect<=physical_tolerance
      if(rank==0)write(*,'(a,i0,6(a,es16.8))')'[HYBRID-SCF] iteration=',iterations,&
        ' density=',density_residual,' band_energy_change=',energy_residual,&
        ' eigensystem=',eigensystem_residual,' electrons=',electron_count_defect,&
        ' symmetry=',symmetry_defect,' band_energy=',band_energy_indicator
      if(converged)exit
      reset_history=.false.;reduce_rate=.false.
      if(iterations>1)then
        local_dot=sum(delta*previous_delta);local_value=sum(delta**2);current_norm=sum(previous_delta**2)
        call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        call MPI_Allreduce(local_value,previous_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,current_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        if(global_dot< -0.25d0*sqrt(max(0d0,previous_norm*current_norm)).or.&
          density_residual>1.2d0*previous_density_residual)then;reset_history=.true.;reduce_rate=.true.;endif
      endif
      if(rank==0.and.reset_history)write(*,'(a,i0)')'[HYBRID-SCF] rejected density history at iteration=',iterations
      call mix_density(iterations,current_density,output_density,reset_history,reduce_rate,mixed_density,callback_ok)
      local_bad=merge(0,1,callback_ok.and.all(ieee_is_finite(mixed_density)).and.all(mixed_density>=0d0))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid density mixing failed';return;endif
      previous_delta=delta;previous_density_residual=density_residual
      previous_band_energy=band_energy_indicator;current_density=mixed_density
    enddo
    if(.not.converged)then;call cleanup();message='hybrid self-consistent ground state did not converge';return;endif
    allocate(converged_density(nlocal),stat=allocation_status);local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate converged hybrid density';return;endif
    converged_density=output_density;quantization_scale=1000d0*density_tolerance
    quantization_limit=0.25d0*real(huge(0_int64),real64)*quantization_scale;local_hash=0_int64
    do i=1,nlocal
      if(converged_density(i)>quantization_limit)then;local_bad=1;cycle;endif
      quantized=nint(converged_density(i)/quantization_scale,int64)
      entry_hash=ieor(point_ids(i),ishftc(quantized,19));local_hash=ieor(local_hash,ishftc(entry_hash,mod(int(point_ids(i)),63)))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid SCF fingerprint range is unsafe';return;endif
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    fingerprint=ieor(global_hash,hybrid_basis_fingerprint);fingerprint=ieor(fingerprint,ishftc(metric_fingerprint,13))
    if(fingerprint==0_int64)fingerprint=1879_int64;ok=.true.;message='';call cleanup(.true.)
  contains
    subroutine agree_integer(value,minimum_value,maximum_value,status)
      integer,intent(in)::value;integer,intent(out)::minimum_value,maximum_value,status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine agree_integer
    subroutine agree_receipt(value,minimum_value,maximum_value,status)
      integer(int64),intent(in)::value;integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
    end subroutine agree_receipt
    subroutine agree_real(value,minimum_value,maximum_value,status)
      real(real64),intent(in)::value;integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::status
      integer(int64)::value_bits;value_bits=transfer(value,value_bits);call agree_receipt(value_bits,minimum_value,maximum_value,status)
    end subroutine agree_real
    subroutine callback_consensus(callback_result,bad,status)
      logical,intent(in)::callback_result;integer,intent(out)::bad,status;integer::local
      local=merge(0,1,callback_result);call MPI_Allreduce(local,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine callback_consensus
    logical function valid_tolerance(value)
      real(real64),intent(in)::value
      valid_tolerance=ieee_is_finite(value).and.value>=1d-15.and.value<=1d-2
    end function valid_tolerance
    subroutine cleanup(keep_output)
      logical,intent(in),optional::keep_output;logical::keep
      keep=.false.;if(present(keep_output))keep=keep_output
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(current_density))deallocate(current_density)
      if(allocated(output_density))deallocate(output_density)
      if(allocated(mixed_density))deallocate(mixed_density)
      if(allocated(delta))deallocate(delta)
      if(allocated(previous_delta))deallocate(previous_delta)
      if(.not.keep.and.allocated(converged_density))deallocate(converged_density)
    end subroutine cleanup
#else
    ok=.false.;message='MPI is required for hybrid self-consistent ground state';iterations=0;fingerprint=0_int64
#endif
  end subroutine run_dg_hybrid_self_consistent_ground_state
end module dg_hybrid_scf
