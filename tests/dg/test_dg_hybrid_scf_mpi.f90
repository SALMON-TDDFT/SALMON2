#include "config.h"
program test_dg_hybrid_scf_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_scf,only:run_dg_hybrid_self_consistent_ground_state
  use dg_hybrid_generalized_eigensystem,only:solve_dg_hybrid_generalized_scalapack
  use dg_hybrid_block_cg,only:solve_dg_hybrid_block_cg
  implicit none
  integer,parameter::nglobal=8
  integer,parameter::nstate=2
  integer::comm,rank,nproc,ierr,nlocal,p,position,iterations,reset_count,backend,inner_iterations
  integer(int64),allocatable::point_ids(:)
  real(real64),allocatable::density(:),converged_density(:),target(:),output_density(:)
  real(real64),allocatable::reference_density(:)
  real(real64)::density_residual,energy_residual,eigensystem_residual,electron_count_defect,symmetry_defect
  integer(int64)::fingerprint,reference_fingerprint
  integer(int64)::inner_workspace,inner_fingerprint
  complex(real64),allocatable::hrows(:,:),srows(:,:),occupied_coefficients(:,:),initial_coefficients(:,:)
  complex(real64)::hfull(nglobal,nglobal),sfull(nglobal,nglobal),global_input(nglobal,2*nstate)
  real(real64)::inner_eigenvalues(nstate),inner_residual,inner_orthogonality,inner_projector,reference_energy,&
    outer_solver_residual,defect
  logical::ok,oscillatory,physical_scf
  character(64)::stop_reason
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(point_ids(nlocal),density(nlocal),target(nlocal),output_density(nlocal));position=0
  do p=nglobal,1,-1
    if(mod(p-1,nproc)/=rank)cycle
    position=position+1;point_ids(position)=p;target(position)=0.5d0+0.02d0*p
  enddo
  density=0.2d0;oscillatory=.false.;physical_scf=.false.;reset_count=0
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    40,1d-9,1d-10,1d-10,1d-10,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(maxval(abs(converged_density-target))<2d-8.and.iterations>1.and.iterations<=40,&
    'hybrid SCF did not converge to the nonlinear fixed point')
  oscillatory=.true.;density=0.2d0;reset_count=0
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    60,1d-8,1d-9,1d-10,1d-10,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok.and.reset_count>0,'oscillating hybrid SCF did not reset/reduce Pulay history')

  ! Both inner solvers must converge the same nonlinear fixed-basis problem.
  physical_scf=.true.;oscillatory=.false.;allocate(hrows(nlocal,nglobal),srows(nlocal,nglobal),&
    initial_coefficients(nlocal,nstate),reference_density(nlocal))
  sfull=(0d0,0d0);do p=1,nglobal;sfull(p,p)=1d0;enddo
  do position=1,nlocal
    srows(position,:)=sfull(int(point_ids(position)),:)
    initial_coefficients(position,1)=cmplx(1d0/(int(point_ids(position))+1),0.02d0*int(point_ids(position)),real64)
    initial_coefficients(position,2)=cmplx(0.1d0*int(point_ids(position)),-1d0/(int(point_ids(position))+2),real64)
  enddo
  backend=1;density=0.25d0;outer_solver_residual=1d-2
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    50,2d-8,2d-9,2d-9,2d-9,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok,trim(message));reference_density=converged_density;reference_energy=inner_eigenvalues(1)+inner_eigenvalues(2)
  backend=2;density=0.25d0;outer_solver_residual=1d-2
  do position=1,nlocal
    initial_coefficients(position,1)=cmplx(1d0/(int(point_ids(position))+1),0.02d0*int(point_ids(position)),real64)
    initial_coefficients(position,2)=cmplx(0.1d0*int(point_ids(position)),-1d0/(int(point_ids(position))+2),real64)
  enddo
  call run_dg_hybrid_self_consistent_ground_state(comm,nglobal,point_ids,density,9911_int64,8822_int64,&
    update_potential,assemble_hamiltonian,solve_occupied_states,reconstruct_density,pulay_mix,&
    50,2d-8,2d-9,2d-9,2d-9,converged_density,iterations,density_residual,energy_residual,&
    eigensystem_residual,electron_count_defect,symmetry_defect,fingerprint,ok,message)
  call require(ok,trim(message));defect=maxval(abs(converged_density-reference_density))
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<2d-7.and.abs(sum(inner_eigenvalues)-reference_energy)<2d-7,&
    'ScaLAPACK and adaptive block-CG SCF fixed points differ')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_SCF ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid SCF on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine update_potential(input_density,callback_ok)
    real(real64),intent(in)::input_density(:);logical,intent(out)::callback_ok
    callback_ok=all(input_density>=0d0)
    if(physical_scf)then
      hfull=(0d0,0d0)
      do p=1,nglobal;hfull(p,p)=0.15d0*p+0.08d0*global_density_value(p,input_density);enddo
      do p=1,nglobal-1
        hfull(p,p+1)=cmplx(-0.018d0,0.006d0,real64);hfull(p+1,p)=conjg(hfull(p,p+1))
      enddo
      do position=1,nlocal;hrows(position,:)=hfull(int(point_ids(position)),:);enddo
    endif
  end subroutine update_potential
  subroutine assemble_hamiltonian(iteration,callback_ok)
    integer,intent(in)::iteration;logical,intent(out)::callback_ok
    callback_ok=iteration>0
  end subroutine assemble_hamiltonian
  subroutine solve_occupied_states(iteration,total_energy,residual,electron_defect,group_defect,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(out)::total_energy,residual,electron_defect,group_defect
    logical,intent(out)::callback_ok
    if(.not.physical_scf)then
      total_energy=sum((density-target)**2);call MPI_Allreduce(MPI_IN_PLACE,total_energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      residual=min(1d-12,1d-13*iteration);electron_defect=0d0;group_defect=0d0;callback_ok=.true.;return
    endif
    if(backend==1)then
      call solve_dg_hybrid_generalized_scalapack(comm,nglobal,nstate,point_ids,hrows,srows,1d-11,&
        occupied_coefficients,inner_eigenvalues,inner_residual,inner_orthogonality,inner_projector,&
        inner_workspace,inner_fingerprint,callback_ok,message)
    else
      call solve_dg_hybrid_block_cg(comm,nglobal,point_ids,initial_coefficients,apply_h,apply_s,&
        outer_solver_residual,1d-11,32,occupied_coefficients,inner_eigenvalues,inner_iterations,&
        inner_residual,stop_reason,inner_workspace,inner_fingerprint,callback_ok,message)
      if(callback_ok)initial_coefficients=occupied_coefficients
    endif
    total_energy=sum(inner_eigenvalues);residual=inner_residual;group_defect=0d0
    electron_defect=0d0;if(callback_ok)electron_defect=abs(global_norm(occupied_coefficients)-2d0)
  end subroutine solve_occupied_states
  subroutine reconstruct_density(reconstructed,callback_ok)
    real(real64),intent(out)::reconstructed(:);logical,intent(out)::callback_ok
    if(physical_scf)then
      reconstructed=sum(abs(occupied_coefficients)**2,dim=2)
    elseif(oscillatory)then;reconstructed=target-1.2d0*(density-target)
    else;reconstructed=target+0.25d0*(density-target);endif
    output_density=reconstructed;callback_ok=.true.
  end subroutine reconstruct_density
  subroutine pulay_mix(iteration,input_density,reconstructed,reset_history,reduce_rate,mixed_density,callback_ok)
    integer,intent(in)::iteration
    real(real64),intent(in)::input_density(:),reconstructed(:)
    logical,intent(in)::reset_history,reduce_rate
    real(real64),intent(out)::mixed_density(:);logical,intent(out)::callback_ok
    real(real64)::beta
    beta=0.7d0;if(reduce_rate)beta=0.2d0
    if(reset_history)reset_count=reset_count+1
    mixed_density=input_density+beta*(reconstructed-input_density);density=mixed_density
    outer_solver_residual=sqrt(sum((reconstructed-input_density)**2)/real(max(1,size(input_density)),real64))
    call MPI_Allreduce(MPI_IN_PLACE,outer_solver_residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);callback_ok=.true.
  end subroutine pulay_mix
  subroutine apply_h(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    call apply_dense(hfull,input,output);callback_ok=.true.
  end subroutine apply_h
  subroutine apply_s(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    call apply_dense(sfull,input,output);callback_ok=.true.
  end subroutine apply_s
  subroutine apply_dense(matrix,input,output)
    complex(real64),intent(in)::matrix(:,:),input(:,:);complex(real64),intent(out)::output(:,:)
    global_input(:,1:size(input,2))=(0d0,0d0)
    do position=1,nlocal;global_input(int(point_ids(position)),1:size(input,2))=input(position,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_input,nglobal*size(input,2),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do position=1,nlocal;output(position,:)=matmul(matrix(int(point_ids(position)),:),global_input(:,1:size(input,2)));enddo
  end subroutine apply_dense
  real(real64) function global_density_value(global_row,local_density)
    integer,intent(in)::global_row;real(real64),intent(in)::local_density(:);real(real64)::local_value
    local_value=0d0
    do position=1,nlocal;if(int(point_ids(position))==global_row)local_value=local_density(position);enddo
    call MPI_Allreduce(local_value,global_density_value,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  end function global_density_value
  real(real64) function global_norm(vectors)
    complex(real64),intent(in)::vectors(:,:);real(real64)::local_value
    local_value=sum(abs(vectors)**2);call MPI_Allreduce(local_value,global_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  end function global_norm
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_scf_mpi
