#include "config.h"
program test_dg_hybrid_low_energy_symmetry_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_ground_state_types,only:s_dg_hybrid_spectral_certification
  use dg_hybrid_low_energy_symmetry,only:select_dg_hybrid_symmetry_target,&
    evaluate_dg_hybrid_low_energy_symmetry,certify_dg_hybrid_energy_window
  implicit none
  integer,parameter::ncert=6
  integer::comm,rank,nproc,ierr,i,target_rank,between_rank,on_level_rank,differential_worst_operation
  real(real64)::occupied_defect,target_defect,energy_defect
  real(real64)::clustered(4),unresolved(3),eigenvalues(3),cert_eigenvalues(ncert)
  complex(real64)::metric(3,3),representation(3,3,2),coefficients(3,3)
  complex(real64)::cert_metric(ncert,ncert),cert_representation(ncert,ncert,2),&
    cert_coefficients(ncert,ncert),cert_action(ncert,ncert)
  type(s_dg_hybrid_spectral_certification)::certification
  integer(int64)::reference_fingerprint
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)

  metric=0d0;coefficients=(0d0,0d0);representation=(0d0,0d0)
  do i=1,3
    metric(i,i)=1d0;coefficients(i,i)=1d0;representation(i,i,1)=1d0
  enddo
  representation(1,1,2)=1d0;representation(2,2,2)=1d0;representation(3,3,2)=0.5d0
  eigenvalues=[-1d0,0d0,1d0]
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(ok,'closed low-energy space was rejected: '//trim(message))
  call require(max(occupied_defect,target_defect,energy_defect)<1d-13,&
    'unused high-energy nonunitarity contaminated target defects')

  representation(2,2,2)=0d0;representation(3,2,2)=1d0
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(.not.ok.and.target_defect>0.5d0,'target leakage was accepted')

  representation=(0d0,0d0)
  do i=1,3
    representation(i,i,1)=1d0
  enddo
  representation(1,2,2)=1d0;representation(2,1,2)=1d0;representation(3,3,2)=1d0
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
    1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(.not.ok.and.target_defect<1d-13.and.energy_defect>0.5d0,&
    'nondegenerate target energy mixing was accepted')

  if(nproc>1)then
    if(rank>0)coefficients(1,1)=2d0
    call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
      1,2,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
    call require(.not.ok.and.index(message,'orthonormal')>0,&
      'rank-local metric orthogonality failure was not rejected collectively')
    coefficients(1,1)=1d0
  endif

  clustered=[-1d0,0d0,0d0,1d0]
  call select_dg_hybrid_symmetry_target(comm,clustered,2,1d-12,target_rank,ok,message)
  call require(ok.and.target_rank==3,'degenerate target boundary was not extended')
  unresolved=[-1d0,0d0,0d0]
  call select_dg_hybrid_symmetry_target(comm,unresolved,2,1d-12,target_rank,ok,message)
  call require(.not.ok,'unresolved terminal degeneracy was accepted')
  call select_dg_hybrid_symmetry_target(comm,clustered,2+merge(1,0,rank>0),1d-12,target_rank,ok,message)
  if(nproc>1)call require(.not.ok,'rank-disagreeing target request was accepted')

  call reset_certification_fixture()
  cert_eigenvalues=[-1d0,-0.2d0,0.4d0,0.8d0,1.2d0,2d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),0d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%requested_rank==1.and.certification%certified_rank==1.and.&
    certification%proof_state_present.and.certification%proof_energy==cert_eigenvalues(2),&
    'zero energy window did not certify the occupied cluster')

  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1.1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok,'cutoff-between-levels certification failed: '//trim(message))
  between_rank=certification%requested_rank
  call require(between_rank==2.and.certification%certified_rank==2,&
    'cutoff between levels selected the wrong rank')

  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1.4d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok,'cutoff-on-level certification failed: '//trim(message))
  on_level_rank=certification%requested_rank
  call require(on_level_rank==3.and.on_level_rank/=between_rank,&
    'cutoff on a level or spectrum-dependent requested rank is incorrect')

  cert_eigenvalues=[-1d0,-0.2d0,0.05d0,0.8d0,1.2d0,2d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1.1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%requested_rank==3.and.certification%requested_rank/=between_rank,&
    'different complete spectra did not select different requested ranks for the same window')

  cert_eigenvalues=[-1d0,0d0,0d0,0.8d0,1.2d0,2d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%requested_rank==3.and.certification%boundary_cluster_rank==3,&
    'exactly degenerate cutoff cluster was not completed')

  cert_eigenvalues=[-1d0,0d0,5d-13,0.8d0,1.2d0,2d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%requested_rank==3.and.certification%boundary_cluster_rank==3,&
    'numerically split cutoff cluster was not completed')

  cert_eigenvalues=[-1d0,0d0,1d-6,0.7d0,1.5d0,2d0]
  cert_representation(2,2,2)=0d0;cert_representation(3,3,2)=0d0
  cert_representation(2,3,2)=1d0;cert_representation(3,2,2)=1d0
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1d0,2,1d-12,1d-5,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%requested_rank==2.and.certification%certified_rank==3.and.&
    certification%extension_states==1.and.abs(certification%extension_energy-1d-6)<1d-14.and.&
    certification%proof_energy==0.7d0.and.certification%worst_operation==2.and.&
    abs(certification%worst_operation_defect-max(certification%occupied_subspace_defect,&
      certification%target_subspace_defect,certification%target_energy_defect))<1d-14,&
    'first passing complete cluster was not selected with complete receipts')

  call reset_certification_fixture()
  cert_eigenvalues=[-1d0,0d0,0.4d0,0.8d0,1.2d0,2d0]
  cert_representation(2,2,2)=0d0;cert_representation(6,6,2)=0d0
  cert_representation(2,6,2)=1d0;cert_representation(6,2,2)=1d0
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(.not.ok.and.index(message,'capacity')>0,&
    'absence of a passing cluster before the basis ceiling was accepted')

  call reset_certification_fixture()
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),10d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(.not.ok.and.index(message,'capacity')>0.and.index(message,'proof')>0,&
    'energy space without a proof state was certified')

  cert_representation(1,1,2)=0d0;cert_representation(2,2,2)=0d0
  cert_representation(1,2,2)=1d0;cert_representation(2,1,2)=1d0
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1d0,2,1d-12,1d-5,0d0,0d0,certification,ok,message)
  call require(.not.ok.and.index(message,'occupied')>0,&
    'occupied symmetry failure was repaired by adding empty states')

  call reset_certification_fixture()
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),0d0,2,1d-12,1d-5,2d-5,0d0,certification,ok,message)
  call require(.not.ok.and.index(message,'occupied-projector')>0,&
    'occupied-projector failure was repaired by adding empty states')
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),0d0,2,1d-12,1d-5,0d0,2d-5,certification,ok,message)
  call require(.not.ok.and.index(message,'density')>0,&
    'density failure was repaired by adding empty states')
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),0d0,2,1d-12,1d-5,5d-6,4d-6,certification,ok,message)
  call require(ok.and.abs(certification%maximum_physical_defect-5d-6)<1d-14.and.&
    certification%worst_operation_defect<1d-12,&
    'external physical defects were conflated with the worst symmetry operation')

  cert_eigenvalues=[-1d0,0d0,0d0,1d0,2d0,3d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),-1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%compatibility_dynamic_rank.and.&
    certification%requested_rank==2.and.certification%certified_rank==3,&
    'exact -1 compatibility path did not retain dynamic-rank selection')
  reference_fingerprint=certification%fingerprint
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),-1d0,ncert,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%compatibility_dynamic_rank.and.&
    certification%certified_rank==ncert.and..not.certification%proof_state_present,&
    'exact -1 compatibility path changed legacy full-rank acceptance')
  cert_eigenvalues=[-1d0,0d0,1d0,2d0,2d0,2d0]
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),-1d0,4,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%compatibility_dynamic_rank.and.certification%requested_rank==4.and.&
    certification%certified_rank==ncert.and..not.certification%proof_state_present,&
    'exact -1 compatibility path lost the legacy unresolved-tail full-rank fallback')
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),-1d0-epsilon(1d0),2,1d-12,1d-10,0d0,0d0,certification,ok,message)
  call require(.not.ok,'a negative non-sentinel energy window was accepted')

  if(nproc>1)then
    cert_eigenvalues=[-1d0,0d0,merge(0.1d0,0d0,rank>0),1d0,2d0,3d0]
    call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
      1,cert_eigenvalues(1),1d0,2,1d-12,1d-10,0d0,0d0,certification,ok,message)
    call require(.not.ok.and.index(message,'rank')>0,&
      'rank-dependent complete LCFO spectrum was accepted')
  endif

  call configure_nontrivial_certification_fixture()
  cert_eigenvalues=[-1d0,-0.2d0,0.4d0,0.8d0,1.2d0,2d0]
  call evaluate_dg_hybrid_low_energy_symmetry(comm,cert_metric,cert_representation,cert_coefficients,&
    cert_eigenvalues,1,3,1d0,occupied_defect,target_defect,energy_defect,ok,message,&
    differential_worst_operation)
  call require(ok.and.differential_worst_operation==2,&
    'nontrivial complex differential evaluator fixture is invalid')
  call certify_dg_hybrid_energy_window(comm,cert_metric,cert_representation,cert_coefficients,cert_eigenvalues,&
    1,cert_eigenvalues(1),1.4d0,2,1d-12,1d0,0d0,0d0,certification,ok,message)
  call require(ok.and.certification%certified_rank==3.and.certification%worst_operation==differential_worst_operation.and.&
    abs(certification%occupied_subspace_defect-occupied_defect)<1d-12.and.&
    abs(certification%target_subspace_defect-target_defect)<1d-12.and.&
    abs(certification%target_energy_defect-energy_defect)<1d-12,&
    'prefix precomputation differs from the public evaluator for a complex metric fixture')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_SPECTRAL_CERTIFICATION ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS low-energy Hybrid symmetry on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine reset_certification_fixture()
    integer::j
    cert_metric=(0d0,0d0);cert_coefficients=(0d0,0d0);cert_representation=(0d0,0d0)
    do j=1,ncert
      cert_metric(j,j)=1d0;cert_coefficients(j,j)=1d0
      cert_representation(j,j,1)=1d0;cert_representation(j,j,2)=1d0
    enddo
    cert_eigenvalues=[-1d0,0d0,0.4d0,0.8d0,1.2d0,2d0]
  end subroutine reset_certification_fixture

  subroutine configure_nontrivial_certification_fixture()
    integer::j
    real(real64)::mix
    mix=1d0/sqrt(2d0)
    cert_metric=(0d0,0d0);cert_coefficients=(0d0,0d0);cert_action=(0d0,0d0)
    do j=1,ncert
      cert_metric(j,j)=real(j,real64);cert_action(j,j)=1d0
    enddo
    cert_coefficients(1,1)=mix;cert_coefficients(1,2)=cmplx(0d0,mix,kind=real64)
    cert_coefficients(2,1)=cmplx(0d0,mix/sqrt(2d0),kind=real64)
    cert_coefficients(2,2)=mix/sqrt(2d0)
    do j=3,ncert;cert_coefficients(j,j)=1d0/sqrt(real(j,real64));enddo
    cert_representation=(0d0,0d0)
    do j=1,ncert;cert_representation(j,j,1)=1d0;enddo
    cert_action(2,2)=0d0;cert_action(3,3)=0d0
    cert_action(2,3)=(0d0,1d0);cert_action(3,2)=(0d0,1d0)
    cert_representation(:,:,2)=matmul(cert_coefficients,&
      matmul(cert_action,matmul(conjg(transpose(cert_coefficients)),cert_metric)))
  end subroutine configure_nontrivial_certification_fixture

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_low_energy_symmetry_mpi
