program test_dg_dc_local_basis_scf_mpi
  use mpi
  use dg_dc_local_basis_ground_state
  implicit none
  type(s_dg_dc_local_basis_layout) :: layout
  complex(8) :: overlap(3,3),hamiltonian(3,3),coefficients(3,2)
  complex(8), allocatable :: h_rows(:,:),s_rows(:,:),coefficient_rows(:,:),distributed_coefficients(:,:)
  complex(8), allocatable :: volume_block(:,:),interface_rows(:,:),composed_rows(:,:)
  complex(8), allocatable :: local_basis(:,:),local_coefficients(:,:),psi(:,:)
  complex(8), allocatable :: raw_basis(:,:),orthonormal_basis(:,:),basis_transform(:,:),basis_gram(:,:)
  complex(8) :: full_raw_basis(5,2),full_transformed_basis(5,2)
  complex(8) :: projected_volume(2,2),volume_basis(3,2),volume_hbasis(3,2)
  real(8) :: eigenvalues(2),occupations(2),density(2),expected_density(2)
  real(8) :: projector_overlap,orthogonality_defect,hermiticity_defect
  complex(8) :: residual(3,2),metric(2,2)
  integer :: ierr,rank,nproc,first,last,i,effective_basis_count
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'

  allocate(raw_basis(3,rank+1),orthonormal_basis(3,rank+1),basis_transform(rank+1,rank+1))
  raw_basis(:,1)=[cmplx(1d0,0.2d0,8),cmplx(0.3d0,-0.1d0,8),cmplx(-0.2d0,0.4d0,8)]
  if(rank==1) raw_basis(:,2)=0.4d0*raw_basis(:,1)+ &
    [cmplx(0.1d0,-0.3d0,8),cmplx(0.7d0,0.2d0,8),cmplx(0.2d0,0.1d0,8)]
  call orthonormalize_dg_dc_fragment_core_basis(raw_basis,0.5d0,1d-10,orthonormal_basis, &
    basis_transform,effective_basis_count,ok,message)
  call require(ok .and. effective_basis_count==rank+1,'fragment core basis orthonormalization')
  allocate(basis_gram(rank+1,rank+1))
  basis_gram=0.5d0*matmul(conjg(transpose(orthonormal_basis(:,1:rank+1))), &
    orthonormal_basis(:,1:rank+1))
  call require(maxval(abs(basis_gram-identity(rank+1)))<1d-11,'fragment core metric is identity')
  call require(maxval(abs(orthonormal_basis-matmul(raw_basis,basis_transform)))<1d-12, &
    'fragment basis transform reconstructs orthonormal basis')
  if(rank==1) then
    raw_basis(:,2)=raw_basis(:,1)*(1d0+1d-14)
    call orthonormalize_dg_dc_fragment_core_basis(raw_basis,0.5d0,1d-10,orthonormal_basis, &
      basis_transform,effective_basis_count,ok,message)
    call require(ok .and. effective_basis_count==1,'fragment core metric rank reduction')
    call require(maxval(abs(basis_transform(:,2:2)))<1d-14 .and. &
      maxval(abs(orthonormal_basis(:,2:2)))<1d-14,'discarded fragment basis columns are zero')
  else
    call require(.true.,'fragment core metric rank reduction')
    call require(.true.,'discarded fragment basis columns are zero')
  end if
  raw_basis=cmplx(0.5d0*huge(1d0),0d0,8)
  call orthonormalize_dg_dc_fragment_core_basis(raw_basis,0.5d0,1d-10,orthonormal_basis, &
    basis_transform,effective_basis_count,ok,message)
  call require(.not.ok .and. effective_basis_count==0,'derived core metric overflow rejection')
  deallocate(raw_basis,orthonormal_basis,basis_transform,basis_gram)

  full_raw_basis=reshape([(cmplx(real(i,8),0d0,8),i=1,10)],[5,2])
  allocate(basis_transform(2,2))
  basis_transform=reshape([cmplx(1d0,0d0,8),cmplx(0.5d0,0d0,8), &
    cmplx(0d0,0d0,8),cmplx(2d0,0d0,8)],[2,2])
  call transform_dg_dc_fragment_buffer_basis(full_raw_basis,basis_transform,2, &
    full_transformed_basis,ok,message)
  call require(ok .and. maxval(abs(full_transformed_basis- &
    matmul(full_raw_basis,basis_transform)))<1d-14, &
    'core metric transform is applied to the full core plus buffer basis')
  deallocate(basis_transform)

  volume_basis=reshape([cmplx(1d0,0d0,8),cmplx(0d0,1d0,8),cmplx(0.5d0,0d0,8), &
    cmplx(0.2d0,-0.1d0,8),cmplx(0.7d0,0d0,8),cmplx(0d0,0.3d0,8)],[3,2])
  volume_hbasis=reshape([cmplx(2d0,0d0,8),cmplx(0d0,2d0,8),cmplx(1d0,0d0,8), &
    cmplx(0.4d0,-0.2d0,8),cmplx(1.4d0,0d0,8),cmplx(0d0,0.6d0,8)],[3,2])
  call project_dg_dc_local_basis_volume(volume_basis,volume_hbasis,0.25d0,projected_volume,ok,message)
  call require(ok .and. maxval(abs(projected_volume-0.25d0* &
    matmul(conjg(transpose(volume_basis)),volume_hbasis)))<1d-14, &
    'DC volume action projects into the fragment local basis')

  call initialize_dg_dc_local_basis_layout(layout,rank+1,rank+1,2,101_8,202_8, &
    MPI_COMM_WORLD,ok,message)
  call require(ok,'1+2 local basis layout')
  overlap=(0d0,0d0)
  do i=1,3
    overlap(i,i)=1d0
  end do
  hamiltonian=reshape([ &
    cmplx(1d0,0d0,8),cmplx(0.2d0,-0.1d0,8),(0d0,0d0), &
    cmplx(0.2d0,0.1d0,8),cmplx(2d0,0d0,8),cmplx(0.3d0,-0.05d0,8), &
    (0d0,0d0),cmplx(0.3d0,0.05d0,8),cmplx(4d0,0d0,8)],[3,3])
  call solve_dg_dc_local_basis_bands_reference(layout,overlap,hamiltonian,MPI_COMM_WORLD,1d-10, &
    eigenvalues,coefficients,ok,message)
  call require(ok,'global local-basis band solve')
  residual=matmul(hamiltonian,coefficients)- &
    matmul(overlap,coefficients)*spread(cmplx(eigenvalues,0d0,8),1,3)
  metric=matmul(conjg(transpose(coefficients)),matmul(overlap,coefficients))
  call require(maxval(abs(residual))<1d-11,'generalized eigen residual')
  call require(maxval(abs(metric-reshape([cmplx(1d0,0d0,8),(0d0,0d0),(0d0,0d0), &
    cmplx(1d0,0d0,8)],[2,2])))<1d-11,'metric orthonormality')
  call require(size(coefficients,2)==layout%global_band_count,'requested global bands retained')

  first=layout%basis_offsets(rank)+1
  last=layout%basis_offsets(rank+1)
  allocate(local_basis(2,layout%local_basis_count),local_coefficients(layout%local_basis_count,2),psi(2,2))
  do i=1,layout%local_basis_count
    local_basis(:,i)=[cmplx(0.4d0*(first+i),0.1d0*i,8),cmplx(-0.2d0*i,0.05d0*first,8)]
  end do
  local_coefficients=coefficients(first:last,:)
  occupations=[2d0,0d0]
  call reconstruct_dg_dc_local_basis_density(local_basis,local_coefficients,occupations,density,ok,message)
  call require(ok,'local coefficient density reconstruction')
  psi=matmul(local_basis,local_coefficients)
  expected_density=occupations(1)*abs(psi(:,1))**2+occupations(2)*abs(psi(:,2))**2
  call require(maxval(abs(density-expected_density))<1d-13,'density uses coefficients and occupations')
  call require(maxval(abs(coefficients(:,2)))>1d-8,'empty band coefficient retained')

  first=layout%basis_offsets(rank)+1
  last=layout%basis_offsets(rank+1)
  allocate(h_rows(layout%local_basis_count,3),s_rows(layout%local_basis_count,3))
  allocate(volume_block(layout%local_basis_count,layout%local_basis_count))
  allocate(interface_rows(layout%local_basis_count,3),composed_rows(layout%local_basis_count,3))
  allocate(coefficient_rows(layout%local_basis_count,2),distributed_coefficients(3,2))
  call initialize_dg_dc_local_basis_coefficients(layout,coefficient_rows,ok,message)
  call require(ok,'distributed coefficient initialization')
  call gather_coefficients(coefficient_rows,distributed_coefficients)
  call require(maxval(abs(distributed_coefficients(:,1)-[cmplx(1d0,0d0,8), &
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8)]))<1d-14 .and. &
    maxval(abs(distributed_coefficients(:,2)-[cmplx(0d0,0d0,8), &
    cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)]))<1d-14, &
    'initial coefficient rows retain occupied and empty global bands')
  call assign_dg_dc_local_basis_occupations(2d0,occupations,ok,message)
  call require(ok .and. maxval(abs(occupations-[2d0,0d0]))<1d-14, &
    'global occupations retain an empty band')
  call assign_dg_dc_local_basis_occupations(5d0,occupations,ok,message)
  call require(.not.ok,'insufficient global bands reject the electron count')
  occupations=[2d0,0d0]
  h_rows=hamiltonian(first:last,:)
  volume_block=hamiltonian(first:last,first:last)
  interface_rows=h_rows
  interface_rows(:,first:last)=(0d0,0d0)
  call compose_dg_dc_distributed_hamiltonian_rows(layout,volume_block,interface_rows,1d0, &
    composed_rows,ok,message)
  call require(ok .and. maxval(abs(composed_rows-h_rows))<1d-14, &
    'distributed rows compose DC volume and SIPG interface')
  s_rows=overlap(first:last,:)
  call solve_dg_dc_local_basis_bands_cg(layout,h_rows,s_rows,MPI_COMM_WORLD,80,1d-11, &
    coefficient_rows,eigenvalues,ok,message)
  call require(ok,'distributed coefficient CG solve')
  call diagnose_dg_dc_local_basis_continuation(layout,h_rows,coefficient_rows,coefficient_rows, &
    occupations,MPI_COMM_WORLD,projector_overlap,orthogonality_defect,hermiticity_defect,ok,message)
  call require(ok .and. projector_overlap>1d0-1d-12 .and. orthogonality_defect<1d-10 .and. &
    hermiticity_defect<1d-12,'continuation diagnostics accept the same occupied projector')
  call gather_coefficients(coefficient_rows,distributed_coefficients)
  residual=matmul(hamiltonian,distributed_coefficients)- &
    matmul(overlap,distributed_coefficients)*spread(cmplx(eigenvalues,0d0,8),1,3)
  metric=matmul(conjg(transpose(distributed_coefficients)),matmul(overlap,distributed_coefficients))
  call require(maxval(abs(residual))<1d-9,'distributed coefficient CG residual')
  call require(maxval(abs(metric-reshape([cmplx(1d0,0d0,8),(0d0,0d0),(0d0,0d0), &
    cmplx(1d0,0d0,8)],[2,2])))<1d-10,'distributed coefficient S orthonormality')
  s_rows(1,1)=s_rows(1,1)+1d-6
  call solve_dg_dc_local_basis_bands_cg(layout,h_rows,s_rows,MPI_COMM_WORLD,80,1d-11, &
    coefficient_rows,eigenvalues,ok,message)
  call require(.not.ok,'coefficient CG rejects non-orthonormal distributed metric')

  overlap=(0d0,0d0)
  overlap(1,1)=1d0
  overlap(2,2)=1d0
  overlap(3,3)=1d-14
  call solve_dg_dc_local_basis_bands_reference(layout,overlap,hamiltonian,MPI_COMM_WORLD,1d-10, &
    eigenvalues,coefficients,ok,message)
  call require(.not.ok,'near-singular overlap rejected collectively')

  overlap=(0d0,0d0)
  do i=1,3
    overlap(i,i)=1d0
  end do
  if(rank==1) hamiltonian(1,1)=hamiltonian(1,1)+1d-3
  call solve_dg_dc_local_basis_bands_reference(layout,overlap,hamiltonian,MPI_COMM_WORLD,1d-10, &
    eigenvalues,coefficients,ok,message)
  call require(.not.ok,'rank-disagreeing matrices rejected collectively')

  occupations=[2d0,0d0]
  density=[0.6d0,1.4d0]
  call validate_dg_dc_local_basis_density(occupations,2d0,2d0,density,[0.5d0,0.5d0], &
    MPI_COMM_WORLD,ok,message)
  call require(ok,'rank-consistent occupations and global charge')
  if(rank==1) occupations(1)=1.9d0
  call validate_dg_dc_local_basis_density(occupations,2d0,2d0,density,[0.5d0,0.5d0], &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreeing occupations rejected collectively')
  occupations=[2d0,0d0]
  call validate_dg_dc_local_basis_density(occupations,2d0,2.1d0,density,[0.5d0,0.5d0], &
    MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'electron-count mismatch rejected collectively')

  if(rank==0) print '(a)','PASS DG-DC local-basis global bands and density'
  call MPI_Finalize(ierr)

contains
  function identity(n) result(matrix)
    integer, intent(in) :: n
    complex(8) :: matrix(n,n)
    integer :: diagonal
    matrix=(0d0,0d0)
    do diagonal=1,n
      matrix(diagonal,diagonal)=1d0
    end do
  end function identity

  subroutine gather_coefficients(local,global)
    complex(8), intent(in) :: local(:,:)
    complex(8), intent(out) :: global(:,:)
    integer :: counts(2),displacements(2),iband,mpi_error
    counts=layout%basis_offsets(1:2)-layout%basis_offsets(0:1)
    displacements=layout%basis_offsets(0:1)
    do iband=1,size(local,2)
      call MPI_Allgatherv(local(:,iband),size(local,1),MPI_DOUBLE_COMPLEX,global(:,iband), &
        counts,displacements,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,mpi_error)
    end do
  end subroutine gather_coefficients

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: global
    integer :: mpi_error
    call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if(.not.global) then
      if(rank==0) write(*,'(a,1x,a)') 'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    end if
  end subroutine require
end program test_dg_dc_local_basis_scf_mpi
