program test_dg_ground_state_checkpoint_mpi
  use mpi
  use dg_dc_local_basis_ground_state, only: s_dg_dc_local_basis_production_state
  use dg_dc_ground_state, only: s_dg_dc_gs_result,s_dg_dc_gs_diagnostics,s_dg_dc_gs_controls, &
    default_dg_dc_gs_controls,DG_DC_ACCEPTED
  use dg_ground_state_checkpoint, only: s_dg_ground_state_checkpoint, &
    populate_dg_ground_state_checkpoint,publish_dg_ground_state_checkpoint, &
    read_dg_ground_state_checkpoint,validate_dg_ground_state_checkpoint
  implicit none
  type(s_dg_dc_local_basis_production_state) :: state
  type(s_dg_dc_gs_result) :: result
  type(s_dg_dc_gs_diagnostics) :: diagnostics
  type(s_dg_dc_gs_controls) :: controls
  type(s_dg_ground_state_checkpoint) :: checkpoint,loaded
  real(8), allocatable :: density(:,:,:,:),hartree(:,:,:),vxc(:,:,:,:),vlocal(:,:,:,:)
  integer :: ierr,rank,nproc,info,unit,ios,file_size,i
  integer :: local_nx
  integer(8) :: expected_geometry
  integer(1) :: corrupt_byte
  logical :: ok,reusable
  logical :: exists
  character(512) :: root,message,path

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call get_command_argument(1,root)
  local_nx=merge(1,3,rank==0)
  state%fragment_id=rank+1
  state%local_basis_count=2
  state%global_basis_count=4
  state%global_band_count=2
  state%raw_basis_count=2
  state%core_size=[local_nx,1,1]
  state%full_spatial_shape=[local_nx+2,3,3]
  state%scf_iterations=9
  state%geometry_fingerprint=101_8
  state%operator_fingerprint=202_8
  state%hamiltonian_operator_fingerprint=201_8
  state%density_residual=1d-9
  state%interface_scale=1d0
  allocate(state%coefficient_rows(2,2),state%full_fragment_basis((local_nx+2)*9,2), &
    state%basis_transform(2,2),state%eigenvalues(2),state%occupations(2), &
    state%basis_offsets(0:2),state%fragment_ids(2))
  state%coefficient_rows=cmplx(real(rank+1,8),0d0,8)
  state%full_fragment_basis=cmplx(real(rank+2,8),0d0,8)
  state%basis_transform=cmplx(real(rank+3,8),0d0,8)
  state%eigenvalues=[-1d0,1d0]
  state%occupations=[2d0,0d0]
  state%basis_offsets=[0,2,4]
  state%fragment_ids=[1,2]
  state%ready=.true.
  allocate(density(2,1,1,1),hartree(2,1,1),vxc(2,1,1,1),vlocal(2,1,1,1))
  density=real(rank+1,8); hartree=3d0+rank; vxc=4d0+rank; vlocal=5d0+rank
  allocate(result%lambda_history(3),result%lambda_steps(3))
  result%phase=DG_DC_ACCEPTED
  result%naccepted_steps=3
  result%nrollbacks=1
  result%mixing_reset_count=4
  result%total_scf_iterations=9
  result%lambda=1d0
  result%maximum_interface_scale=1d0
  result%minimum_projector_overlap=0.95d0
  result%accepted=.true.
  result%unmixed_verified=.true.
  result%lambda_history=[0d0,0.5d0,1d0]
  result%lambda_steps=[0d0,0.5d0,0.5d0]
  diagnostics%orbital_residual=1d-9
  diagnostics%density_residual=2d-9
  diagnostics%subspace_residual=3d-9
  diagnostics%projector_overlap=0.99d0
  diagnostics%hermiticity_defect=4d-12
  diagnostics%orthogonality_defect=5d-12
  diagnostics%face_balance_defect=6d-12
  diagnostics%electron_number=2d0
  diagnostics%expected_electron_number=2d0
  diagnostics%interface_scale=1d0
  diagnostics%eigensolver_iterations=7
  diagnostics%hamiltonian_potential_epoch=11_8
  diagnostics%updated_potential_epoch=12_8
  diagnostics%eigensolver_converged=.true.
  diagnostics%finite=.true.
  result%final_diagnostics=diagnostics
  controls=default_dg_dc_gs_controls()

  call populate_dg_ground_state_checkpoint(checkpoint,state,result,controls,density,hartree,vxc,vlocal, &
    [4,1,1],reshape([0,0,0,1,0,0],[3,2]),reshape([1,1,1,3,1,1],[3,2]), &
    [merge(1,0,rank==0),merge(1,0,rank==0),rank,rank,rank,rank], &
    'DG_DC_GS',ok,message)
  call require(ok,'populate complete accepted checkpoint')
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(ok,'validate complete accepted checkpoint')
  checkpoint%face_neighbors(1)=mod(checkpoint%face_neighbors(1)+1,nproc)
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'geometrically incorrect periodic face neighbor rejected')
  checkpoint%face_neighbors(1)=mod(checkpoint%state%fragment_id-2+nproc,nproc)
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  if(.not.ok)write(*,'(a,i0,2a)')'checkpoint publication rank ',rank,': ',trim(message)
  call require(ok,'transactional rank payload and manifest publication')
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. reusable,'collective round trip')
  call require(all(loaded%state%coefficient_rows==checkpoint%state%coefficient_rows) .and. &
    all(loaded%density==density) .and. all(loaded%hartree==hartree) .and. &
    all(loaded%vxc==vxc) .and. all(loaded%vlocal==vlocal),'payload round trip')
  call require(size(loaded%face_neighbors)==6,'six physical faces retained')
  call require(loaded%result%unmixed_verified .and. loaded%result%lambda==1d0 .and. &
    loaded%result%final_diagnostics%face_balance_defect==diagnostics%face_balance_defect, &
    'acceptance and every final diagnostic round trip')
  if(rank==0)checkpoint%fragment_origins(1,2)=2
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreed checkpoint payload rejected before publication')
  if(rank==0)checkpoint%fragment_origins(1,2)=1
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. reusable,'failed republication preserves previous manifest generation')
  if(rank==0)checkpoint%result%maximum_interface_scale=0.9d0
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreed result metadata rejected before publication')
  if(rank==0)checkpoint%result%maximum_interface_scale=1d0
  if(rank==0)checkpoint%state%eigenvalues(1)=checkpoint%state%eigenvalues(1)+1d0
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreed eigenvalues rejected before publication')
  if(rank==0)checkpoint%state%eigenvalues(1)=checkpoint%state%eigenvalues(1)-1d0
  if(rank==0)checkpoint%state%occupations(1)=1d0
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(.not.ok,'rank-disagreed occupations rejected before publication')
  if(rank==0)checkpoint%state%occupations(1)=2d0
  checkpoint%state%basis_offsets(1)=checkpoint%state%basis_offsets(1)+1000000
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'basis offsets must match global basis geometry')
  checkpoint%state%basis_offsets(1)=2
  checkpoint%state%fragment_ids(2)=1
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'fragment IDs must uniquely cover the distributed layout')
  checkpoint%state%fragment_ids(2)=2
  checkpoint%result%lambda_steps(2)=0.4d0
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'continuation step must equal the lambda-history increment')
  checkpoint%result%lambda_steps(2)=0.5d0
  checkpoint%result%minimum_projector_overlap=1d0
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'trajectory projector minimum cannot exceed the final projector')
  checkpoint%result%minimum_projector_overlap=0.95d0
  write(path,'(a,".g",i0,".rank",i0,".dg_gs.tmp")')trim(root),checkpoint%publication_generation,rank
  inquire(file=trim(path),exist=exists)
  call require(.not.exists,'transactional rank temporary is not published')
  if(rank==0)then
    write(path,'(a,".g",i0,".rank0.dg_gs")')trim(root),loaded%publication_generation
    call corrupt_payload_byte(trim(path),157)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(.not.ok .and. .not.reusable,'in-record metadata corruption rejected collectively')
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(ok,'republish after metadata corruption')
  if(rank==0)then
    write(path,'(a,".g",i0,".rank0.dg_gs")')trim(root),checkpoint%publication_generation
    inquire(file=trim(path),size=file_size)
    call corrupt_payload_byte(trim(path),file_size-7)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(.not.ok .and. .not.reusable,'in-record array corruption rejected collectively')
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(ok,'republish after array corruption')
  if(rank==0)then
    write(path,'(a,".g",i0,".rank0.dg_gs")')trim(root),checkpoint%publication_generation
    open(newunit=unit,file=trim(path),access='stream',form='unformatted',status='old', &
      position='append',action='write',iostat=ios)
    corrupt_byte=1
    if(ios==0)write(unit,iostat=ios)corrupt_byte
    if(ios==0)close(unit,iostat=ios)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(.not.ok .and. .not.reusable,'corrupt rank payload rejected collectively')
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(ok,'republish after corruption')
  expected_geometry=101_8
  if(rank==0)expected_geometry=999_8
  call read_dg_ground_state_checkpoint(trim(root),expected_geometry,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. .not.reusable,'rank-disagreed reuse rejected collectively')
  if(rank==0)then
    open(newunit=unit,file=trim(root)//'.manifest',status='old',iostat=ios)
    if(ios==0)close(unit,status='delete',iostat=ios)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(.not.ok .and. .not.reusable,'incomplete publication without manifest rejected')
  call publish_dg_ground_state_checkpoint(trim(root),checkpoint,MPI_COMM_WORLD,ok,message)
  call require(ok,'republish after incomplete publication')

  call read_dg_ground_state_checkpoint(trim(root),999_8,202_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. .not.reusable,'stale geometry rejected')
  call read_dg_ground_state_checkpoint(trim(root),101_8,999_8,'DG_DC_GS',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. .not.reusable,'stale operator rejected')
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'DG_WPW',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. .not.reusable,'WPW-mislabeled checkpoint rejected')
  call read_dg_ground_state_checkpoint(trim(root),101_8,202_8,'LCFO',MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  call require(ok .and. .not.reusable,'LCFO-mislabeled checkpoint rejected')
  checkpoint%result%accepted=.false.
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'unaccepted checkpoint rejected')
  checkpoint%result%accepted=.true.
  checkpoint%result%unmixed_verified=.false.
  call validate_dg_ground_state_checkpoint(checkpoint,101_8,202_8,'DG_DC_GS',ok,message)
  call require(.not.ok,'unverified fixed point rejected')

  if(rank==0) print '(a)','PASS verified DG ground-state checkpoint'
  call MPI_Finalize(ierr)
contains
  subroutine corrupt_payload_byte(filename,position)
    character(*),intent(in)::filename
    integer,intent(in)::position
    integer(1)::value
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='old', &
      action='readwrite',iostat=ios)
    if(ios==0)read(unit,pos=position,iostat=ios)value
    if(ios==0)then
      value=ieor(value,1_1)
      write(unit,pos=position,iostat=ios)value
    end if
    close(unit)
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    if(.not.condition)then
      write(*,'(a,i0,2a)')'FAIL rank ',rank,': ',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,info)
    end if
  end subroutine require
end program test_dg_ground_state_checkpoint_mpi
