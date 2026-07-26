program test_dg_dc_direct_checkpoint_mpi
  use mpi
  use dg_ground_state_checkpoint, only: s_dg_dc_direct_checkpoint_state, &
    publish_dg_dc_direct_checkpoint,read_dg_dc_direct_checkpoint
  implicit none
  type(s_dg_dc_direct_checkpoint_state) :: state,loaded
  integer :: ierr,rank,nproc,unit,ios,version,file_nproc
  logical :: ok,reusable
  character(512) :: root,message,path
  character(18) :: magic
  integer(1) :: byte
  integer(8) :: generation,manifest_geometry,manifest_operator
  integer(8),allocatable :: checksums(:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'fixture requires two ranks'
  call get_command_argument(1,root)
  state%ready=.true.
  state%accepted=.true.
  state%unmixed_verified=.true.
  state%fragment_id=rank+1
  state%nstate=3
  state%state_start=1
  state%state_count=3
  state%nspin=1
  state%scf_iterations=17
  state%core_size=[merge(1,3,rank==0),1,1]
  state%full_spatial_shape=state%core_size+2
  state%geometry_fingerprint=101_8
  state%operator_fingerprint=202_8
  state%hamiltonian_operator_fingerprint=202_8
  state%density_residual=1d-10
  state%interface_scale=1d0
  state%global_size=[4,1,1]
  state%fragment_origin=[merge(0,1,rank==0),0,0]
  state%fragment_size=state%core_size
  state%density_origin=state%fragment_origin
  state%density_size=state%core_size
  state%orbital_spatial_lower_bound=state%fragment_origin-1
  state%orbital_core_local_origin=state%fragment_origin
  state%orbital_core_origin=state%density_origin
  state%orbital_core_size=state%density_size
  state%face_neighbors=[mod(rank+1,nproc)+1,mod(rank-1+nproc,nproc)+1,rank+1,rank+1,rank+1,rank+1]
  state%continuation_rollbacks=1
  state%projector_generation=7
  state%projector_retained_rank=3
  state%projector_required_rank=3
  state%projector_fingerprint=303_8
  state%orbital_residual=dble(rank+1)*2d-10
  state%orthogonality_defect=3d-11
  state%hermiticity_defect=4d-12
  state%charge_error=5d-13
  state%face_balance_defect=6d-12
  state%density_tolerance=1d-8
  state%orbital_tolerance=1d-8
  state%orthogonality_tolerance=1d-8
  state%hermiticity_tolerance=1d-8
  state%charge_tolerance=1d-8
  state%face_balance_tolerance=1d-8
  state%projector_projection_residual=dble(rank+1)*1d-12
  state%projector_escape_norm=dble(rank+1)*2d-12
  state%projector_residual_limit=1d-8
  state%projector_escape_limit=1d-8
  state%face_action_norm=[1d0,1d0,2d0,2d0,3d0,3d0]
  state%face_pair_balance=0d0
  allocate(state%fragment_orbitals(state%full_spatial_shape(1),state%full_spatial_shape(2), &
    state%full_spatial_shape(3),1,3,1,1),state%occupations(3,1,1))
  allocate(state%continuation_scale_history(3),state%continuation_step_history(3))
  allocate(state%density(state%density_size(1),1,1,1),state%hartree(state%density_size(1),1,1), &
    state%vxc(state%density_size(1),1,1,1),state%vlocal(state%density_size(1),1,1,1))
  state%fragment_orbitals=dble(rank+1)
  state%fragment_orbitals(:,:,:,:,2,:,:)=2d0*dble(rank+1)
  state%fragment_orbitals(:,:,:,:,3,:,:)=3d0*dble(rank+1)
  state%occupations=reshape([2d0,2d0,0d0],[3,1,1])
  state%continuation_scale_history=[0d0,0.5d0,1d0]
  state%continuation_step_history=[0d0,0.5d0,0.5d0]
  state%density=dble(rank+1)
  state%hartree=2d0*dble(rank+1)
  state%vxc=3d0*dble(rank+1)
  state%vlocal=4d0*dble(rank+1)
  call publish_dg_dc_direct_checkpoint(trim(root),state,MPI_COMM_WORLD,ok,message)
  if(.not.ok) error stop trim(message)
  call read_dg_dc_direct_checkpoint(trim(root),101_8,202_8,MPI_COMM_WORLD,loaded,reusable,ok,message)
  if(.not.ok .or. .not.reusable) error stop trim(message)
  if(any(shape(loaded%fragment_orbitals)/=shape(state%fragment_orbitals))) error stop 'orbital shape changed'
  if(maxval(abs(loaded%fragment_orbitals-state%fragment_orbitals))>0d0) error stop 'orbitals changed'
  if(any(loaded%occupations/=state%occupations)) error stop 'occupations changed'
  if(any(loaded%continuation_scale_history/=state%continuation_scale_history)) error stop 'history changed'
  if(any(loaded%density/=state%density) .or. any(loaded%hartree/=state%hartree) .or. &
     any(loaded%vxc/=state%vxc) .or. any(loaded%vlocal/=state%vlocal)) error stop 'fields changed'
  if(any(loaded%face_neighbors/=state%face_neighbors)) error stop 'topology changed'
  if(loaded%projector_generation/=state%projector_generation .or. &
     loaded%projector_fingerprint/=state%projector_fingerprint) error stop 'projector identity changed'

  open(newunit=unit,file=trim(root)//'.manifest',access='stream',form='unformatted', &
    status='old',action='read',iostat=ios)
  if(ios/=0)error stop 'cannot open direct manifest'
  read(unit,iostat=ios)magic,version,generation,file_nproc,manifest_geometry,manifest_operator
  if(ios/=0.or.file_nproc/=nproc)error stop 'cannot read direct manifest'
  allocate(checksums(nproc))
  read(unit,iostat=ios)checksums
  close(unit)
  write(path,'(a,".g",i0,".rank",i0,".dg_direct")')trim(root),generation,rank
  if(rank==0)then
    open(newunit=unit,file=trim(path),access='stream',form='unformatted',status='old',action='readwrite',iostat=ios)
    if(ios/=0) error stop 'cannot open payload for corruption'
    read(unit,pos=1,iostat=ios)byte
    byte=ieor(byte,int(1,1))
    write(unit,pos=1,iostat=ios)byte
    close(unit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call read_dg_dc_direct_checkpoint(trim(root),101_8,202_8,MPI_COMM_WORLD,loaded,reusable,ok,message)
  if(.not.ok .or. reusable) error stop 'corrupt checkpoint was reusable'

  state%geometry_fingerprint=101_8+rank
  call publish_dg_dc_direct_checkpoint(trim(root)//'_disagree',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'rank metadata disagreement was published'
  state%geometry_fingerprint=101_8
  state%scf_iterations=17+rank
  call publish_dg_dc_direct_checkpoint(trim(root)//'_scf_disagree',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'common continuation metadata disagreement was published'
  state%scf_iterations=17
  if(rank==0)state%face_neighbors(1)=999
  call publish_dg_dc_direct_checkpoint(trim(root)//'_bad_topology',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'invalid direct topology was published'
  state%face_neighbors=[mod(rank+1,nproc)+1,mod(rank-1+nproc,nproc)+1,rank+1,rank+1,rank+1,rank+1]
  state%orbital_residual=2d0*state%orbital_tolerance
  call publish_dg_dc_direct_checkpoint(trim(root)//'_bad_gate',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'unaccepted direct diagnostics were published'
  state%orbital_residual=2d-10*dble(rank+1)

  state%fragment_id=1
  state%nstate=4
  state%state_start=2*rank+1
  state%state_count=2
  state%core_size=[4,1,1]
  state%full_spatial_shape=[6,3,3]
  state%orbital_spatial_lower_bound=[-1,0,0]
  state%orbital_core_local_origin=[0,0,0]
  state%orbital_core_origin=[0,0,0]
  state%orbital_core_size=[4,1,1]
  state%global_size=[4,1,1]
  state%fragment_origin=[0,0,0]
  state%fragment_size=state%core_size
  state%density_origin=[0,0,0]
  state%density_size=[4,1,1]
  state%face_neighbors=1
  deallocate(state%fragment_orbitals,state%occupations,state%density,state%hartree,state%vxc,state%vlocal)
  allocate(state%fragment_orbitals(6,3,3,1,2,1,1),state%occupations(2,1,1))
  allocate(state%density(4,1,1,1),state%hartree(4,1,1),state%vxc(4,1,1,1),state%vlocal(4,1,1,1))
  state%fragment_orbitals=dble(rank+1)
  state%occupations=reshape(merge([2d0,0d0],[2d0,2d0],rank==1),[2,1,1])
  state%density=1d0
  state%hartree=2d0
  state%vxc=3d0
  state%vlocal=4d0
  call publish_dg_dc_direct_checkpoint(trim(root)//'_state_split',state,MPI_COMM_WORLD,ok,message)
  if(.not.ok) error stop 'state-axis distributed direct checkpoint was rejected'
  call read_dg_dc_direct_checkpoint(trim(root)//'_state_split',101_8,202_8,MPI_COMM_WORLD, &
    loaded,reusable,ok,message)
  if(.not.ok .or. .not.reusable) error stop 'state-axis distributed direct checkpoint was not reusable'
  if(loaded%state_start/=state%state_start .or. loaded%state_count/=state%state_count) &
    error stop 'state-axis ownership changed'
  if(any(loaded%orbital_spatial_lower_bound/=state%orbital_spatial_lower_bound)) &
    error stop 'orbital spatial lower bounds changed'
  if(any(loaded%orbital_core_origin/=state%orbital_core_origin) .or. &
     any(loaded%orbital_core_local_origin/=state%orbital_core_local_origin) .or. &
     any(loaded%orbital_core_size/=state%orbital_core_size)) error stop 'orbital owned core changed'
  if(rank==1)state%state_start=1
  call publish_dg_dc_direct_checkpoint(trim(root)//'_duplicate_state_shard',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'duplicate state shard and matching state hole were published'
  state%state_start=2*rank+1
  if(rank==1)state%orbital_spatial_lower_bound(1)=0
  call publish_dg_dc_direct_checkpoint(trim(root)//'_inconsistent_spatial_shards',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'state intervals with inconsistent spatial shard sets were published'
  state%orbital_spatial_lower_bound=[-1,0,0]
  if(rank==1)state%density=9d0
  call publish_dg_dc_direct_checkpoint(trim(root)//'_inconsistent_density_replica',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'inconsistent replicated density payload was published'
  state%density=1d0
  state%orbital_spatial_lower_bound=[2,0,0]
  call publish_dg_dc_direct_checkpoint(trim(root)//'_missing_orbital_core',state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'orbital storage that misses the owned core was published'
  state%nstate=2
  state%state_start=1
  state%state_count=2
  state%orbital_spatial_lower_bound=[rank-1,0,0]
  state%orbital_core_local_origin=[0,0,0]
  state%orbital_core_origin=[0,0,0]
  state%orbital_core_size=[4,1,1]
  state%occupations=1d0
  call publish_dg_dc_direct_checkpoint(trim(root)//'_duplicate_orbital_core_owner', &
    state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'duplicate nonzero orbital core/state ownership was published'
  state%nstate=2
  state%state_start=1
  state%state_count=2
  state%full_spatial_shape=[4,3,3]
  state%orbital_spatial_lower_bound=[2*rank-1,0,0]
  state%orbital_core_local_origin=[2*rank,0,0]
  state%orbital_core_origin=[2*rank,0,0]
  state%orbital_core_size=[2,1,1]
  state%density_origin=state%orbital_core_origin
  state%density_size=state%orbital_core_size
  deallocate(state%fragment_orbitals,state%density,state%hartree,state%vxc,state%vlocal)
  allocate(state%fragment_orbitals(4,3,3,1,2,1,1),state%density(2,1,1,1), &
    state%hartree(2,1,1),state%vxc(2,1,1,1),state%vlocal(2,1,1,1))
  state%fragment_orbitals=dble(rank+1)
  state%density=dble(rank+1)
  state%hartree=2d0*dble(rank+1)
  state%vxc=3d0*dble(rank+1)
  state%vlocal=4d0*dble(rank+1)
  state%occupations=dble(rank+1)
  call publish_dg_dc_direct_checkpoint(trim(root)//'_inconsistent_occupation_replica', &
    state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'inconsistent replicated occupation payload was published'
  state%nstate=4
  state%state_start=2*rank+1
  state%full_spatial_shape=[6,3,3]
  state%orbital_spatial_lower_bound=[-1,0,0]
  deallocate(state%fragment_orbitals)
  allocate(state%fragment_orbitals(6,3,3,1,2,1,1))
  state%fragment_orbitals=dble(rank+1)
  call publish_dg_dc_direct_checkpoint(trim(root)//'_missing_state_core_cross_product', &
    state,MPI_COMM_WORLD,ok,message)
  if(ok) error stop 'state intervals missing part of the fragment core were published'
  if(rank==0)print '(a)','PASS direct DC checkpoint MPI fixture'
  call MPI_Finalize(ierr)
end program test_dg_dc_direct_checkpoint_mpi
