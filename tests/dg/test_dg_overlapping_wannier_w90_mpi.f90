program test_dg_overlapping_wannier_w90_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_w90,only:estimate_dg_w90_coordinator_bytes,&
    validate_dg_w90_result,setup_dg_w90_gamma_library,run_dg_w90_gamma_library,&
    assemble_dg_w90_gamma_matrices,apply_dg_w90_gamma_transform,&
    validate_dg_w90_convergence_log
  use dg_overlapping_wannier_w90,only:validate_dg_w90_generator_covariance
  use dg_overlapping_wannier_w90,only:align_dg_w90_character_sector_gauge
  use dg_overlapping_wannier_w90,only:validate_dg_w90_localization_cluster
  use dg_overlapping_wannier_w90,only:align_dg_w90_cross_character_sector_gauge
  use dg_overlapping_wannier_w90,only:align_dg_w90_character_sectors_by_periodic_phase
  use dg_overlapping_wannier_w90,only:sew_dg_w90_periodic_phase_conjugate_sector
  use dg_overlapping_wannier_w90,only:anchor_dg_w90_reference_character_sector
  use dg_overlapping_wannier_w90,only:project_dg_w90_reference_sector_operators
  use dg_overlapping_wannier_w90,only:inherit_dg_w90_affine_receipts
  use dg_overlapping_wannier_w90,only:build_dg_sector_periodic_position_tuple
  use dg_overlapping_wannier_w90,only:canonicalize_dg_sector_periodic_position_gauge
  use dg_overlapping_wannier_w90,only:jointly_canonicalize_dg_sector_periodic_position_gauge
  use dg_overlapping_wannier_w90,only:build_dg_orbital_major_periodic_position_tuple,&
    apply_dg_orbital_rotation_tiled
  use dg_overlapping_wannier_w90,only:export_dg_w90_replay_bundle
  use dg_overlapping_wannier_w90,only:convert_dg_w90_library_geometry
  implicit none
  integer::ierr,rank,nproc,b,i,m,n,p,nlocal
  integer::convergence_iterations,log_unit,win_unit,win_io
  complex(8)::transform(2,2)
  real(8)::centers(3,2),spreads(2),spread(3)
  integer(8)::bytes
  logical::ok,matrix_matches,win_has_random_projection,replay_exists
  character(256)::message,win_line
  complex(8),allocatable::local_values(:,:),local_anchors(:,:)
  complex(8),allocatable::assembled_m(:,:,:),assembled_a(:,:)
  real(8),allocatable::local_weights(:),local_fractional(:,:)
  integer::test_nncell(3,2),global_point
  integer(8)::matrix_peak,matrix_estimate
  complex(8)::local_m_reference(2,2,2),local_a_reference(2,2),m_reference(2,2,2),a_reference(2,2),phase
  real(8)::angle
  real(8)::test_atomic_lattice(3,3),test_atomic_reciprocal(3,3),test_atomic_atoms(3,1),&
    library_lattice_units(3,3),library_reciprocal_units(3,3),library_atom_units(3,1)
  real(8)::inherited_identity,inherited_unitarity,inherited_closure
  complex(8)::covariance_band(2,2),covariance_wann(2,2),covariance_transform(2,2)
  real(8)::covariance_defect
  integer(8)::covariance_workspace
  complex(8),allocatable::gauge_values(:,:),gauge_gradients(:,:,:)
  complex(8)::gauge_transform(2,2)
  real(8)::gauge_centers(3,2)
  real(8)::gauge_spreads(2)
  integer(8),allocatable::gauge_ids(:)
  integer(8),allocatable::sector_ids(:),anchor_duplicate_ids(:)
  complex(8),allocatable::sector_frame(:,:),sector_reference(:,:),sector_gamma(:,:),&
    sector_conjugate(:,:),sector_aligned(:,:),sector_conjugate_aligned(:,:)
  complex(8),allocatable::sector_rotated(:,:),sector_reference_permuted(:,:),sector_trial_aligned(:,:),&
    sector_trial_conjugate(:,:)
  complex(8),allocatable::sector_gamma_trial(:,:)
  complex(8),allocatable::cross_reference_sector(:,:),cross_target_sector(:,:),cross_w90_reference(:,:),&
    cross_aligned_target(:,:),cross_rotated_target(:,:)
  complex(8),allocatable::cross_periodic_phase(:)
  complex(8),allocatable::anchor_sector(:,:),anchor_w90(:,:),anchor_lcfo(:,:),anchor_result(:,:),anchor_trial(:,:)
  complex(8)::anchor_w90_operator(2,2),anchor_lcfo_operator(2,2),anchor_rotation(2,2)
  real(8)::anchor_defect
  integer(8)::anchor_fingerprint,anchor_trial_fingerprint,anchor_workspace
  complex(8),allocatable::anchor_projected_w90(:,:),anchor_projected_lcfo(:,:)
  complex(8),allocatable::cross_gamma_rows(:,:),cross_conjugate_sector(:,:),cross_aligned_conjugate(:,:)
  complex(8),allocatable::implicit_gamma_rows(:,:)
  complex(8),allocatable::weighted_spatial_sector(:,:)
  real(8),allocatable::spatial_integration_weights(:)
  real(8)::cross_localization_weights(4)
  real(8),allocatable::cross_singular_values(:)
  real(8)::cross_polar_defect
  real(8)::cross_gamma_defect
  integer(8)::cross_fingerprint,cross_workspace
  integer(8)::cross_phase_payload_fingerprint,cross_phase_bits
  integer(8)::cross_gamma_fingerprint
  integer(8)::implicit_gamma_fingerprint
  real(8),allocatable::sector_singular_values(:)
  real(8)::sector_polar_defect,sector_gamma_defect
  real(8)::sector_local_cost,sector_global_cost,sector_swapped_local_cost,sector_swapped_global_cost
  real(8)::localization_cluster_spectrum(4)
  integer(8)::sector_alignment_fingerprint,sector_alignment_workspace
  complex(8),allocatable::sector_position_phases(:,:),sector_position_tuple(:,:,:)
  complex(8),allocatable::orbital_major_values(:,:),orbital_major_tuple(:,:,:)
  complex(8)::sector_position_reference(2,2,3),sector_position_local(2,2,3)
  integer(8)::sector_position_fingerprint,sector_position_workspace
  real(8)::sector_position_gram_defect
  complex(8),allocatable::position_canonical_rows(:,:),position_trial_rows(:,:),position_trial_tuple(:,:,:),&
    position_rotated_sector(:,:),position_weighted_sector(:,:)
  real(8),allocatable::position_weights(:)
  complex(8)::position_lcfo_operator(2,2),position_trial_operator(2,2),position_input_rotation(2,2),&
    position_canonical_rotation(2,2),position_trial_rotation(2,2)
  real(8)::position_canonical_defect,position_trial_defect
  integer(8)::position_canonical_fingerprint,position_trial_canonical_fingerprint,&
    position_canonical_workspace,position_trial_position_fingerprint
  complex(8),allocatable::joint_canonical_rows(:,:),joint_trial_rows(:,:)
  complex(8)::joint_rotation(2,2),joint_trial_rotation(2,2)
  complex(8)::joint_point_representations(2,2,2)
  real(8),allocatable::joint_centers(:,:),joint_trial_centers(:,:)
  real(8)::joint_objective,joint_trial_objective,joint_update,joint_trial_update,joint_defect,joint_trial_defect
  integer::joint_sweeps,joint_trial_sweeps
  integer(8)::joint_fingerprint,joint_trial_fingerprint,joint_workspace
  integer(8)::sector_trial_fingerprint,sector_trial_workspace
  complex(8),allocatable::replay_m(:,:,:),replay_a(:,:)
  real(8),allocatable::replay_eigenvalues(:)
  integer::replay_gvec(3,2),replay_unit,replay_ios
  integer(8),allocatable::sector_reference_keys(:),sector_permuted_keys(:)
#ifdef USE_WANNIER90
  integer::nntot
  integer,allocatable::nncell(:,:)
  complex(8),allocatable::m_matrix(:,:,:),a_matrix(:,:),library_transform(:,:)
  real(8),allocatable::library_centers(:,:),library_spreads(:)
  real(8)::lattice(3,3),reciprocal(3,3),atoms_cart(3,1),library_spread(3),eigenvalues(1)
  character(2)::atom_symbols(1)
#endif
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0
  centers=reshape([0.1d0,0.2d0,0.3d0,0.6d0,0.2d0,0.3d0],[3,2])
  spreads=[0.4d0,0.5d0];spread=[0.9d0,0.2d0,0.7d0]
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(ok,trim(message))
  call inherit_dg_w90_affine_receipts(transform,3d-13,1d-12,inherited_identity,&
    inherited_unitarity,inherited_closure,bytes,ok,message)
  call require(ok.and.inherited_identity<=3d-13.and.inherited_unitarity<1d-12.and.&
    inherited_closure<=3d-13.and.bytes>0_8,trim(message))
  covariance_band=reshape([(0d0,0d0),(1d0,0d0),(1d0,0d0),(0d0,0d0)],[2,2])
  covariance_wann=covariance_band;covariance_transform=covariance_band
  call validate_dg_w90_generator_covariance(MPI_COMM_WORLD,covariance_transform,&
    covariance_band,covariance_wann,1d-12,covariance_defect,covariance_workspace,ok,message)
  call require(ok.and.covariance_defect<1d-12.and.covariance_workspace>0_8,&
    'dense internal Wannier representation satisfies measured covariance')
  covariance_transform=(0d0,0d0);covariance_transform(1,1)=1d0;covariance_transform(2,2)=-1d0
  call validate_dg_w90_generator_covariance(MPI_COMM_WORLD,covariance_transform,&
    covariance_band,covariance_wann,1d-12,covariance_defect,covariance_workspace,ok,message)
  call require(.not.ok.and.covariance_defect>1d0,&
    'noncommuting Wannier transform fails measured covariance')
  call estimate_dg_w90_coordinator_bytes(384,384,12,1,bytes,ok,message)
  call require(ok.and.bytes>0_8,'finite Si64 Wannier90 byte estimate')
  call estimate_dg_w90_coordinator_bytes(huge(0),huge(0),12,1,bytes,ok,message)
  call require(.not.ok,'Wannier90 byte estimate rejects integer overflow')
  test_atomic_lattice=0d0;test_atomic_reciprocal=0d0;test_atomic_atoms=0d0
  test_atomic_lattice(1,1)=10d0;test_atomic_lattice(2,2)=10d0;test_atomic_lattice(3,3)=10d0
  test_atomic_reciprocal(1,1)=2d0*acos(-1d0)/10d0
  test_atomic_reciprocal(2,2)=test_atomic_reciprocal(1,1)
  test_atomic_reciprocal(3,3)=test_atomic_reciprocal(1,1)
  test_atomic_atoms(:,1)=[1d0,2d0,3d0]
  call convert_dg_w90_library_geometry(test_atomic_lattice,test_atomic_reciprocal,test_atomic_atoms,&
    library_lattice_units,library_reciprocal_units,library_atom_units,ok,message)
  call require(ok.and.abs(library_lattice_units(1,1)-5.2917721067d0)<1d-12.and.&
    abs(library_reciprocal_units(1,1)-test_atomic_reciprocal(1,1)/0.52917721067d0)<1d-12.and.&
    maxval(abs(library_atom_units(:,1)-[0.52917721067d0,1.05835442134d0,1.58753163201d0]))<1d-12,&
    'Wannier90 library geometry converts SALMON atomic units to Angstrom')
  replay_gvec=reshape([1,0,0,0,1,0],[3,2])
  if(rank==0)then
    allocate(replay_m(2,2,2),replay_a(2,2),replay_eigenvalues(2))
    replay_eigenvalues=[-1d0,2d0]
    replay_a=reshape([cmplx(11d0,-11d0,8),cmplx(21d0,-21d0,8),&
      cmplx(12d0,-12d0,8),cmplx(22d0,-22d0,8)],[2,2])
    replay_m(:,:,1)=replay_a;replay_m(:,:,2)=2d0*replay_a
    open(newunit=replay_unit,file='replay_source.win',status='replace',iostat=replay_ios)
    if(replay_ios==0)write(replay_unit,'(a)')'num_wann = 2'
    if(replay_ios==0)close(replay_unit)
    open(newunit=replay_unit,file='replay_source.dmn',status='replace',iostat=replay_ios)
    if(replay_ios==0)write(replay_unit,'(a)')'replay dmn fixture'
    if(replay_ios==0)close(replay_unit)
  else
    allocate(replay_m(0,0,0),replay_a(0,0),replay_eigenvalues(0))
  endif
  call export_dg_w90_replay_bundle(MPI_COMM_WORLD,'.','replay_source','.',&
    'replay_bundle',replay_eigenvalues,replay_a,replay_m,replay_gvec,ok,message)
  call require(ok,trim(message))
  replay_exists=.false.
  if(rank==0)inquire(file='replay_bundle.win',exist=replay_exists)
  call MPI_Bcast(replay_exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call require(replay_exists,'replay bundle contains .win')
  if(rank==0)then
    open(newunit=log_unit,file='w90_converged_fixture.wout',status='replace')
    write(log_unit,'(a)')'      7  -0.100E-13  0.0  1.0  0.0 <-- CONV'
    write(log_unit,'(a)')'             <<< Wannierisation convergence criteria satisfied >>>'
    write(log_unit,'(a)')' Final State';close(log_unit)
    call validate_dg_w90_convergence_log('w90_converged_fixture.wout',200,&
      convergence_iterations,ok,message)
  endif
  call MPI_Bcast(ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(convergence_iterations,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call require(ok.and.convergence_iterations==7,trim(message))
  if(rank==0)then
    open(newunit=log_unit,file='w90_exhausted_fixture.wout',status='replace')
    write(log_unit,'(a)')'    200  -0.100E-02  0.1  1.0  0.0 <-- CONV'
    write(log_unit,'(a)')' Final State';close(log_unit)
    call validate_dg_w90_convergence_log('w90_exhausted_fixture.wout',200,&
      convergence_iterations,ok,message)
  endif
  call MPI_Bcast(ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call require(.not.ok,'Wannier90 iteration-limit exhaustion rejected')
  nlocal=8/nproc
  allocate(local_values(2,nlocal),local_anchors(2,nlocal),local_weights(nlocal),local_fractional(3,nlocal))
  do p=1,nlocal
    global_point=rank*nlocal+p
    local_values(:,p)=[cmplx(global_point,0d0,8),cmplx((-1d0)**global_point,0d0,8)]
    local_anchors(:,p)=[cmplx(1d0,0d0,8),cmplx(global_point,0d0,8)]
    local_fractional(:,p)=[real(global_point-1,8)/8d0,0.25d0*mod(global_point-1,4),0d0]
  enddo
  local_weights=0.125d0;test_nncell=reshape([1,0,0,0,1,0],[3,2])
  local_m_reference=(0d0,0d0);local_a_reference=(0d0,0d0)
  do b=1,2;do p=1,nlocal
    angle=-2d0*acos(-1d0)*dot_product(real(test_nncell(:,b),8),local_fractional(:,p))
    phase=cmplx(cos(angle),sin(angle),8)
    do n=1,2;do m=1,2
      local_m_reference(m,n,b)=local_m_reference(m,n,b)+&
        local_weights(p)*conjg(local_values(m,p))*phase*local_values(n,p)
    enddo;enddo
  enddo;enddo
  do p=1,nlocal;do n=1,2;do m=1,2
    local_a_reference(m,n)=local_a_reference(m,n)+&
      local_weights(p)*conjg(local_values(m,p))*local_anchors(n,p)
  enddo;enddo;enddo
  call MPI_Reduce(local_m_reference,m_reference,size(m_reference),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
    MPI_COMM_WORLD,ierr)
  call MPI_Reduce(local_a_reference,a_reference,size(a_reference),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
    MPI_COMM_WORLD,ierr)
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,huge(0_8),assembled_m,assembled_a,matrix_estimate,&
    matrix_peak,ok,message)
  call require(ok.and.matrix_estimate>0_8.and.matrix_peak<=matrix_estimate,trim(message))
  if(rank==0)then
    call require(all(shape(assembled_m)==[2,2,2]).and.all(shape(assembled_a)==[2,2]),&
      'coordinator owns complete Wannier90 M/A matrices')
    matrix_matches=maxval(abs(assembled_m-m_reference))<1d-12.and.&
      maxval(abs(assembled_a-a_reference))<1d-12
    write(*,'(a,1x,es24.16)')'W90_MATRIX_FINGERPRINT',sum(abs(assembled_m))+sum(abs(assembled_a))
  else
    call require(size(assembled_m)==0.and.size(assembled_a)==0,&
      'noncoordinator does not own complete Wannier90 M/A matrices')
    matrix_matches=.true.
  endif
  call require(matrix_matches,'distributed Wannier90 M/A match dense references')
  allocate(sector_ids(nlocal),sector_frame(nlocal,2),sector_reference(nlocal,2),sector_gamma(nlocal,8),&
    sector_conjugate(nlocal,2))
  do p=1,nlocal
    global_point=rank*nlocal+p
    sector_ids(p)=global_point
    sector_frame(p,1)=exp(cmplx(0d0,2d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    sector_frame(p,2)=exp(cmplx(0d0,6d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
  enddo
  allocate(sector_position_phases(nlocal,3))
  sector_position_local=(0d0,0d0)
  do p=1,nlocal
    global_point=rank*nlocal+p
    do b=1,3
      angle=2d0*acos(-1d0)*real(b*global_point,8)/17d0
      sector_position_phases(p,b)=cmplx(cos(angle),sin(angle),8)
      do i=1,2;do m=1,2
        sector_position_local(i,m,b)=sector_position_local(i,m,b)+&
          conjg(sector_frame(p,i))*sector_position_phases(p,b)*sector_frame(p,m)
      enddo;enddo
    enddo
  enddo
  call MPI_Allreduce(sector_position_local,sector_position_reference,size(sector_position_reference),&
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call build_dg_sector_periodic_position_tuple(MPI_COMM_WORLD,sector_ids,8,sector_frame,&
    sector_position_phases,1d-12,7781_8,sector_position_tuple,sector_position_gram_defect,&
    sector_position_fingerprint,sector_position_workspace,ok,message)
  call require(ok.and.maxval(abs(sector_position_tuple-sector_position_reference))<1d-12.and.&
    sector_position_gram_defect<1d-12.and.sector_position_fingerprint/=0_8.and.&
    sector_position_workspace>0_8,'distributed sector periodic-position tuple matches direct sum')
  allocate(orbital_major_values(2,nlocal),source=transpose(sector_frame))
  call build_dg_orbital_major_periodic_position_tuple(MPI_COMM_WORLD,orbital_major_values,&
    [(1d0,i=1,nlocal)],transpose(sector_position_phases),1d-12,7781_8,orbital_major_tuple,&
    sector_position_gram_defect,position_trial_position_fingerprint,sector_position_workspace,ok,message)
  call require(ok.and.maxval(abs(orbital_major_tuple-sector_position_reference))<1d-12.and.&
    sector_position_gram_defect<1d-12.and.sector_position_workspace>0_8,&
    'orbital-major periodic-position tuple matches row-major direct sum without a transpose copy')
  position_input_rotation=reshape([cmplx(1d0,0d0,8),cmplx(0d0,1d0,8),&
    cmplx(0d0,1d0,8),cmplx(1d0,0d0,8)],[2,2])/sqrt(2d0)
  call apply_dg_orbital_rotation_tiled(MPI_COMM_WORLD,orbital_major_values,position_input_rotation,ok,message)
  call require(ok.and.maxval(abs(orbital_major_values-transpose(matmul(sector_frame,position_input_rotation))))<1d-12,&
    'tiled orbital-major rotation matches the row-major channel rotation')
  deallocate(orbital_major_values,orbital_major_tuple)
  allocate(position_weighted_sector(nlocal,2),position_weights(nlocal))
  position_weights=2d0;position_weighted_sector=sector_frame/sqrt(2d0)
  sector_position_local=(0d0,0d0)
  do p=1,nlocal;do b=1,3;do i=1,2;do m=1,2
    sector_position_local(i,m,b)=sector_position_local(i,m,b)+position_weights(p)*&
      conjg(position_weighted_sector(p,i))*sector_position_phases(p,b)*position_weighted_sector(p,m)
  enddo;enddo;enddo;enddo
  call MPI_Allreduce(sector_position_local,sector_position_reference,size(sector_position_reference),&
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call build_dg_sector_periodic_position_tuple(MPI_COMM_WORLD,sector_ids,8,position_weighted_sector,&
    sector_position_phases,1d-12,7781_8,position_trial_tuple,sector_position_gram_defect,&
    position_trial_position_fingerprint,sector_position_workspace,ok,message,position_weights)
  call require(ok.and.maxval(abs(position_trial_tuple-sector_position_reference))<1d-12.and.&
    sector_position_gram_defect<1d-12,'sector periodic-position tuple uses the spatial integration measure')
  if(rank==0)write(*,'(a,1x,i0)')'W90_POSITION_TUPLE_FINGERPRINT',sector_position_fingerprint
  position_lcfo_operator=(0d0,0d0);position_lcfo_operator(1,1)=0.2d0;position_lcfo_operator(2,2)=0.7d0
  position_input_rotation=reshape([cmplx(1d0,0d0,8),cmplx(0d0,1d0,8),&
    cmplx(0d0,1d0,8),cmplx(1d0,0d0,8)],[2,2])/sqrt(2d0)
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,sector_position_fingerprint,&
    position_canonical_rows,position_canonical_rotation,position_canonical_defect,&
    position_canonical_fingerprint,position_canonical_workspace,ok,message)
  call require(ok.and.position_canonical_defect<1d-10.and.position_canonical_workspace>0_8,&
    'periodic-position tuple defines a canonical internal frame')
  allocate(position_rotated_sector(nlocal,2))
  position_rotated_sector=matmul(sector_frame,position_input_rotation)
  call build_dg_sector_periodic_position_tuple(MPI_COMM_WORLD,sector_ids,8,position_rotated_sector,&
    sector_position_phases,1d-12,7781_8,position_trial_tuple,sector_position_gram_defect,&
    position_trial_position_fingerprint,sector_position_workspace,ok,message)
  position_trial_operator=matmul(conjg(transpose(position_input_rotation)),&
    matmul(position_lcfo_operator,position_input_rotation))
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,position_trial_position_fingerprint,&
    position_trial_rows,position_trial_rotation,position_trial_defect,&
    position_trial_canonical_fingerprint,position_canonical_workspace,ok,message)
  call require(ok.and.position_trial_canonical_fingerprint==position_canonical_fingerprint.and.&
    maxval(abs(position_trial_rows-position_canonical_rows))<1d-10,&
    'canonical periodic-position gauge is invariant under input sector rotation')
  sector_position_tuple=(0d0,0d0)
  position_trial_tuple=(0d0,0d0)
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,811_8,&
    position_canonical_rows,position_canonical_rotation,position_canonical_defect,&
    position_canonical_fingerprint,position_canonical_workspace,ok,message)
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,823_8,&
    position_trial_rows,position_trial_rotation,position_trial_defect,&
    position_trial_canonical_fingerprint,position_canonical_workspace,ok,message)
  call require(ok.and.maxval(abs(position_trial_rows-position_canonical_rows))<1d-10.and.&
    position_trial_canonical_fingerprint==position_canonical_fingerprint,&
    'LCFO operator resolves a periodic-position-degenerate internal block')
  position_lcfo_operator=(0d0,0d0);position_trial_operator=(0d0,0d0)
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,827_8,&
    position_canonical_rows,position_canonical_rotation,position_canonical_defect,&
    position_canonical_fingerprint,position_canonical_workspace,ok,message)
  call canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,829_8,&
    position_trial_rows,position_trial_rotation,position_trial_defect,&
    position_trial_canonical_fingerprint,position_canonical_workspace,ok,message)
  call require(ok.and.position_trial_canonical_fingerprint==position_canonical_fingerprint,&
    'an exact internal multiplet preserves one common projector fingerprint')
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,1,1)=exp(cmplx(0d0,0.3d0,8))
  sector_position_tuple(2,2,1)=exp(cmplx(0d0,1.1d0,8))
  sector_position_tuple(1,1,2)=exp(cmplx(0d0,0.5d0,8))
  sector_position_tuple(2,2,2)=exp(cmplx(0d0,1.4d0,8))
  sector_position_tuple(1,1,3)=exp(cmplx(0d0,0.7d0,8))
  sector_position_tuple(2,2,3)=exp(cmplx(0d0,1.8d0,8))
  joint_point_representations=(0d0,0d0)
  joint_point_representations(1,1,1)=1d0;joint_point_representations(2,2,1)=1d0
  joint_point_representations(1,2,2)=1d0;joint_point_representations(2,1,2)=1d0
  do b=1,3
    position_trial_tuple(:,:,b)=matmul(conjg(transpose(position_input_rotation)),&
      matmul(sector_position_tuple(:,:,b),position_input_rotation))
  enddo
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,901_8,joint_canonical_rows,joint_rotation,&
    joint_centers,joint_objective,joint_update,joint_sweeps,joint_defect,joint_fingerprint,joint_workspace,ok,message,&
    point_representations=joint_point_representations)
  call require(ok.and.joint_objective<1d-20.and.joint_defect<1d-10,'joint center gauge resolves distinct centers')
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,897_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message,&
    point_representations=joint_point_representations(:,:,1:1))
  call require(ok.and.joint_trial_defect<1d-10,&
    'multiple seeds complete an internal space for the trivial point group')
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,1,1)=exp(cmplx(0d0,0.4d0*acos(-1d0),8))
  sector_position_tuple(2,2,1)=exp(cmplx(0d0,1.4d0*acos(-1d0),8))
  sector_position_tuple(1,1,2)=1d0;sector_position_tuple(2,2,2)=1d0
  sector_position_tuple(1,1,3)=1d0;sector_position_tuple(2,2,3)=1d0
  joint_point_representations(:,:,2)=reshape([1d0,1d0,1d0,-1d0],[2,2])/sqrt(2d0)
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,899_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message,point_representations=joint_point_representations)
  call require(ok.and.joint_trial_defect<1d-10,&
    'point-orbit blocks close a noncommuting projected-position tuple')
  joint_point_representations=(0d0,0d0)
  joint_point_representations(1,1,1)=1d0;joint_point_representations(2,2,1)=1d0
  joint_point_representations(1,2,2)=1d0;joint_point_representations(2,1,2)=1d0
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,1,1)=exp(cmplx(0d0,0.3d0,8))
  sector_position_tuple(2,2,1)=exp(cmplx(0d0,1.1d0,8))
  sector_position_tuple(1,1,2)=exp(cmplx(0d0,0.5d0,8))
  sector_position_tuple(2,2,2)=exp(cmplx(0d0,1.4d0,8))
  sector_position_tuple(1,1,3)=exp(cmplx(0d0,0.7d0,8))
  sector_position_tuple(2,2,3)=exp(cmplx(0d0,1.8d0,8))
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,1,1)=0.8d0;sector_position_tuple(2,2,1)=-0.2d0
  sector_position_tuple(1,2,1)=cmplx(0.3d0,0.4d0,8)
  sector_position_tuple(2,1,1)=conjg(sector_position_tuple(1,2,1))
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,9001_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok.and.joint_trial_objective<1d-20,&
    'complex Jacobi rotation monotonically diagonalizes one Hermitian component')
  ! Projected periodic-position components need not commute exactly.  For an
  ! isotropic Pauli triple the joint objective is stationary under every
  ! two-column rotation, so convergence must be based on objective progress
  ! rather than requiring the (non-unique) Jacobi angle to vanish.
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,2,1)=1d0;sector_position_tuple(2,1,1)=1d0
  sector_position_tuple(1,2,2)=cmplx(0d0,-1d0,8)
  sector_position_tuple(2,1,2)=cmplx(0d0,1d0,8)
  sector_position_tuple(1,1,3)=1d0;sector_position_tuple(2,2,3)=-1d0
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,902_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok.and.joint_trial_objective<huge(1d0).and.joint_trial_sweeps<100,&
    'stationary noncommuting periodic-position tuple terminates deterministically')
  sector_position_tuple=(0d0,0d0)
  sector_position_tuple(1,1,1)=exp(cmplx(0d0,0.3d0,8))
  sector_position_tuple(2,2,1)=exp(cmplx(0d0,1.1d0,8))
  sector_position_tuple(1,1,2)=exp(cmplx(0d0,0.5d0,8))
  sector_position_tuple(2,2,2)=exp(cmplx(0d0,1.4d0,8))
  sector_position_tuple(1,1,3)=exp(cmplx(0d0,0.7d0,8))
  sector_position_tuple(2,2,3)=exp(cmplx(0d0,1.8d0,8))
  joint_point_representations(1,2,2)=2d0
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,903_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message,point_representations=joint_point_representations)
  call require(.not.ok,'joint center gauge rejects a nonunitary point representation')
  joint_point_representations(1,2,2)=1d0
  if(nproc>1)then
    if(rank==0)joint_point_representations(1,2,2)=-1d0
    call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
      sector_position_tuple,position_lcfo_operator,1d-12,905_8,joint_trial_rows,joint_trial_rotation,&
      joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
      joint_trial_fingerprint,joint_workspace,ok,message,point_representations=joint_point_representations)
    call require(.not.ok,'joint center gauge rejects rank-disagreeing point representations')
    if(rank==0)joint_point_representations(1,2,2)=1d0
  endif
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,907_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok.and.maxval(abs(joint_trial_rows-joint_canonical_rows))<1d-10.and.&
    maxval(abs(joint_trial_centers-joint_centers))<1d-10.and.joint_trial_fingerprint==joint_fingerprint,&
    'joint center gauge is invariant under an input-sector unitary rotation')
  position_trial_tuple(:,:,1)=sector_position_tuple(:,:,2)
  position_trial_tuple(:,:,2)=sector_position_tuple(:,:,1)
  position_trial_tuple(:,:,3)=sector_position_tuple(:,:,3)
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    position_trial_tuple,position_lcfo_operator,1d-12,909_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok.and.maxval(abs(joint_trial_rows-joint_canonical_rows))<1d-10.and.&
    maxval(abs(joint_trial_centers(1,:)-joint_centers(2,:)))<1d-10.and.&
    maxval(abs(joint_trial_centers(2,:)-joint_centers(1,:)))<1d-10.and.&
    maxval(abs(joint_trial_centers(3,:)-joint_centers(3,:)))<1d-10.and.&
    joint_trial_fingerprint==joint_fingerprint,&
    'joint center gauge projector receipt is invariant under a Cartesian point rotation')
  do b=1,3
    sector_position_tuple(:,:,b)=(0d0,0d0)
    sector_position_tuple(1,1,b)=exp(cmplx(0d0,0.4d0*real(b,8),8))
    sector_position_tuple(2,2,b)=sector_position_tuple(1,1,b)
    position_trial_tuple(:,:,b)=matmul(conjg(transpose(position_input_rotation)),&
      matmul(sector_position_tuple(:,:,b),position_input_rotation))
  enddo
  position_lcfo_operator=(0d0,0d0);position_lcfo_operator(1,1)=0.2d0;position_lcfo_operator(2,2)=0.7d0
  position_trial_operator=matmul(conjg(transpose(position_input_rotation)),&
    matmul(position_lcfo_operator,position_input_rotation))
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,911_8,joint_canonical_rows,joint_rotation,&
    joint_centers,joint_objective,joint_update,joint_sweeps,joint_defect,joint_fingerprint,joint_workspace,ok,message)
  call require(ok,trim(message))
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,919_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok,trim(message))
  call require(ok.and.maxval(abs(joint_trial_rows-joint_canonical_rows))<1d-10.and.&
    joint_trial_fingerprint==joint_fingerprint,&
    'LCFO operator resolves a repeated-center block in the joint gauge')
  position_lcfo_operator=(0d0,0d0);position_trial_operator=(0d0,0d0)
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_position_tuple,position_lcfo_operator,1d-12,923_8,joint_canonical_rows,joint_rotation,&
    joint_centers,joint_objective,joint_update,joint_sweeps,joint_defect,joint_fingerprint,joint_workspace,ok,message)
  call require(ok,trim(message))
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,position_rotated_sector,&
    position_trial_tuple,position_trial_operator,1d-12,929_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(ok.and.maxval(abs(matmul(joint_trial_rows,conjg(transpose(joint_trial_rows)))-&
    matmul(joint_canonical_rows,conjg(transpose(joint_canonical_rows)))))<1d-10.and.&
    joint_trial_fingerprint==joint_fingerprint,&
    'an exact repeated-center multiplet preserves its projector receipt')
  if(rank==0)write(*,'(a,1x,i0)')'W90_JOINT_CENTER_FINGERPRINT',joint_fingerprint
  if(nlocal>1)then
    allocate(anchor_duplicate_ids,source=sector_ids)
    if(rank==0)anchor_duplicate_ids(1)=anchor_duplicate_ids(2)
    call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,anchor_duplicate_ids,sector_frame,&
      sector_position_tuple,position_lcfo_operator,1d-12,937_8,joint_trial_rows,joint_trial_rotation,&
      joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
      joint_trial_fingerprint,joint_workspace,ok,message)
    call require(.not.ok,'joint center gauge rejects duplicate row ownership')
    deallocate(anchor_duplicate_ids)
  endif
  position_trial_tuple=sector_position_tuple
  if(rank==0)position_trial_tuple(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    position_trial_tuple,position_lcfo_operator,1d-12,941_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(.not.ok,'joint center gauge rejects a nonfinite tuple collectively')
  position_trial_tuple=sector_position_tuple
  position_trial_tuple(1,1,1)=cmplx(0.25d0*huge(1d0),0d0,8)
  call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    position_trial_tuple,position_lcfo_operator,1d-12,943_8,joint_trial_rows,joint_trial_rotation,&
    joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
    joint_trial_fingerprint,joint_workspace,ok,message)
  call require(.not.ok,'joint center gauge rejects finite-huge tuple entries before squaring')
  if(nproc>1)then
    call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
      sector_position_tuple,position_lcfo_operator,merge(1d-11,1d-12,rank==0),947_8,&
      joint_trial_rows,joint_trial_rotation,joint_trial_centers,joint_trial_objective,joint_trial_update,&
      joint_trial_sweeps,joint_trial_defect,joint_trial_fingerprint,joint_workspace,ok,message)
    call require(.not.ok,'joint center gauge rejects rank-disagreeing tolerance')
    position_trial_tuple=sector_position_tuple
    if(rank==0)position_trial_tuple(1,1,1)=position_trial_tuple(1,1,1)+cmplx(1d-4,0d0,8)
    call jointly_canonicalize_dg_sector_periodic_position_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
      position_trial_tuple,position_lcfo_operator,1d-12,953_8,joint_trial_rows,joint_trial_rotation,&
      joint_trial_centers,joint_trial_objective,joint_trial_update,joint_trial_sweeps,joint_trial_defect,&
      joint_trial_fingerprint,joint_workspace,ok,message)
    call require(.not.ok,'joint center gauge rejects rank-disagreeing tuple payload')
  endif
  sector_reference(:,1)=(sector_frame(:,1)+cmplx(0.3d0,0.4d0,8)*sector_frame(:,2))/sqrt(1.25d0)
  sector_reference(:,2)=(-cmplx(0.3d0,-0.4d0,8)*sector_frame(:,1)+sector_frame(:,2))/sqrt(1.25d0)
  sector_gamma=(0d0,0d0)
  do p=1,nlocal;sector_gamma(p,rank*nlocal+p)=1d0;enddo
  sector_conjugate=conjg(sector_frame)
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,sector_reference,&
    sector_gamma,sector_conjugate,0d0,1d-12,sector_aligned,sector_conjugate_aligned,&
    sector_singular_values,sector_reference_keys,sector_polar_defect,sector_gamma_defect,sector_alignment_fingerprint,&
    sector_alignment_workspace,ok,message)
  call require(ok.and.minval(sector_singular_values)>0.9d0.and.sector_polar_defect<1d-10.and.&
    sector_gamma_defect<1d-10.and.sector_alignment_workspace>0_8.and.&
    maxval(abs(sector_conjugate_aligned-conjg(sector_aligned)))<1d-10,trim(message))
  if(rank==0)write(*,'(a,1x,i0)')'W90_SECTOR_FINGERPRINT',sector_alignment_fingerprint
  if(nproc>1)then
    call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,sector_reference,&
      sector_gamma,sector_conjugate,0d0,merge(1d-11,1d-12,rank==0),sector_trial_aligned,&
      sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
      sector_trial_fingerprint,sector_trial_workspace,ok,message)
    call require(.not.ok.and.index(trim(message),'metadata')>0,&
      'Wannier90 alignment rejects rank-inconsistent canonical metadata')
  endif
  allocate(sector_gamma_trial(nlocal,8));sector_gamma_trial=sector_gamma
  do p=1,nlocal
    global_point=rank*nlocal+p
    do n=1,8
      sector_gamma_trial(p,n)=sector_gamma_trial(p,n)+0.1d0*&
        exp(cmplx(0d0,4d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)*&
        (exp(cmplx(0d0,2d0*acos(-1d0)*real(n-1,8)/8d0,8))/sqrt(8d0)+&
        cmplx(0.3d0,0.4d0,8)*exp(cmplx(0d0,6d0*acos(-1d0)*real(n-1,8)/8d0,8))/sqrt(8d0))/sqrt(1.25d0)
    enddo
  enddo
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,sector_reference,&
    sector_gamma_trial,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(.not.ok.and.index(trim(message),'Gamma')>0,&
    'Wannier90 alignment rejects Gamma leakage outside the conjugate sector')
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,sector_reference,&
    sector_gamma,sector_conjugate,1d-4,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(.not.ok,'Wannier90 alignment rejects an unaccepted Gamma sewing receipt')
  sector_conjugate=2d0*sector_conjugate
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,sector_reference,&
    sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(.not.ok.and.index(trim(message),'conjugate sector frame')>0,&
    'Wannier90 alignment rejects a nonorthonormal conjugate sector')
  sector_conjugate=0.5d0*sector_conjugate
  allocate(sector_rotated(nlocal,2),sector_reference_permuted(nlocal,2))
  sector_rotated(:,1)=(sector_frame(:,1)+cmplx(0d0,1d0,8)*sector_frame(:,2))/sqrt(2d0)
  sector_rotated(:,2)=(cmplx(0d0,1d0,8)*sector_frame(:,1)+sector_frame(:,2))/sqrt(2d0)
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_rotated,sector_reference,&
    sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,sector_trial_conjugate,&
    sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,sector_trial_fingerprint,&
    sector_trial_workspace,ok,message)
  call require(ok.and.sector_trial_fingerprint==sector_alignment_fingerprint.and.&
    maxval(abs(sector_trial_aligned-sector_aligned))<1d-10,&
    'Wannier90 sector alignment is invariant under input sector-frame rotation')
  sector_reference_permuted(:,1)=exp(cmplx(0d0,0.7d0,8))*sector_reference(:,2)
  sector_reference_permuted(:,2)=exp(cmplx(0d0,-0.4d0,8))*sector_reference(:,1)
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_reference_permuted,sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(ok.and.sector_trial_fingerprint==sector_alignment_fingerprint.and.&
    maxval(abs(sector_trial_aligned-sector_aligned))<1d-10.and.&
    maxval(abs(sector_trial_conjugate-sector_conjugate_aligned))<1d-10,&
    'Wannier90 sector alignment is invariant under reference phases and Wannier numbering')
  do p=1,nlocal
    global_point=rank*nlocal+p
    sector_reference_permuted(p,2)=exp(cmplx(0d0,4d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
  enddo
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_reference_permuted,sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(.not.ok.and.index(trim(message),'singular')>0,&
    'Wannier90 sector alignment rejects a singular localization link')
  do p=1,nlocal
    global_point=rank*nlocal+p
    sector_reference_permuted(p,1)=1.5d-12*sector_frame(p,1)+sqrt(1d0-(1.5d-12)**2)*&
      exp(cmplx(0d0,4d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    sector_reference_permuted(p,2)=0.5d-12*sector_frame(p,2)+sqrt(1d0-(0.5d-12)**2)*&
      exp(cmplx(0d0,10d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
  enddo
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_reference_permuted,sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  call require(.not.ok.and.index(trim(message),'splits')>0,&
    'Wannier90 alignment rejects a threshold that splits a singular-value cluster')
  do p=1,nlocal
    global_point=rank*nlocal+p
    sector_reference_permuted(p,1)=cos(0.4d0)*sector_frame(p,1)+sin(0.4d0)*&
      exp(cmplx(0d0,4d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    sector_reference_permuted(p,2)=sector_frame(p,2)
  enddo
  call align_dg_w90_character_sector_gauge(MPI_COMM_WORLD,sector_ids,sector_frame,&
    sector_reference_permuted,sector_gamma,sector_conjugate,0d0,1d-12,sector_trial_aligned,&
    sector_trial_conjugate,sector_singular_values,sector_permuted_keys,sector_polar_defect,sector_gamma_defect,&
    sector_trial_fingerprint,sector_trial_workspace,ok,message)
  sector_local_cost=sum(abs(sector_trial_aligned-sector_reference_permuted)**2)
  sector_swapped_local_cost=sum(abs(sector_trial_aligned(:,1)-sector_reference_permuted(:,2))**2)+&
    sum(abs(sector_trial_aligned(:,2)-sector_reference_permuted(:,1))**2)
  call MPI_Allreduce(sector_local_cost,sector_global_cost,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(sector_swapped_local_cost,sector_swapped_global_cost,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,ierr)
  sector_global_cost=min(sector_global_cost,sector_swapped_global_cost)
  call require(ok.and.abs(sector_singular_values(1)-1d0)<1d-10.and.&
    abs(sector_singular_values(2)-cos(0.4d0))<1d-10.and.&
    abs(sector_global_cost-(4d0-2d0*sum(sector_singular_values)))<1d-10,&
    'Wannier90 SVD returns the closest polar frame for a nonunitary full-rank link')
  allocate(cross_reference_sector(nlocal,2),cross_target_sector(nlocal,2),cross_w90_reference(nlocal,4),&
    cross_rotated_target(nlocal,2))
  do p=1,nlocal
    global_point=rank*nlocal+p
    cross_reference_sector(p,1)=exp(cmplx(0d0,2d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    cross_reference_sector(p,2)=exp(cmplx(0d0,6d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    cross_target_sector(p,1)=exp(cmplx(0d0,4d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
    cross_target_sector(p,2)=exp(cmplx(0d0,8d0*acos(-1d0)*real(global_point-1,8)/8d0,8))/sqrt(8d0)
  enddo
  cross_w90_reference(:,1)=(cross_reference_sector(:,1)+exp(cmplx(0d0,0.4d0,8))*&
    cross_target_sector(:,1))/sqrt(2d0)
  cross_w90_reference(:,2)=(cross_reference_sector(:,2)+exp(cmplx(0d0,-0.7d0,8))*&
    cross_target_sector(:,2))/sqrt(2d0)
  cross_w90_reference(:,3)=(cross_reference_sector(:,1)-exp(cmplx(0d0,0.4d0,8))*&
    cross_target_sector(:,1))/sqrt(2d0)
  cross_w90_reference(:,4)=(cross_reference_sector(:,2)-exp(cmplx(0d0,-0.7d0,8))*&
    cross_target_sector(:,2))/sqrt(2d0)
  cross_localization_weights=[1d0,2d0,-1d0,-2d0]
  cross_rotated_target(:,1)=(cross_target_sector(:,1)+cmplx(0d0,1d0,8)*cross_target_sector(:,2))/sqrt(2d0)
  cross_rotated_target(:,2)=(cmplx(0d0,1d0,8)*cross_target_sector(:,1)+cross_target_sector(:,2))/sqrt(2d0)
  allocate(cross_periodic_phase(nlocal))
  do p=1,nlocal
    global_point=rank*nlocal+p
    cross_periodic_phase(p)=exp(cmplx(0d0,2d0*acos(-1d0)*real(global_point-1,8)/8d0,8))
  enddo
  cross_phase_payload_fingerprint=int(z'243F6A8885A308D3',8)
  do global_point=1,8
    phase=exp(cmplx(0d0,2d0*acos(-1d0)*real(global_point-1,8)/8d0,8))
    cross_phase_bits=transfer(real(phase,8),cross_phase_bits)
    cross_phase_payload_fingerprint=ieor(ishftc(cross_phase_payload_fingerprint,11),cross_phase_bits)
    cross_phase_bits=transfer(aimag(phase),cross_phase_bits)
    cross_phase_payload_fingerprint=ieor(ishftc(cross_phase_payload_fingerprint,11),cross_phase_bits)
  enddo
  if(cross_phase_payload_fingerprint==0_8)cross_phase_payload_fingerprint=1_8
  call align_dg_w90_character_sectors_by_periodic_phase(MPI_COMM_WORLD,sector_ids,&
    cross_reference_sector,cross_target_sector,cross_periodic_phase,8,777_8,cross_phase_payload_fingerprint,1d-12,&
    cross_aligned_target,cross_singular_values,cross_polar_defect,cross_fingerprint,&
    cross_workspace,ok,message)
  call require(ok.and.maxval(abs(cross_aligned_target(:,1)-cross_target_sector(:,1)))<1d-10.and.&
    maxval(abs(cross_aligned_target(:,2)-cross_target_sector(:,2)))<1d-10,&
    'periodic spatial phase maps the reference character sector into the target sector')
  call align_dg_w90_character_sectors_by_periodic_phase(MPI_COMM_WORLD,sector_ids,&
    cross_reference_sector,cross_rotated_target,cross_periodic_phase,8,777_8,cross_phase_payload_fingerprint,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(ok.and.maxval(abs(sector_trial_aligned-cross_aligned_target))<1d-10,&
    'periodic-phase alignment is invariant under target-sector unitary rotations')
  allocate(weighted_spatial_sector(nlocal,2),spatial_integration_weights(nlocal))
  weighted_spatial_sector=2d0*cross_reference_sector
  spatial_integration_weights=0.25d0
  call align_dg_w90_character_sectors_by_periodic_phase(MPI_COMM_WORLD,sector_ids,&
    weighted_spatial_sector,2d0*cross_target_sector,cross_periodic_phase,8,777_8,cross_phase_payload_fingerprint,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message,spatial_integration_weights)
  call require(ok.and.maxval(abs(sector_trial_aligned-2d0*cross_target_sector))<1d-10,&
    'periodic-phase alignment uses the real-space integration measure')
  deallocate(weighted_spatial_sector,spatial_integration_weights)
  call require(sector_trial_fingerprint==cross_fingerprint,&
    'periodic-phase projector fingerprint is invariant under target-sector rotations')
  allocate(anchor_sector(nlocal,2),anchor_w90(nlocal,4),anchor_lcfo(nlocal,4))
  anchor_sector=(0d0,0d0);anchor_w90=(0d0,0d0);anchor_lcfo=(0d0,0d0)
  do p=1,nlocal
    global_point=rank*nlocal+p
    if(global_point<=4)then
      anchor_w90(p,global_point)=1d0;anchor_lcfo(p,global_point)=1d0
      if(global_point<=2)anchor_sector(p,global_point)=1d0
    endif
  enddo
  anchor_w90_operator=(0d0,0d0);anchor_lcfo_operator=(0d0,0d0);anchor_rotation=(0d0,0d0)
  anchor_w90_operator(1,1)=1d0;anchor_w90_operator(2,2)=1d0
  anchor_lcfo_operator(1,1)=1d0;anchor_lcfo_operator(2,2)=2d0
  anchor_rotation=reshape([cmplx(sqrt(0.5d0),0d0,8),cmplx(0d0,sqrt(0.5d0),8),&
    cmplx(0d0,sqrt(0.5d0),8),cmplx(sqrt(0.5d0),0d0,8)],[2,2])
  call project_dg_w90_reference_sector_operators(MPI_COMM_WORLD,sector_ids,anchor_sector,anchor_w90,&
    [1d0,1d0,3d0,4d0],anchor_lcfo,[1d0,2d0,3d0,4d0],8,321_8,654_8,0d0,0d0,1d-12,&
    anchor_projected_w90,anchor_projected_lcfo,anchor_defect,anchor_fingerprint,anchor_workspace,ok,message)
  call require(ok.and.maxval(abs(anchor_projected_w90-anchor_w90_operator))<1d-10.and.&
    maxval(abs(anchor_projected_lcfo-anchor_lcfo_operator))<1d-10,&
    'full W90 and LCFO physical operators stream-project into the reference sector')
  call anchor_dg_w90_reference_character_sector(MPI_COMM_WORLD,sector_ids,anchor_sector,&
    anchor_w90_operator,anchor_lcfo_operator,8,321_8,654_8,0d0,0d0,1d-12,anchor_result,&
    anchor_defect,anchor_fingerprint,anchor_workspace,ok,message)
  call require(ok.and.anchor_defect<1d-10,'W90 reference sector resolves repeated-center channels with LCFO anchors')
  anchor_w90(:,1)=(anchor_lcfo(:,1)+cmplx(0d0,1d0,8)*anchor_lcfo(:,2))/sqrt(2d0)
  anchor_w90(:,2)=(cmplx(0d0,1d0,8)*anchor_lcfo(:,1)+anchor_lcfo(:,2))/sqrt(2d0)
  call anchor_dg_w90_reference_character_sector(MPI_COMM_WORLD,sector_ids,matmul(anchor_sector,anchor_rotation),&
    matmul(conjg(transpose(anchor_rotation)),matmul(anchor_w90_operator,anchor_rotation)),&
    matmul(conjg(transpose(anchor_rotation)),matmul(anchor_lcfo_operator,anchor_rotation)),&
    8,321_8,654_8,0d0,0d0,1d-12,anchor_trial,&
    anchor_defect,anchor_trial_fingerprint,anchor_workspace,ok,message)
  call require(ok.and.anchor_trial_fingerprint==anchor_fingerprint.and.&
    maxval(abs(matmul(anchor_trial,conjg(transpose(anchor_trial)))-&
      matmul(anchor_result,conjg(transpose(anchor_result)))))<1d-10,&
    'reference-sector anchor is invariant under W90 rotations inside a repeated-center block')
  anchor_w90_operator=(0d0,0d0);anchor_lcfo_operator=(0d0,0d0)
  anchor_w90_operator(1,1)=1d0;anchor_w90_operator(2,2)=1d0
  anchor_lcfo_operator(1,1)=2d0;anchor_lcfo_operator(2,2)=2d0
  call anchor_dg_w90_reference_character_sector(MPI_COMM_WORLD,sector_ids,anchor_sector,&
    anchor_w90_operator,anchor_lcfo_operator,8,321_8,654_8,0d0,0d0,1d-12,anchor_result,&
    anchor_defect,anchor_fingerprint,anchor_workspace,ok,message)
  call require(ok.and.anchor_defect<1d-10,&
    'symmetry-enforced occupied multiplet remains a complete unresolved internal block')
  call anchor_dg_w90_reference_character_sector(MPI_COMM_WORLD,sector_ids,matmul(anchor_sector,anchor_rotation),&
    anchor_w90_operator,anchor_lcfo_operator,8,321_8,654_8,0d0,0d0,1d-12,anchor_trial,&
    anchor_defect,anchor_trial_fingerprint,anchor_workspace,ok,message)
  call require(ok.and.anchor_trial_fingerprint==anchor_fingerprint.and.&
    maxval(abs(matmul(anchor_trial,conjg(transpose(anchor_trial)))-&
      matmul(anchor_result,conjg(transpose(anchor_result)))))<1d-10,&
    'unresolved symmetry multiplet receipt is internal-gauge invariant')
  if(nlocal>1)then
    allocate(anchor_duplicate_ids,source=sector_ids)
    if(rank==0)anchor_duplicate_ids(1)=anchor_duplicate_ids(2)
    call anchor_dg_w90_reference_character_sector(MPI_COMM_WORLD,anchor_duplicate_ids,anchor_sector,&
      anchor_w90_operator,anchor_lcfo_operator,8,321_8,654_8,0d0,0d0,1d-12,anchor_trial,&
      anchor_defect,anchor_trial_fingerprint,anchor_workspace,ok,message)
    call require(.not.ok,'reference-sector anchor rejects same-rank duplicate row ownership')
    deallocate(anchor_duplicate_ids)
  endif
  allocate(cross_gamma_rows(nlocal,8),cross_conjugate_sector(nlocal,2))
  cross_gamma_rows=(0d0,0d0)
  do p=1,nlocal;cross_gamma_rows(p,rank*nlocal+p)=1d0;enddo
  cross_gamma_fingerprint=int(z'6A09E667F3BCC909',8)
  do global_point=1,8
    do p=1,8
      phase=merge((1d0,0d0),(0d0,0d0),p==global_point)
      cross_phase_bits=transfer(real(phase,8),cross_phase_bits)
      cross_gamma_fingerprint=ieor(ishftc(cross_gamma_fingerprint,9),cross_phase_bits)
      cross_phase_bits=transfer(aimag(phase),cross_phase_bits)
      cross_gamma_fingerprint=ieor(ishftc(cross_gamma_fingerprint,9),cross_phase_bits)
    enddo
  enddo
  if(cross_gamma_fingerprint==0_8)cross_gamma_fingerprint=1_8
  cross_conjugate_sector=conjg(cross_target_sector)
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,cross_aligned_target,&
    cross_gamma_rows,cross_conjugate_sector,8,cross_gamma_fingerprint,.false.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,&
    cross_workspace,ok,message)
  call require(ok.and.cross_gamma_defect<1d-10.and.&
    maxval(abs(cross_aligned_conjugate-conjg(cross_aligned_target)))<1d-10,&
    'periodic-phase target and conjugate sectors share one Gamma sewing gauge')
  cross_conjugate_sector=2d0*cross_conjugate_sector
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,cross_aligned_target,&
    cross_gamma_rows,cross_conjugate_sector,8,cross_gamma_fingerprint,.false.,0d0,1d-12,&
    sector_trial_conjugate,cross_gamma_defect,&
    cross_workspace,ok,message)
  call require(.not.ok,'periodic-phase Gamma sewing rejects a nonunitary conjugate sector')
  cross_conjugate_sector=0.5d0*cross_conjugate_sector
  sector_trial_aligned=(0d0,0d0)
  do p=1,nlocal
    global_point=rank*nlocal+p
    if(global_point==1)sector_trial_aligned(p,:)=[cmplx(1d0,0d0,8),cmplx(0d0,1d0,8)]/sqrt(2d0)
    if(global_point==2)sector_trial_aligned(p,:)=[cmplx(0d0,1d0,8),cmplx(1d0,0d0,8)]/sqrt(2d0)
  enddo
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
    cross_gamma_rows,sector_trial_aligned,8,cross_gamma_fingerprint,.true.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message)
  call require(ok.and.cross_gamma_defect<1d-10.and.maxval(abs(aimag(cross_aligned_conjugate)))<1d-10,&
    'self-conjugate periodic-phase sector is fixed to one Gamma-real frame')
  cross_conjugate_sector(:,1)=sector_trial_aligned(:,2)
  cross_conjugate_sector(:,2)=sector_trial_aligned(:,1)
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
    cross_gamma_rows,cross_conjugate_sector,8,cross_gamma_fingerprint,.true.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message)
  call require(ok.and.maxval(abs(aimag(cross_aligned_conjugate)))<1d-10,&
    'self-conjugate Gamma fixing is invariant to the validation-frame gauge')
  allocate(implicit_gamma_rows(nlocal,0))
  implicit_gamma_fingerprint=ieor(ishftc(int(z'6A09E667F3BCC909',8),9),8_8)
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
    implicit_gamma_rows,cross_conjugate_sector,8,implicit_gamma_fingerprint,.true.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message,.true.)
  call require(ok.and.cross_gamma_defect<1d-10.and.maxval(abs(aimag(cross_aligned_conjugate)))<1d-10,&
    'implicit identity Gamma fixes a self-conjugate spatial sector without a dense grid operator')
  allocate(weighted_spatial_sector(nlocal,2),spatial_integration_weights(nlocal))
  weighted_spatial_sector=2d0*sector_trial_aligned
  spatial_integration_weights=0.25d0
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,weighted_spatial_sector,&
    implicit_gamma_rows,weighted_spatial_sector,8,implicit_gamma_fingerprint,.true.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message,.true.,spatial_integration_weights)
  call require(ok.and.cross_gamma_defect<1d-10,&
    'Gamma sewing uses the real-space integration measure for materialized sectors')
  deallocate(weighted_spatial_sector,spatial_integration_weights)
  if(nproc>1)then
    call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
      cross_gamma_rows,cross_conjugate_sector,8,implicit_gamma_fingerprint,.true.,0d0,1d-12,&
      cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message,rank==0)
    call require(.not.ok,'implicit identity Gamma rejects a rank-disagreeing operator branch')
  endif
  deallocate(implicit_gamma_rows)
  if(nlocal>1)then
    global_point=int(sector_ids(2));sector_ids(2)=sector_ids(1)
  endif
  call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
    cross_gamma_rows,cross_conjugate_sector,8,cross_gamma_fingerprint,.true.,0d0,1d-12,&
    cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message)
  call require(merge(.not.ok,ok,nlocal>1),'Gamma sewing rejects duplicate rows owned by one rank')
  if(nlocal>1)sector_ids(2)=int(global_point,8)
  if(nproc>1)then
    call sew_dg_w90_periodic_phase_conjugate_sector(MPI_COMM_WORLD,sector_ids,sector_trial_aligned,&
      cross_gamma_rows,cross_conjugate_sector,8,cross_gamma_fingerprint,rank==0,0d0,1d-12,&
      cross_aligned_conjugate,cross_gamma_defect,cross_workspace,ok,message)
    call require(.not.ok,'Gamma sewing rejects rank-disagreeing self-conjugacy metadata')
  endif
  cross_conjugate_sector=conjg(cross_target_sector)
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,0_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(.not.ok,'cross-character alignment rejects a zero W90 frame receipt')
  sector_ids(1)=9_8
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(.not.ok,'cross-character alignment rejects a row outside the global extent')
  sector_ids(1)=int(rank*nlocal+1,8)
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    cross_aligned_target,cross_singular_values,&
    cross_polar_defect,cross_fingerprint,cross_workspace,ok,message)
  call require(ok.and.cross_polar_defect<1d-10.and.minval(cross_singular_values)>0.4d0.and.&
    maxval(abs(cross_aligned_target(:,1)-exp(cmplx(0d0,0.4d0,8))*cross_target_sector(:,1)))<1d-10.and.&
    maxval(abs(cross_aligned_target(:,2)-exp(cmplx(0d0,-0.7d0,8))*cross_target_sector(:,2)))<1d-10,&
    'Wannier90 reference links orthogonal character sectors without center matching')
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_rotated_target,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,&
    cross_polar_defect,sector_trial_fingerprint,cross_workspace,ok,message)
  call require(ok.and.sector_trial_fingerprint==cross_fingerprint.and.&
    maxval(abs(sector_trial_aligned-cross_aligned_target))<1d-10,&
    'cross-character gauge is invariant under target-sector frame rotation')
  cross_localization_weights=0d0
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(.not.ok,'cross-character alignment rejects a singular weighted localization link')
  cross_localization_weights=[1d0,2d0,-1d0,-2d0]
  if(rank==0.and.nproc>1)cross_localization_weights(1)=1.5d0
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(nproc==1.or..not.ok,'cross-character alignment rejects rank-disagreeing localization weights')
  cross_localization_weights=[1d0,2d0,-1d0,-2d0]
  cross_localization_weights(1)=1d300
  call align_dg_w90_cross_character_sector_gauge(MPI_COMM_WORLD,sector_ids,cross_reference_sector,&
    cross_target_sector,cross_w90_reference,cross_localization_weights,8,4,123_8,0d0,1d-12,&
    sector_trial_aligned,cross_singular_values,cross_polar_defect,sector_trial_fingerprint,&
    cross_workspace,ok,message)
  call require(.not.ok,'cross-character alignment rejects finite-huge localization weights before multiplication')
  cross_localization_weights=[1d0,2d0,-1d0,-2d0]
  localization_cluster_spectrum=[0d0,1d0,1d0,2d0]
  call validate_dg_w90_localization_cluster(localization_cluster_spectrum,2,1d-12,ok,message)
  call require(.not.ok.and.index(trim(message),'splits')>0,&
    'Wannier90 alignment rejects a split degenerate localization cluster')
  call validate_dg_w90_localization_cluster(localization_cluster_spectrum,3,1d-12,ok,message)
  call require(ok,'Wannier90 alignment retains a complete degenerate localization cluster')
  p=merge(nlocal,max(0,nlocal-1),rank==0)
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values(:,1:p),local_anchors(:,1:p),&
    local_weights(1:p),local_fractional(:,1:p),test_nncell,huge(0_8),assembled_m,assembled_a,&
    matrix_estimate,matrix_peak,ok,message)
  call require(ok,'Wannier90 matrix assembly accepts uneven local point ownership')
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,1_8,assembled_m,assembled_a,matrix_estimate,matrix_peak,&
    ok,message)
  call require(.not.ok,'Wannier90 matrix assembly rejects a low coordinator byte limit')
  call assemble_dg_w90_gamma_matrices(MPI_COMM_WORLD,local_values,local_anchors,local_weights,&
    local_fractional,test_nncell,merge(1_8,huge(0_8),rank==0),assembled_m,assembled_a,&
    matrix_estimate,matrix_peak,ok,message)
  call require(.not.ok,'Wannier90 matrix assembly rejects rank-inconsistent byte limits')
  allocate(gauge_values,source=local_values)
  allocate(gauge_gradients(3,2,nlocal));gauge_gradients=(0d0,0d0)
  allocate(gauge_ids(nlocal));gauge_ids=[(int(rank*nlocal+p,8),p=1,nlocal)]
  gauge_transform=reshape([cmplx(0d0,0d0,8),cmplx(-1d0,0d0,8),&
    cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)],[2,2])
  gauge_centers=reshape([0.2d0,0d0,0d0,0.2d0,0d0,0d0],[3,2])
  gauge_spreads=[20d0,10d0]
  call apply_dg_w90_gamma_transform(MPI_COMM_WORLD,gauge_ids,gauge_values,gauge_gradients,&
    gauge_transform,gauge_centers,1d-12,ok,message,gauge_spreads)
  call require(ok.and.maxval(abs(gauge_values(1,:)-local_values(1,:)))<1d-12.and.&
    maxval(abs(gauge_values(2,:)+local_values(2,:)))<1d-12.and.&
    abs(gauge_transform(1,1)-1d0)<1d-12.and.abs(gauge_transform(2,2)+1d0)<1d-12.and.&
    maxval(abs(gauge_transform-reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
      cmplx(0d0,0d0,8),cmplx(-1d0,0d0,8)],[2,2])))<1d-12.and.&
    maxval(abs(gauge_spreads-[10d0,20d0]))<1d-12,trim(message))
  transform(1,1)=2d0
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'nonunitary Wannier90 transform rejection')
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=cmplx(1d0,1d-4,8)
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'complex Gamma Wannier90 gauge rejection')
  transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0;spread(3)=0.9d0
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'increased Wannier90 gauge-dependent spread rejection')
  spread(3)=0.7d0;centers(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call validate_dg_w90_result(transform,centers,spreads,spread,0.8d0,1d-12,ok,message)
  call require(.not.ok,'nonfinite Wannier90 center rejection')
#ifdef USE_WANNIER90
  lattice=0d0;reciprocal=0d0
  lattice(1,1)=10d0;lattice(2,2)=10d0;lattice(3,3)=10d0
  reciprocal(1,1)=2d0*acos(-1d0)/10d0
  reciprocal(2,2)=reciprocal(1,1);reciprocal(3,3)=reciprocal(1,1)
  atoms_cart=0d0;atom_symbols(1)='H ';eigenvalues=0d0
  if(rank==0)then
    open(newunit=log_unit,file='ow_w90_one_band.dmn',status='old',iostat=p)
    if(p==0)close(log_unit,status='delete')
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call setup_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,1,1,10,nntot,nncell,ok,message)
  call require(.not.ok,'Wannier90 setup rejects missing DMN')
  if(rank==0)then
    open(newunit=log_unit,file='ow_w90_one_band.dmn',status='replace')
    write(log_unit,'(a)')'SALMON SAWF Gamma-only symmetry data'
    write(log_unit,'(4i9)')1,1,1,1
    write(log_unit,*);write(log_unit,*)1;write(log_unit,*);write(log_unit,*)1
    write(log_unit,*);write(log_unit,*)1;write(log_unit,*)
    write(log_unit,*)cmplx(1d0,0d0,8);write(log_unit,*);write(log_unit,*)cmplx(1d0,0d0,8)
    close(log_unit)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call setup_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,1,1,10,nntot,nncell,ok,message)
  call require(ok.and.nntot>0,trim(message))
  win_has_random_projection=.false.
  if(rank==0)then
    open(newunit=win_unit,file='ow_w90_one_band.win',status='old',action='read',iostat=win_io)
    if(win_io==0)then
      do
        read(win_unit,'(a)',iostat=win_io)win_line
        if(win_io/=0)exit
        if(index(adjustl(win_line),'random')==1)win_has_random_projection=.true.
      enddo
      close(win_unit)
    endif
  endif
  call MPI_Bcast(win_io,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call require(win_io==0.or.win_io<0,'Wannier90 setup writes its input file')
  call MPI_Bcast(win_has_random_projection,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call require(.not.win_has_random_projection,&
    'externally supplied Wannier90 A matrices must not retain random projections')
  if(rank==0)then
    allocate(m_matrix(1,1,nntot),a_matrix(1,1));m_matrix=(1d0,0d0);a_matrix=(1d0,0d0)
  else
    allocate(m_matrix(0,0,0),a_matrix(0,0))
  endif
  call run_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d6,1d-10,10,library_transform,&
    library_centers,library_spreads,library_spread,ok,message)
  call require(ok,trim(message))
  call require(abs(abs(library_transform(1,1))-1d0)<1d-10,&
    'one-band Wannier90 library returns a unitary Gamma transform')
  if(rank==0)m_matrix(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call run_dg_w90_gamma_library(MPI_COMM_WORLD,'ow_w90_one_band',lattice,reciprocal,&
    atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d6,1d-10,10,library_transform,&
    library_centers,library_spreads,library_spread,ok,message)
  call require(.not.ok,'nonfinite Wannier90 M matrix is rejected before library entry')
#endif
  if(rank==0)write(*,'(a)')'PASS Wannier90 MLWF adapter validation'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_bad/=0)error stop label
  end subroutine
end program
