#include "config.h"
program test_dg_overlapping_wannier_fragment_symmetry_mpi
  use mpi
  use dg_overlapping_wannier_symmetry, only: select_dg_exact_fragment_subgroup, &
    promote_dg_exact_global_subgroup,build_dg_fragment_site_stabilizer, &
    fingerprint_dg_exact_fragment_symmetry,build_dg_fragment_permuted_representation, &
    build_dg_fragment_symmetry_orbits,build_dg_symmetry_constrained_pair_generator, &
    factor_dg_affine_translation_cocycle,measure_dg_hamiltonian_density_commutators,&
    symmetrize_dg_distributed_pencil_rows
  use iso_fortran_env,only:int64
  implicit none
  integer :: ierr,rank,nproc,i,j
  integer(int64) :: c4_fingerprint,c1_fingerprint,tolerance_fingerprint,translation_fingerprint
  integer :: product_table(4,4)
  integer :: invalid_product_table(4,4)
  integer :: affine_rotation(3,3,6)
  integer,allocatable :: affine_product(:,:)
  integer,allocatable :: translation_subgroup(:),point_representatives(:),point_product(:,:),&
    translation_cocycle(:,:)
  integer,allocatable :: subgroup(:)
  logical :: fragment_exact(2,4)
  logical :: affine_allowed(6)
  real(8) :: scalar_block_residual(3,4),vector_block_residual(3,4)
  real(8) :: affine_translation(3,6),fragment_center(3),site_residual
  real(8) :: pair_centers(3,2),inversion_cartesian(3,3,1)
  complex(8) :: local_pair_representation(1,1,1)
  complex(8),allocatable :: global_pair_representation(:,:,:)
  complex(8) :: dense_representation(4,4,2),broken_representation(4,4,2)
  complex(8) :: gate_representation(2,2,2),gate_hamiltonian(2,2),gate_density(2,2)
  complex(8) :: gate_overlap(2,2)
  complex(8),allocatable :: constrained_generator(:,:)
  integer,allocatable :: generator_active_indices(:)
  integer :: generator_product(2,2)
  real(8) :: antihermiticity_defect,commutator_defect
  real(8),allocatable :: hamiltonian_commutator(:),density_commutator(:)
  integer(int64),allocatable :: pencil_row_ids(:)
  complex(8),allocatable :: pencil_h_rows(:,:),pencil_s_rows(:,:),pencil_rho_rows(:,:),&
    pencil_artifact_rows(:,:),sym_h_rows(:,:),sym_s_rows(:,:),sym_rho_rows(:,:),&
    pencil_component_rows(:,:,:)
  complex(8) :: inversion_generator(2,2,1)
  integer :: inversion_product(2,2),inversion_generators(1)
  real(8) :: pencil_before(3),pencil_after(3),artifact_change,artifact_magnitude
  real(8) :: pencil_component_residual(2)
  real(8) :: pencil_dense_error
  integer(int64) :: pencil_workspace_peak
  integer,allocatable :: fragment_permutation(:,:)
  integer :: multi_orbit_map(4,2),identity_orbit_map(4,1),invalid_orbit_map(4,2)
  integer,allocatable :: fragment_orbit(:),orbit_representative(:)
  real(8) :: atom(4),boundary(4),grid(4),center(4)
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  do i=1,4
    do j=1,4
      product_table(i,j)=modulo(i+j-2,4)+1
    end do
  end do

  atom=1d-13; boundary=2d-13; grid=3d-13; center=4d-13
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,2,3,4]),'roundoff-scale exact C4')

  atom=[1d-13,2d-4,1d-13,2d-4]
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,3]),'physical displacement reduces C4 to C2')

  atom=[1d-13,2d-4,3d-4,4d-4]
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.size(subgroup)==1.and.subgroup(1)==1,'physical displacement permits C1')

  atom=1d-13; grid(2)=2d-4
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,3]),'grid-incompatible operations form closed subgroup')

  grid=1d-13; boundary(1)=2d-4
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(.not.ok.and.index(message,'identity')>0,'invalid identity is fatal')

  invalid_product_table=product_table
  invalid_product_table(2,4)=2
  invalid_product_table(4,2)=2
  boundary=1d-13
  call select_dg_exact_fragment_subgroup(invalid_product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(.not.ok.and.index(message,'group')>0,'non-group product table is fatal')

  fragment_exact=.true.;scalar_block_residual=1d-13;vector_block_residual=1d-13
  scalar_block_residual(:,2)=2d-4;scalar_block_residual(:,4)=2d-4
  call promote_dg_exact_global_subgroup(product_table,fragment_exact,scalar_block_residual, &
    vector_block_residual,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,3]),'cross-block covariance promotes only exact C2')
  fragment_exact(2,3)=.false.
  call promote_dg_exact_global_subgroup(product_table,fragment_exact,scalar_block_residual, &
    vector_block_residual,1d-10,subgroup,ok,message)
  call require(ok.and.size(subgroup)==1.and.subgroup(1)==1, &
    'locally broken neighboring fragment reduces promotion to C1')

  affine_rotation=0;affine_translation=0d0;affine_allowed=.true.
  affine_rotation(:,:,1)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
  affine_rotation(:,:,2)=reshape([0,1,0,-1,0,0,0,0,1],[3,3])
  affine_rotation(:,:,3)=reshape([-1,0,0,0,-1,0,0,0,1],[3,3])
  affine_rotation(:,:,4)=reshape([0,-1,0,1,0,0,0,0,1],[3,3])
  affine_rotation(:,:,5)=affine_rotation(:,:,2);affine_translation(1,5)=0.5d0
  affine_rotation(:,:,6)=reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])
  affine_translation(:,6)=0.5d0
  fragment_center=0d0
  call build_dg_fragment_site_stabilizer(affine_rotation,affine_translation,fragment_center, &
    affine_allowed,1d-10,subgroup,affine_product,site_residual,ok,message)
  call require(ok.and.all(subgroup==[1,2,3,4]),'origin-centered C4 site stabilizer')
  call require(site_residual<1d-13,'C4 site residual')
  affine_allowed=.false.;affine_allowed([1,6])=.true.;fragment_center=0.25d0
  call build_dg_fragment_site_stabilizer(affine_rotation,affine_translation,fragment_center, &
    affine_allowed,1d-10,subgroup,affine_product,site_residual,ok,message)
  call require(ok.and.all(subgroup==[1,6]),'fractional-center inversion site stabilizer')
  call require(all(affine_product==reshape([1,2,2,1],[2,2])),'inversion affine closure')

  affine_rotation=0;affine_translation=0d0
  affine_rotation(:,:,1)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
  affine_rotation(:,:,2)=affine_rotation(:,:,1);affine_translation(3,2)=0.5d0
  affine_rotation(:,:,3)=reshape([-1,0,0,0,-1,0,0,0,1],[3,3])
  affine_translation(3,3)=0.25d0
  affine_rotation(:,:,4)=affine_rotation(:,:,3);affine_translation(3,4)=0.75d0
  product_table=reshape([1,2,3,4,2,1,4,3,3,4,2,1,4,3,1,2],[4,4])
  call factor_dg_affine_translation_cocycle(affine_rotation(:,:,1:4),&
    affine_translation(:,1:4),product_table,1d-12,translation_subgroup,&
    point_representatives,point_product,translation_cocycle,ok,message)
  call require(ok.and.all(translation_subgroup==[1,2]).and.&
    all(point_representatives==[1,3]),'screw group factors into translations and point co-group')
  call require(point_product(2,2)==1.and.translation_cocycle(2,2)==2,&
    'screw-square translation is retained as a nontrivial cocycle')
  affine_rotation(:,:,3)=reshape([1,0,0,0,-1,0,0,0,1],[3,3])
  affine_rotation(:,:,4)=affine_rotation(:,:,3)
  call factor_dg_affine_translation_cocycle(affine_rotation(:,:,1:4),&
    affine_translation(:,1:4),product_table,1d-12,translation_subgroup,&
    point_representatives,point_product,translation_cocycle,ok,message)
  call require(ok.and.translation_cocycle(2,2)==2,&
    'glide-square translation is retained as a nontrivial cocycle')

  gate_representation=(0d0,0d0);gate_representation(1,1,1)=1d0
  gate_representation(2,2,1)=1d0;gate_representation(1,1,2)=1d0
  gate_representation(2,2,2)=-1d0
  gate_hamiltonian=(0d0,0d0);gate_hamiltonian(1,1)=-1d0;gate_hamiltonian(2,2)=2d0
  gate_density=(0d0,0d0);gate_density(1,1)=1d0
  call measure_dg_hamiltonian_density_commutators(gate_representation,gate_hamiltonian,&
    gate_density,1d-12,hamiltonian_commutator,density_commutator,ok,message)
  call require(ok.and.maxval(hamiltonian_commutator)<1d-14.and.&
    maxval(density_commutator)<1d-14,'H and density projector commute with inversion')
  gate_density(1,2)=0.1d0;gate_density(2,1)=0.1d0
  call measure_dg_hamiltonian_density_commutators(gate_representation,gate_hamiltonian,&
    gate_density,1d-12,hamiltonian_commutator,density_commutator,ok,message)
  call require(ok.and.maxval(hamiltonian_commutator)<1d-14.and.&
    density_commutator(2)>0.1d0,'density-projector symmetry breaking is measured independently of H')

  allocate(pencil_row_ids(count([(mod(i-1,nproc)==rank,i=1,2)])))
  j=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    j=j+1;pencil_row_ids(j)=i
  enddo
  allocate(pencil_h_rows(size(pencil_row_ids),2),pencil_s_rows(size(pencil_row_ids),2),&
    pencil_rho_rows(size(pencil_row_ids),2),pencil_artifact_rows(size(pencil_row_ids),2))
  gate_hamiltonian=reshape([cmplx(-1d0,0d0,8),cmplx(0.2d0,0d0,8),&
    cmplx(0.2d0,0d0,8),cmplx(2d0,0d0,8)],[2,2])
  gate_density=reshape([cmplx(1d0,0d0,8),cmplx(-0.15d0,0d0,8),&
    cmplx(-0.15d0,0d0,8),cmplx(0.4d0,0d0,8)],[2,2])
  pencil_h_rows=gate_hamiltonian(int(pencil_row_ids),:)
  gate_overlap=reshape([cmplx(1d0,0d0,8),cmplx(0.1d0,0d0,8),&
    cmplx(0.1d0,0d0,8),cmplx(1.5d0,0d0,8)],[2,2])
  pencil_s_rows=gate_overlap(int(pencil_row_ids),:)
  pencil_rho_rows=gate_density(int(pencil_row_ids),:)
  pencil_artifact_rows=0d0
  do i=1,size(pencil_row_ids)
    if(pencil_row_ids(i)==1)pencil_artifact_rows(i,1)=0.3d0
    if(pencil_row_ids(i)==2)pencil_artifact_rows(i,2)=-0.3d0
  enddo
  inversion_generator=0d0;inversion_generator(1,1,1)=1d0;inversion_generator(2,2,1)=-1d0
  inversion_product=reshape([1,2,2,1],[2,2]);inversion_generators=2
  allocate(pencil_component_rows(size(pencil_row_ids),2,2))
  pencil_component_rows(:,:,1)=pencil_h_rows;pencil_component_rows(:,:,2)=pencil_artifact_rows
  call symmetrize_dg_distributed_pencil_rows(MPI_COMM_WORLD,pencil_row_ids,pencil_h_rows,&
    pencil_s_rows,pencil_rho_rows,pencil_artifact_rows,inversion_generator,inversion_generators,&
    inversion_product,[1],[1,2],1d-12,sym_h_rows,sym_s_rows,sym_rho_rows,pencil_before,pencil_after,&
    artifact_change,artifact_magnitude,pencil_workspace_peak,ok,message,&
    pencil_component_rows,pencil_component_residual)
  call require(ok,trim(message))
  call require(maxval(pencil_before)>0.1d0.and.maxval(pencil_after)<1d-12,&
    'full-group average removes nonsymmetric pencil error')
  call require(abs(pencil_component_residual(1)-pencil_before(1))<1d-14.and.&
    pencil_component_residual(2)<1d-14,&
    'component covariance diagnostic distinguishes broken and invariant operator terms')
  pencil_dense_error=0d0
  do i=1,size(pencil_row_ids)
    if(pencil_row_ids(i)==1)then
      pencil_dense_error=max(pencil_dense_error,abs(sym_h_rows(i,1)+1d0),abs(sym_s_rows(i,1)-1d0),&
        abs(sym_rho_rows(i,1)-1d0),abs(sym_h_rows(i,2)),abs(sym_s_rows(i,2)),abs(sym_rho_rows(i,2)))
    else
      pencil_dense_error=max(pencil_dense_error,abs(sym_h_rows(i,2)-2d0),abs(sym_s_rows(i,2)-1.5d0),&
        abs(sym_rho_rows(i,2)-0.4d0),abs(sym_h_rows(i,1)),abs(sym_s_rows(i,1)),abs(sym_rho_rows(i,1)))
    endif
  enddo
  call require(pencil_dense_error<1d-12,'factorized affine average matches the full dense reference')
  call require(artifact_change<1d-12.and.artifact_magnitude>0.29d0,&
    'symmetric fragment artifact survives and is reported independently')
  call require(pencil_workspace_peak>0_int64,'pencil symmetry workspace is measured')
  inversion_generator(2,2,1)=cmplx(0d0,1d0,8)
  call symmetrize_dg_distributed_pencil_rows(MPI_COMM_WORLD,pencil_row_ids,pencil_h_rows,&
    pencil_s_rows,pencil_rho_rows,pencil_artifact_rows,inversion_generator,inversion_generators,&
    inversion_product,[1],[1,2],1d-12,sym_h_rows,sym_s_rows,sym_rho_rows,pencil_before,pencil_after,&
    artifact_change,artifact_magnitude,pencil_workspace_peak,ok,message)
  call require(.not.ok,'generator representation inconsistent with the affine product is rejected')

  c4_fingerprint=fingerprint_dg_exact_fragment_symmetry(affine_rotation(:,:,1:4),product_table,1d-10)
  c1_fingerprint=fingerprint_dg_exact_fragment_symmetry(affine_rotation(:,:,1:1),reshape([1],[1,1]),1d-10)
  tolerance_fingerprint=fingerprint_dg_exact_fragment_symmetry(&
    affine_rotation(:,:,1:4),product_table,2d-10)
  translation_fingerprint=fingerprint_dg_exact_fragment_symmetry(&
    affine_rotation(:,:,1:4),product_table,1d-10,affine_translation(:,1:4)+0.125d0)
  call require(c4_fingerprint/=0_int64,'exact fragment symmetry fingerprint is nonzero')
  call require(c4_fingerprint/=c1_fingerprint,'C4 and displaced C1 checkpoint evidence differ')
  call require(c4_fingerprint/=tolerance_fingerprint,'symmetry tolerance is checkpoint evidence')
  call require(c4_fingerprint/=translation_fingerprint,'affine translation is checkpoint evidence')
  pair_centers=reshape([0.25d0,0.5d0,0.5d0,0.75d0,0.5d0,0.5d0],[3,2])
  inversion_cartesian(:,:,1)=0d0
  inversion_cartesian(1,1,1)=-1d0;inversion_cartesian(2,2,1)=1d0
  inversion_cartesian(3,3,1)=1d0;local_pair_representation=1d0
  call build_dg_fragment_permuted_representation(local_pair_representation,inversion_cartesian,&
    pair_centers,1d-12,global_pair_representation,fragment_permutation,ok,message)
  call require(ok.and.all(fragment_permutation(:,1)==[2,1]),'point operation permutes fragment centers')
  call require(abs(global_pair_representation(2,1,1)-1d0)<1d-14.and.&
    abs(global_pair_representation(1,2,1)-1d0)<1d-14,&
    'global point representation contains source-to-target fragment blocks')

  multi_orbit_map(:,1)=[1,2,3,4]
  multi_orbit_map(:,2)=[2,1,4,3]
  call build_dg_fragment_symmetry_orbits(multi_orbit_map,fragment_orbit,&
    orbit_representative,ok,message)
  call require(ok.and.all(fragment_orbit==[1,1,2,2]),&
    'disconnected fragment symmetry components form separate orbits')
  call require(all(orbit_representative==[1,1,3,3]),&
    'each fragment orbit has its own deterministic representative')
  identity_orbit_map(:,1)=[1,2,3,4]
  call build_dg_fragment_symmetry_orbits(identity_orbit_map,fragment_orbit,&
    orbit_representative,ok,message)
  call require(ok.and.all(fragment_orbit==[1,2,3,4]).and.&
    all(orbit_representative==[1,2,3,4]),&
    'fully broken symmetry retains independent fragment-local construction')
  invalid_orbit_map=multi_orbit_map;invalid_orbit_map(:,2)=[2,2,4,3]
  call build_dg_fragment_symmetry_orbits(invalid_orbit_map,fragment_orbit,&
    orbit_representative,ok,message)
  call require(.not.ok.and.index(message,'permutation')>0,&
    'non-bijective fragment symmetry operation is rejected')

  dense_representation=(0d0,0d0)
  do i=1,4;dense_representation(i,i,1)=1d0;end do
  dense_representation(1,1,2)=1d0/sqrt(2d0);dense_representation(1,2,2)=1d0/sqrt(2d0)
  dense_representation(2,1,2)=1d0/sqrt(2d0);dense_representation(2,2,2)=-1d0/sqrt(2d0)
  dense_representation(3:4,3:4,2)=dense_representation(1:2,1:2,2)
  generator_product=reshape([1,2,2,1],[2,2])
  call build_dg_symmetry_constrained_pair_generator(1,3,(1d0,0d0),dense_representation,&
    generator_product,1d-12,constrained_generator,antihermiticity_defect,&
    commutator_defect,generator_active_indices,ok,message)
  call require(ok,'dense-representation generator group average')
  call require(antihermiticity_defect<1d-12.and.commutator_defect<1d-12,&
    'constrained generator is anti-Hermitian and symmetry commuting')
  call require(all(generator_active_indices==[1,2,3,4]).and.&
    abs(constrained_generator(2,4))>1d-3,&
    'dense symmetry representation expands a sparse pair seed')
  call build_dg_symmetry_constrained_pair_generator(1,3,(1d0,0d0),&
    dense_representation(:,:,1:1),&
    reshape([1],[1,1]),1d-12,constrained_generator,antihermiticity_defect,&
    commutator_defect,generator_active_indices,ok,message)
  call require(ok.and.all(generator_active_indices==[1,3]).and.&
    maxval(abs(constrained_generator-reshape([(0d0,0d0),(-1d0,0d0),&
      (1d0,0d0),(0d0,0d0)],[2,2])))<1d-14,&
    'identity-only symmetry retains the local pair generator')
  broken_representation=dense_representation;broken_representation(1,1,2)=2d0
  call build_dg_symmetry_constrained_pair_generator(1,3,(1d0,0d0),broken_representation,&
    generator_product,1d-12,constrained_generator,antihermiticity_defect,&
    commutator_defect,generator_active_indices,ok,message)
  call require(.not.ok.and.index(message,'unitary')>0,&
    'nonunitary symmetry representation is rejected')

  if(rank==0)write(*,'(a,i0,a)')'PASS exact buffered-fragment symmetry on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_failure/=0)error stop label
  end subroutine require
end program test_dg_overlapping_wannier_fragment_symmetry_mpi
