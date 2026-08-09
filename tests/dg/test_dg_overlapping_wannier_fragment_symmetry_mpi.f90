#include "config.h"
program test_dg_overlapping_wannier_fragment_symmetry_mpi
  use mpi
  use dg_overlapping_wannier_symmetry, only: select_dg_exact_fragment_subgroup, &
    promote_dg_exact_global_subgroup,build_dg_fragment_site_stabilizer, &
    fingerprint_dg_exact_fragment_symmetry,build_dg_fragment_permuted_representation, &
    build_dg_fragment_symmetry_orbits,build_dg_symmetry_constrained_pair_generator
  use iso_fortran_env,only:int64
  implicit none
  integer :: ierr,rank,nproc,i,j
  integer(int64) :: c4_fingerprint,c1_fingerprint,tolerance_fingerprint,translation_fingerprint
  integer :: product_table(4,4)
  integer :: invalid_product_table(4,4)
  integer :: affine_rotation(3,3,6)
  integer,allocatable :: affine_product(:,:)
  integer,allocatable :: subgroup(:)
  logical :: fragment_exact(2,4)
  logical :: affine_allowed(6)
  real(8) :: scalar_block_residual(3,4),vector_block_residual(3,4)
  real(8) :: affine_translation(3,6),fragment_center(3),site_residual
  real(8) :: pair_centers(3,2),inversion_cartesian(3,3,1)
  complex(8) :: local_pair_representation(1,1,1)
  complex(8),allocatable :: global_pair_representation(:,:,:)
  complex(8) :: dense_representation(4,4,2),broken_representation(4,4,2)
  complex(8),allocatable :: constrained_generator(:,:)
  integer,allocatable :: generator_active_indices(:)
  integer :: generator_product(2,2)
  real(8) :: antihermiticity_defect,commutator_defect
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
