#include "config.h"
module dg_hybrid_fragment_wannier
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_overlapping_wannier_w90,only:DG_W90_UNCONSTRAINED,&
    setup_dg_w90_gamma_library,assemble_dg_w90_gamma_matrices,&
    run_dg_w90_gamma_library,apply_dg_w90_gamma_transform
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  integer,parameter::message_length=1024

  type,public::s_dg_hybrid_fragment_wannier_receipt
    integer::fragment_id=0
    integer::basis_generation=0
    integer::candidate_rank=0
    integer::retained_rank=0
    integer::setup_count=0
    integer::run_count=0
    integer(int64)::seed_fingerprint=0_int64
    integer(int64)::basis_fingerprint=0_int64
    integer(int64)::transform_fingerprint=0_int64
    integer(int64)::replicated_payload_fingerprint=0_int64
    integer(int64)::distributed_wannier_fingerprint=0_int64
    real(real64)::seed_reconstruction_defect=huge(1d0)
  end type s_dg_hybrid_fragment_wannier_receipt

  type,public::s_dg_hybrid_fragment_wannier_cache
    logical::valid=.false.
    type(s_dg_hybrid_fragment_wannier_receipt)::receipt
    integer::fragment_comm_rank=-1
    integer::fragment_comm_size=0
    integer(int64)::local_row_layout_fingerprint=0_int64
    integer(int64),allocatable::local_grid_ids(:)
    complex(real64),allocatable::wannier_values(:,:)
    complex(real64),allocatable::candidate_compression(:,:)
    complex(real64),allocatable::wannier_transform(:,:)
    complex(real64),allocatable::dc_seed_coefficients_in_wannier(:,:)
    real(real64),allocatable::physical_dc_seed_energies(:)
    real(real64),allocatable::physical_dc_seed_occupations(:)
  end type s_dg_hybrid_fragment_wannier_cache

  interface hash_real_array
    module procedure hash_real_vector
    module procedure hash_real_matrix
  end interface hash_real_array

  public::build_dg_hybrid_fragment_wannier

contains

  subroutine build_dg_hybrid_fragment_wannier(comm_total,comm_fragment,fragment_id,&
      basis_generation,seed_directory,grid_ids,grid_weights,dc_seed_values,&
      dc_seed_energies,dc_seed_occupations,buffer_candidate_values,&
      projector_candidate_values,metric_tolerance,fragment_real_lattice,&
      fragment_reciprocal_lattice,atom_symbols,atoms_cart,fractional_coordinates,&
      num_iter,localization_tolerance,coordinator_byte_limit,cache,ok,message)
    integer,intent(in)::comm_total,comm_fragment,fragment_id,basis_generation,num_iter
    character(*),intent(in)::seed_directory
    integer(int64),intent(in)::grid_ids(:)
    real(real64),intent(in)::grid_weights(:),dc_seed_energies(:),dc_seed_occupations(:)
    complex(real64),intent(in)::dc_seed_values(:,:),buffer_candidate_values(:,:),&
      projector_candidate_values(:,:)
    real(real64),intent(in)::metric_tolerance,fragment_real_lattice(3,3),&
      fragment_reciprocal_lattice(3,3),atoms_cart(:,:),fractional_coordinates(:,:),&
      localization_tolerance
    character(*),intent(in)::atom_symbols(:)
    integer(int64),intent(in)::coordinator_byte_limit
    type(s_dg_hybrid_fragment_wannier_cache),intent(inout)::cache
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_fragment_wannier_cache)::working_cache
    real(real64)::inverse_lattice(3,3)
    integer(int64)::seed_fingerprint,basis_fingerprint,local_layout_fingerprint
    character(message_length)::local_message,seed_name
    logical::local_ok,build_required

    ok=.false.;message='';local_message='';seed_name=''
    call validate_fragment_partition(comm_total,comm_fragment,fragment_id,local_ok,local_message)
    call canonical_total_status(comm_total,local_ok,local_message,ok,message)
    if(.not.ok)return

    call validate_fragment_contract(comm_fragment,fragment_id,basis_generation,seed_directory,&
      grid_ids,grid_weights,dc_seed_values,dc_seed_energies,dc_seed_occupations,&
      buffer_candidate_values,projector_candidate_values,metric_tolerance,&
      fragment_real_lattice,fragment_reciprocal_lattice,atom_symbols,atoms_cart,&
      fractional_coordinates,num_iter,localization_tolerance,coordinator_byte_limit,&
      inverse_lattice,seed_fingerprint,basis_fingerprint,local_layout_fingerprint,&
      local_ok,local_message)
    call canonical_total_status(comm_total,local_ok,local_message,ok,message)
    if(.not.ok)return

    call classify_fragment_cache(comm_fragment,cache,fragment_id,basis_generation,&
      size(dc_seed_values,1),size(dc_seed_values,1)+size(buffer_candidate_values,1)+&
      size(projector_candidate_values,1),grid_ids,seed_fingerprint,basis_fingerprint,&
      local_layout_fingerprint,dc_seed_energies,dc_seed_occupations,build_required,&
      local_ok,local_message)
    call canonical_total_status(comm_total,local_ok,local_message,ok,message)
    if(.not.ok)return

    local_ok=.true.;local_message=''
    if(build_required)call create_fragment_namespace(comm_fragment,seed_directory,&
      fragment_id,basis_generation,seed_name,local_ok,local_message)
    call canonical_total_status(comm_total,local_ok,local_message,ok,message)
    if(.not.ok)return

    local_ok=.true.;local_message=''
    if(build_required)then
      call construct_fragment_wannier(comm_fragment,fragment_id,basis_generation,seed_name,&
        grid_ids,grid_weights,dc_seed_values,dc_seed_energies,dc_seed_occupations,&
        buffer_candidate_values,projector_candidate_values,metric_tolerance,&
        fragment_real_lattice,fragment_reciprocal_lattice,inverse_lattice,atom_symbols,&
        atoms_cart,fractional_coordinates,num_iter,localization_tolerance,&
        coordinator_byte_limit,seed_fingerprint,basis_fingerprint,local_layout_fingerprint,&
        working_cache,local_ok,local_message)
    endif
    call canonical_total_status(comm_total,local_ok,local_message,ok,message)
    if(.not.ok)return

    if(build_required)call move_fragment_cache(working_cache,cache)
    ok=.true.;message=''
#else
    ok=.false.;message='fragment Wannier construction requires MPI'
#endif
  end subroutine build_dg_hybrid_fragment_wannier

#ifdef USE_MPI
  subroutine validate_fragment_partition(comm_total,comm_fragment,fragment_id,ok,message)
    integer,intent(in)::comm_total,comm_fragment,fragment_id
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),allocatable::root_fragment_ids(:)
    integer::fragment_rank,total_rank,total_size,minimum_fragment_id,maximum_fragment_id
    integer(int64)::root_fragment_id
    integer::i,ierr,status,global_status,allocation_status
    logical::allocation_ok
    character(message_length)::allocation_message

    ok=.false.;message='';status=0
    call MPI_Comm_rank(comm_fragment,fragment_rank,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Allreduce(fragment_id,minimum_fragment_id,1,MPI_INTEGER,MPI_MIN,comm_fragment,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Allreduce(fragment_id,maximum_fragment_id,1,MPI_INTEGER,MPI_MAX,comm_fragment,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fragment_id/=maximum_fragment_id)status=1
    call MPI_Comm_rank(comm_total,total_rank,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Comm_size(comm_total,total_size,ierr)
    if(ierr/=MPI_SUCCESS.or.total_size<1)status=1
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm_total,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment communicator membership validation failed';return
    endif
    allocation_status=0
    if(total_rank==0)then
      if(.not.extent_product_fits([total_size]))then
        allocation_status=1
      else
        allocate(root_fragment_ids(max(1,total_size)),stat=allocation_status)
      endif
    else
      allocate(root_fragment_ids(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm_total,allocation_status,&
      'fragment communicator root catalog allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    root_fragment_ids=0
    root_fragment_id=merge(int(fragment_id,int64),0_int64,fragment_rank==0)
    call MPI_Gather(root_fragment_id,1,MPI_INTEGER8,root_fragment_ids,1,MPI_INTEGER8,0,&
      comm_total,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    if(total_rank==0.and.status==0)then
      call sort_int64(root_fragment_ids(1:total_size))
      do i=2,total_size
        if(root_fragment_ids(i)>0.and.&
            root_fragment_ids(i)==root_fragment_ids(i-1))status=1
      enddo
    endif
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm_total,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment communicator root is not unique for its fragment ID';return
    endif
    ok=.true.;message=''
  end subroutine validate_fragment_partition

  subroutine validate_fragment_contract(comm,fragment_id,basis_generation,seed_directory,&
      grid_ids,grid_weights,dc_seed_values,dc_seed_energies,dc_seed_occupations,&
      buffer_values,projector_values,metric_tolerance,real_lattice,reciprocal_lattice,&
      atom_symbols,atoms_cart,fractional,num_iter,localization_tolerance,&
      coordinator_byte_limit,inverse_lattice,seed_fingerprint,basis_fingerprint,&
      local_layout_fingerprint,ok,message)
    integer,intent(in)::comm,fragment_id,basis_generation,num_iter
    character(*),intent(in)::seed_directory
    integer(int64),intent(in)::grid_ids(:)
    real(real64),intent(in)::grid_weights(:),dc_seed_energies(:),dc_seed_occupations(:)
    complex(real64),intent(in)::dc_seed_values(:,:),buffer_values(:,:),projector_values(:,:)
    real(real64),intent(in)::metric_tolerance,real_lattice(3,3),reciprocal_lattice(3,3),&
      atoms_cart(:,:),fractional(:,:),localization_tolerance
    character(*),intent(in)::atom_symbols(:)
    integer(int64),intent(in)::coordinator_byte_limit
    real(real64),intent(out)::inverse_lattice(3,3)
    integer(int64),intent(out)::seed_fingerprint,basis_fingerprint,local_layout_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nlocal,nseed,nbuffer,nprojector,natom,ierr,bad,global_bad,fragment_rank,fragment_size
    integer::dimensions(5),minimum_dimensions(5),maximum_dimensions(5)
    integer::total_rows
    integer(int64)::local_rows_64,total_rows_64,minimum_rows_64,maximum_rows_64,candidate_rank_64
    integer(int64)::replica_hash,minimum_hash,maximum_hash,row_seed,row_basis
    real(real64)::determinant

    ok=.false.;message='';inverse_lattice=0d0
    seed_fingerprint=0_int64;basis_fingerprint=0_int64;local_layout_fingerprint=0_int64
    nlocal=size(grid_ids);nseed=size(dc_seed_values,1)
    nbuffer=size(buffer_values,1);nprojector=size(projector_values,1)
    candidate_rank_64=int(nseed,int64)+int(nbuffer,int64)+int(nprojector,int64)
    natom=size(atom_symbols);dimensions=[nseed,nbuffer,nprojector,natom,len(atom_symbols)]
    bad=0
    if(fragment_id<1.or.fragment_id>999999.or.basis_generation<0.or.&
        basis_generation>99999999.or.num_iter<=0.or.coordinator_byte_limit<=0_int64)bad=10
    if(len_trim(seed_directory)<1.or.len_trim(seed_directory)>800.or.&
        index(seed_directory,achar(0))>0)bad=max(bad,11)
    if(nseed<1.or.nlocal<0.or.size(grid_weights)/=nlocal.or.&
        size(dc_seed_values,2)/=nlocal.or.size(buffer_values,2)/=nlocal.or.&
        size(projector_values,2)/=nlocal.or.any(shape(fractional)/=[3,nlocal]).or.&
        size(dc_seed_energies)/=nseed.or.size(dc_seed_occupations)/=nseed.or.&
        natom<1.or.any(shape(atoms_cart)/=[3,natom]))bad=max(bad,12)
    if(candidate_rank_64>int(huge(0),int64))bad=max(bad,12)
    if(.not.ieee_is_finite(metric_tolerance).or.metric_tolerance<=0d0.or.&
        .not.ieee_is_finite(localization_tolerance).or.localization_tolerance<=0d0)bad=max(bad,13)
    if(any(grid_ids<=0_int64).or..not.all(ieee_is_finite(grid_weights)).or.&
        any(grid_weights<=0d0).or..not.all(ieee_is_finite(fractional)))bad=max(bad,14)
    if(.not.finite_complex(dc_seed_values).or..not.finite_complex(buffer_values).or.&
        .not.finite_complex(projector_values))bad=max(bad,15)
    if(.not.all(ieee_is_finite(dc_seed_energies)).or.&
        .not.all(ieee_is_finite(dc_seed_occupations)))bad=max(bad,16)
    if(.not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or.&
        .not.all(ieee_is_finite(atoms_cart)))bad=max(bad,17)
    do ierr=1,natom
      if(len_trim(atom_symbols(ierr))==0)bad=max(bad,18)
    enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call validation_message(global_bad,message);return
    endif

    call MPI_Allreduce(dimensions,minimum_dimensions,5,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(dimensions,maximum_dimensions,5,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions))then
      message='fragment Wannier replicated dimensions disagree';return
    endif
    local_rows_64=int(nlocal,int64)
    call MPI_Allreduce(local_rows_64,total_rows_64,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call MPI_Allreduce(total_rows_64,minimum_rows_64,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(total_rows_64,maximum_rows_64,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.total_rows_64<1_int64.or.&
        total_rows_64>int(huge(0),int64).or.minimum_rows_64/=maximum_rows_64)then
      message='fragment Wannier distributed grid is empty or inconsistent';return
    endif
    total_rows=int(total_rows_64)

    replica_hash=initial_hash(101_int64)
    call hash_integer(replica_hash,fragment_id);call hash_integer(replica_hash,basis_generation)
    call hash_integer(replica_hash,nseed);call hash_integer(replica_hash,nbuffer)
    call hash_integer(replica_hash,nprojector);call hash_integer(replica_hash,natom)
    call hash_integer(replica_hash,num_iter);call hash_int64(replica_hash,coordinator_byte_limit)
    call hash_real(replica_hash,metric_tolerance);call hash_real(replica_hash,localization_tolerance)
    call hash_character(replica_hash,trimmed_directory(seed_directory))
    call hash_real_array(replica_hash,real_lattice);call hash_real_array(replica_hash,reciprocal_lattice)
    call hash_real_array(replica_hash,atoms_cart);call hash_character_array(replica_hash,atom_symbols)
    call hash_real_array(replica_hash,dc_seed_energies);call hash_real_array(replica_hash,dc_seed_occupations)
    call MPI_Allreduce(replica_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(replica_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='fragment Wannier replicated controls, geometry, or physical metadata disagree';return
    endif

    call validate_unique_grid_ids(comm,grid_ids,ok,message)
    if(.not.ok)return
    call validate_physical_seed_rank(comm,dc_seed_values,grid_weights,metric_tolerance,ok,message)
    if(.not.ok)return
    call invert_lattice(real_lattice,inverse_lattice,determinant,ok)
    if(.not.ok)then;message='fragment Wannier real lattice is singular';return;endif

    call distributed_row_hashes(comm,grid_ids,grid_weights,dc_seed_values,buffer_values,&
      projector_values,fractional,row_seed,row_basis,ok)
    if(.not.ok)then;message='fragment Wannier distributed fingerprint reduction failed';return;endif
    seed_fingerprint=initial_hash(211_int64)
    call hash_integer(seed_fingerprint,fragment_id);call hash_integer(seed_fingerprint,basis_generation)
    call hash_integer(seed_fingerprint,total_rows);call hash_integer(seed_fingerprint,nseed)
    call hash_int64(seed_fingerprint,row_seed)
    call hash_real_array(seed_fingerprint,dc_seed_energies)
    call hash_real_array(seed_fingerprint,dc_seed_occupations)
    if(seed_fingerprint==0_int64)seed_fingerprint=1_int64
    basis_fingerprint=initial_hash(307_int64)
    call hash_integer(basis_fingerprint,fragment_id);call hash_integer(basis_fingerprint,basis_generation)
    call hash_integer(basis_fingerprint,total_rows);call hash_integer(basis_fingerprint,nseed)
    call hash_integer(basis_fingerprint,nbuffer);call hash_integer(basis_fingerprint,nprojector)
    call hash_integer(basis_fingerprint,natom);call hash_integer(basis_fingerprint,num_iter)
    call hash_int64(basis_fingerprint,coordinator_byte_limit);call hash_int64(basis_fingerprint,row_basis)
    call hash_real(basis_fingerprint,metric_tolerance);call hash_real(basis_fingerprint,localization_tolerance)
    call hash_real_array(basis_fingerprint,real_lattice)
    call hash_real_array(basis_fingerprint,reciprocal_lattice)
    call hash_real_array(basis_fingerprint,atoms_cart)
    call hash_character_array(basis_fingerprint,atom_symbols)
    if(basis_fingerprint==0_int64)basis_fingerprint=2_int64
    call MPI_Comm_rank(comm,fragment_rank,ierr);call MPI_Comm_size(comm,fragment_size,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment row layout communicator query failed';return;endif
    local_layout_fingerprint=initial_hash(353_int64)
    call hash_integer(local_layout_fingerprint,fragment_rank)
    call hash_integer(local_layout_fingerprint,fragment_size)
    call hash_integer(local_layout_fingerprint,nlocal)
    do ierr=1,nlocal;call hash_int64(local_layout_fingerprint,grid_ids(ierr));enddo
    if(local_layout_fingerprint==0_int64)local_layout_fingerprint=4_int64
    ok=.true.;message=''
  end subroutine validate_fragment_contract

  subroutine validation_message(code,message)
    integer,intent(in)::code
    character(*),intent(out)::message
    select case(code)
    case(10);message='invalid fragment Wannier fragment/generation or iteration controls'
    case(11);message='invalid fragment Wannier seed directory'
    case(12);message='invalid fragment Wannier state-major array dimensions'
    case(13);message='invalid fragment Wannier metric or localization tolerance'
    case(14);message='fragment Wannier grid IDs or integration weights are invalid'
    case(15);message='fragment Wannier distributed candidate values are nonfinite'
    case(16);message='fragment Wannier physical seed energies or occupations are nonfinite'
    case(17);message='fragment Wannier geometry is nonfinite'
    case(18);message='fragment Wannier atom symbol is empty'
    case default;message='fragment Wannier collective validation failed'
    end select
  end subroutine validation_message

  subroutine validate_physical_seed_rank(comm,seeds,weights,tolerance,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::seeds(:,:)
    real(real64),intent(in)::weights(:),tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:),work(:)
    complex(real64)::symmetrized_value
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::threshold,scale
    integer::nseed,i,j,rank,ierr,info,status,lwork,allocation_status
    logical::allocation_ok
    character(message_length)::allocation_message
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character,intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine zheev
    end interface
    ok=.false.;message='';nseed=size(seeds,1);status=0
    allocation_status=0
    if(.not.extent_product_fits([nseed,nseed]))then
      allocation_status=1
    else
      allocate(gram(nseed,nseed),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'physical seed Gram allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    gram=(0d0,0d0)
    do j=1,nseed;do i=1,nseed
      gram(i,j)=sum(weights*conjg(seeds(i,:))*seeds(j,:))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,size(gram),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.finite_complex(gram))then
      message='physical seed metric reduction failed';return
    endif
    do j=1,nseed
      gram(j,j)=cmplx(real(gram(j,j),real64),0d0,real64)
      do i=1,j-1
        symmetrized_value=0.5d0*(gram(i,j)+conjg(gram(j,i)))
        gram(i,j)=symmetrized_value;gram(j,i)=conjg(symmetrized_value)
      enddo
    enddo
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='physical seed communicator rank query failed';return;endif
    allocation_status=0
    if(.not.extent_product_fits([nseed]))then
      allocation_status=1
    else
      allocate(eigenvalues(nseed),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'physical seed eigenvalue allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    eigenvalues=0d0
    allocation_status=0
    if(rank==0)then
      lwork=max(1,2*nseed-1)
      if(.not.extent_product_fits([lwork]).or.&
          .not.extent_product_fits([max(1,3*nseed-2)]))then
        allocation_status=1
      else
        allocate(work(lwork),rwork(max(1,3*nseed-2)),stat=allocation_status)
      endif
    endif
    call collective_allocation_status(comm,allocation_status,&
      'physical seed eigensolver workspace allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    if(rank==0)then
      call zheev('N','U',nseed,gram,nseed,eigenvalues,work,lwork,rwork,info)
      if(info/=0.or..not.all(ieee_is_finite(eigenvalues)))status=1
      if(status==0)then
        scale=max(1d0,maxval(abs(eigenvalues)));threshold=tolerance*scale
        if(any(eigenvalues< -threshold))status=2
        if(count(eigenvalues>threshold)/=nseed)status=3
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      select case(status)
      case(1);message='physical seed metric eigensolver failed'
      case(2);message='physical seed metric has a negative mode'
      case default;message='physical seed metric rank is smaller than the physical seed count'
      end select
      return
    endif
    ok=.true.;message=''
  end subroutine validate_physical_seed_rank

  subroutine classify_fragment_cache(comm,cache,fragment_id,basis_generation,nseed,&
      candidate_rank,grid_ids,seed_fingerprint,basis_fingerprint,local_layout_fingerprint,&
      energies,occupations,&
      build_required,ok,message)
    integer,intent(in)::comm,fragment_id,basis_generation,nseed,candidate_rank
    integer(int64),intent(in)::grid_ids(:),local_layout_fingerprint
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::cache
    integer(int64),intent(in)::seed_fingerprint,basis_fingerprint
    real(real64),intent(in)::energies(:),occupations(:)
    logical,intent(out)::build_required,ok
    character(*),intent(out)::message
    integer::valid_flag,minimum_flag,maximum_flag,code,global_code,ierr,retained
    integer::fragment_rank,fragment_size,nlocal
    integer(int64)::payload_hash,minimum_hash,maximum_hash,transform_hash,wannier_hash
    logical::hash_ok
    character(message_length)::hash_message

    ok=.false.;message='';valid_flag=merge(1,0,cache%valid);nlocal=size(grid_ids)
    call MPI_Comm_rank(comm,fragment_rank,ierr);call MPI_Comm_size(comm,fragment_size,ierr)
    call MPI_Allreduce(valid_flag,minimum_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(valid_flag,maximum_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_flag/=maximum_flag)then
      build_required=.false.;message='fragment Wannier cache validity disagrees across ranks';return
    endif
    build_required=valid_flag==0
    if(build_required)then;ok=.true.;return;endif

    code=0
    if(cache%receipt%fragment_id/=fragment_id)code=10
    if(code==0.and.cache%receipt%basis_generation/=basis_generation)code=20
    if(code==0.and.cache%receipt%seed_fingerprint/=seed_fingerprint)code=30
    if(code==0.and.cache%receipt%basis_fingerprint/=basis_fingerprint)code=40
    if(code==0.and.(cache%fragment_comm_rank/=fragment_rank.or.&
        cache%fragment_comm_size/=fragment_size.or.&
        cache%local_row_layout_fingerprint/=local_layout_fingerprint))code=45
    if(code==0.and..not.allocated(cache%local_grid_ids))code=45
    if(code==0)then
      if(size(cache%local_grid_ids)/=nlocal)then
        code=45
      else if(any(cache%local_grid_ids/=grid_ids))then
        code=45
      endif
    endif
    retained=cache%receipt%retained_rank
    if(code==0.and.(cache%receipt%candidate_rank/=candidate_rank.or.retained<1.or.&
        cache%receipt%setup_count/=1.or.cache%receipt%run_count/=1.or.&
        cache%receipt%transform_fingerprint==0_int64.or.&
        cache%receipt%replicated_payload_fingerprint==0_int64.or.&
        cache%receipt%distributed_wannier_fingerprint==0_int64.or.&
        .not.ieee_is_finite(cache%receipt%seed_reconstruction_defect)))code=50
    if(code==0.and.(.not.allocated(cache%wannier_values).or.&
        .not.allocated(cache%candidate_compression).or.&
        .not.allocated(cache%wannier_transform).or.&
        .not.allocated(cache%dc_seed_coefficients_in_wannier).or.&
        .not.allocated(cache%physical_dc_seed_energies).or.&
        .not.allocated(cache%physical_dc_seed_occupations)))code=50
    if(code==0)then
      if(any(shape(cache%wannier_values)/=[retained,nlocal]).or.&
          any(shape(cache%candidate_compression)/=[candidate_rank,retained]).or.&
          any(shape(cache%wannier_transform)/=[retained,retained]).or.&
          any(shape(cache%dc_seed_coefficients_in_wannier)/=[retained,nseed]).or.&
          size(cache%physical_dc_seed_energies)/=nseed.or.&
          size(cache%physical_dc_seed_occupations)/=nseed)code=50
      if(code==0)then
        if(.not.finite_complex(cache%wannier_values).or.&
            .not.finite_complex(cache%candidate_compression).or.&
            .not.finite_complex(cache%wannier_transform).or.&
            .not.finite_complex(cache%dc_seed_coefficients_in_wannier).or.&
            .not.all(ieee_is_finite(cache%physical_dc_seed_energies)).or.&
            .not.all(ieee_is_finite(cache%physical_dc_seed_occupations)))code=50
        if(.not.bitwise_real_equal(cache%physical_dc_seed_energies,energies).or.&
            .not.bitwise_real_equal(cache%physical_dc_seed_occupations,occupations))code=50
      endif
    endif
    if(code==0)then
      transform_hash=hash_complex_matrix(cache%wannier_transform,401_int64)
      if(transform_hash/=cache%receipt%transform_fingerprint)code=60
    endif
    if(code==0)then
      payload_hash=replicated_cache_integrity_hash(cache)
      if(payload_hash/=cache%receipt%replicated_payload_fingerprint)code=70
    endif
    call MPI_Allreduce(code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      select case(global_code)
      case(10);message='stale fragment-Wannier cache fragment identity'
      case(20);message='stale fragment-Wannier cache generation'
      case(30);message='stale fragment-Wannier cache seed fingerprint'
      case(40);message='stale fragment-Wannier cache basis fingerprint'
      case(45);message='stale fragment-Wannier cache local row layout'
      case(60);message='invalid fragment-Wannier cache transform integrity receipt'
      case(70);message='invalid fragment-Wannier replicated cache payload integrity receipt'
      case default;message='invalid fragment-Wannier cache payload'
      end select
      return
    endif
    call MPI_Allreduce(payload_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='fragment-Wannier replicated cache payload disagrees across ranks';return
    endif
    call distributed_wannier_integrity_hash(comm,cache%fragment_comm_rank,&
      cache%fragment_comm_size,cache%local_row_layout_fingerprint,cache%local_grid_ids,&
      cache%wannier_values,wannier_hash,hash_ok,hash_message)
    if(.not.hash_ok)then;message=hash_message;return;endif
    code=merge(0,80,wannier_hash==cache%receipt%distributed_wannier_fingerprint)
    call MPI_Allreduce(code,global_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_code/=0)then
      message='invalid fragment-Wannier distributed local cache payload integrity receipt';return
    endif
    ok=.true.;message=''
  end subroutine classify_fragment_cache

  subroutine create_fragment_namespace(comm,base_directory,fragment_id,basis_generation,&
      seed_name,ok,message)
    integer,intent(in)::comm,fragment_id,basis_generation
    character(*),intent(in)::base_directory
    character(*),intent(out)::seed_name
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(message_length)::generation_directory,base
    character(15)::fragment_component
    character(19)::generation_component
    integer::rank,ierr,status
    integer(c_int)::retcode
    interface
      subroutine posix_mkdir(dirpath,retcode)bind(C,name='posix_mkdir')
        import c_char,c_int
        character(kind=c_char),intent(in)::dirpath
        integer(c_int),intent(out)::retcode
      end subroutine posix_mkdir
    end interface
    ok=.false.;message='';seed_name='';status=0
    base=trimmed_directory(base_directory)
    write(fragment_component,'("fragment-",i6.6)')fragment_id
    write(generation_component,'("generation-",i8.8)')basis_generation
    generation_directory=trim(base)//'/'//trim(fragment_component)//'/'//trim(generation_component)
    seed_name=trim(generation_directory)//'/construction_wannier'
    if(len_trim(generation_directory)>=message_length-2.or.len_trim(seed_name)>=len(seed_name))status=1
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    if(rank==0.and.status==0)then
      call posix_mkdir(trim(generation_directory)//c_null_char,retcode)
      if(retcode/=0_c_int)status=2
      if(status==0)then
        call posix_mkdir(trim(seed_name)//c_null_char,retcode)
        if(retcode/=0_c_int)status=3
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      if(status==2)then
        message='fragment Wannier generation directory already exists or cannot be created'
      else if(status==3)then
        message='fragment Wannier construction namespace cannot be created'
      else
        message='fragment Wannier artifact namespace is too long'
      endif
      return
    endif
    ok=.true.;message=''
  end subroutine create_fragment_namespace

  subroutine construct_fragment_wannier(comm,fragment_id,basis_generation,seed_name,&
      grid_ids,grid_weights,dc_seed_values,dc_seed_energies,dc_seed_occupations,&
      buffer_values,projector_values,metric_tolerance,real_lattice,reciprocal_lattice,&
      inverse_lattice,atom_symbols,atoms_cart,fractional,num_iter,localization_tolerance,&
      coordinator_byte_limit,seed_fingerprint,basis_fingerprint,local_layout_fingerprint,&
      working_cache,ok,message)
    integer,intent(in)::comm,fragment_id,basis_generation,num_iter
    character(*),intent(in)::seed_name
    integer(int64),intent(in)::grid_ids(:),coordinator_byte_limit,seed_fingerprint,&
      basis_fingerprint,local_layout_fingerprint
    real(real64),intent(in)::grid_weights(:),dc_seed_energies(:),dc_seed_occupations(:),&
      metric_tolerance,real_lattice(3,3),reciprocal_lattice(3,3),inverse_lattice(3,3),&
      atoms_cart(:,:),fractional(:,:),localization_tolerance
    complex(real64),intent(in)::dc_seed_values(:,:),buffer_values(:,:),projector_values(:,:)
    character(*),intent(in)::atom_symbols(:)
    type(s_dg_hybrid_fragment_wannier_cache),intent(out)::working_cache
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::raw(:,:),compression(:,:),retained_values(:,:),&
      m_matrix(:,:,:),a_matrix(:,:),transform(:,:),seed_coefficients(:,:)
    real(real64),allocatable::eigenvalues(:),centers(:,:),spreads(:)
    integer,allocatable::nncell(:,:)
    real(real64)::spread(3),seed_defect,density_defect,orthogonality_defect,&
      seed_certificate_tolerance
    integer(int64)::coordinator_bytes,workspace_peak_bytes,transform_fingerprint,&
      minimum_fingerprint,maximum_fingerprint,payload_fingerprint,wannier_fingerprint
    integer::nseed,nbuffer,nprojector,ncandidate,nlocal,retained_rank,nntot
    integer::ierr,rank,fragment_size,allocation_status
    character(message_length)::adapter_message
    logical::step_ok,allocation_ok

    ok=.false.;message='';adapter_message='';working_cache%valid=.false.
    nseed=size(dc_seed_values,1);nbuffer=size(buffer_values,1)
    nprojector=size(projector_values,1);ncandidate=nseed+nbuffer+nprojector
    nlocal=size(grid_ids)
    allocation_status=0
    if(.not.extent_product_fits([ncandidate,nlocal]))then
      allocation_status=1
    else
      allocate(raw(ncandidate,nlocal),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment raw candidate allocation failed',allocation_ok,adapter_message)
    if(.not.allocation_ok)then;message=adapter_message;return;endif
    raw(1:nseed,:)=dc_seed_values
    if(nbuffer>0)raw(nseed+1:nseed+nbuffer,:)=buffer_values
    if(nprojector>0)raw(nseed+nbuffer+1:ncandidate,:)=projector_values
    call metric_compress_candidates(comm,raw,grid_weights,metric_tolerance,compression,&
      retained_values,retained_rank,orthogonality_defect,step_ok,adapter_message)
    if(.not.step_ok)then;message='fragment candidate metric failure: '//trim(adapter_message);return;endif
    deallocate(raw)
    seed_certificate_tolerance=max(10d0*metric_tolerance,&
      256d0*epsilon(1d0)*real(max(1,ncandidate),real64))
    call certify_seed_span(comm,dc_seed_values,retained_values,grid_weights,&
      seed_certificate_tolerance,seed_coefficients,seed_defect,&
      density_defect,dc_seed_occupations,step_ok,adapter_message)
    if(.not.step_ok)then;message='fragment DC seed span failure: '//trim(adapter_message);return;endif

    call setup_dg_w90_gamma_library(comm,trim(seed_name),real_lattice,reciprocal_lattice,&
      atom_symbols,atoms_cart,retained_rank,retained_rank,num_iter,'spectral',&
      DG_W90_UNCONSTRAINED,nntot,nncell,step_ok,adapter_message)
    if(.not.step_ok)then
      call fragment_error(fragment_id,'Wannier90 setup failed',adapter_message,message);return
    endif
    call assemble_dg_w90_gamma_matrices(comm,retained_values,retained_values,grid_weights,&
      fractional,nncell,coordinator_byte_limit,m_matrix,a_matrix,coordinator_bytes,&
      workspace_peak_bytes,step_ok,adapter_message)
    if(.not.step_ok)then
      call fragment_error(fragment_id,'Wannier90 matrix assembly failed',adapter_message,message);return
    endif
    allocation_status=0
    if(.not.extent_product_fits([retained_rank]))then
      allocation_status=1
    else
      allocate(eigenvalues(retained_rank),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment auxiliary eigenvalue allocation failed',allocation_ok,adapter_message)
    if(.not.allocation_ok)then;message=adapter_message;return;endif
    eigenvalues=0d0
    call run_dg_w90_gamma_library(comm,trim(seed_name),real_lattice,reciprocal_lattice,&
      atom_symbols,atoms_cart,m_matrix,a_matrix,eigenvalues,1d100,localization_tolerance,&
      num_iter,transform,centers,spreads,spread,step_ok,adapter_message,&
      require_nonincreasing_spread=.false.)
    if(.not.step_ok)then
      call fragment_error(fragment_id,'Wannier90 run failed',adapter_message,message);return
    endif
    centers=matmul(inverse_lattice,centers)
    call apply_dg_w90_gamma_transform(comm,grid_ids,retained_values,transform=transform,&
      centers=centers,tolerance=localization_tolerance,ok=step_ok,message=adapter_message,&
      spreads=spreads)
    if(.not.step_ok)then
      call fragment_error(fragment_id,'Wannier90 transform application failed',adapter_message,message);return
    endif

    call certify_seed_span(comm,dc_seed_values,retained_values,grid_weights,&
      seed_certificate_tolerance,seed_coefficients,seed_defect,&
      density_defect,dc_seed_occupations,step_ok,adapter_message)
    if(.not.step_ok)then
      message='localized fragment seed certificate failed: '//trim(adapter_message);return
    endif
    transform_fingerprint=hash_complex_matrix(transform,401_int64)
    call MPI_Allreduce(transform_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(transform_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='fragment Wannier transform differs across fragment ranks';return
    endif
    if(transform_fingerprint==0_int64)transform_fingerprint=3_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then
      message='fragment Wannier communicator rank query failed';return
    endif
    call MPI_Comm_size(comm,fragment_size,ierr);if(ierr/=MPI_SUCCESS)then
      message='fragment Wannier communicator size query failed';return
    endif

    allocation_status=0
    if(.not.extent_product_fits([nlocal]).or..not.extent_product_fits([nseed]))then
      allocation_status=1
    else
      allocate(working_cache%local_grid_ids(nlocal),&
        working_cache%physical_dc_seed_energies(nseed),&
        working_cache%physical_dc_seed_occupations(nseed),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment published cache metadata allocation failed',allocation_ok,adapter_message)
    if(.not.allocation_ok)then;message=adapter_message;return;endif

    working_cache%receipt%fragment_id=fragment_id
    working_cache%receipt%basis_generation=basis_generation
    working_cache%receipt%candidate_rank=ncandidate
    working_cache%receipt%retained_rank=retained_rank
    working_cache%receipt%setup_count=1
    working_cache%receipt%run_count=1
    working_cache%receipt%seed_fingerprint=seed_fingerprint
    working_cache%receipt%basis_fingerprint=basis_fingerprint
    working_cache%receipt%transform_fingerprint=transform_fingerprint
    working_cache%receipt%seed_reconstruction_defect=seed_defect
    working_cache%fragment_comm_rank=rank
    working_cache%fragment_comm_size=fragment_size
    working_cache%local_row_layout_fingerprint=local_layout_fingerprint
    working_cache%local_grid_ids=grid_ids
    call move_alloc(retained_values,working_cache%wannier_values)
    call move_alloc(compression,working_cache%candidate_compression)
    call move_alloc(transform,working_cache%wannier_transform)
    call move_alloc(seed_coefficients,working_cache%dc_seed_coefficients_in_wannier)
    working_cache%physical_dc_seed_energies=dc_seed_energies
    working_cache%physical_dc_seed_occupations=dc_seed_occupations
    payload_fingerprint=replicated_cache_integrity_hash(working_cache)
    call MPI_Allreduce(payload_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(payload_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='fragment published replicated cache payload differs across ranks';return
    endif
    if(payload_fingerprint==0_int64)payload_fingerprint=5_int64
    working_cache%receipt%replicated_payload_fingerprint=payload_fingerprint
    call distributed_wannier_integrity_hash(comm,rank,fragment_size,&
      local_layout_fingerprint,working_cache%local_grid_ids,working_cache%wannier_values,&
      wannier_fingerprint,step_ok,adapter_message)
    if(.not.step_ok)then;message=adapter_message;return;endif
    if(wannier_fingerprint==0_int64)wannier_fingerprint=6_int64
    working_cache%receipt%distributed_wannier_fingerprint=wannier_fingerprint
    working_cache%valid=.true.
    ok=.true.;message=''
  end subroutine construct_fragment_wannier

  subroutine metric_compress_candidates(comm,raw,weights,tolerance,compression,retained,&
      retained_rank,orthogonality_defect,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::raw(:,:)
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),allocatable,intent(out)::compression(:,:),retained(:,:)
    integer,intent(out)::retained_rank
    real(real64),intent(out)::orthogonality_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:),work(:),retained_gram(:,:)
    complex(real64)::symmetrized_value
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::threshold,scale
    integer::n,i,j,k,p,rank,ierr,info,status,lwork,allocation_status
    logical::allocation_ok
    character(message_length)::allocation_message
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character,intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine zheev
    end interface
    ok=.false.;message='';retained_rank=0;orthogonality_defect=huge(1d0)
    n=size(raw,1);allocation_status=0
    if(.not.extent_product_fits([n,n]))then
      allocation_status=1
    else
      allocate(gram(n,n),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'candidate Gram allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    gram=(0d0,0d0)
    do j=1,n;do i=1,n
      gram(i,j)=sum(weights*conjg(raw(i,:))*raw(j,:))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.finite_complex(gram))then
      message='distributed candidate Gram reduction failed';return
    endif
    do j=1,n
      gram(j,j)=cmplx(real(gram(j,j),real64),0d0,real64)
      do i=1,j-1
        symmetrized_value=0.5d0*(gram(i,j)+conjg(gram(j,i)))
        gram(i,j)=symmetrized_value;gram(j,i)=conjg(symmetrized_value)
      enddo
    enddo
    call MPI_Comm_rank(comm,rank,ierr);status=0
    if(ierr/=MPI_SUCCESS)then;message='candidate communicator rank query failed';return;endif
    allocation_status=0
    if(.not.extent_product_fits([n]))then
      allocation_status=1
    else
      allocate(eigenvalues(n),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'candidate Gram eigenvalue allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    eigenvalues=0d0;allocation_status=0
    if(rank==0)then
      lwork=max(1,2*n-1)
      if(.not.extent_product_fits([lwork]).or.&
          .not.extent_product_fits([max(1,3*n-2)]))then
        allocation_status=1
      else
        allocate(work(lwork),rwork(max(1,3*n-2)),stat=allocation_status)
      endif
    endif
    call collective_allocation_status(comm,allocation_status,&
      'candidate Gram eigensolver workspace allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    if(rank==0)then
      call zheev('V','U',n,gram,n,eigenvalues,work,lwork,rwork,info)
      if(info/=0.or..not.all(ieee_is_finite(eigenvalues)).or..not.finite_complex(gram))status=1
      if(status==0)then
        scale=max(1d0,maxval(abs(eigenvalues)));threshold=tolerance*scale
        if(any(eigenvalues< -threshold))status=2
        retained_rank=count(eigenvalues>threshold)
        if(retained_rank<1)status=3
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(retained_rank,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      select case(status)
      case(1);message='candidate Gram eigensolver failed'
      case(2);message='candidate Gram has a negative metric mode'
      case default;message='candidate Gram has no retained metric modes'
      end select
      return
    endif
    allocation_status=0
    if(.not.extent_product_fits([n,retained_rank]))then
      allocation_status=1
    else
      allocate(compression(n,retained_rank),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'candidate compression allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    compression=(0d0,0d0)
    if(rank==0)then
      k=0
      do i=1,n
        if(eigenvalues(i)<=tolerance*max(1d0,maxval(abs(eigenvalues))))cycle
        k=k+1;compression(:,k)=gram(:,i)/sqrt(eigenvalues(i))
      enddo
    endif
    call MPI_Bcast(compression,size(compression),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='candidate compression broadcast failed';return;endif
    allocation_status=0
    if(.not.extent_product_fits([retained_rank,size(raw,2)]))then
      allocation_status=1
    else
      allocate(retained(retained_rank,size(raw,2)),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'retained candidate allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    retained=(0d0,0d0)
    do p=1,size(raw,2);do k=1,retained_rank;do i=1,n
      retained(k,p)=retained(k,p)+compression(i,k)*raw(i,p)
    enddo;enddo;enddo
    allocation_status=0
    if(.not.extent_product_fits([retained_rank,retained_rank]))then
      allocation_status=1
    else
      allocate(retained_gram(retained_rank,retained_rank),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'retained candidate Gram allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    retained_gram=(0d0,0d0)
    do j=1,retained_rank;do i=1,retained_rank
      retained_gram(i,j)=sum(weights*conjg(retained(i,:))*retained(j,:))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,retained_gram,size(retained_gram),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.finite_complex(retained_gram))then
      message='retained candidate metric reduction failed';return
    endif
    do i=1,retained_rank;retained_gram(i,i)=retained_gram(i,i)-1d0;enddo
    orthogonality_defect=maxval(abs(retained_gram))
    if(.not.ieee_is_finite(orthogonality_defect).or.&
        orthogonality_defect>max(100d0*epsilon(1d0)*real(n,real64),10d0*tolerance))then
      message='retained candidates are not metric orthonormal';return
    endif
    ok=.true.;message=''
  end subroutine metric_compress_candidates

  subroutine certify_seed_span(comm,seeds,basis,weights,tolerance,coefficients,defect,&
      density_defect,occupations,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::seeds(:,:),basis(:,:)
    real(real64),intent(in)::weights(:),tolerance,occupations(:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::defect,density_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::reconstructed(:,:)
    real(real64),allocatable::residual_norm2(:)
    real(real64)::local_density_scale,density_scale
    integer::i,j,p,ierr,allocation_status
    logical::allocation_ok
    character(message_length)::allocation_message
    ok=.false.;message='';defect=huge(1d0);density_defect=huge(1d0)
    allocation_status=0
    if(.not.extent_product_fits([size(basis,1),size(seeds,1)]))then
      allocation_status=1
    else
      allocate(coefficients(size(basis,1),size(seeds,1)),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'seed coefficient allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    coefficients=(0d0,0d0)
    do j=1,size(seeds,1);do i=1,size(basis,1)
      coefficients(i,j)=sum(weights*conjg(basis(i,:))*seeds(j,:))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,coefficients,size(coefficients),MPI_DOUBLE_COMPLEX,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.finite_complex(coefficients))then
      message='seed coefficient reduction failed';return
    endif
    allocation_status=0
    if(.not.extent_product_fits([size(seeds,1),size(seeds,2)]))then
      allocation_status=1
    else
      allocate(reconstructed(size(seeds,1),size(seeds,2)),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'seed reconstruction allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    reconstructed=(0d0,0d0)
    do p=1,size(seeds,2);do j=1,size(seeds,1);do i=1,size(basis,1)
      reconstructed(j,p)=reconstructed(j,p)+coefficients(i,j)*basis(i,p)
    enddo;enddo;enddo
    allocation_status=0
    if(.not.extent_product_fits([size(seeds,1)]))then
      allocation_status=1
    else
      allocate(residual_norm2(size(seeds,1)),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'seed residual allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    residual_norm2=0d0
    do i=1,size(seeds,1)
      do p=1,size(seeds,2)
        residual_norm2(i)=residual_norm2(i)+weights(p)*&
          abs(reconstructed(i,p)-seeds(i,p))**2
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,residual_norm2,size(residual_norm2),MPI_DOUBLE_PRECISION,&
      MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(residual_norm2)))then
      message='seed reconstruction reduction failed';return
    endif
    defect=sqrt(max(0d0,maxval(residual_norm2)))
    density_defect=0d0;local_density_scale=0d0
    do j=1,size(seeds,2)
      density_defect=max(density_defect,abs(sum(occupations*abs(seeds(:,j))**2)-&
        sum(occupations*abs(reconstructed(:,j))**2)))
      local_density_scale=max(local_density_scale,abs(sum(occupations*abs(seeds(:,j))**2)))
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,density_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_density_scale,density_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(density_defect))then
      message='occupation-weighted density certificate reduction failed';return
    endif
    if(defect>tolerance.or.density_defect>tolerance*max(1d0,density_scale))then
      message='DC seeds, occupied projector, or density are not preserved';return
    endif
    ok=.true.;message=''
  end subroutine certify_seed_span

  subroutine validate_unique_grid_ids(comm,grid_ids,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::grid_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::total_64
    integer::rank,nrank,nlocal,total,i,ierr,status,allocation_status,local_failure,global_failure
    logical::allocation_ok
    character(message_length)::allocation_message
    ok=.false.;message='';status=0;nlocal=size(grid_ids)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then
      message='fragment grid communicator rank query failed';return
    endif
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS.or.nrank<1)then
      message='fragment grid communicator size query failed';return
    endif
    allocation_status=0
    if(rank==0)then
      if(.not.extent_product_fits([nrank]))then
        allocation_status=1
      else
        allocate(counts(nrank),displacements(nrank),stat=allocation_status)
      endif
    else
      allocate(counts(1),displacements(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment grid ownership metadata allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    counts=0;displacements=0
    call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm,ierr)
    local_failure=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)then
      message='fragment grid ownership count gather failed';return
    endif
    total=0
    if(rank==0)then
      total_64=0_int64
      do i=1,nrank
        total_64=total_64+int(counts(i),int64)
      enddo
      if(total_64<1_int64.or.total_64>int(huge(0),int64))then
        status=1
      else
        total=int(total_64)
        do i=2,nrank;displacements(i)=displacements(i-1)+counts(i-1);enddo
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      message='fragment grid ID ownership extent is invalid';return
    endif
    allocation_status=0
    if(rank==0)then
      if(.not.extent_product_fits([total]))then
        allocation_status=1
      else
        allocate(all_ids(max(1,total)),stat=allocation_status)
      endif
    else
      allocate(all_ids(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment grid ID catalog allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    all_ids=0_int64
    call MPI_Gatherv(grid_ids,nlocal,MPI_INTEGER8,all_ids,counts,displacements,MPI_INTEGER8,&
      0,comm,ierr)
    local_failure=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)then
      message='fragment grid ID gather failed';return
    endif
    if(rank==0)then
      call sort_int64(all_ids(1:total))
      if(any(all_ids(1:total)<=0_int64))status=1
      do i=2,total
        if(all_ids(i)==all_ids(i-1))status=2
      enddo
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      if(status==2)then
        message='fragment grid IDs are duplicated inside the fragment communicator'
      else
        message='fragment grid ID ownership validation failed'
      endif
      return
    endif
    ok=.true.;message=''
  end subroutine validate_unique_grid_ids

  subroutine distributed_row_hashes(comm,ids,weights,seeds,buffers,projectors,fractional,&
      seed_hash,basis_hash,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::ids(:)
    real(real64),intent(in)::weights(:),fractional(:,:)
    complex(real64),intent(in)::seeds(:,:),buffers(:,:),projectors(:,:)
    integer(int64),intent(out)::seed_hash,basis_hash
    logical,intent(out)::ok
    integer(int64),allocatable::local_seed_rows(:),local_basis_rows(:),all_ids(:),&
      all_seed_rows(:),all_basis_rows(:)
    integer,allocatable::counts(:),displacements(:)
    integer(int64)::row_hash,total_64
    integer::p,i,ierr,rank,nrank,total,status,allocation_status,local_failure,global_failure
    logical::allocation_ok
    character(message_length)::allocation_message
    ok=.false.;seed_hash=0_int64;basis_hash=0_int64;status=0;total=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS.or.nrank<1)return
    allocation_status=0
    if(.not.extent_product_fits([size(ids)]))then
      allocation_status=1
    else
      allocate(local_seed_rows(size(ids)),local_basis_rows(size(ids)),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'distributed row hash allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)return
    do p=1,size(ids)
      row_hash=initial_hash(503_int64);call hash_int64(row_hash,ids(p));call hash_real(row_hash,weights(p))
      do i=1,size(seeds,1);call hash_complex(row_hash,seeds(i,p));enddo
      local_seed_rows(p)=row_hash
      row_hash=initial_hash(601_int64);call hash_int64(row_hash,ids(p));call hash_real(row_hash,weights(p))
      do i=1,3;call hash_real(row_hash,fractional(i,p));enddo
      do i=1,size(buffers,1);call hash_complex(row_hash,buffers(i,p));enddo
      do i=1,size(projectors,1);call hash_complex(row_hash,projectors(i,p));enddo
      local_basis_rows(p)=row_hash
    enddo
    allocation_status=0
    if(rank==0)then
      if(.not.extent_product_fits([nrank]))then
        allocation_status=1
      else
        allocate(counts(nrank),displacements(nrank),stat=allocation_status)
      endif
    else
      allocate(counts(1),displacements(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'distributed row hash metadata allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)return
    counts=0;displacements=0
    call MPI_Gather(size(ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm,ierr)
    local_failure=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)return
    if(rank==0)then
      total_64=0_int64
      do i=1,nrank;total_64=total_64+int(counts(i),int64);enddo
      if(total_64<1_int64.or.total_64>int(huge(0),int64))then
        status=1
      else
        total=int(total_64)
        do i=2,nrank;displacements(i)=displacements(i-1)+counts(i-1);enddo
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)return
    allocation_status=0
    if(rank==0)then
      if(.not.extent_product_fits([total]))then
        allocation_status=1
      else
        allocate(all_ids(max(1,total)),all_seed_rows(max(1,total)),&
          all_basis_rows(max(1,total)),stat=allocation_status)
      endif
    else
      allocate(all_ids(1),all_seed_rows(1),all_basis_rows(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'distributed row hash catalog allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)return
    all_ids=0_int64;all_seed_rows=0_int64;all_basis_rows=0_int64;status=0
    call MPI_Gatherv(ids,size(ids),MPI_INTEGER8,all_ids,counts,displacements,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Gatherv(local_seed_rows,size(ids),MPI_INTEGER8,all_seed_rows,counts,&
      displacements,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Gatherv(local_basis_rows,size(ids),MPI_INTEGER8,all_basis_rows,counts,&
      displacements,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)status=1
    call MPI_Allreduce(status,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)return
    if(rank==0.and.status==0)then
      call sort_row_hashes(all_ids(1:total),all_seed_rows(1:total),all_basis_rows(1:total))
      seed_hash=initial_hash(619_int64);basis_hash=initial_hash(631_int64)
      do i=1,total
        call hash_int64(seed_hash,all_ids(i));call hash_int64(seed_hash,all_seed_rows(i))
        call hash_int64(basis_hash,all_ids(i));call hash_int64(basis_hash,all_basis_rows(i))
      enddo
    endif
    call MPI_Bcast(seed_hash,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(basis_hash,1,MPI_INTEGER8,0,comm,ierr)
    ok=ierr==MPI_SUCCESS
  end subroutine distributed_row_hashes

  subroutine distributed_wannier_integrity_hash(comm,cached_rank,cached_size,&
      local_layout_fingerprint,ids,values,fingerprint,ok,message)
    integer,intent(in)::comm,cached_rank,cached_size
    integer(int64),intent(in)::local_layout_fingerprint,ids(:)
    complex(real64),intent(in)::values(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),allocatable::rank_hashes(:)
    integer(int64)::local_hash
    integer::rank,nrank,i,p,ierr,status,global_status,allocation_status
    logical::allocation_ok
    character(message_length)::allocation_message

    ok=.false.;message='';fingerprint=0_int64;status=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then
      message='fragment-Wannier integrity communicator rank query failed';return
    endif
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS.or.nrank<1)then
      message='fragment-Wannier integrity communicator size query failed';return
    endif
    if(cached_rank/=rank.or.cached_size/=nrank.or.size(values,2)/=size(ids).or.&
        local_layout_fingerprint==0_int64.or.any(ids<=0_int64).or.&
        .not.finite_complex(values))status=1
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='invalid fragment-Wannier local payload for integrity hashing';return
    endif

    local_hash=initial_hash(733_int64)
    call hash_integer(local_hash,rank);call hash_integer(local_hash,nrank)
    call hash_int64(local_hash,local_layout_fingerprint)
    call hash_integer(local_hash,size(ids));call hash_integer(local_hash,size(values,1))
    do p=1,size(ids)
      call hash_int64(local_hash,ids(p))
      do i=1,size(values,1);call hash_complex(local_hash,values(i,p));enddo
    enddo

    allocation_status=0
    if(rank==0)then
      if(.not.extent_product_fits([nrank]))then
        allocation_status=1
      else
        allocate(rank_hashes(nrank),stat=allocation_status)
      endif
    else
      allocate(rank_hashes(1),stat=allocation_status)
    endif
    call collective_allocation_status(comm,allocation_status,&
      'fragment-Wannier integrity catalog allocation failed',allocation_ok,allocation_message)
    if(.not.allocation_ok)then;message=allocation_message;return;endif
    rank_hashes=0_int64
    call MPI_Gather(local_hash,1,MPI_INTEGER8,rank_hashes,1,MPI_INTEGER8,0,comm,ierr)
    status=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment-Wannier integrity hash gather failed';return
    endif
    if(rank==0)then
      fingerprint=initial_hash(739_int64);call hash_integer(fingerprint,nrank)
      do i=1,nrank;call hash_int64(fingerprint,rank_hashes(i));enddo
      if(fingerprint==0_int64)fingerprint=6_int64
    endif
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='fragment-Wannier integrity hash broadcast failed';return
    endif
    ok=.true.;message=''
  end subroutine distributed_wannier_integrity_hash

  subroutine collective_allocation_status(comm,local_status,context,ok,message)
    integer,intent(in)::comm,local_status
    character(*),intent(in)::context
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::local_failure,global_failure,ierr
    local_failure=merge(0,1,local_status==0)
    global_failure=1
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_failure==0
    if(ok)then
      message=''
    else if(ierr/=MPI_SUCCESS)then
      message='fragment Wannier allocation-status reduction failed'
    else
      message=context
    endif
  end subroutine collective_allocation_status

  logical function extent_product_fits(extents)
    integer,intent(in)::extents(:)
    integer(int64)::product,limit,extent
    integer::i
    extent_product_fits=.false.;product=1_int64;limit=int(huge(0),int64)
    do i=1,size(extents)
      extent=int(extents(i),int64)
      if(extent<0_int64)return
      if(extent==0_int64)then
        product=0_int64
      else if(product>0_int64)then
        if(extent>limit/product)return
        product=product*extent
      endif
    enddo
    extent_product_fits=.true.
  end function extent_product_fits

  subroutine canonical_total_status(comm,local_ok,local_message,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    character(*),intent(in)::local_message
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(message_length)::canonical_message
    integer::rank,ierr,failing_rank,local_failing_rank
    call MPI_Comm_rank(comm,rank,ierr)
    local_failing_rank=merge(huge(0),rank,local_ok)
    call MPI_Allreduce(local_failing_rank,failing_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='fragment Wannier total communicator reduction failed';return;endif
    if(failing_rank==huge(0))then;ok=.true.;message='';return;endif
    canonical_message='';if(rank==failing_rank)canonical_message=local_message
    call MPI_Bcast(canonical_message,len(canonical_message),MPI_CHARACTER,failing_rank,comm,ierr)
    ok=.false.
    if(ierr/=MPI_SUCCESS)then
      message='fragment Wannier total failure-message broadcast failed'
    else
      message=canonical_message
    endif
  end subroutine canonical_total_status

  subroutine move_fragment_cache(source,destination)
    type(s_dg_hybrid_fragment_wannier_cache),intent(inout)::source,destination
    logical::source_valid
    source_valid=source%valid;destination%valid=.false.
    destination%receipt=source%receipt
    destination%fragment_comm_rank=source%fragment_comm_rank
    destination%fragment_comm_size=source%fragment_comm_size
    destination%local_row_layout_fingerprint=source%local_row_layout_fingerprint
    call move_alloc(source%local_grid_ids,destination%local_grid_ids)
    call move_alloc(source%wannier_values,destination%wannier_values)
    call move_alloc(source%candidate_compression,destination%candidate_compression)
    call move_alloc(source%wannier_transform,destination%wannier_transform)
    call move_alloc(source%dc_seed_coefficients_in_wannier,&
      destination%dc_seed_coefficients_in_wannier)
    call move_alloc(source%physical_dc_seed_energies,destination%physical_dc_seed_energies)
    call move_alloc(source%physical_dc_seed_occupations,destination%physical_dc_seed_occupations)
    destination%valid=source_valid;source%valid=.false.
  end subroutine move_fragment_cache

  subroutine fragment_error(fragment_id,context,detail,message)
    integer,intent(in)::fragment_id
    character(*),intent(in)::context,detail
    character(*),intent(out)::message
    character(24)::fragment_label
    write(fragment_label,'("fragment ",i0)')fragment_id
    message=trim(fragment_label)//' '//trim(context)//': '//trim(detail)
  end subroutine fragment_error

  subroutine invert_lattice(lattice,inverse,determinant,ok)
    real(real64),intent(in)::lattice(3,3)
    real(real64),intent(out)::inverse(3,3),determinant
    logical,intent(out)::ok
    real(real64)::scale
    determinant=lattice(1,1)*(lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2))-&
      lattice(1,2)*(lattice(2,1)*lattice(3,3)-lattice(2,3)*lattice(3,1))+&
      lattice(1,3)*(lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1))
    scale=max(1d0,maxval(abs(lattice)))
    ok=ieee_is_finite(determinant).and.abs(determinant)>epsilon(1d0)*scale**3
    inverse=0d0;if(.not.ok)return
    inverse(1,1)=(lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2))/determinant
    inverse(1,2)=(lattice(1,3)*lattice(3,2)-lattice(1,2)*lattice(3,3))/determinant
    inverse(1,3)=(lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2))/determinant
    inverse(2,1)=(lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3))/determinant
    inverse(2,2)=(lattice(1,1)*lattice(3,3)-lattice(1,3)*lattice(3,1))/determinant
    inverse(2,3)=(lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3))/determinant
    inverse(3,1)=(lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1))/determinant
    inverse(3,2)=(lattice(1,2)*lattice(3,1)-lattice(1,1)*lattice(3,2))/determinant
    inverse(3,3)=(lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1))/determinant
    ok=all(ieee_is_finite(inverse))
  end subroutine invert_lattice

  subroutine sort_int64(values)
    integer(int64),intent(inout)::values(:)
    integer::start,finish
    integer(int64)::temporary
    do start=size(values)/2,1,-1
      call sift_down_int64(values,start,size(values))
    enddo
    do finish=size(values),2,-1
      temporary=values(1);values(1)=values(finish);values(finish)=temporary
      call sift_down_int64(values,1,finish-1)
    enddo
  end subroutine sort_int64

  subroutine sift_down_int64(values,start,finish)
    integer(int64),intent(inout)::values(:)
    integer,intent(in)::start,finish
    integer::root,child
    integer(int64)::temporary
    root=start
    do while(2*root<=finish)
      child=2*root
      if(child<finish)then
        if(values(child)<values(child+1))child=child+1
      endif
      if(values(root)>=values(child))exit
      temporary=values(root);values(root)=values(child);values(child)=temporary
      root=child
    enddo
  end subroutine sift_down_int64

  subroutine sort_row_hashes(ids,seed_rows,basis_rows)
    integer(int64),intent(inout)::ids(:),seed_rows(:),basis_rows(:)
    integer::start,finish
    do start=size(ids)/2,1,-1
      call sift_down_row_hashes(ids,seed_rows,basis_rows,start,size(ids))
    enddo
    do finish=size(ids),2,-1
      call swap_row_hashes(ids,seed_rows,basis_rows,1,finish)
      call sift_down_row_hashes(ids,seed_rows,basis_rows,1,finish-1)
    enddo
  end subroutine sort_row_hashes

  subroutine sift_down_row_hashes(ids,seed_rows,basis_rows,start,finish)
    integer(int64),intent(inout)::ids(:),seed_rows(:),basis_rows(:)
    integer,intent(in)::start,finish
    integer::root,child
    root=start
    do while(2*root<=finish)
      child=2*root
      if(child<finish)then
        if(ids(child)<ids(child+1))child=child+1
      endif
      if(ids(root)>=ids(child))exit
      call swap_row_hashes(ids,seed_rows,basis_rows,root,child);root=child
    enddo
  end subroutine sift_down_row_hashes

  subroutine swap_row_hashes(ids,seed_rows,basis_rows,left,right)
    integer(int64),intent(inout)::ids(:),seed_rows(:),basis_rows(:)
    integer,intent(in)::left,right
    integer(int64)::temporary
    temporary=ids(left);ids(left)=ids(right);ids(right)=temporary
    temporary=seed_rows(left);seed_rows(left)=seed_rows(right);seed_rows(right)=temporary
    temporary=basis_rows(left);basis_rows(left)=basis_rows(right);basis_rows(right)=temporary
  end subroutine swap_row_hashes

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    integer::i,j
    finite_complex=.false.
    do j=1,size(values,2);do i=1,size(values,1)
      if(.not.ieee_is_finite(real(values(i,j),real64)).or.&
          .not.ieee_is_finite(aimag(values(i,j))))return
    enddo;enddo
    finite_complex=.true.
  end function finite_complex

  logical function bitwise_real_equal(left,right)
    real(real64),intent(in)::left(:),right(:)
    integer(int64)::left_bits,right_bits
    integer::i
    bitwise_real_equal=.false.;if(size(left)/=size(right))return
    do i=1,size(left)
      left_bits=transfer(left(i),left_bits);right_bits=transfer(right(i),right_bits)
      if(left_bits/=right_bits)return
    enddo
    bitwise_real_equal=.true.
  end function bitwise_real_equal

  character(message_length) function trimmed_directory(value)
    character(*),intent(in)::value
    integer::last
    trimmed_directory='';last=len_trim(value)
    do while(last>1.and.value(last:last)=='/');last=last-1;enddo
    if(last>0)trimmed_directory=value(1:last)
  end function trimmed_directory

  integer(int64) function initial_hash(tag)
    integer(int64),intent(in)::tag
    initial_hash=fnv1a_word(int(z'CBF29CE484222325',int64),tag)
  end function initial_hash

  subroutine hash_int64(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    hash=fnv1a_word(hash,value)
  end subroutine hash_int64

  integer(int64) function fnv1a_word(hash,value)
    integer(int64),intent(in)::hash,value
    integer(int64)::byte_value
    integer::byte
    fnv1a_word=hash
    do byte=0,7
      byte_value=iand(shiftr(value,8*byte),int(z'FF',int64))
      fnv1a_word=fnv_prime_multiply(ieor(fnv1a_word,byte_value))
    enddo
  end function fnv1a_word

  integer(int64) function fnv_prime_multiply(value)
    integer(int64),intent(in)::value
    fnv_prime_multiply=0_int64
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,value)
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,1))
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,4))
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,5))
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,7))
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,8))
    fnv_prime_multiply=add_mod64(fnv_prime_multiply,shiftl(value,40))
  end function fnv_prime_multiply

  integer(int64) function add_mod64(left,right)
    integer(int64),intent(in)::left,right
    integer(int64)::sum_bits,carry_bits,next_sum
    integer::iteration
    sum_bits=left;carry_bits=right
    do iteration=1,64
      if(carry_bits==0_int64)exit
      next_sum=ieor(sum_bits,carry_bits)
      carry_bits=shiftl(iand(sum_bits,carry_bits),1)
      sum_bits=next_sum
    enddo
    add_mod64=sum_bits
  end function add_mod64

  subroutine hash_integer(hash,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::value
    call hash_int64(hash,int(value,int64))
  end subroutine hash_integer

  subroutine hash_real(hash,value)
    integer(int64),intent(inout)::hash
    real(real64),intent(in)::value
    integer(int64)::bits
    bits=transfer(value,bits);call hash_int64(hash,bits)
  end subroutine hash_real

  subroutine hash_complex(hash,value)
    integer(int64),intent(inout)::hash
    complex(real64),intent(in)::value
    call hash_real(hash,real(value,real64));call hash_real(hash,aimag(value))
  end subroutine hash_complex

  subroutine hash_real_vector(hash,values)
    integer(int64),intent(inout)::hash
    real(real64),intent(in)::values(:)
    integer::i
    do i=1,size(values);call hash_real(hash,values(i));enddo
  end subroutine hash_real_vector

  subroutine hash_real_matrix(hash,values)
    integer(int64),intent(inout)::hash
    real(real64),intent(in)::values(:,:)
    integer::i,j
    do j=1,size(values,2);do i=1,size(values,1)
      call hash_real(hash,values(i,j))
    enddo;enddo
  end subroutine hash_real_matrix

  subroutine hash_character(hash,value)
    integer(int64),intent(inout)::hash
    character(*),intent(in)::value
    integer::i
    call hash_integer(hash,len(value))
    do i=1,len(value);call hash_integer(hash,iachar(value(i:i)));enddo
  end subroutine hash_character

  subroutine hash_character_array(hash,values)
    integer(int64),intent(inout)::hash
    character(*),intent(in)::values(:)
    integer::i
    do i=1,size(values);call hash_character(hash,values(i));enddo
  end subroutine hash_character_array

  integer(int64) function hash_complex_matrix(values,tag)
    complex(real64),intent(in)::values(:,:)
    integer(int64),intent(in)::tag
    integer::i,j
    hash_complex_matrix=initial_hash(tag)
    do j=1,size(values,2);do i=1,size(values,1)
      call hash_complex(hash_complex_matrix,values(i,j))
    enddo;enddo
    if(hash_complex_matrix==0_int64)hash_complex_matrix=tag
  end function hash_complex_matrix

  integer(int64) function replicated_cache_integrity_hash(cache)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::cache
    integer::i,j
    replicated_cache_integrity_hash=initial_hash(709_int64)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%fragment_id)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%basis_generation)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%candidate_rank)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%retained_rank)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%setup_count)
    call hash_integer(replicated_cache_integrity_hash,cache%receipt%run_count)
    call hash_int64(replicated_cache_integrity_hash,cache%receipt%seed_fingerprint)
    call hash_int64(replicated_cache_integrity_hash,cache%receipt%basis_fingerprint)
    call hash_int64(replicated_cache_integrity_hash,cache%receipt%transform_fingerprint)
    call hash_real(replicated_cache_integrity_hash,cache%receipt%seed_reconstruction_defect)
    call hash_integer(replicated_cache_integrity_hash,cache%fragment_comm_size)
    do j=1,size(cache%candidate_compression,2);do i=1,size(cache%candidate_compression,1)
      call hash_complex(replicated_cache_integrity_hash,cache%candidate_compression(i,j))
    enddo;enddo
    do j=1,size(cache%wannier_transform,2);do i=1,size(cache%wannier_transform,1)
      call hash_complex(replicated_cache_integrity_hash,cache%wannier_transform(i,j))
    enddo;enddo
    do j=1,size(cache%dc_seed_coefficients_in_wannier,2);do i=1,size(cache%dc_seed_coefficients_in_wannier,1)
      call hash_complex(replicated_cache_integrity_hash,&
        cache%dc_seed_coefficients_in_wannier(i,j))
    enddo;enddo
    call hash_real_array(replicated_cache_integrity_hash,cache%physical_dc_seed_energies)
    call hash_real_array(replicated_cache_integrity_hash,cache%physical_dc_seed_occupations)
    if(replicated_cache_integrity_hash==0_int64)replicated_cache_integrity_hash=5_int64
  end function replicated_cache_integrity_hash
#endif
end module dg_hybrid_fragment_wannier
