#include "config.h"
module rt_dg_hybrid_initialization
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_checkpoint,only:s_rt_dg_hybrid_ground_state_payload,&
    read_rt_dg_hybrid_ground_state_checkpoint_coalesced,&
    authenticate_rt_dg_hybrid_ground_state_payload,rt_dg_hybrid_ground_state_checkpoint_version,&
    fingerprint_rt_dg_hybrid_component
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter::salmon_xctype_none=0,salmon_xctype_pz=1,salmon_xctype_pzm=2,&
    salmon_xctype_pbe=3,salmon_xctype_pw=7
  integer(int64),parameter::cell_wrapped_position_convention_fingerprint=&
    int(z'43454C4C57524150',int64)
  type,public::s_rt_dg_hybrid_state
    logical::valid=.false.,initial_invariants_valid=.false.,density_freshly_reconstructed=.false.
    real(real64)::startup_operator_covariance=huge(1d0),startup_projector_covariance=huge(1d0)
    integer::certified_rank=0,global_count=0,noccupied=0,operation_count=0,&
      nonidentity_operation_count=0
    integer(int64)::payload_fingerprint=0_int64,operator_structure_fingerprint=0_int64,&
      operator_value_fingerprint=0_int64,scope_fingerprint=0_int64
    type(s_dg_hybrid_sparse_metric)::metric
    type(s_dg_hybrid_sparse_operators)::operators
    integer(int64),allocatable::owned_row_ids(:),grid_ids(:)
    complex(real64),allocatable::coefficients(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      local_rows(:),sipg_rows(:,:),basis_values(:,:)
    real(real64),allocatable::density(:),grid_weights(:),occupations(:),eigenvalues(:),energy_receipt(:)
  end type s_rt_dg_hybrid_state
  type,public::s_rt_dg_hybrid_v3_startup_receipt
    logical::valid=.false.
    integer::certified_rank=0
    integer(int64)::payload_fingerprint=0_int64
    real(real64)::orbital_residual=huge(0d0),metric_defect=huge(0d0),&
      embedding_defect=huge(0d0),unitarity_defect=huge(0d0),basis_defect=huge(0d0),density_defect=huge(0d0),&
      electron_defect=huge(0d0),target_closure_defect=huge(0d0),&
      energy_covariance_defect=huge(0d0),projector_defect=huge(0d0),&
      operator_component_defect=huge(0d0),fixed_operator_covariance_defect=huge(0d0)
  end type s_rt_dg_hybrid_v3_startup_receipt
  type::s_rt_dg_hybrid_v3_named_fingerprints
    integer(int64)::top(11)=0_int64,catalog(6)=0_int64,certified(7)=0_int64,&
      rt(14)=0_int64,receipts(2)=0_int64,handoff(6)=0_int64
  end type s_rt_dg_hybrid_v3_named_fingerprints
  public::initialize_rt_dg_hybrid_from_checkpoint,validate_rt_dg_hybrid_v3_startup,&
    fingerprint_rt_dg_hybrid_scope,stamp_rt_dg_hybrid_v3_fingerprints
contains
  subroutine initialize_rt_dg_hybrid_from_checkpoint(comm,path,theory,periodic,nspin,spinorbit,plus_u,hse,&
      fix_func,jm,xctype,tolerances,state,ok,message)
    integer,intent(in)::comm,nspin,xctype(:)
    character(*),intent(in)::path,theory
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    real(real64),intent(in)::tolerances(4)
    type(s_rt_dg_hybrid_state),intent(out)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_rt_dg_hybrid_ground_state_payload)::payload
    type(s_rt_dg_hybrid_v3_startup_receipt)::receipt
    integer::i,ierr,local_bad,global_bad,file_nproc,active_xc
    integer(int64)::payload_fingerprint
    complex(real64),allocatable::startup_projected_position(:,:,:),startup_projected_basis(:,:)
    ok=.false.;message='';state=s_rt_dg_hybrid_state()
    local_bad=0
    if((trim(theory)/='tddft_response'.and.trim(theory)/='tddft_pulse').or..not.periodic.or.nspin/=1.or.&
      spinorbit.or.plus_u.or.hse.or.fix_func.or.jm.or.size(xctype)<1) local_bad=1
    if(any(.not.ieee_is_finite(tolerances)).or.any(tolerances<=0d0))local_bad=1
    active_xc=0
    do i=1,size(xctype)
      if(xctype(i)/=salmon_xctype_none)active_xc=active_xc+1
      if(.not.(xctype(i)==salmon_xctype_none.or.xctype(i)==salmon_xctype_pz.or.&
        xctype(i)==salmon_xctype_pzm.or.xctype(i)==salmon_xctype_pbe.or.xctype(i)==salmon_xctype_pw))local_bad=1
    enddo
    if(active_xc==0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='unsupported local hybrid RT scope';return;endif
    call probe_ground_state_rank_count(comm,path,file_nproc,ok,message);if(.not.ok)return
    call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,path,payload,payload_fingerprint,ok,message)
    if(.not.ok)return
    local_bad=merge(0,1,ok);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;return;endif
    call validate_rt_dg_hybrid_v3_startup(comm,payload,payload_fingerprint,tolerances,receipt,ok,message,&
      startup_projected_position,startup_projected_basis)
    if(.not.ok)return
    local_bad=merge(0,1,payload%scope_fingerprint/=0_int64.and.allocated(payload%scope_selectors).and.&
      allocated(payload%xc_types))
    if(local_bad==0)then
      if(size(payload%scope_selectors)/=8.or.any(payload%scope_selectors/=[1,1,1,0,0,0,0,0]).or.&
        size(payload%xc_types)/=size(xctype).or.any(payload%xc_types/=xctype))local_bad=1
      if(local_bad==0.and.payload%scope_fingerprint/=&
        fingerprint_rt_dg_hybrid_scope(payload%scope_selectors,payload%xc_types))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='checkpoint/local hybrid RT scope mismatch';return;endif
    call build_certified_rt_state(comm,payload,startup_projected_position,startup_projected_basis,state,ok,message)
    if(.not.ok)return
    state%payload_fingerprint=payload_fingerprint
    state%startup_operator_covariance=receipt%fixed_operator_covariance_defect
    state%startup_projector_covariance=receipt%projector_defect
    state%initial_invariants_valid=.true.;state%valid=.true.;ok=.true.;message=''
#else
    ok=.false.;message='hybrid RT initialization requires MPI'
#endif
  end subroutine initialize_rt_dg_hybrid_from_checkpoint
  integer(int64) function fingerprint_rt_dg_hybrid_scope(selectors,xctype) result(hash)
    integer,intent(in)::selectors(:),xctype(:)
    integer::i
    hash=int(z'A4093822299F31D0',int64)
    if(size(selectors)/=8)then;hash=0_int64;return;endif
    hash=ieor(ishftc(hash,7),int(selectors(1),int64))
    hash=ieor(ishftc(hash,7),int(selectors(3),int64))
    hash=ieor(ishftc(hash,7),int(selectors(2),int64))
    do i=4,8;hash=ieor(ishftc(hash,7),int(selectors(i),int64));enddo
    hash=ieor(ishftc(hash,7),int(size(xctype),int64))
    do i=1,size(xctype);hash=ieor(ishftc(hash,7),int(xctype(i),int64));enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_rt_dg_hybrid_scope

  subroutine stamp_rt_dg_hybrid_v3_fingerprints(comm,payload,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_rt_dg_hybrid_v3_named_fingerprints)::fingerprints
    integer::nproc,ierr
    call compute_rt_dg_hybrid_v3_named_fingerprints(comm,payload,fingerprints,ok,message)
    if(.not.ok)return
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='certified v3 fingerprint communicator query failed';return;endif
    payload%construction_catalog%ids_fingerprint=fingerprints%catalog(1)
    payload%construction_catalog%generation_fingerprint=fingerprints%catalog(2)
    payload%construction_catalog%ordering_fingerprint=fingerprints%catalog(3)
    payload%construction_catalog%ownership_fingerprint=fingerprints%catalog(4)
    payload%construction_catalog%provenance_fingerprint=fingerprints%catalog(5)
    payload%construction_catalog%catalog_fingerprint=fingerprints%catalog(6)
    payload%catalog_fingerprint=fingerprints%catalog(6)
    payload%state_fingerprint=fingerprints%top(1)
    payload%operator_structure_fingerprint=fingerprints%top(2+mod(nproc,2))
    payload%operator_value_fingerprint=fingerprints%top(4)
    payload%kinetic_fingerprint=fingerprints%top(5)
    payload%nonlocal_fingerprint=fingerprints%top(6)
    payload%local_fingerprint=fingerprints%top(7)
    payload%sipg_fingerprint=fingerprints%top(8)
    payload%scope_fingerprint=fingerprints%top(9)
    payload%analysis_fingerprint=fingerprints%top(10)
    payload%energy_fingerprint=fingerprints%top(11)
    payload%certified_basis%c_cert_fingerprint=fingerprints%certified(1)
    payload%certified_basis%u_rt_fingerprint=fingerprints%certified(2)
    payload%certified_basis%b_rt_fingerprint=fingerprints%certified(3)
    payload%certified_basis%initial_state_fingerprint=fingerprints%certified(4)
    payload%certified_basis%transformation_fingerprint=fingerprints%certified(5)
    payload%certified_basis%operator_fingerprint=fingerprints%certified(6)
    payload%certified_basis%fingerprint=fingerprints%certified(7)
    payload%rt_space%metric_fingerprint=fingerprints%rt(1)
    payload%rt_space%kinetic_fingerprint=fingerprints%rt(2)
    payload%rt_space%nonlocal_fingerprint=fingerprints%rt(3)
    payload%rt_space%local_fingerprint=fingerprints%rt(4)
    payload%rt_space%sipg_fingerprint=fingerprints%rt(5)
    payload%rt_space%hamiltonian_fingerprint=fingerprints%rt(6)
    payload%rt_space%basis_fingerprint=fingerprints%rt(7)
    payload%rt_space%density_fingerprint=fingerprints%rt(8)
    payload%rt_space%ownership_fingerprint=fingerprints%rt(9)
    payload%rt_space%scalar_fingerprint=fingerprints%rt(10)
    payload%rt_space%vector_fingerprint=fingerprints%rt(11)
    payload%rt_space%tensor_fingerprint=fingerprints%rt(12)
    payload%rt_space%representation_fingerprint=fingerprints%rt(13)
    payload%rt_space%fingerprint=fingerprints%rt(14)
    payload%energy_window%fingerprint=fingerprints%receipts(1)
    payload%symmetry_receipt%fingerprint=fingerprints%receipts(2)
    payload%handoff_receipts%position_fingerprint=fingerprints%handoff(1)
    payload%handoff_receipts%nonlocal_fingerprint=fingerprints%handoff(2)
    payload%handoff_receipts%face_fingerprint=fingerprints%handoff(3)
    payload%handoff_receipts%pseudopotential_fingerprint=fingerprints%handoff(4)
    payload%handoff_receipts%transformation_fingerprint=fingerprints%handoff(5)
    payload%handoff_receipts%fingerprint=fingerprints%handoff(6)
    ok=.true.;message=''
#else
    ok=.false.;message='certified Hybrid RT fingerprint stamping requires MPI'
#endif
  end subroutine stamp_rt_dg_hybrid_v3_fingerprints

  subroutine validate_rt_dg_hybrid_v3_startup(comm,payload,expected_fingerprint,tolerances,receipt,ok,message,&
      cached_projected_position,cached_projected_basis)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::expected_fingerprint
    real(real64),intent(in)::tolerances(4)
    type(s_rt_dg_hybrid_v3_startup_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable,intent(out),optional::cached_projected_position(:,:,:),cached_projected_basis(:,:)
    receipt=s_rt_dg_hybrid_v3_startup_receipt()
#ifdef USE_MPI
    call validate_rt_dg_hybrid_v3_startup_mpi(comm,payload,expected_fingerprint,tolerances,receipt,ok,message,&
      cached_projected_position,cached_projected_basis)
#else
    ok=.false.;message='certified Hybrid RT startup validation requires MPI'
#endif
  end subroutine validate_rt_dg_hybrid_v3_startup
#ifdef USE_MPI
  subroutine probe_ground_state_rank_count(comm,path,file_nproc,ok,message)
    integer,intent(in)::comm;character(*),intent(in)::path;integer,intent(out)::file_nproc
    logical,intent(out)::ok;character(*),intent(out)::message
    integer::rank,unit,status,ierr,version;character(16)::magic
    call MPI_Comm_rank(comm,rank,ierr);status=0;file_nproc=0
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=status)
      if(status==0)read(unit,iostat=status)magic,version,file_nproc
      if(status==0)close(unit)
      if(status==0.and.(magic/='SALMON_DG_GS001 '.or.&
        version/=rt_dg_hybrid_ground_state_checkpoint_version))status=1
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr);call MPI_Bcast(file_nproc,1,MPI_INTEGER,0,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.status==0.and.file_nproc>0
    if(ok)then;message='';else;message='cannot probe complete DG ground-state checkpoint';endif
  end subroutine probe_ground_state_rank_count

  subroutine compute_rt_dg_hybrid_v3_named_fingerprints(comm,payload,fingerprints,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    type(s_rt_dg_hybrid_v3_named_fingerprints),intent(out)::fingerprints
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::full_u(:,:),full_s(:,:),full_h(:,:),full_scalar(:,:,:),&
      full_vector(:,:,:,:),full_tensor(:,:,:,:,:)
    integer::n,r,nrtrow,nscalar,nvector,ntensor,i,row,ierr
    integer(int64)::grid_ownership_fingerprint,rotation_fingerprint
    logical::fingerprint_ok
    fingerprints=s_rt_dg_hybrid_v3_named_fingerprints();ok=.false.;message=''
    n=payload%global_count;r=payload%certified_basis%certified_count
    nrtrow=size(payload%rt_space%row_ids);nscalar=payload%rt_space%scalar_count
    nvector=payload%rt_space%vector_count;ntensor=payload%rt_space%tensor_count
    allocate(full_u(r,r),full_s(r,r),full_h(r,r),full_scalar(r,r,nscalar),&
      full_vector(r,r,3,nvector),full_tensor(r,r,3,3,ntensor))
    call collect_complex_rows(comm,payload%certified_basis%transformation_row_ids,&
      payload%certified_basis%u_rt,r,full_u,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    call collect_complex_rows(comm,payload%rt_space%row_ids,payload%rt_space%metric_rows,r,full_s,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    call collect_complex_rows(comm,payload%rt_space%row_ids,payload%rt_space%hamiltonian_rows,r,full_h,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    full_scalar=(0d0,0d0);full_vector=(0d0,0d0);full_tensor=(0d0,0d0)
    do i=1,nrtrow
      row=int(payload%rt_space%row_ids(i))
      full_scalar(row,:,:)=payload%rt_space%scalar_operator_rows(i,:,:)
      full_vector(row,:,:,:)=payload%rt_space%vector_operator_rows(i,:,:,:)
      full_tensor(row,:,:,:,:)=payload%rt_space%tensor_operator_rows(i,:,:,:,:)
    enddo
    if(size(full_scalar)>0)then
      call allreduce_complex_bits(comm,full_scalar,size(full_scalar),ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif
    if(size(full_vector)>0)then
      call allreduce_complex_bits(comm,full_vector,size(full_vector),ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif
    if(size(full_tensor)>0)then
      call allreduce_complex_bits(comm,full_tensor,size(full_tensor),ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif

    fingerprints%catalog(1)=fingerprint_checkpoint_integer64(payload%construction_catalog%ids)
    fingerprints%catalog(2)=fingerprint_checkpoint_integer(payload%construction_catalog%generations)
    fingerprints%catalog(3)=fingerprint_checkpoint_integer(payload%construction_catalog%ordering)
    fingerprints%catalog(4)=fingerprint_checkpoint_integer(payload%construction_catalog%ownership)
    fingerprints%catalog(5)=payload%basis_fingerprint
    fingerprints%catalog(6)=int(z'6A09E667F3BCC909',int64)
    do i=1,5
      fingerprints%catalog(6)=ieor(ishftc(fingerprints%catalog(6),7),fingerprints%catalog(i))
    enddo
    if(fingerprints%catalog(6)==0_int64)fingerprints%catalog(6)=1_int64
    call fingerprint_variational_operator(comm,payload,fingerprints%top(2),fingerprints%top(3),fingerprint_ok)
    if(.not.fingerprint_ok)goto 900
    call fingerprint_distributed_matrix(comm,payload%row_ids,payload%hamiltonian_rows,&
      fingerprints%top(4),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%kinetic_rows,&
      fingerprints%top(5),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%nonlocal_rows,&
      fingerprints%top(6),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%local_rows,&
      fingerprints%top(7),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%row_ids,payload%sipg_rows,&
      fingerprints%top(8),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    fingerprints%top(9)=fingerprint_rt_dg_hybrid_scope(payload%scope_selectors,payload%xc_types)
    fingerprints%top(10)=fingerprint_checkpoint_complex_rank3(payload%symmetry_representation)
    fingerprints%top(11)=fingerprint_checkpoint_real(payload%energy_receipt)
    call fingerprint_ground_state(comm,payload,fingerprints%top(4),fingerprints%top(1),fingerprint_ok)
    if(.not.fingerprint_ok)goto 900

    call fingerprint_certified_rows(comm,payload%certified_basis%construction_row_ids,&
      payload%certified_basis%c_cert,n,fingerprints%certified(1),fingerprint_ok)
    if(.not.fingerprint_ok)goto 900
    fingerprints%certified(2)=fingerprint_certified_localization(full_u,payload)
    call fingerprint_certified_rows(comm,payload%certified_basis%construction_row_ids,&
      payload%certified_basis%b_rt,n,fingerprints%certified(3),fingerprint_ok)
    if(.not.fingerprint_ok)goto 900
    fingerprints%certified(4)=fingerprint_checkpoint_complex_matrix(&
      payload%certified_basis%initial_occupied_amplitudes)
    fingerprints%certified(5)=fingerprints%certified(2)
    fingerprints%certified(6)=fingerprint_certified_operator(payload,full_s,full_h,&
      full_scalar,full_vector,full_tensor)
    fingerprints%certified(7)=fingerprint_certified_receipt(payload,fingerprints%certified)

    fingerprints%receipts(1)=fingerprint_energy_window_receipt(payload)
    fingerprints%receipts(2)=fingerprint_symmetry_receipt(payload,fingerprints%receipts(1),&
      fingerprints%certified(7))
    fingerprints%handoff(1)=payload%position_convention_fingerprint
    fingerprints%handoff(2)=fingerprints%top(6)
    fingerprints%handoff(3)=payload%face_fingerprint
    fingerprints%handoff(4)=payload%pseudopotential_fingerprint
    fingerprints%handoff(5)=fingerprints%certified(5)
    fingerprints%handoff(6)=ieor(fingerprints%handoff(1),ishftc(fingerprints%handoff(2),7))
    fingerprints%handoff(6)=ieor(fingerprints%handoff(6),ishftc(fingerprints%handoff(3),13))
    fingerprints%handoff(6)=ieor(fingerprints%handoff(6),ishftc(fingerprints%handoff(4),17))
    fingerprints%handoff(6)=ieor(fingerprints%handoff(6),ishftc(fingerprints%handoff(5),19))
    if(fingerprints%handoff(6)==0_int64)fingerprints%handoff(6)=1_int64

    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%metric_rows,fingerprints%rt(1),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%kinetic_rows,fingerprints%rt(2),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%nonlocal_rows,fingerprints%rt(3),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%local_rows,fingerprints%rt(4),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%sipg_rows,fingerprints%rt(5),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_rt_dg_hybrid_component(comm,payload%rt_space%row_ids,&
      payload%rt_space%hamiltonian_rows,fingerprints%rt(6),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_grid_complex(comm,payload%grid_ids,payload%rt_space%basis_values,&
      fingerprints%rt(7),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_grid_real(comm,payload%grid_ids,payload%rt_space%density,&
      fingerprints%rt(8),fingerprint_ok);if(.not.fingerprint_ok)goto 900
    call fingerprint_grid_integer(comm,payload%grid_ids,payload%rt_space%grid_owner_keys,&
      grid_ownership_fingerprint,fingerprint_ok);if(.not.fingerprint_ok)goto 900
    fingerprints%rt(9)=fingerprint_checkpoint_integer(payload%rt_space%row_owner_keys)
    fingerprints%rt(9)=ieor(ishftc(fingerprints%rt(9),7),grid_ownership_fingerprint)
    if(fingerprints%rt(9)==0_int64)fingerprints%rt(9)=1_int64
    fingerprints%rt(10)=fingerprint_checkpoint_complex_rank3(full_scalar)
    fingerprints%rt(11)=fingerprint_checkpoint_complex_rank4(full_vector)
    fingerprints%rt(12)=fingerprint_checkpoint_complex_rank5(full_tensor)
    fingerprints%rt(13)=fingerprint_checkpoint_complex_rank3(payload%rt_space%representation)
    rotation_fingerprint=fingerprint_checkpoint_real_rank3(payload%rt_space%cartesian_rotations)
    fingerprints%rt(14)=int(z'5BE0CD19137E2179',int64)
    do i=1,13;fingerprints%rt(14)=ieor(ishftc(fingerprints%rt(14),7),fingerprints%rt(i));enddo
    fingerprints%rt(14)=ieor(ishftc(fingerprints%rt(14),7),rotation_fingerprint)
    if(fingerprints%rt(14)==0_int64)fingerprints%rt(14)=1_int64
    ok=.true.;message='';return
900 message='certified Hybrid RT named fingerprint recomputation failed'
  end subroutine compute_rt_dg_hybrid_v3_named_fingerprints

  subroutine validate_rt_dg_hybrid_v3_startup_mpi(comm,payload,expected_fingerprint,tolerances,receipt,ok,message,&
      cached_projected_position,cached_projected_basis)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::expected_fingerprint
    real(real64),intent(in)::tolerances(4)
    type(s_rt_dg_hybrid_v3_startup_receipt),intent(inout)::receipt
    type(s_rt_dg_hybrid_v3_named_fingerprints)::observed_fingerprints
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable,intent(out),optional::cached_projected_position(:,:,:),cached_projected_basis(:,:)
    integer::n,r,nocc,nop,nscalar,nvector,ntensor,nrow,nrtrow,npoint
    integer::i,j,op,item,a,b,c,d,ierr,local_bad,global_bad
    integer(int64)::bits,minimum_bits,maximum_bits,minimum_fingerprint,maximum_fingerprint
    logical::authenticated
    character(256)::authentication_message
    real(real64)::local_value,global_value,scale,local_scale,global_scale,local_charge,global_charge
    real(real64),allocatable::reconstructed_density(:)
    complex(real64),allocatable::full_u(:,:),full_rt_s(:,:),full_rt_h(:,:),sc_local(:,:),hc_local(:,:),&
      projected_s(:,:),projected_h(:,:),projected_position(:,:,:),projected_components(:,:,:),projected_basis(:,:),&
      identity(:,:),energy(:,:),gram(:,:),work(:,:),projector(:,:),scalar_operators(:,:,:),&
      vector_operators(:,:,:,:),tensor_operators(:,:,:,:,:),transformed(:,:),expected(:,:)

    ok=.false.;message='';receipt=s_rt_dg_hybrid_v3_startup_receipt()
    call MPI_Allreduce(expected_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(expected_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='rank-disagreeing certified Hybrid RT payload fingerprint';return
    endif
    local_bad=0
    do i=1,4
      bits=transfer(tolerances(i),bits)
      call MPI_Allreduce(bits,minimum_bits,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Allreduce(bits,maximum_bits,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)local_bad=1
    enddo
    if(any(.not.ieee_is_finite(tolerances)).or.any(tolerances<=0d0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid or rank-disagreeing certified Hybrid RT startup tolerances';return
    endif
    call authenticate_rt_dg_hybrid_ground_state_payload(comm,payload,expected_fingerprint,&
      authenticated,authentication_message)
    if(.not.authenticated)then
      message='certified Hybrid RT payload authentication failed: '//trim(authentication_message);return
    endif
    call compute_rt_dg_hybrid_v3_named_fingerprints(comm,payload,observed_fingerprints,authenticated,&
      authentication_message)
    if(.not.authenticated)then
      message='certified Hybrid RT named fingerprint recomputation failed: '//trim(authentication_message);return
    endif
    local_bad=0
    if(any(observed_fingerprints%catalog/=[payload%construction_catalog%ids_fingerprint,&
      payload%construction_catalog%generation_fingerprint,payload%construction_catalog%ordering_fingerprint,&
      payload%construction_catalog%ownership_fingerprint,payload%construction_catalog%provenance_fingerprint,&
      payload%construction_catalog%catalog_fingerprint]))local_bad=1
    if(payload%catalog_fingerprint/=observed_fingerprints%catalog(6))local_bad=1
    if(payload%state_fingerprint/=observed_fingerprints%top(1).or.&
      (payload%operator_structure_fingerprint/=observed_fingerprints%top(2).and.&
      payload%operator_structure_fingerprint/=observed_fingerprints%top(3)).or.&
      any([payload%operator_value_fingerprint,payload%kinetic_fingerprint,payload%nonlocal_fingerprint,&
      payload%local_fingerprint,payload%sipg_fingerprint,payload%scope_fingerprint,&
      payload%analysis_fingerprint,payload%energy_fingerprint]/=observed_fingerprints%top(4:11)))local_bad=1
    if(.not.all(observed_fingerprints%certified==[&
      payload%certified_basis%c_cert_fingerprint,payload%certified_basis%u_rt_fingerprint,&
      payload%certified_basis%b_rt_fingerprint,payload%certified_basis%initial_state_fingerprint,&
      payload%certified_basis%transformation_fingerprint,payload%certified_basis%operator_fingerprint,&
      payload%certified_basis%fingerprint]))local_bad=1
    if(.not.all(observed_fingerprints%rt==[&
      payload%rt_space%metric_fingerprint,payload%rt_space%kinetic_fingerprint,&
      payload%rt_space%nonlocal_fingerprint,payload%rt_space%local_fingerprint,&
      payload%rt_space%sipg_fingerprint,payload%rt_space%hamiltonian_fingerprint,&
      payload%rt_space%basis_fingerprint,payload%rt_space%density_fingerprint,&
      payload%rt_space%ownership_fingerprint,payload%rt_space%scalar_fingerprint,&
      payload%rt_space%vector_fingerprint,payload%rt_space%tensor_fingerprint,&
      payload%rt_space%representation_fingerprint,payload%rt_space%fingerprint]))local_bad=1
    if(any(observed_fingerprints%receipts/=[payload%energy_window%fingerprint,&
      payload%symmetry_receipt%fingerprint]))local_bad=1
    if(any(observed_fingerprints%handoff/=[payload%handoff_receipts%position_fingerprint,&
      payload%handoff_receipts%nonlocal_fingerprint,payload%handoff_receipts%face_fingerprint,&
      payload%handoff_receipts%pseudopotential_fingerprint,&
      payload%handoff_receipts%transformation_fingerprint,payload%handoff_receipts%fingerprint]))local_bad=1
    if(payload%position_convention_fingerprint/=&
      cell_wrapped_position_convention_fingerprint)local_bad=1
    ! Complete-v3 does not serialize the producer preimages for these upstream
    ! provenance tokens (the electron token, for example, includes the full
    ! construction spectrum).  The aggregate payload authenticator binds the
    ! nonzero tokens; startup below independently recomputes every serialized
    ! physical consequence instead of pretending to regenerate those hashes.
    if(any([payload%metric_fingerprint,payload%basis_fingerprint,payload%face_fingerprint,&
      payload%dc_seed_fingerprint,payload%continuation_fingerprint,payload%selection_fingerprint,&
      payload%pseudopotential_fingerprint,payload%electron_count%fingerprint]==0_int64))local_bad=1
    if(size(payload%construction_catalog%ids)/=size(payload%effective_ids))then
      local_bad=1
    else if(any(payload%construction_catalog%ids/=int(payload%effective_ids,int64)))then
      local_bad=1
    endif
    if(payload%symmetry_receipt%scalar_covariance_defect/=&
      payload%certified_basis%scalar_covariance_defect.or.&
      payload%symmetry_receipt%vector_covariance_defect/=&
      payload%certified_basis%vector_covariance_defect.or.&
      payload%symmetry_receipt%tensor_covariance_defect/=&
      payload%certified_basis%tensor_covariance_defect.or.&
      payload%symmetry_receipt%final_basis_defect/=max(payload%certified_basis%rt_metric_defect,&
      payload%certified_basis%embedding_defect,payload%certified_basis%projector_invariance_defect).or.&
      payload%symmetry_receipt%maximum_physical_defect/=max(&
      payload%symmetry_receipt%worst_operation_defect,payload%symmetry_receipt%occupied_projector_defect,&
      payload%symmetry_receipt%density_defect,payload%certified_basis%scalar_covariance_defect,&
      payload%certified_basis%vector_covariance_defect,payload%certified_basis%tensor_covariance_defect,&
      payload%symmetry_receipt%final_basis_defect))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='certified Hybrid RT named fingerprint mismatch';return
    endif

    n=payload%global_count;r=payload%certified_basis%certified_count
    nocc=payload%certified_basis%occupied_count;nop=payload%rt_space%operation_count
    nscalar=payload%rt_space%scalar_count;nvector=payload%rt_space%vector_count
    ntensor=payload%rt_space%tensor_count;nrow=size(payload%row_ids)
    nrtrow=size(payload%rt_space%row_ids);npoint=size(payload%grid_ids)
    allocate(full_u(r,r),full_rt_s(r,r),full_rt_h(r,r),projected_s(r,r),projected_h(r,r),&
      projected_position(3,r,r),projected_components(r,r,4),identity(r,r),energy(r,r),gram(r,r),&
      work(r,r),projector(r,r),scalar_operators(r,r,nscalar),&
      vector_operators(r,r,3,nvector),tensor_operators(r,r,3,3,ntensor),transformed(r,r),expected(r,r))
    call collect_complex_rows(comm,payload%certified_basis%transformation_row_ids,&
      payload%certified_basis%u_rt,r,full_u,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call collect_complex_rows(comm,payload%rt_space%row_ids,payload%rt_space%metric_rows,r,full_rt_s,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    call collect_complex_rows(comm,payload%rt_space%row_ids,payload%rt_space%hamiltonian_rows,r,full_rt_h,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    call compute_construction_projection(comm,payload,full_u,sc_local,hc_local,projected_s,projected_h,&
      projected_position,projected_components,projected_basis,receipt%embedding_defect,ierr)
    if(ierr/=MPI_SUCCESS)goto 900

    local_value=0d0;local_scale=1d0
    if(npoint>0)then
      local_value=maximum_complex_matrix(projected_basis-payload%rt_space%basis_values)
      local_scale=max(local_scale,maximum_complex_matrix(projected_basis),&
        maximum_complex_matrix(payload%rt_space%basis_values))
    endif
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    receipt%basis_defect=global_value/global_scale

    identity=(0d0,0d0);energy=(0d0,0d0)
    do i=1,r
      identity(i,i)=(1d0,0d0);energy(i,i)=cmplx(payload%certified_basis%certified_eigenvalues(i),0d0,real64)
    enddo
    local_scale=1d0;if(nrow>0)local_scale=max(local_scale,maximum_complex_matrix(hc_local))
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    do j=1,r;hc_local(:,j)=hc_local(:,j)-payload%certified_basis%certified_eigenvalues(j)*sc_local(:,j);enddo
    local_value=0d0;if(nrow>0)local_value=maximum_complex_matrix(hc_local)/global_scale
    call MPI_Allreduce(local_value,receipt%orbital_residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    scale=max(1d0,maximum_complex_matrix(full_rt_h),maximum_complex_matrix(projected_h))
    receipt%orbital_residual=max(receipt%orbital_residual,maximum_complex_matrix(projected_h-full_rt_h)/scale,&
      maximum_complex_matrix(full_rt_h-matmul(conjg(transpose(full_u)),matmul(energy,full_u)))/scale)
    gram=matmul(conjg(transpose(payload%certified_basis%c_cert)),sc_local)
    call MPI_Allreduce(MPI_IN_PLACE,gram,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    receipt%metric_defect=max(maximum_complex_matrix(gram-identity),&
      maximum_complex_matrix(projected_s-full_rt_s),maximum_complex_matrix(full_rt_s-identity))
    receipt%embedding_defect=max(receipt%embedding_defect,maximum_complex_matrix(&
      payload%certified_basis%initial_occupied_amplitudes-conjg(transpose(full_u(1:nocc,:)))))
    receipt%unitarity_defect=maximum_complex_matrix(matmul(conjg(transpose(full_u)),full_u)-identity)

    receipt%target_closure_defect=0d0;receipt%energy_covariance_defect=0d0
    local_scale=1d0;if(nrow>0)local_scale=max(local_scale,maximum_complex_matrix(payload%certified_basis%c_cert))
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call compute_construction_target_closure(comm,payload,full_u,global_scale,&
      receipt%target_closure_defect,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    do op=1,nop
      work=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
        payload%rt_space%representation(:,:,op))-identity
      receipt%target_closure_defect=max(receipt%target_closure_defect,maximum_complex_matrix(work))
      work=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
        matmul(full_rt_s,payload%rt_space%representation(:,:,op)))-full_rt_s
      receipt%target_closure_defect=max(receipt%target_closure_defect,maximum_complex_matrix(work))
      work=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
        matmul(full_rt_h,payload%rt_space%representation(:,:,op)))-full_rt_h
      receipt%energy_covariance_defect=max(receipt%energy_covariance_defect,&
        maximum_complex_matrix(work)/max(1d0,maximum_complex_matrix(full_rt_h)))
    enddo

    projector=(0d0,0d0)
    do j=1,nocc
      do i=1,r
        projector(i,:)=projector(i,:)+payload%certified_basis%occupations(j)*&
          payload%certified_basis%initial_occupied_amplitudes(i,j)*&
          conjg(payload%certified_basis%initial_occupied_amplitudes(:,j))
      enddo
    enddo
    receipt%projector_defect=0d0
    do op=1,nop
      work=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
        matmul(projector,payload%rt_space%representation(:,:,op)))-projector
      receipt%projector_defect=max(receipt%projector_defect,&
        maximum_complex_matrix(work)/max(1d0,maximum_complex_matrix(projector)))
    enddo

    scalar_operators=(0d0,0d0);vector_operators=(0d0,0d0);tensor_operators=(0d0,0d0)
    do i=1,nrtrow
      scalar_operators(int(payload%rt_space%row_ids(i)),:,:)=payload%rt_space%scalar_operator_rows(i,:,:)
      vector_operators(int(payload%rt_space%row_ids(i)),:,:,:)=payload%rt_space%vector_operator_rows(i,:,:,:)
      tensor_operators(int(payload%rt_space%row_ids(i)),:,:,:,:)=payload%rt_space%tensor_operator_rows(i,:,:,:,:)
    enddo
    local_bad=merge(0,1,nscalar>=4)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    if(global_bad/=0)then;message='certified Hybrid RT scalar component receipt is incomplete';return;endif
    local_value=0d0;local_scale=1d0
    if(nrtrow>0)then
      local_value=max(maximum_complex_matrix(payload%rt_space%kinetic_rows-&
        payload%rt_space%scalar_operator_rows(:,:,1)),&
        maximum_complex_matrix(payload%rt_space%nonlocal_rows-&
        payload%rt_space%scalar_operator_rows(:,:,2)),&
        maximum_complex_matrix(payload%rt_space%local_rows-&
        payload%rt_space%scalar_operator_rows(:,:,3)),&
        maximum_complex_matrix(payload%rt_space%sipg_rows-&
        payload%rt_space%scalar_operator_rows(:,:,4)),&
        maximum_complex_matrix(payload%rt_space%hamiltonian_rows-&
        (payload%rt_space%kinetic_rows+payload%rt_space%nonlocal_rows+&
        payload%rt_space%local_rows+payload%rt_space%sipg_rows)))
      if(nscalar>=5)local_value=max(local_value,maximum_complex_matrix(&
        payload%rt_space%hamiltonian_rows-payload%rt_space%scalar_operator_rows(:,:,5)))
      local_scale=max(local_scale,maximum_complex_matrix(payload%rt_space%hamiltonian_rows),&
        maximum_complex_rank3(payload%rt_space%scalar_operator_rows))
    endif
    call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    receipt%operator_component_defect=global_value/global_scale
    if(size(scalar_operators)>0)then
      call MPI_Allreduce(MPI_IN_PLACE,scalar_operators,size(scalar_operators),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif
    scale=max(1d0,maximum_complex_rank3(projected_components),maximum_complex_rank3(scalar_operators(:,:,1:4)))
    receipt%operator_component_defect=max(receipt%operator_component_defect,&
      maximum_complex_rank3(projected_components-scalar_operators(:,:,1:4))/scale)
    if(size(vector_operators)>0)then
      call MPI_Allreduce(MPI_IN_PLACE,vector_operators,size(vector_operators),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif
    if(size(tensor_operators)>0)then
      call MPI_Allreduce(MPI_IN_PLACE,tensor_operators,size(tensor_operators),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)goto 900
    endif
    ! Cell-wrapped position transforms affinely under translations and
    ! nonsymmorphic operations.  Its translation term is not serialized, so
    ! validate projection provenance and Hermiticity here; the homogeneous
    ! vector certificate below is the stored canonical momentum.
    receipt%fixed_operator_covariance_defect=0d0
    scale=max(1d0,maximum_complex_rank3(projected_position))
    do a=1,3
      receipt%fixed_operator_covariance_defect=max(receipt%fixed_operator_covariance_defect,&
        maximum_complex_matrix(projected_position(a,:,:)-conjg(transpose(projected_position(a,:,:))))/scale)
    enddo
    do op=1,nop
      call update_scalar_covariance(full_rt_s,payload%rt_space%representation(:,:,op),&
        receipt%fixed_operator_covariance_defect)
      call update_scalar_covariance(full_rt_h,payload%rt_space%representation(:,:,op),&
        receipt%fixed_operator_covariance_defect)
      do item=1,nscalar
        call update_scalar_covariance(scalar_operators(:,:,item),payload%rt_space%representation(:,:,op),&
          receipt%fixed_operator_covariance_defect)
      enddo
      do item=1,nvector
        scale=max(1d0,maximum_complex_rank3(vector_operators(:,:,:,item)))
        do a=1,3
          transformed=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
            matmul(vector_operators(:,:,a,item),payload%rt_space%representation(:,:,op)))
          expected=(0d0,0d0)
          do b=1,3;expected=expected+payload%rt_space%cartesian_rotations(a,b,op)*vector_operators(:,:,b,item);enddo
          receipt%fixed_operator_covariance_defect=max(receipt%fixed_operator_covariance_defect,&
            maximum_complex_matrix(transformed-expected)/scale)
        enddo
      enddo
      do item=1,ntensor
        scale=max(1d0,maximum_complex_rank4(tensor_operators(:,:,:,:,item)))
        do b=1,3;do a=1,3
          transformed=matmul(conjg(transpose(payload%rt_space%representation(:,:,op))),&
            matmul(tensor_operators(:,:,a,b,item),payload%rt_space%representation(:,:,op)))
          expected=(0d0,0d0)
          do d=1,3;do c=1,3
            expected=expected+payload%rt_space%cartesian_rotations(a,c,op)*&
              payload%rt_space%cartesian_rotations(b,d,op)*tensor_operators(:,:,c,d,item)
          enddo;enddo
          receipt%fixed_operator_covariance_defect=max(receipt%fixed_operator_covariance_defect,&
            maximum_complex_matrix(transformed-expected)/scale)
        enddo;enddo
      enddo
    enddo

    call reconstruct_certified_density(payload,reconstructed_density)
    local_value=0d0
    if(npoint>0)local_value=max(maxval(abs(reconstructed_density-payload%rt_space%density)),&
      maxval(abs(reconstructed_density-payload%density)))
    call MPI_Allreduce(local_value,receipt%density_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    local_charge=sum(reconstructed_density*payload%grid_weights)
    call MPI_Allreduce(local_charge,global_charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    receipt%electron_defect=max(abs(payload%electron_count%expected_count-payload%electron_count%actual_count),&
      abs(sum(payload%certified_basis%occupations)-payload%electron_count%actual_count),&
      abs(global_charge-payload%electron_count%actual_count),payload%electron_count%defect,&
      payload%electron_count%omitted_tail)

    local_bad=merge(0,1,all(ieee_is_finite([receipt%orbital_residual,receipt%metric_defect,&
      receipt%embedding_defect,receipt%unitarity_defect,receipt%basis_defect,receipt%density_defect,receipt%electron_defect,&
      receipt%target_closure_defect,receipt%energy_covariance_defect,receipt%projector_defect,&
      receipt%operator_component_defect,receipt%fixed_operator_covariance_defect])))
    if(receipt%orbital_residual>tolerances(1).or.receipt%metric_defect>tolerances(1).or.&
      receipt%embedding_defect>tolerances(1).or.receipt%unitarity_defect>tolerances(1).or.&
      receipt%basis_defect>tolerances(1).or.&
      receipt%density_defect>tolerances(2).or.receipt%electron_defect>tolerances(3).or.&
      receipt%target_closure_defect>tolerances(4).or.receipt%energy_covariance_defect>tolerances(4).or.&
      receipt%projector_defect>tolerances(4).or.&
      receipt%operator_component_defect>tolerances(4).or.&
      receipt%fixed_operator_covariance_defect>tolerances(4))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)goto 900
    if(global_bad/=0)then;message='certified Hybrid RT startup invariant tolerance exceeded';return;endif
    receipt%certified_rank=r;receipt%payload_fingerprint=expected_fingerprint
    if(present(cached_projected_position))allocate(cached_projected_position,source=projected_position)
    if(present(cached_projected_basis))allocate(cached_projected_basis,source=projected_basis)
    receipt%valid=.true.;ok=.true.;message='';return
900 message='certified Hybrid RT startup collective validation failed'
  contains
    subroutine update_scalar_covariance(operator,representation,defect)
      complex(real64),intent(in)::operator(:,:),representation(:,:)
      real(real64),intent(inout)::defect
      transformed=matmul(conjg(transpose(representation)),matmul(operator,representation))
      defect=max(defect,maximum_complex_matrix(transformed-operator)/max(1d0,maximum_complex_matrix(operator)))
    end subroutine update_scalar_covariance
  end subroutine validate_rt_dg_hybrid_v3_startup_mpi

  subroutine build_certified_rt_state(comm,payload,projected_position,projected_basis,state,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    complex(real64),intent(in)::projected_position(:,:,:),projected_basis(:,:)
    type(s_rt_dg_hybrid_state),intent(out)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::r,nocc,nowned,i,j,metric_edge,operator_edge,metric_nnz,operator_nnz,ierr
    integer(int64)::structure_fingerprint
    real(real64)::local_maximum,global_maximum,entry_scale,entry_tolerance
    real(real64),allocatable::reconstructed_density(:)
    r=payload%rt_space%rank;nocc=payload%certified_basis%occupied_count
    nowned=size(payload%rt_space%row_ids);state=s_rt_dg_hybrid_state();ok=.false.;message=''

    allocate(state%owned_row_ids(nowned),state%coefficients(nowned,nocc),state%kinetic_rows(nowned,r),&
      state%nonlocal_rows(nowned,r),state%sipg_rows(nowned,r))
    state%owned_row_ids=payload%rt_space%row_ids
    do i=1,nowned
      state%coefficients(i,:)=payload%certified_basis%initial_occupied_amplitudes(&
        int(state%owned_row_ids(i)),:)
    enddo
    state%kinetic_rows=payload%rt_space%kinetic_rows
    state%nonlocal_rows=payload%rt_space%nonlocal_rows
    state%sipg_rows=payload%rt_space%sipg_rows

    entry_scale=max(1d0,maximum_complex_matrix(payload%rt_space%metric_rows),&
      maximum_complex_matrix(payload%rt_space%hamiltonian_rows),&
      maximum_complex_matrix(payload%rt_space%kinetic_rows),&
      maximum_complex_matrix(payload%rt_space%nonlocal_rows),&
      maximum_complex_matrix(payload%rt_space%local_rows),&
      maximum_complex_matrix(payload%rt_space%sipg_rows),maximum_complex_rank3(projected_position))
    entry_tolerance=100d0*epsilon(1d0)*entry_scale
    metric_nnz=0;operator_nnz=0
    do i=1,nowned
      do j=1,r
        if(payload%rt_space%metric_rows(i,j)/=(0d0,0d0))metric_nnz=metric_nnz+1
        if(abs(payload%rt_space%hamiltonian_rows(i,j))>entry_tolerance.or.&
          payload%rt_space%metric_rows(i,j)/=(0d0,0d0).or.&
          abs(payload%rt_space%kinetic_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%nonlocal_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%local_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%sipg_rows(i,j))>entry_tolerance.or.&
          any(abs(projected_position(:,int(state%owned_row_ids(i)),j))>entry_tolerance))operator_nnz=operator_nnz+1
      enddo
    enddo
    allocate(state%metric%owned_row_ids(nowned),state%metric%row_offsets(nowned+1),&
      state%metric%column_ids(metric_nnz),state%metric%values(metric_nnz),state%metric%active_rows(r),&
      state%metric%packet_ids(r),state%operators%owned_row_ids(nowned),state%operators%row_offsets(nowned+1),&
      state%operators%column_ids(operator_nnz),state%operators%metric_values(operator_nnz),&
      state%operators%hamiltonian_values(operator_nnz),state%operators%position_values(3,operator_nnz),&
      state%local_rows(operator_nnz))
    state%metric%row_offsets(1)=1;state%operators%row_offsets(1)=1
    metric_edge=0;operator_edge=0
    do i=1,nowned
      do j=1,r
        if(payload%rt_space%metric_rows(i,j)/=(0d0,0d0))then
          metric_edge=metric_edge+1;state%metric%column_ids(metric_edge)=j
          state%metric%values(metric_edge)=payload%rt_space%metric_rows(i,j)
        endif
        if(abs(payload%rt_space%hamiltonian_rows(i,j))>entry_tolerance.or.&
          payload%rt_space%metric_rows(i,j)/=(0d0,0d0).or.&
          abs(payload%rt_space%kinetic_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%nonlocal_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%local_rows(i,j))>entry_tolerance.or.&
          abs(payload%rt_space%sipg_rows(i,j))>entry_tolerance.or.&
          any(abs(projected_position(:,int(state%owned_row_ids(i)),j))>entry_tolerance))then
          operator_edge=operator_edge+1;state%operators%column_ids(operator_edge)=j
          state%operators%metric_values(operator_edge)=payload%rt_space%metric_rows(i,j)
          state%operators%hamiltonian_values(operator_edge)=payload%rt_space%hamiltonian_rows(i,j)
          state%operators%position_values(:,operator_edge)=projected_position(:,int(state%owned_row_ids(i)),j)
          state%local_rows(operator_edge)=payload%rt_space%local_rows(i,j)
        endif
      enddo
      state%metric%row_offsets(i+1)=metric_edge+1;state%operators%row_offsets(i+1)=operator_edge+1
    enddo
    state%metric%owned_row_ids=state%owned_row_ids;state%metric%global_count=r
    state%metric%numerical_rank=r;state%metric%max_row_nnz=0;state%metric%active_rows=.true.
    do i=1,nowned
      state%metric%max_row_nnz=max(state%metric%max_row_nnz,&
        state%metric%row_offsets(i+1)-state%metric%row_offsets(i))
    enddo
    state%metric%packet_ids=1;state%metric%valid=.true.
    state%metric%fingerprint=payload%rt_space%metric_fingerprint;state%metric%condition_estimate=1d0
    local_maximum=0d0
    if(nowned>0)local_maximum=maximum_complex_matrix(payload%rt_space%metric_rows)
    call MPI_Allreduce(local_maximum,global_maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='certified Hybrid RT metric reduction failed';return;endif
    state%metric%maximum_value=global_maximum

    structure_fingerprint=ieor(ishftc(payload%rt_space%ownership_fingerprint,7),&
      payload%position_convention_fingerprint)
    structure_fingerprint=ieor(ishftc(structure_fingerprint,11),int(r,int64))
    if(structure_fingerprint==0_int64)structure_fingerprint=1_int64
    state%operators%owned_row_ids=state%owned_row_ids;state%operators%global_count=r
    state%operators%metric_fingerprint=payload%rt_space%metric_fingerprint
    state%operators%selection_fingerprint=payload%selection_fingerprint
    state%operators%window_fingerprint=payload%energy_window%fingerprint
    state%operators%complement_fingerprint=payload%certified_basis%fingerprint
    state%operators%position_convention_fingerprint=payload%position_convention_fingerprint
    state%operators%fingerprint=structure_fingerprint;state%operators%valid=.true.

    call reconstruct_density_from_basis(projected_basis,payload%certified_basis%initial_occupied_amplitudes,&
      payload%certified_basis%occupations,reconstructed_density)
    allocate(state%grid_ids,source=payload%grid_ids)
    allocate(state%grid_weights,source=payload%grid_weights)
    allocate(state%density,source=reconstructed_density)
    allocate(state%basis_values,source=projected_basis)
    allocate(state%occupations,source=payload%certified_basis%occupations)
    allocate(state%eigenvalues,source=payload%certified_basis%certified_eigenvalues)
    allocate(state%energy_receipt,source=payload%energy_receipt)
    state%certified_rank=r;state%global_count=r;state%noccupied=nocc
    state%operation_count=payload%operation_count
    state%nonidentity_operation_count=payload%nonidentity_operation_count
    state%operator_structure_fingerprint=structure_fingerprint
    state%operator_value_fingerprint=payload%rt_space%hamiltonian_fingerprint
    state%scope_fingerprint=payload%scope_fingerprint
    ok=.true.;message=''
  end subroutine build_certified_rt_state

  subroutine compute_construction_projection(comm,payload,full_u,sc_local,hc_local,projected_s,projected_h,&
      projected_position,projected_components,projected_basis,embedding_defect,ierr)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    complex(real64),intent(in)::full_u(:,:)
    complex(real64),allocatable,intent(out)::sc_local(:,:),hc_local(:,:)
    complex(real64),intent(out)::projected_s(:,:),projected_h(:,:),projected_position(:,:,:),&
      projected_components(:,:,:)
    complex(real64),allocatable,intent(out)::projected_basis(:,:)
    real(real64),intent(out)::embedding_defect
    integer,intent(out)::ierr
    integer::n,r,nrow,npoint,rank,column,row_position,i,component,point
    integer,allocatable::owners(:)
    complex(real64),allocatable::c_buffer(:),b_buffer(:),sb_local(:,:),hb_local(:,:),xb_local(:,:,:),&
      component_b_local(:,:,:),local_matrix(:,:)
    real(real64)::local_defect
    n=payload%global_count;r=payload%certified_basis%certified_count;nrow=size(payload%row_ids)
    npoint=size(payload%grid_ids)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(owners(n),c_buffer(r),b_buffer(r),sc_local(nrow,r),hc_local(nrow,r),&
      sb_local(nrow,r),hb_local(nrow,r),xb_local(3,nrow,r),component_b_local(4,nrow,r),&
      local_matrix(r,r),projected_basis(r,npoint))
    call construction_owner_directory(comm,payload%row_ids,n,owners,ierr);if(ierr/=MPI_SUCCESS)return
    sc_local=(0d0,0d0);hc_local=(0d0,0d0);sb_local=(0d0,0d0);hb_local=(0d0,0d0)
    xb_local=(0d0,0d0);component_b_local=(0d0,0d0);projected_basis=(0d0,0d0)
    do column=1,n
      c_buffer=(0d0,0d0);b_buffer=(0d0,0d0)
      if(rank==owners(column))then
        row_position=find_row_position(payload%row_ids,column)
        c_buffer=payload%certified_basis%c_cert(row_position,:)
        b_buffer=payload%certified_basis%b_rt(row_position,:)
      endif
      call MPI_Bcast(c_buffer,r,MPI_DOUBLE_COMPLEX,owners(column),comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Bcast(b_buffer,r,MPI_DOUBLE_COMPLEX,owners(column),comm,ierr);if(ierr/=MPI_SUCCESS)return
      do point=1,npoint
        projected_basis(:,point)=projected_basis(:,point)+b_buffer*payload%basis_values(column,point)
      enddo
      do i=1,nrow
        sc_local(i,:)=sc_local(i,:)+payload%metric_rows(i,column)*c_buffer
        hc_local(i,:)=hc_local(i,:)+payload%hamiltonian_rows(i,column)*c_buffer
        sb_local(i,:)=sb_local(i,:)+payload%metric_rows(i,column)*b_buffer
        hb_local(i,:)=hb_local(i,:)+payload%hamiltonian_rows(i,column)*b_buffer
        component_b_local(1,i,:)=component_b_local(1,i,:)+payload%kinetic_rows(i,column)*b_buffer
        component_b_local(2,i,:)=component_b_local(2,i,:)+payload%nonlocal_rows(i,column)*b_buffer
        component_b_local(3,i,:)=component_b_local(3,i,:)+payload%local_rows(i,column)*b_buffer
        component_b_local(4,i,:)=component_b_local(4,i,:)+payload%sipg_rows(i,column)*b_buffer
        do component=1,3
          xb_local(component,i,:)=xb_local(component,i,:)+payload%position_rows(component,i,column)*b_buffer
        enddo
      enddo
    enddo
    local_matrix=matmul(conjg(transpose(payload%certified_basis%b_rt)),sb_local)
    call MPI_Allreduce(local_matrix,projected_s,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_matrix=matmul(conjg(transpose(payload%certified_basis%b_rt)),hb_local)
    call MPI_Allreduce(local_matrix,projected_h,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do component=1,3
      local_matrix=matmul(conjg(transpose(payload%certified_basis%b_rt)),xb_local(component,:,:))
      call MPI_Allreduce(MPI_IN_PLACE,local_matrix,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      projected_position(component,:,:)=local_matrix
    enddo
    do component=1,4
      local_matrix=matmul(conjg(transpose(payload%certified_basis%b_rt)),component_b_local(component,:,:))
      call MPI_Allreduce(MPI_IN_PLACE,local_matrix,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      projected_components(:,:,component)=local_matrix
    enddo
    local_defect=0d0
    if(nrow>0)local_defect=maximum_complex_matrix(payload%certified_basis%b_rt-&
      matmul(payload%certified_basis%c_cert,full_u))
    call MPI_Allreduce(local_defect,embedding_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  end subroutine compute_construction_projection

  subroutine compute_construction_target_closure(comm,payload,full_u,scale,defect,ierr)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    complex(real64),intent(in)::full_u(:,:)
    real(real64),intent(in)::scale
    real(real64),intent(out)::defect
    integer,intent(out)::ierr
    integer(int64),parameter::target_action_elements=4194304_int64
    integer::n,r,nrow,nop,nproc,rank,owner,count,max_count,first_column,last_column
    integer::first_operation,operation_count,tile_operation,operation,tile_size,i,local_bad,global_bad
    integer(int64)::elements_per_operation
    integer,allocatable::owner_counts(:),owner_offsets(:)
    integer(int64),allocatable::owner_ids(:)
    real(real64),allocatable::local_defects(:),global_defects(:)
    complex(real64),allocatable::c_block(:,:),action(:,:,:),representation_c(:,:)
    n=payload%global_count;r=payload%certified_basis%certified_count
    nrow=size(payload%row_ids);nop=payload%rt_space%operation_count
    defect=0d0;ierr=MPI_SUCCESS
    if(nop<=0)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,size(payload%certified_basis%c_cert,1)==nrow.and.&
      size(payload%certified_basis%c_cert,2)==r)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)then;ierr=-1;return;endif
    allocate(owner_counts(nproc),owner_offsets(nproc))
    call MPI_Allgather(nrow,1,MPI_INTEGER,owner_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    owner_offsets(1)=0
    do owner=2,nproc;owner_offsets(owner)=owner_offsets(owner-1)+owner_counts(owner-1);enddo
    if(sum(owner_counts)/=n)then;ierr=-1;return;endif
    max_count=maxval(owner_counts)
    if(max_count<=0)then;ierr=-1;return;endif
    allocate(owner_ids(n))
    call MPI_Allgatherv(payload%row_ids,nrow,MPI_INTEGER8,owner_ids,owner_counts,owner_offsets,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    elements_per_operation=max(1_int64,int(max_count,int64)*int(r,int64))
    tile_size=int(max(1_int64,min(int(nop,int64),target_action_elements/elements_per_operation)))
    allocate(c_block(r,max_count),action(nrow,r,tile_size),representation_c(r,r),&
      local_defects(tile_size),global_defects(tile_size))
    do first_operation=1,nop,tile_size
      operation_count=min(tile_size,nop-first_operation+1)
      action=(0d0,0d0)
      do owner=1,nproc
        count=owner_counts(owner);if(count==0)cycle
        first_column=owner_offsets(owner)+1;last_column=first_column+count-1
        if(rank==owner-1)c_block(:,1:count)=transpose(payload%certified_basis%c_cert)
        call MPI_Bcast(c_block,count*r,MPI_DOUBLE_COMPLEX,owner-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)return
        do tile_operation=1,operation_count
          operation=first_operation+tile_operation-1
          do i=1,nrow
            action(i,:,tile_operation)=action(i,:,tile_operation)+matmul(c_block(:,1:count),&
              payload%symmetry_representation(int(payload%row_ids(i)),&
              int(owner_ids(first_column:last_column)),operation))
          enddo
        enddo
      enddo
      local_defects=0d0
      do tile_operation=1,operation_count
        operation=first_operation+tile_operation-1
        representation_c=matmul(full_u,matmul(payload%rt_space%representation(:,:,operation),&
          conjg(transpose(full_u))))
        if(nrow>0)local_defects(tile_operation)=maximum_complex_matrix(action(:,:,tile_operation)-&
          matmul(payload%certified_basis%c_cert,representation_c))/scale
      enddo
      call MPI_Allreduce(local_defects,global_defects,operation_count,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      defect=max(defect,maxval(global_defects(1:operation_count)))
    enddo
  end subroutine compute_construction_target_closure

  subroutine construction_owner_directory(comm,row_ids,global_count,owners,ierr)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    integer,intent(out)::owners(global_count),ierr
    integer::rank,i,counts(global_count)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    owners=-1;counts=0
    do i=1,size(row_ids)
      counts(int(row_ids(i)))=counts(int(row_ids(i)))+1;owners(int(row_ids(i)))=rank
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(counts/=1).or.any(owners<0))ierr=1
  end subroutine construction_owner_directory

  subroutine collect_complex_rows(comm,row_ids,local_rows,global_count,full_rows,ierr)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::local_rows(:,:)
    complex(real64),intent(out)::full_rows(:,:)
    integer,intent(out)::ierr
    integer::i,local_bad,global_bad
    integer,allocatable::local_counts(:),global_counts(:)
    local_bad=merge(0,1,size(local_rows,1)==size(row_ids).and.size(full_rows,1)==global_count.and.&
      size(full_rows,2)==size(local_rows,2))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)then;ierr=1;return;endif
    local_bad=merge(0,1,all(row_ids>=1_int64).and.all(row_ids<=int(global_count,int64)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)then;ierr=1;return;endif
    allocate(local_counts(global_count),global_counts(global_count));local_counts=0
    do i=1,size(row_ids);local_counts(int(row_ids(i)))=local_counts(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(local_counts,global_counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(any(global_counts/=1))then;ierr=1;return;endif
    full_rows=(0d0,0d0)
    do i=1,size(row_ids);full_rows(int(row_ids(i)),:)=local_rows(i,:);enddo
    call allreduce_complex_bits(comm,full_rows,size(full_rows),ierr)
  end subroutine collect_complex_rows

  subroutine allreduce_complex_bits(comm,values,value_count,ierr)
    integer,intent(in)::comm,value_count
    complex(real64),intent(inout)::values(*)
    integer,intent(out)::ierr
    integer(int64),allocatable::bits(:)
    if(value_count==0)then;ierr=MPI_SUCCESS;return;endif
    allocate(bits(2*value_count));bits=transfer(values(1:value_count),bits)
    call MPI_Allreduce(MPI_IN_PLACE,bits,size(bits),MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr==MPI_SUCCESS)values(1:value_count)=transfer(bits,values(1:value_count))
  end subroutine allreduce_complex_bits

  integer function find_row_position(row_ids,target) result(position)
    integer(int64),intent(in)::row_ids(:)
    integer,intent(in)::target
    integer::i
    position=0
    do i=1,size(row_ids);if(row_ids(i)==int(target,int64))then;position=i;return;endif;enddo
  end function find_row_position

  subroutine fingerprint_distributed_matrix(comm,row_ids,matrix,fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::matrix(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr,local_bad,global_bad
    integer(int64)::local_hash,real_bits,imaginary_bits,entry_hash
    local_bad=merge(0,1,size(matrix,1)==size(row_ids).and.&
      all(ieee_is_finite(real(matrix))).and.all(ieee_is_finite(aimag(matrix))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;fingerprint=0_int64;ok=.false.;return;endif
    local_hash=0_int64
    do j=1,size(matrix,2);do i=1,size(row_ids)
      real_bits=transfer(real(matrix(i,j),real64),real_bits)
      imaginary_bits=transfer(aimag(matrix(i,j)),imaginary_bits)
      entry_hash=ieor(ishftc(real_bits,modulo(int(row_ids(i)),63)),&
        ishftc(imaginary_bits,modulo(j+11,63)))
      entry_hash=ieor(entry_hash,ishftc(ieor(row_ids(i),ishft(int(j,int64),21)),17))
      local_hash=ieor(local_hash,entry_hash)
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_int64;ok=.false.;return;endif
    fingerprint=ieor(fingerprint,ishftc(int(size(matrix,2),int64),29))
    fingerprint=ieor(fingerprint,int(z'6A09E667F3BCC909',int64))
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
  end subroutine fingerprint_distributed_matrix

  subroutine fingerprint_variational_operator(comm,payload,even_fingerprint,odd_fingerprint,ok)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(out)::even_fingerprint,odd_fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr
    integer(int64)::common_hash,local_hash,global_hash,bits
    common_hash=ieor(payload%basis_fingerprint,ishftc(payload%metric_fingerprint,11))
    common_hash=ieor(common_hash,ishftc(payload%face_fingerprint,23));local_hash=0_int64
    do i=1,size(payload%row_ids);do j=1,payload%global_count
      call hash_variational_complex(local_hash,payload%row_ids(i),j,payload%metric_rows(i,j))
      call hash_variational_complex(local_hash,payload%row_ids(i),j+payload%global_count,&
        payload%kinetic_rows(i,j))
      call hash_variational_complex(local_hash,payload%row_ids(i),j+2*payload%global_count,&
        payload%nonlocal_rows(i,j))
      call hash_variational_complex(local_hash,payload%row_ids(i),j+3*payload%global_count,&
        payload%sipg_rows(i,j))
    enddo;enddo
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;even_fingerprint=0_int64;odd_fingerprint=0_int64;ok=.false.;return;endif
    bits=int(payload%global_count,int64)
    even_fingerprint=ieor(global_hash,ishftc(bits,37))
    odd_fingerprint=ieor(ieor(global_hash,common_hash),ishftc(bits,37))
    if(even_fingerprint==0_int64)even_fingerprint=1543_int64
    if(odd_fingerprint==0_int64)odd_fingerprint=1543_int64
    ok=.true.
  contains
    subroutine hash_variational_complex(hash,row,column,value)
      integer(int64),intent(inout)::hash
      integer(int64),intent(in)::row
      integer,intent(in)::column
      complex(real64),intent(in)::value
      integer(int64)::value_bits
      value_bits=transfer(real(value,real64),value_bits)
      hash=ieor(hash,ishftc(ieor(value_bits,row),mod(7*column,63)))
      value_bits=transfer(aimag(value),value_bits)
      hash=ieor(hash,ishftc(ieor(value_bits,ishftc(row,9)),mod(13*column,63)))
    end subroutine hash_variational_complex
  end subroutine fingerprint_variational_operator

  subroutine fingerprint_ground_state(comm,payload,operator_fingerprint,fingerprint,ok)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::operator_fingerprint
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr
    integer(int64)::local_hash,global_hash,entry_hash,bits
    local_hash=0_int64
    do i=1,size(payload%row_ids);do j=1,payload%noccupied
      entry_hash=ieor(payload%row_ids(i),ishftc(int(j,int64),11))
      bits=transfer(real(payload%coefficients(i,j),real64),bits);entry_hash=ieor(entry_hash,ishftc(bits,19))
      bits=transfer(aimag(payload%coefficients(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,37))
      local_hash=ieor(local_hash,entry_hash)
    enddo;enddo
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_int64;ok=.false.;return;endif
    fingerprint=ieor(global_hash,payload%basis_fingerprint)
    fingerprint=ieor(fingerprint,ishftc(payload%metric_fingerprint,7))
    fingerprint=ieor(fingerprint,ishftc(operator_fingerprint,13))
    fingerprint=ieor(fingerprint,ishftc(payload%position_convention_fingerprint,23))
    do i=1,payload%noccupied
      bits=transfer(payload%occupations(i),bits);fingerprint=ieor(fingerprint,ishftc(bits,mod(5*i,63)))
      bits=transfer(payload%eigenvalues(i),bits);fingerprint=ieor(fingerprint,ishftc(bits,mod(9*i,63)))
    enddo
    if(fingerprint==0_int64)fingerprint=ieor(global_hash,719_int64)
    ok=.true.
  end subroutine fingerprint_ground_state

  pure integer(int64) function fingerprint_energy_window_receipt(payload) result(hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    hash=1469598103934665603_int64
    call mix_integer(payload%energy_window%solved_rank)
    call mix_integer(payload%energy_window%occupied_rank)
    call mix_integer(payload%energy_window%requested_rank)
    call mix_integer(payload%energy_window%boundary_cluster_rank)
    call mix_integer(payload%energy_window%certified_rank)
    call mix_integer(payload%energy_window%extension_states)
    call mix_integer(payload%symmetry_receipt%worst_operation)
    call mix_real(payload%energy_window%window_size)
    call mix_real(payload%energy_window%e_homo)
    call mix_real(payload%energy_window%requested_cutoff)
    call mix_real(payload%energy_window%certified_cutoff)
    call mix_real(payload%energy_window%extension_energy)
    call mix_real(payload%energy_window%proof_energy)
    call mix_real(payload%symmetry_receipt%occupied_subspace_defect)
    call mix_real(payload%symmetry_receipt%occupied_projector_defect)
    call mix_real(payload%symmetry_receipt%target_subspace_defect)
    call mix_real(payload%symmetry_receipt%target_energy_defect)
    call mix_real(payload%symmetry_receipt%density_defect)
    call mix_real(payload%symmetry_receipt%worst_operation_defect)
    call mix_real(max(payload%symmetry_receipt%worst_operation_defect,&
      payload%symmetry_receipt%occupied_projector_defect,payload%symmetry_receipt%density_defect))
    if(payload%energy_window%compatibility_dynamic_rank)hash=ieor(hash,97_int64)
    if(payload%energy_window%proof_state_present)hash=ieor(hash,193_int64)
    if(hash==0_int64)hash=389_int64
  contains
    pure subroutine mix_integer(value)
      integer,intent(in)::value
      hash=ieor(ishftc(hash,7),int(value,int64))
    end subroutine mix_integer
    pure subroutine mix_real(value)
      real(real64),intent(in)::value
      hash=ieor(ishftc(hash,7),transfer(value,0_int64))
    end subroutine mix_real
  end function fingerprint_energy_window_receipt

  pure integer(int64) function fingerprint_symmetry_receipt(payload,window_fingerprint,&
      certified_fingerprint) result(hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::window_fingerprint,certified_fingerprint
    hash=fingerprint_checkpoint_real([payload%symmetry_receipt%occupied_subspace_defect,&
      payload%symmetry_receipt%occupied_projector_defect,payload%symmetry_receipt%target_subspace_defect,&
      payload%symmetry_receipt%target_energy_defect,payload%symmetry_receipt%density_defect,&
      payload%symmetry_receipt%scalar_covariance_defect,payload%symmetry_receipt%vector_covariance_defect,&
      payload%symmetry_receipt%tensor_covariance_defect,payload%symmetry_receipt%final_basis_defect,&
      payload%symmetry_receipt%worst_operation_defect,payload%symmetry_receipt%maximum_physical_defect])
    hash=ieor(ishftc(hash,7),window_fingerprint)
    hash=ieor(ishftc(hash,7),certified_fingerprint)
    hash=ieor(ishftc(hash,7),int(payload%symmetry_receipt%worst_operation,int64))
    if(hash==0_int64)hash=1_int64
  end function fingerprint_symmetry_receipt

  subroutine reconstruct_certified_density(payload,density)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    real(real64),allocatable,intent(out)::density(:)
    call reconstruct_density_from_basis(payload%rt_space%basis_values,&
      payload%certified_basis%initial_occupied_amplitudes,payload%certified_basis%occupations,density)
  end subroutine reconstruct_certified_density

  subroutine reconstruct_density_from_basis(basis,amplitudes,occupations,density)
    complex(real64),intent(in)::basis(:,:),amplitudes(:,:)
    real(real64),intent(in)::occupations(:)
    real(real64),allocatable,intent(out)::density(:)
    complex(real64),allocatable::orbital_values(:)
    integer::point
    allocate(density(size(basis,2)),orbital_values(size(amplitudes,2)))
    do point=1,size(basis,2)
      orbital_values=matmul(basis(:,point),amplitudes)
      density(point)=sum(occupations*abs(orbital_values)**2)
    enddo
  end subroutine reconstruct_density_from_basis

  subroutine fingerprint_certified_rows(comm,row_ids,rows,global_count,fingerprint,ok)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::rows(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr_hash,ierr_count,local_count,global_rows
    integer(int64)::entry,local_hash,global_hash
    local_hash=0_int64
    do i=1,size(row_ids)
      entry=certified_mix_hash(int(z'243F6A8885A308D3',int64),row_ids(i))
      entry=certified_mix_hash(entry,int(size(rows,2),int64))
      do j=1,size(rows,2);entry=certified_mix_complex(entry,rows(i,j));enddo
      local_hash=ieor(local_hash,entry)
    enddo
    local_count=size(row_ids);global_hash=0_int64;global_rows=0
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr_hash)
    call MPI_Allreduce(local_count,global_rows,1,MPI_INTEGER,MPI_SUM,comm,ierr_count)
    if(ierr_hash/=MPI_SUCCESS.or.ierr_count/=MPI_SUCCESS.or.global_rows/=global_count)then
      fingerprint=0_int64;ok=.false.;return
    endif
    fingerprint=certified_mix_hash(int(z'13198A2E03707344',int64),int(global_count,int64))
    fingerprint=certified_mix_hash(fingerprint,int(size(rows,2),int64))
    fingerprint=certified_mix_hash(fingerprint,global_hash)
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
  end subroutine fingerprint_certified_rows

  pure integer(int64) function fingerprint_certified_localization(transform,payload) result(hash)
    complex(real64),intent(in)::transform(:,:)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    hash=certified_hash_complex_2(int(z'BB67AE8584CAA73B',int64),transform)
    hash=certified_hash_real_2(hash,payload%certified_basis%centers)
    hash=certified_hash_real_1(hash,payload%certified_basis%spreads_before)
    hash=certified_hash_real_1(hash,payload%certified_basis%spreads_after)
    hash=certified_mix_hash(hash,int(payload%certified_basis%localization_iterations,int64))
    hash=certified_mix_hash(hash,merge(1_int64,0_int64,&
      payload%certified_basis%localization_symmetry_constrained))
    hash=certified_mix_hash(hash,merge(1_int64,0_int64,payload%certified_basis%localization_converged))
    if(hash==0_int64)hash=1_int64
  end function fingerprint_certified_localization

  pure integer(int64) function fingerprint_certified_operator(payload,metric,hamiltonian,&
      scalar_operators,vector_operators,tensor_operators) result(hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    complex(real64),intent(in)::metric(:,:),hamiltonian(:,:),scalar_operators(:,:,:),&
      vector_operators(:,:,:,:),tensor_operators(:,:,:,:,:)
    hash=certified_hash_real_1(int(z'3C6EF372FE94F82B',int64),&
      payload%certified_basis%certified_eigenvalues)
    hash=certified_hash_complex_2(hash,payload%certified_basis%initial_occupied_amplitudes)
    hash=certified_hash_complex_2(hash,metric);hash=certified_hash_complex_2(hash,hamiltonian)
    hash=certified_hash_complex_3(hash,payload%rt_space%representation)
    hash=certified_hash_real_3(hash,payload%rt_space%cartesian_rotations)
    hash=certified_hash_complex_3(hash,scalar_operators)
    hash=certified_hash_complex_4(hash,vector_operators)
    hash=certified_hash_complex_5(hash,tensor_operators)
    if(hash==0_int64)hash=1_int64
  end function fingerprint_certified_operator

  pure integer(int64) function fingerprint_certified_receipt(payload,fingerprints) result(hash)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    integer(int64),intent(in)::fingerprints(7)
    hash=certified_mix_hash(int(z'A54FF53A5F1D36F1',int64),int(payload%global_count,int64))
    hash=certified_mix_hash(hash,int(payload%certified_basis%certified_count,int64))
    hash=certified_mix_hash(hash,int(payload%certified_basis%occupied_count,int64))
    hash=certified_mix_hash(hash,fingerprints(1));hash=certified_mix_hash(hash,fingerprints(2))
    hash=certified_mix_hash(hash,fingerprints(3));hash=certified_mix_hash(hash,fingerprints(6))
    hash=certified_mix_hash(hash,int(payload%certified_basis%localization_iterations,int64))
    hash=certified_mix_real(hash,payload%certified_basis%spread_before_total)
    hash=certified_mix_real(hash,payload%certified_basis%spread_after_total)
    hash=certified_mix_real(hash,payload%certified_basis%spread_improvement)
    hash=certified_mix_real(hash,payload%certified_basis%transform_unitarity_defect)
    hash=certified_mix_real(hash,payload%certified_basis%certified_metric_defect)
    hash=certified_mix_real(hash,payload%certified_basis%rt_metric_defect)
    hash=certified_mix_real(hash,payload%certified_basis%embedding_defect)
    hash=certified_mix_real(hash,payload%certified_basis%projector_invariance_defect)
    hash=certified_mix_real(hash,payload%certified_basis%target_symmetry_defect_before)
    hash=certified_mix_real(hash,payload%certified_basis%target_symmetry_defect_after)
    hash=certified_mix_real(hash,payload%certified_basis%energy_symmetry_defect_before)
    hash=certified_mix_real(hash,payload%certified_basis%energy_symmetry_defect_after)
    hash=certified_mix_real(hash,payload%certified_basis%symmetry_defect_invariance)
    hash=certified_mix_real(hash,payload%certified_basis%scalar_covariance_defect)
    hash=certified_mix_real(hash,payload%certified_basis%vector_covariance_defect)
    hash=certified_mix_real(hash,payload%certified_basis%tensor_covariance_defect)
    if(hash==0_int64)hash=1_int64
  end function fingerprint_certified_receipt

  pure integer(int64) function certified_mix_hash(left,right) result(mixed)
    integer(int64),intent(in)::left,right
    mixed=certified_add_modulo_64(ishftc(left,11),right)
    mixed=ieor(mixed,ishftc(mixed,25))
    mixed=certified_add_modulo_64(mixed,int(z'9E3779B97F4A7C15',int64))
    mixed=ieor(mixed,ishft(mixed,-27))
    mixed=certified_add_modulo_64(mixed,ishftc(right,17))
    mixed=ieor(mixed,ishftc(mixed,42))
  end function certified_mix_hash

  pure integer(int64) function certified_add_modulo_64(left,right) result(sum_value)
    integer(int64),intent(in)::left,right
    integer(int64),parameter::low_mask=int(z'00000000FFFFFFFF',int64)
    integer(int64)::low_value,high_value,carry
    low_value=iand(left,low_mask)+iand(right,low_mask);carry=ishft(low_value,-32)
    high_value=iand(ishft(left,-32),low_mask)+iand(ishft(right,-32),low_mask)+carry
    sum_value=ior(iand(low_value,low_mask),ishft(iand(high_value,low_mask),32))
  end function certified_add_modulo_64

  pure integer(int64) function certified_mix_real(hash,value) result(mixed)
    integer(int64),intent(in)::hash
    real(real64),intent(in)::value
    mixed=certified_mix_hash(hash,transfer(value,0_int64))
  end function certified_mix_real

  pure integer(int64) function certified_mix_complex(hash,value) result(mixed)
    integer(int64),intent(in)::hash
    complex(real64),intent(in)::value
    mixed=certified_mix_real(hash,real(value,real64));mixed=certified_mix_real(mixed,aimag(value))
  end function certified_mix_complex

  pure integer(int64) function certified_hash_real_1(seed,values) result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:)
    integer::i
    hash=certified_mix_hash(seed,int(size(values),int64))
    do i=1,size(values);hash=certified_mix_real(hash,values(i));enddo
  end function certified_hash_real_1

  pure integer(int64) function certified_hash_real_2(seed,values) result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:,:)
    integer::i,j
    hash=certified_mix_hash(certified_mix_hash(seed,int(size(values,1),int64)),int(size(values,2),int64))
    do j=1,size(values,2);do i=1,size(values,1);hash=certified_mix_real(hash,values(i,j));enddo;enddo
  end function certified_hash_real_2

  pure integer(int64) function certified_hash_real_3(seed,values) result(hash)
    integer(int64),intent(in)::seed
    real(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    hash=certified_mix_hash(certified_mix_hash(certified_mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64))
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=certified_mix_real(hash,values(i,j,k))
    enddo;enddo;enddo
  end function certified_hash_real_3

  pure integer(int64) function certified_hash_complex_2(seed,values) result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:)
    integer::i,j
    hash=certified_mix_hash(certified_mix_hash(seed,int(size(values,1),int64)),int(size(values,2),int64))
    do j=1,size(values,2);do i=1,size(values,1);hash=certified_mix_complex(hash,values(i,j));enddo;enddo
  end function certified_hash_complex_2

  pure integer(int64) function certified_hash_complex_3(seed,values) result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    hash=certified_mix_hash(certified_mix_hash(certified_mix_hash(seed,int(size(values,1),int64)),&
      int(size(values,2),int64)),int(size(values,3),int64))
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=certified_mix_complex(hash,values(i,j,k))
    enddo;enddo;enddo
  end function certified_hash_complex_3

  pure integer(int64) function certified_hash_complex_4(seed,values) result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:,:)
    integer::i,j,k,l
    hash=certified_mix_hash(certified_mix_hash(certified_mix_hash(certified_mix_hash(seed,&
      int(size(values,1),int64)),int(size(values,2),int64)),int(size(values,3),int64)),&
      int(size(values,4),int64))
    do l=1,size(values,4);do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      hash=certified_mix_complex(hash,values(i,j,k,l))
    enddo;enddo;enddo;enddo
  end function certified_hash_complex_4

  pure integer(int64) function certified_hash_complex_5(seed,values) result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::values(:,:,:,:,:)
    integer::i,j,k,l,m
    hash=certified_mix_hash(certified_mix_hash(certified_mix_hash(certified_mix_hash(&
      certified_mix_hash(seed,int(size(values,1),int64)),int(size(values,2),int64)),&
      int(size(values,3),int64)),int(size(values,4),int64)),int(size(values,5),int64))
    do m=1,size(values,5);do l=1,size(values,4);do k=1,size(values,3)
      do j=1,size(values,2);do i=1,size(values,1)
        hash=certified_mix_complex(hash,values(i,j,k,l,m))
      enddo;enddo
    enddo;enddo;enddo
  end function certified_hash_complex_5

  pure integer(int64) function fingerprint_checkpoint_integer(values) result(hash)
    integer,intent(in)::values(:)
    integer::i
    hash=int(z'BB67AE8584CAA73B',int64);hash=ieor(ishftc(hash,7),int(size(values),int64))
    do i=1,size(values);hash=ieor(ishftc(hash,11),ieor(int(values(i),int64),int(i,int64)));enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_integer

  pure integer(int64) function fingerprint_checkpoint_integer64(values) result(hash)
    integer(int64),intent(in)::values(:)
    integer::i
    hash=int(z'6A09E667F3BCC909',int64);hash=ieor(ishftc(hash,7),int(size(values),int64))
    do i=1,size(values);hash=ieor(ishftc(hash,11),ieor(values(i),int(i,int64)));enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_integer64

  pure integer(int64) function fingerprint_checkpoint_real(values) result(hash)
    real(real64),intent(in)::values(:)
    integer::i
    hash=int(z'9E3779B97F4A7C15',int64)
    do i=1,size(values);hash=ieor(ishftc(hash,11),ieor(transfer(values(i),0_int64),int(i,int64)));enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_real

  pure integer(int64) function fingerprint_checkpoint_complex_matrix(values) result(hash)
    complex(real64),intent(in)::values(:,:)
    integer::i,j
    integer(int64)::bits
    hash=int(z'3C6EF372FE94F82B',int64)
    hash=ieor(ishftc(hash,7),int(size(values,1),int64))
    hash=ieor(ishftc(hash,7),int(size(values,2),int64))
    do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(real(values(i,j)),bits);hash=ieor(ishftc(hash,11),bits)
      bits=transfer(aimag(values(i,j)),bits);hash=ieor(ishftc(hash,13),bits)
    enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_complex_matrix

  pure integer(int64) function fingerprint_checkpoint_complex_rank3(values) result(hash)
    complex(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    integer(int64)::bits
    hash=int(size(values,3),int64)
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(real(values(i,j,k),real64),bits);hash=ieor(ishftc(hash,9),bits)
      bits=transfer(aimag(values(i,j,k)),bits);hash=ieor(ishftc(hash,9),bits)
    enddo;enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_complex_rank3

  pure integer(int64) function fingerprint_checkpoint_complex_rank4(values) result(hash)
    complex(real64),intent(in)::values(:,:,:,:)
    integer::i,j,k,l
    integer(int64)::bits
    hash=int(z'A54FF53A5F1D36F1',int64)
    hash=ieor(hash,int(size(values),int64));hash=ieor(hash,ishftc(int(size(values,4),int64),17))
    do l=1,size(values,4);do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(real(values(i,j,k,l)),bits);hash=ieor(ishftc(hash,7),bits)
      bits=transfer(aimag(values(i,j,k,l)),bits);hash=ieor(ishftc(hash,11),bits)
    enddo;enddo;enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_complex_rank4

  pure integer(int64) function fingerprint_checkpoint_complex_rank5(values) result(hash)
    complex(real64),intent(in)::values(:,:,:,:,:)
    integer::i,j,k,l,m
    integer(int64)::bits
    hash=int(z'510E527FADE682D1',int64)
    hash=ieor(hash,int(size(values),int64));hash=ieor(hash,ishftc(int(size(values,5),int64),17))
    do m=1,size(values,5);do l=1,size(values,4);do k=1,size(values,3)
      do j=1,size(values,2);do i=1,size(values,1)
        bits=transfer(real(values(i,j,k,l,m)),bits);hash=ieor(ishftc(hash,7),bits)
        bits=transfer(aimag(values(i,j,k,l,m)),bits);hash=ieor(ishftc(hash,11),bits)
      enddo;enddo
    enddo;enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_complex_rank5

  pure integer(int64) function fingerprint_checkpoint_real_rank3(values) result(hash)
    real(real64),intent(in)::values(:,:,:)
    integer::i,j,k
    integer(int64)::bits
    hash=int(z'1F83D9ABFB41BD6B',int64)
    hash=ieor(ishftc(hash,7),int(size(values,1),int64))
    hash=ieor(ishftc(hash,7),int(size(values,2),int64))
    hash=ieor(ishftc(hash,7),int(size(values,3),int64))
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(values(i,j,k),bits);hash=ieor(ishftc(hash,11),bits)
    enddo;enddo;enddo
    if(hash==0_int64)hash=1_int64
  end function fingerprint_checkpoint_real_rank3

  subroutine fingerprint_grid_complex(comm,grid_ids,values,fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::grid_ids(:)
    complex(real64),intent(in)::values(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,j,ierr,local_bad,global_bad,local_rows,minimum_rows,maximum_rows
    integer(int64)::local_hash,bits,local_count,global_count
    local_bad=merge(0,1,size(values,2)==size(grid_ids))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;fingerprint=0_int64;ok=.false.;return;endif
    local_rows=size(values,1)
    call MPI_Allreduce(local_rows,minimum_rows,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_rows,maximum_rows,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_rows/=maximum_rows)then;fingerprint=0_int64;ok=.false.;return;endif
    local_hash=0_int64
    do i=1,size(grid_ids);do j=1,size(values,1)
      bits=transfer(real(values(j,i)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,grid_ids(i)),mod(7*j,63)))
      bits=transfer(aimag(values(j,i)),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,ishftc(grid_ids(i),17)),mod(11*j,63)))
    enddo;enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    local_count=int(size(grid_ids),int64);global_count=0_int64
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    fingerprint=ieor(fingerprint,ishftc(int(size(values,1),int64),31))
    fingerprint=ieor(fingerprint,ishftc(global_count,23));if(fingerprint==0_int64)fingerprint=1_int64
    ok=ierr==MPI_SUCCESS
  end subroutine fingerprint_grid_complex

  subroutine fingerprint_grid_real(comm,grid_ids,values,fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::grid_ids(:)
    real(real64),intent(in)::values(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,ierr,local_bad,global_bad
    integer(int64)::local_hash,bits,local_count,global_count
    local_bad=merge(0,1,size(values)==size(grid_ids))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;fingerprint=0_int64;ok=.false.;return;endif
    local_hash=0_int64
    do i=1,size(values)
      bits=transfer(values(i),bits)
      local_hash=ieor(local_hash,ishftc(ieor(bits,ishftc(grid_ids(i),17)),&
        int(modulo(grid_ids(i),63_int64))))
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    local_count=int(size(values),int64);global_count=0_int64
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    fingerprint=ieor(fingerprint,ishftc(global_count,31));if(fingerprint==0_int64)fingerprint=1_int64
    ok=ierr==MPI_SUCCESS
  end subroutine fingerprint_grid_real

  subroutine fingerprint_grid_integer(comm,grid_ids,values,fingerprint,ok)
    integer,intent(in)::comm,values(:)
    integer(int64),intent(in)::grid_ids(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,ierr,local_bad,global_bad
    integer(int64)::local_hash,entry,local_count,global_count
    local_bad=merge(0,1,size(values)==size(grid_ids).and.all(values>0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;fingerprint=0_int64;ok=.false.;return;endif
    local_hash=0_int64
    do i=1,size(values)
      entry=ieor(grid_ids(i),ishftc(int(values(i),int64),17))
      local_hash=ieor(local_hash,ishftc(entry,int(modulo(grid_ids(i),63_int64))))
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    local_count=int(size(values),int64);global_count=0_int64
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    fingerprint=ieor(fingerprint,ishftc(global_count,31));if(fingerprint==0_int64)fingerprint=1_int64
    ok=ierr==MPI_SUCCESS
  end subroutine fingerprint_grid_integer

  pure real(real64) function maximum_complex_matrix(values) result(maximum_value)
    complex(real64),intent(in)::values(:,:)
    maximum_value=0d0;if(size(values)>0)maximum_value=maxval(abs(values))
  end function maximum_complex_matrix

  pure real(real64) function maximum_complex_rank3(values) result(maximum_value)
    complex(real64),intent(in)::values(:,:,:)
    maximum_value=0d0;if(size(values)>0)maximum_value=maxval(abs(values))
  end function maximum_complex_rank3

  pure real(real64) function maximum_complex_rank4(values) result(maximum_value)
    complex(real64),intent(in)::values(:,:,:,:)
    maximum_value=0d0;if(size(values)>0)maximum_value=maxval(abs(values))
  end function maximum_complex_rank4
#endif
end module rt_dg_hybrid_initialization
