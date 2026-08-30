#include "config.h"
module rt_dg_hybrid_initialization
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_checkpoint,only:s_rt_dg_hybrid_ground_state_payload,&
    read_rt_dg_hybrid_ground_state_checkpoint,read_rt_dg_hybrid_ground_state_checkpoint_coalesced
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter::salmon_xctype_none=0,salmon_xctype_pz=1,salmon_xctype_pzm=2,&
    salmon_xctype_pbe=3,salmon_xctype_pw=7
  type,public::s_rt_dg_hybrid_state
    logical::valid=.false.,initial_invariants_valid=.false.
    integer::global_count=0,noccupied=0
    integer(int64)::payload_fingerprint=0_int64,operator_structure_fingerprint=0_int64,&
      operator_value_fingerprint=0_int64,scope_fingerprint=0_int64
    type(s_dg_hybrid_sparse_metric)::metric
    type(s_dg_hybrid_sparse_operators)::operators
    integer(int64),allocatable::owned_row_ids(:),grid_ids(:)
    complex(real64),allocatable::coefficients(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      local_rows(:,:),sipg_rows(:,:),basis_values(:,:)
    real(real64),allocatable::density(:),grid_weights(:),occupations(:),eigenvalues(:)
  end type s_rt_dg_hybrid_state
  public::initialize_rt_dg_hybrid_from_checkpoint,fingerprint_rt_dg_hybrid_scope
contains
  subroutine initialize_rt_dg_hybrid_from_checkpoint(comm,path,theory,periodic,nspin,spinorbit,plus_u,hse,&
      fix_func,jm,xctype,state,ok,message)
    integer,intent(in)::comm,nspin,xctype(:)
    character(*),intent(in)::path,theory
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    type(s_rt_dg_hybrid_state),intent(out)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_rt_dg_hybrid_ground_state_payload)::payload
    integer::rank,nproc,ierr,n,nocc,nowned,row,i,j,position,local_bad,global_bad,file_nproc,read_comm,active_xc
    integer(int64)::payload_fingerprint
    complex(real64),allocatable::metric(:,:),kinetic(:,:),nonlocal(:,:),local(:,:),sipg(:,:),hamiltonian(:,:),coeff(:,:),&
      position_rows(:,:,:)
    integer,allocatable::metric_graph(:,:),operator_graph(:,:)
    real(real64)::residual,orthogonality,electron_defect,hermiticity,covariance,projector_covariance,scale,&
      density_charge
    ok=.false.;message='';state%valid=.false.
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='hybrid RT communicator query failed';return;endif
    local_bad=0
    if((trim(theory)/='tddft_response'.and.trim(theory)/='tddft_pulse').or..not.periodic.or.nspin/=1.or.&
      spinorbit.or.plus_u.or.hse.or.fix_func.or.jm.or.size(xctype)<1) local_bad=1
    active_xc=0
    do i=1,size(xctype)
      if(xctype(i)/=salmon_xctype_none)active_xc=active_xc+1
      if(.not.(xctype(i)==salmon_xctype_none.or.xctype(i)==salmon_xctype_pz.or.&
        xctype(i)==salmon_xctype_pzm.or.xctype(i)==salmon_xctype_pbe.or.xctype(i)==salmon_xctype_pw))local_bad=1
    enddo
    if(active_xc==0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='unsupported local hybrid RT scope';return;endif
    call probe_ground_state_rank_count(comm,path,file_nproc,ok,message)
    if(.not.ok)return
    if(file_nproc>nproc)then
      call read_rt_dg_hybrid_ground_state_checkpoint_coalesced(comm,path,payload,payload_fingerprint,ok,message)
    else
      call MPI_Comm_split(comm,merge(0,MPI_UNDEFINED,rank<file_nproc),rank,read_comm,ierr)
      if(rank<file_nproc)then
        call read_rt_dg_hybrid_ground_state_checkpoint(read_comm,path,payload,payload_fingerprint,ok,message)
        call MPI_Comm_free(read_comm,ierr)
      else
        ok=.true.;payload_fingerprint=0_int64
      endif
      local_bad=merge(0,1,ok);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;return;endif
      call broadcast_ground_state_for_expansion(comm,rank,file_nproc,payload,payload_fingerprint,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT checkpoint rank expansion failed';return;endif
    endif
    local_bad=merge(0,1,ok);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;return;endif
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
    n=payload%global_count;nocc=payload%noccupied
    allocate(metric(n,n),kinetic(n,n),nonlocal(n,n),local(n,n),sipg(n,n),hamiltonian(n,n),coeff(n,nocc),&
      position_rows(3,n,n),metric_graph(n,n),operator_graph(n,n));metric=(0d0,0d0);kinetic=(0d0,0d0)
    nonlocal=(0d0,0d0);local=(0d0,0d0)
    sipg=(0d0,0d0);hamiltonian=(0d0,0d0);coeff=(0d0,0d0);position_rows=(0d0,0d0)
    metric_graph=0;operator_graph=0
    do i=1,size(payload%row_ids)
      row=int(payload%row_ids(i));metric(row,:)=payload%metric_rows(i,:);kinetic(row,:)=payload%kinetic_rows(i,:)
      nonlocal(row,:)=payload%nonlocal_rows(i,:);local(row,:)=payload%local_rows(i,:);sipg(row,:)=payload%sipg_rows(i,:)
      hamiltonian(row,:)=payload%hamiltonian_rows(i,:);coeff(row,:)=payload%coefficients(i,:)
      position_rows(:,row,:)=payload%position_rows(:,i,:)
      do j=payload%metric_row_offsets(i),payload%metric_row_offsets(i+1)-1
        metric_graph(row,payload%metric_column_ids(j))=1
      enddo
      do j=payload%operator_row_offsets(i),payload%operator_row_offsets(i+1)-1
        operator_graph(row,payload%operator_column_ids(j))=1
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,metric,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,kinetic,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,nonlocal,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,local,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,sipg,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,hamiltonian,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,coeff,n*nocc,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,position_rows,3*n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,metric_graph,n*n,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,operator_graph,n*n,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT global payload reconstruction failed';return;endif
    density_charge=sum(payload%density*payload%grid_weights)
    call MPI_Allreduce(MPI_IN_PLACE,density_charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call validate_invariants(metric,hamiltonian,coeff,payload%occupations,payload%eigenvalues,density_charge,&
      payload%symmetry_representation,residual,orthogonality,electron_defect,hermiticity,covariance,projector_covariance)
    scale=max(1d0,maxval(abs(hamiltonian)))
    local_bad=merge(0,1,residual<=1d-11*scale.and.orthogonality<=1d-11.and.electron_defect<=1d-11.and.&
      hermiticity<=1d-11*scale.and.covariance<=1d-11*scale.and.projector_covariance<=1d-11)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='invalid hybrid RT startup invariants';return;endif
    nowned=count([(mod(n-row,nproc)==rank,row=1,n)])
    allocate(state%owned_row_ids(nowned),state%coefficients(nowned,nocc),state%kinetic_rows(nowned,n),&
      state%nonlocal_rows(nowned,n),state%local_rows(nowned,n),state%sipg_rows(nowned,n))
    position=0
    do row=1,n
      if(mod(n-row,nproc)/=rank)cycle
      position=position+1;state%owned_row_ids(position)=row;state%coefficients(position,:)=coeff(row,:)
      state%kinetic_rows(position,:)=kinetic(row,:);state%nonlocal_rows(position,:)=nonlocal(row,:)
      state%local_rows(position,:)=local(row,:);state%sipg_rows(position,:)=sipg(row,:)
    enddo
    call build_sparse_state(n,state%owned_row_ids,metric,hamiltonian,position_rows,metric_graph,operator_graph,payload,state)
    allocate(state%grid_ids,source=payload%grid_ids);allocate(state%density,source=payload%density)
    allocate(state%grid_weights,source=payload%grid_weights);allocate(state%basis_values,source=payload%basis_values)
    allocate(state%occupations,source=payload%occupations);allocate(state%eigenvalues,source=payload%eigenvalues)
    state%global_count=n;state%noccupied=nocc;state%payload_fingerprint=payload_fingerprint
    state%operator_structure_fingerprint=payload%operator_structure_fingerprint
    state%operator_value_fingerprint=payload%operator_value_fingerprint;state%scope_fingerprint=payload%scope_fingerprint
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
      if(status==0.and.(magic/='SALMON_DG_GS001 '.or.version/=1))status=1
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr);call MPI_Bcast(file_nproc,1,MPI_INTEGER,0,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.status==0.and.file_nproc>0
    if(ok)then;message='';else;message='cannot probe complete DG ground-state checkpoint';endif
  end subroutine probe_ground_state_rank_count

  subroutine broadcast_ground_state_for_expansion(comm,rank,file_nproc,p,fp,ierr)
    integer,intent(in)::comm,rank,file_nproc
    type(s_rt_dg_hybrid_ground_state_payload),intent(inout)::p
    integer(int64),intent(inout)::fp
    integer,intent(out)::ierr
    integer::header(4);logical::flags(4)
    integer(int64)::fingerprints(18)
    if(rank==0)then
      header=[p%global_count,p%noccupied,p%operation_count,p%nonidentity_operation_count]
      flags=[p%valid,p%final_refresh_complete,p%analysis_complete,p%identity_only]
      fingerprints=[p%catalog_fingerprint,p%state_fingerprint,p%metric_fingerprint,p%operator_structure_fingerprint,&
        p%operator_value_fingerprint,p%kinetic_fingerprint,p%nonlocal_fingerprint,p%local_fingerprint,p%sipg_fingerprint,&
        p%basis_fingerprint,p%face_fingerprint,p%dc_seed_fingerprint,p%continuation_fingerprint,p%scope_fingerprint,&
        p%analysis_fingerprint,p%selection_fingerprint,p%pseudopotential_fingerprint,p%energy_fingerprint]
    endif
    call MPI_Bcast(header,4,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(flags,4,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(fingerprints,18,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(fp,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank>=file_nproc)then
      p%global_count=header(1);p%noccupied=header(2);p%operation_count=header(3);p%nonidentity_operation_count=header(4)
      p%valid=flags(1);p%final_refresh_complete=flags(2);p%analysis_complete=flags(3);p%identity_only=flags(4)
      p%catalog_fingerprint=fingerprints(1);p%state_fingerprint=fingerprints(2);p%metric_fingerprint=fingerprints(3)
      p%operator_structure_fingerprint=fingerprints(4);p%operator_value_fingerprint=fingerprints(5)
      p%kinetic_fingerprint=fingerprints(6);p%nonlocal_fingerprint=fingerprints(7);p%local_fingerprint=fingerprints(8)
      p%sipg_fingerprint=fingerprints(9);p%basis_fingerprint=fingerprints(10);p%face_fingerprint=fingerprints(11)
      p%dc_seed_fingerprint=fingerprints(12);p%continuation_fingerprint=fingerprints(13);p%scope_fingerprint=fingerprints(14)
      p%analysis_fingerprint=fingerprints(15);p%selection_fingerprint=fingerprints(16)
      p%pseudopotential_fingerprint=fingerprints(17);p%energy_fingerprint=fingerprints(18)
      allocate(p%row_ids(0),p%metric_rows(0,header(1)),p%kinetic_rows(0,header(1)),p%nonlocal_rows(0,header(1)),&
        p%local_rows(0,header(1)),p%sipg_rows(0,header(1)),p%hamiltonian_rows(0,header(1)),p%coefficients(0,header(2)),&
        p%position_rows(3,0,header(1)),p%metric_row_offsets(1),p%metric_column_ids(0),p%operator_row_offsets(1),&
        p%operator_column_ids(0),p%grid_ids(0),p%grid_weights(0),p%density(0),p%basis_values(header(1),0))
      p%metric_row_offsets=1;p%operator_row_offsets=1
    endif
    call bcast_i1(p%scope_selectors);call bcast_i1(p%xc_types);call bcast_r1(p%occupations);call bcast_r1(p%eigenvalues)
    call bcast_z3(p%symmetry_representation)
  contains
    subroutine bcast_i1(a)
      integer,allocatable,intent(inout)::a(:);integer::count
      if(rank==0)count=size(a);call MPI_Bcast(count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(.not.allocated(a))allocate(a(count));call MPI_Bcast(a,count,MPI_INTEGER,0,comm,ierr)
    end subroutine
    subroutine bcast_r1(a)
      real(real64),allocatable,intent(inout)::a(:);integer::count
      if(rank==0)count=size(a);call MPI_Bcast(count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(.not.allocated(a))allocate(a(count));call MPI_Bcast(a,count,MPI_DOUBLE_PRECISION,0,comm,ierr)
    end subroutine
    subroutine bcast_z3(a)
      complex(real64),allocatable,intent(inout)::a(:,:,:);integer::dims(3)
      if(rank==0)dims=shape(a);call MPI_Bcast(dims,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(.not.allocated(a))allocate(a(dims(1),dims(2),dims(3)))
      call MPI_Bcast(a,product(dims),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    end subroutine
  end subroutine broadcast_ground_state_for_expansion

  subroutine validate_invariants(s,h,c,occupations,eigenvalues,density_charge,representation,residual,orthogonality,&
      electron_defect,hermiticity,covariance,projector_covariance)
    complex(real64),intent(in)::s(:,:),h(:,:),c(:,:)
    real(real64),intent(in)::occupations(:),eigenvalues(:)
    real(real64),intent(in)::density_charge
    complex(real64),intent(in)::representation(:,:,:)
    real(real64),intent(out)::residual,orthogonality,electron_defect,hermiticity,covariance,projector_covariance
    complex(real64),allocatable::r(:,:),gram(:,:),projector(:,:),transformed(:,:)
    integer::i,op
    allocate(r(size(c,1),size(c,2)),gram(size(c,2),size(c,2)))
    r=matmul(h,c)-matmul(s,c)
    do i=1,size(c,2);r(:,i)=r(:,i)*eigenvalues(i);enddo
    ! Restore H*C - S*C*epsilon after column scaling above.
    r=matmul(h,c)
    do i=1,size(c,2);r(:,i)=r(:,i)-eigenvalues(i)*matmul(s,c(:,i));enddo
    gram=matmul(conjg(transpose(c)),matmul(s,c))
    do i=1,size(gram,1);gram(i,i)=gram(i,i)-(1d0,0d0);enddo
    residual=maxval(abs(r));orthogonality=maxval(abs(gram));electron_defect=abs(sum(occupations)-density_charge)
    hermiticity=max(maxval(abs(h-conjg(transpose(h)))),maxval(abs(s-conjg(transpose(s)))))
    allocate(projector(size(c,1),size(c,1)),transformed(size(c,1),size(c,1)))
    projector=(0d0,0d0)
    do i=1,size(c,2);projector=projector+occupations(i)*matmul(reshape(c(:,i),[size(c,1),1]),&
      reshape(conjg(c(:,i)),[1,size(c,1)]));enddo
    covariance=0d0;projector_covariance=0d0
    do op=1,size(representation,3)
      transformed=matmul(conjg(transpose(representation(:,:,op))),matmul(h,representation(:,:,op)))
      covariance=max(covariance,maxval(abs(transformed-h)))
      transformed=matmul(conjg(transpose(representation(:,:,op))),matmul(s,representation(:,:,op)))
      covariance=max(covariance,maxval(abs(transformed-s)))
      transformed=matmul(representation(:,:,op),matmul(projector,conjg(transpose(representation(:,:,op)))))
      projector_covariance=max(projector_covariance,maxval(abs(transformed-projector)))
    enddo
  end subroutine validate_invariants
  subroutine build_sparse_state(n,owned,s,h,position,metric_graph,operator_graph,payload,state)
    integer,intent(in)::n
    integer(int64),intent(in)::owned(:)
    complex(real64),intent(in)::s(:,:),h(:,:),position(:,:,:)
    integer,intent(in)::metric_graph(:,:),operator_graph(:,:)
    type(s_rt_dg_hybrid_ground_state_payload),intent(in)::payload
    type(s_rt_dg_hybrid_state),intent(inout)::state
    integer::i,j,edge,row,nowned,metric_edges,operator_edges
    nowned=size(owned)
    metric_edges=0;operator_edges=0
    do i=1,nowned
      row=int(owned(i));metric_edges=metric_edges+count(metric_graph(row,:)==1)
      operator_edges=operator_edges+count(operator_graph(row,:)==1)
    enddo
    allocate(state%metric%owned_row_ids(nowned),state%metric%row_offsets(nowned+1),&
      state%metric%column_ids(metric_edges),state%metric%values(metric_edges),state%metric%active_rows(n),&
      state%metric%packet_ids(n))
    allocate(state%operators%owned_row_ids(nowned),state%operators%row_offsets(nowned+1),&
      state%operators%column_ids(operator_edges),state%operators%metric_values(operator_edges),&
      state%operators%hamiltonian_values(operator_edges),state%operators%position_values(3,operator_edges))
    edge=0;state%metric%row_offsets(1)=1
    do i=1,nowned
      row=int(owned(i))
      do j=1,n
        if(metric_graph(row,j)==0)cycle
        edge=edge+1;state%metric%column_ids(edge)=j;state%metric%values(edge)=s(row,j)
      enddo
      state%metric%row_offsets(i+1)=edge+1
    enddo
    edge=0;state%operators%row_offsets(1)=1
    do i=1,nowned
      row=int(owned(i))
      do j=1,n
        if(operator_graph(row,j)==0)cycle
        edge=edge+1;state%operators%column_ids(edge)=j;state%operators%metric_values(edge)=s(row,j)
        state%operators%hamiltonian_values(edge)=h(row,j);state%operators%position_values(:,edge)=position(:,row,j)
      enddo
      state%operators%row_offsets(i+1)=edge+1
    enddo
    state%metric%owned_row_ids=owned;state%metric%global_count=n;state%metric%numerical_rank=n
    state%metric%active_rows=.true.;state%metric%packet_ids=1;state%metric%valid=.true.
    state%metric%fingerprint=payload%metric_fingerprint;state%metric%condition_estimate=1d0
    state%metric%maximum_value=maxval(abs(s));state%metric%max_row_nnz=n
    state%operators%owned_row_ids=owned;state%operators%global_count=n;state%operators%valid=.true.
    state%operators%metric_fingerprint=payload%metric_fingerprint
    state%operators%fingerprint=payload%operator_structure_fingerprint
  end subroutine build_sparse_state
#endif
end module rt_dg_hybrid_initialization
