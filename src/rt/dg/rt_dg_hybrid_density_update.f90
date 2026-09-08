#include "config.h"
module rt_dg_hybrid_density_update
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state
  use rt_dg_hybrid_sparse_projection,only:validate_rt_dg_hybrid_sparse_hermiticity
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine project_rt_density(row_ids,row_offsets,column_ids,grid_ids,density,local_values,ok,message)
      import::int64,real64
      integer(int64),intent(in)::row_ids(:)
      integer,intent(in)::row_offsets(:),column_ids(:)
      integer(int64),intent(in)::grid_ids(:)
      real(real64),intent(in)::density(:)
      complex(real64),intent(out)::local_values(:)
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine project_rt_density
  end interface
  public::update_rt_dg_hybrid_density,reconstruct_rt_dg_hybrid_density
contains
  subroutine reconstruct_rt_dg_hybrid_density(comm,state,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::basis_batch(:,:),orbital_values(:,:),reduced_orbital_values(:,:)
    integer,allocatable::point_counts(:)
    integer::i,p,owner,rank,nproc,point_count,max_point_count,ierr,local_bad,global_bad
    state%density_freshly_reconstructed=.false.
    call validate_certified_rt_state(comm,state,ok,message);if(.not.ok)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(point_counts(nproc));point_count=size(state%grid_ids)
    call MPI_Allgather(point_count,1,MPI_INTEGER,point_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT grid-count exchange failed';return;endif
    max_point_count=max(1,maxval(point_counts))
    allocate(basis_batch(state%certified_rank,max_point_count),&
      orbital_values(state%noccupied,max_point_count),reduced_orbital_values(state%noccupied,max_point_count))
    do owner=0,nproc-1
      point_count=point_counts(owner+1);if(point_count==0)cycle
      basis_batch(:,1:point_count)=(0d0,0d0)
      if(rank==owner)basis_batch(:,1:point_count)=state%basis_values
      call MPI_Bcast(basis_batch,state%certified_rank*point_count,MPI_DOUBLE_COMPLEX,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT basis-slab exchange failed';return;endif
      orbital_values(:,1:point_count)=(0d0,0d0)
      do i=1,size(state%owned_row_ids)
        do p=1,point_count
          orbital_values(:,p)=orbital_values(:,p)+&
            basis_batch(int(state%owned_row_ids(i)),p)*state%coefficients(i,:)
        enddo
      enddo
      call MPI_Reduce(orbital_values,reduced_orbital_values,state%noccupied*point_count,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT grid-local orbital reduction failed';return;endif
      if(rank==owner)then
        do p=1,point_count
          state%density(p)=sum(state%occupations*abs(reduced_orbital_values(:,p))**2)
        enddo
      endif
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(state%density)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0
    if(ok)then
      state%density_freshly_reconstructed=.true.;message=''
    else
      message='nonfinite hybrid RT reconstructed density'
    endif
#else
    ok=.false.;message='hybrid RT density reconstruction requires MPI'
#endif
  end subroutine reconstruct_rt_dg_hybrid_density

  subroutine update_rt_dg_hybrid_density(comm,state,density,project_local,ok,message,establish_fixed_density_reference)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_state),intent(inout)::state
    real(real64),intent(in)::density(:)
    procedure(project_rt_density)::project_local
    logical,intent(in),optional::establish_fixed_density_reference
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::new_local(:),new_h(:)
    integer::i,j,edge,ierr,local_bad,global_bad
    integer(int64)::local_hash,global_xor,global_sum,bits,pair_hash
    logical::callback_ok,establish_reference
    character(256)::callback_message
    ok=.false.;message='';local_bad=0;establish_reference=.false.
    if(present(establish_fixed_density_reference))establish_reference=establish_fixed_density_reference
    if(establish_reference.and.state%fixed_density_reference_valid)local_bad=1
    if(.not.state%density_freshly_reconstructed)then
      call validate_certified_rt_state(comm,state,ok,message);if(.not.ok)return
    endif
    state%density_freshly_reconstructed=.false.
    ok=.false.
    if(size(density)/=size(state%grid_ids).or..not.all(ieee_is_finite(density)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT density update';return;endif
    allocate(new_local(size(state%operators%column_ids)),new_h(size(state%operators%column_ids)))
    new_local=(0d0,0d0);callback_ok=.false.;callback_message=''
    call project_local(state%owned_row_ids,state%operators%row_offsets,state%operators%column_ids,&
      state%grid_ids,density,new_local,callback_ok,callback_message)
    local_bad=merge(0,1,callback_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid RT local projection failed';return;endif
    local_bad=merge(0,1,all(ieee_is_finite(real(new_local))).and.all(ieee_is_finite(aimag(new_local))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='nonfinite hybrid RT local projection';return;endif
    do i=1,size(state%owned_row_ids)
      do edge=state%operators%row_offsets(i),state%operators%row_offsets(i+1)-1
        j=state%operators%column_ids(edge)
        new_h(edge)=state%kinetic_rows(i,j)+state%nonlocal_rows(i,j)+new_local(edge)+state%sipg_rows(i,j)
      enddo
    enddo
    if(establish_reference)then
      allocate(state%local_reference_correction(size(new_local)),&
        state%hamiltonian_reference_correction(size(new_h)))
      state%local_reference_correction=state%local_rows-new_local
      state%hamiltonian_reference_correction=state%operators%hamiltonian_values-new_h
      state%reference_refresh_defect=0d0;state%reference_refresh_scale=1d0
      if(size(new_h)>0)then
        state%reference_refresh_defect=maxval(abs(state%hamiltonian_reference_correction))
        state%reference_refresh_scale=max(1d0,maxval(abs(state%operators%hamiltonian_values)))
      endif
      state%fixed_density_reference_valid=.true.
    endif
    if(state%fixed_density_reference_valid)then
      if(size(state%local_reference_correction)/=size(new_local).or.&
          size(state%hamiltonian_reference_correction)/=size(new_h))then
        message='invalid hybrid RT fixed-density reference correction';return
      endif
      new_local=new_local+state%local_reference_correction
      new_h=new_h+state%hamiltonian_reference_correction
    endif
    local_bad=merge(0,1,all(ieee_is_finite(real(new_h))).and.all(ieee_is_finite(aimag(new_h))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='nonfinite hybrid RT updated Hamiltonian';return;endif
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,state%certified_rank,state%owned_row_ids,&
      state%operators%row_offsets,state%operators%column_ids,new_h,100d0*epsilon(1d0),ok,message)
    if(.not.ok)then;message='hybrid RT Hamiltonian Hermiticity failed: '//trim(message);return;endif
    local_hash=0_int64;global_sum=0_int64
    do i=1,size(state%owned_row_ids)
      do edge=state%operators%row_offsets(i),state%operators%row_offsets(i+1)-1
        j=state%operators%column_ids(edge)
        bits=transfer(real(new_h(edge)),bits);pair_hash=ieor(ishftc(state%owned_row_ids(i),17),int(j,int64))
        pair_hash=ieor(ishftc(pair_hash,9),bits)
        bits=transfer(aimag(new_h(edge)),bits);pair_hash=ieor(ishftc(pair_hash,9),bits)
        local_hash=ieor(local_hash,pair_hash);global_sum=global_sum+pair_hash
      enddo
    enddo
    call MPI_Allreduce(local_hash,global_xor,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,global_sum,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='hybrid RT value fingerprint reduction failed';return;endif
    state%local_rows=new_local;state%density=density
    do i=1,size(state%owned_row_ids)
      do edge=state%operators%row_offsets(i),state%operators%row_offsets(i+1)-1
        state%operators%hamiltonian_values(edge)=new_h(edge)
      enddo
    enddo
    state%operator_value_fingerprint=ieor(global_xor,ishftc(global_sum,13))
    if(state%operator_value_fingerprint==0_int64)state%operator_value_fingerprint=1_int64
    state%operators%fingerprint=state%operator_structure_fingerprint
    ok=.true.;message=''
#else
    ok=.false.;message='hybrid RT density update requires MPI'
#endif
  end subroutine update_rt_dg_hybrid_density

  subroutine validate_certified_rt_state(comm,state,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_state),intent(in)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::dimensions(3),minimum_dimensions(3),maximum_dimensions(3)
    integer::r,nocc,nowned,npoint,i,j,edge,ierr,local_bad,global_bad
    integer,allocatable::local_counts(:),global_counts(:)
    ok=.false.;message=''
    dimensions=[state%certified_rank,state%global_count,state%noccupied]
    call MPI_Allreduce(dimensions,minimum_dimensions,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(dimensions,maximum_dimensions,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions))then
      message='rank-disagreeing certified hybrid RT dimensions';return
    endif
    r=minimum_dimensions(1);nocc=minimum_dimensions(3)
    local_bad=merge(0,1,state%valid.and.r>0.and.minimum_dimensions(2)==r.and.nocc>0.and.nocc<=r)
    if(.not.allocated(state%owned_row_ids).or..not.allocated(state%coefficients).or.&
      .not.allocated(state%kinetic_rows).or..not.allocated(state%nonlocal_rows).or.&
      .not.allocated(state%local_rows).or..not.allocated(state%sipg_rows).or.&
      .not.allocated(state%basis_values).or..not.allocated(state%density).or.&
      .not.allocated(state%grid_ids).or..not.allocated(state%grid_weights).or.&
      .not.allocated(state%occupations).or..not.allocated(state%eigenvalues).or.&
      .not.allocated(state%metric%owned_row_ids).or..not.allocated(state%metric%row_offsets).or.&
      .not.allocated(state%metric%column_ids).or..not.allocated(state%metric%values).or.&
      .not.allocated(state%metric%active_rows).or..not.allocated(state%metric%packet_ids).or.&
      .not.allocated(state%operators%owned_row_ids).or..not.allocated(state%operators%row_offsets).or.&
      .not.allocated(state%operators%column_ids).or..not.allocated(state%operators%metric_values).or.&
      .not.allocated(state%operators%hamiltonian_values).or.&
      .not.allocated(state%operators%position_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid certified hybrid RT state allocation';return;endif

    nowned=size(state%owned_row_ids);npoint=size(state%grid_ids);local_bad=0
    if(size(state%coefficients,1)/=nowned.or.size(state%coefficients,2)/=nocc)local_bad=1
    if(size(state%kinetic_rows,1)/=nowned.or.size(state%kinetic_rows,2)/=r.or.&
      size(state%nonlocal_rows,1)/=nowned.or.size(state%nonlocal_rows,2)/=r.or.&
      size(state%local_rows)/=size(state%operators%column_ids).or.&
      size(state%sipg_rows,1)/=nowned.or.size(state%sipg_rows,2)/=r)local_bad=1
    if(size(state%basis_values,1)/=r.or.size(state%basis_values,2)/=npoint.or.&
      size(state%density)/=npoint.or.size(state%grid_weights)/=npoint.or.&
      size(state%occupations)/=nocc.or.size(state%eigenvalues)/=r)local_bad=1
    if(.not.state%metric%valid.or.state%metric%global_count/=r.or.&
      size(state%metric%owned_row_ids)/=nowned.or.size(state%metric%row_offsets)/=nowned+1.or.&
      size(state%metric%values)/=size(state%metric%column_ids).or.&
      size(state%metric%active_rows)/=r.or.size(state%metric%packet_ids)/=r)local_bad=1
    if(.not.state%operators%valid.or.state%operators%global_count/=r.or.&
      state%operators%metric_fingerprint/=state%metric%fingerprint.or.&
      size(state%operators%owned_row_ids)/=nowned.or.size(state%operators%row_offsets)/=nowned+1.or.&
      size(state%operators%metric_values)/=size(state%operators%column_ids).or.&
      size(state%operators%hamiltonian_values)/=size(state%operators%column_ids).or.&
      size(state%operators%position_values,1)/=3.or.&
      size(state%operators%position_values,2)/=size(state%operators%column_ids))local_bad=1
    if(size(state%metric%owned_row_ids)==nowned)then
      if(any(state%metric%owned_row_ids/=state%owned_row_ids))local_bad=1
    endif
    if(size(state%operators%owned_row_ids)==nowned)then
      if(any(state%operators%owned_row_ids/=state%owned_row_ids))local_bad=1
    endif
    if(.not.all(ieee_is_finite(real(state%coefficients))).or.&
      .not.all(ieee_is_finite(aimag(state%coefficients))).or.&
      .not.all(ieee_is_finite(real(state%basis_values))).or.&
      .not.all(ieee_is_finite(aimag(state%basis_values))).or.&
      .not.all(ieee_is_finite(state%density)).or..not.all(ieee_is_finite(state%grid_weights)).or.&
      .not.all(ieee_is_finite(state%occupations)).or..not.all(ieee_is_finite(state%eigenvalues)))local_bad=1
    if(.not.all(ieee_is_finite(real(state%kinetic_rows))).or.&
      .not.all(ieee_is_finite(aimag(state%kinetic_rows))).or.&
      .not.all(ieee_is_finite(real(state%nonlocal_rows))).or.&
      .not.all(ieee_is_finite(aimag(state%nonlocal_rows))).or.&
      .not.all(ieee_is_finite(real(state%local_rows))).or.&
      .not.all(ieee_is_finite(aimag(state%local_rows))).or.&
      .not.all(ieee_is_finite(real(state%sipg_rows))).or.&
      .not.all(ieee_is_finite(aimag(state%sipg_rows))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid certified hybrid RT state extent';return;endif

    allocate(local_counts(r),global_counts(r));local_counts=0;local_bad=0
    do i=1,nowned
      if(state%owned_row_ids(i)<1_int64.or.state%owned_row_ids(i)>int(r,int64))then
        local_bad=1
      else
        local_counts(int(state%owned_row_ids(i)))=local_counts(int(state%owned_row_ids(i)))+1
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='out-of-range certified hybrid RT row';return;endif
    call MPI_Allreduce(local_counts,global_counts,r,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_counts/=1))then;message='invalid certified hybrid RT row ownership';return;endif

    local_bad=0
    if(state%metric%row_offsets(1)/=1.or.state%operators%row_offsets(1)/=1.or.&
      state%metric%row_offsets(nowned+1)/=size(state%metric%column_ids)+1.or.&
      state%operators%row_offsets(nowned+1)/=size(state%operators%column_ids)+1)local_bad=1
    if(any(state%metric%column_ids<1).or.any(state%metric%column_ids>r).or.&
      any(state%operators%column_ids<1).or.any(state%operators%column_ids>r))local_bad=1
    do i=1,nowned
      if(state%metric%row_offsets(i)<1.or.state%metric%row_offsets(i+1)<state%metric%row_offsets(i).or.&
        state%metric%row_offsets(i+1)>size(state%metric%column_ids)+1.or.&
        state%operators%row_offsets(i)<1.or.&
        state%operators%row_offsets(i+1)<state%operators%row_offsets(i).or.&
        state%operators%row_offsets(i+1)>size(state%operators%column_ids)+1)local_bad=1
      do edge=state%metric%row_offsets(i)+1,state%metric%row_offsets(i+1)-1
        if(state%metric%column_ids(edge)<=state%metric%column_ids(edge-1))local_bad=1
      enddo
      do edge=state%operators%row_offsets(i)+1,state%operators%row_offsets(i+1)-1
        if(state%operators%column_ids(edge)<=state%operators%column_ids(edge-1))local_bad=1
      enddo
    enddo
    if(.not.all(ieee_is_finite(real(state%metric%values))).or.&
      .not.all(ieee_is_finite(aimag(state%metric%values))).or.&
      .not.all(ieee_is_finite(real(state%operators%metric_values))).or.&
      .not.all(ieee_is_finite(aimag(state%operators%metric_values))).or.&
      .not.all(ieee_is_finite(real(state%operators%hamiltonian_values))).or.&
      .not.all(ieee_is_finite(aimag(state%operators%hamiltonian_values))).or.&
      .not.all(ieee_is_finite(real(state%operators%position_values))).or.&
      .not.all(ieee_is_finite(aimag(state%operators%position_values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid certified hybrid RT sparse layout';return;endif
    ok=.true.;message=''
#else
    ok=.false.;message='certified hybrid RT state validation requires MPI'
#endif
  end subroutine validate_certified_rt_state
end module rt_dg_hybrid_density_update
