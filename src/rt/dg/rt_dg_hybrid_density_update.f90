#include "config.h"
module rt_dg_hybrid_density_update
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_hybrid_initialization,only:s_rt_dg_hybrid_state
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine project_rt_density(row_ids,grid_ids,density,local_rows,ok,message)
      import::int64,real64
      integer(int64),intent(in)::row_ids(:)
      integer(int64),intent(in)::grid_ids(:)
      real(real64),intent(in)::density(:)
      complex(real64),intent(out)::local_rows(:,:)
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
    complex(real64),allocatable::global_coefficients(:,:),orbital_values(:)
    integer::i,p,ierr,local_bad,global_bad
    state%density_freshly_reconstructed=.false.
    call validate_certified_rt_state(comm,state,ok,message);if(.not.ok)return
    allocate(global_coefficients(state%certified_rank,state%noccupied),orbital_values(state%noccupied))
    global_coefficients=(0d0,0d0)
    do i=1,size(state%owned_row_ids)
      global_coefficients(int(state%owned_row_ids(i)),:)=state%coefficients(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,state%certified_rank*state%noccupied,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT coefficient redistribution failed';return;endif
    do p=1,size(state%grid_ids)
      orbital_values=matmul(state%basis_values(:,p),global_coefficients)
      state%density(p)=sum(state%occupations*abs(orbital_values)**2)
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

  subroutine update_rt_dg_hybrid_density(comm,state,density,project_local,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_state),intent(inout)::state
    real(real64),intent(in)::density(:)
    procedure(project_rt_density)::project_local
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::new_local(:,:),new_h(:,:)
    integer::i,j,edge,row,ierr,local_bad,global_bad
    integer(int64)::local_hash,global_xor,global_sum,bits,pair_hash
    logical::callback_ok
    character(256)::callback_message
    ok=.false.;message='';local_bad=0
    if(.not.state%density_freshly_reconstructed)then
      call validate_certified_rt_state(comm,state,ok,message);if(.not.ok)return
    endif
    state%density_freshly_reconstructed=.false.
    ok=.false.
    if(size(density)/=size(state%grid_ids).or..not.all(ieee_is_finite(density)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT density update';return;endif
    allocate(new_local(size(state%owned_row_ids),state%certified_rank),&
      new_h(size(state%owned_row_ids),state%certified_rank))
    new_local=(0d0,0d0);callback_ok=.false.;callback_message=''
    call project_local(state%owned_row_ids,state%grid_ids,density,new_local,callback_ok,callback_message)
    local_bad=merge(0,1,callback_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid RT local projection failed';return;endif
    local_bad=merge(0,1,all(ieee_is_finite(real(new_local))).and.all(ieee_is_finite(aimag(new_local))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='nonfinite hybrid RT local projection';return;endif
    new_h=state%kinetic_rows+state%nonlocal_rows+new_local+state%sipg_rows
    local_bad=merge(0,1,all(ieee_is_finite(real(new_h))).and.all(ieee_is_finite(aimag(new_h))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='nonfinite hybrid RT updated Hamiltonian';return;endif
    local_bad=0
    do i=1,size(state%owned_row_ids)
      do j=1,state%certified_rank
        if(new_h(i,j)==(0d0,0d0))cycle
        if(.not.graph_contains(state,i,j))local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='density update escaped the frozen operator-union envelope';return;endif
    local_hash=0_int64;global_sum=0_int64
    do i=1,size(state%owned_row_ids)
      do j=1,state%certified_rank
        bits=transfer(real(new_h(i,j)),bits);pair_hash=ieor(ishftc(state%owned_row_ids(i),17),int(j,int64))
        pair_hash=ieor(ishftc(pair_hash,9),bits)
        bits=transfer(aimag(new_h(i,j)),bits);pair_hash=ieor(ishftc(pair_hash,9),bits)
        local_hash=ieor(local_hash,pair_hash);global_sum=global_sum+pair_hash
      enddo
    enddo
    call MPI_Allreduce(local_hash,global_xor,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,global_sum,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='hybrid RT value fingerprint reduction failed';return;endif
    state%local_rows=new_local;state%density=density
    do i=1,size(state%owned_row_ids)
      row=int(state%owned_row_ids(i))
      do edge=state%operators%row_offsets(i),state%operators%row_offsets(i+1)-1
        j=state%operators%column_ids(edge);state%operators%hamiltonian_values(edge)=new_h(i,j)
      enddo
    enddo
    state%operator_value_fingerprint=ieor(global_xor,ishftc(global_sum,13))
    if(state%operator_value_fingerprint==0_int64)state%operator_value_fingerprint=1_int64
    state%operators%fingerprint=state%operator_structure_fingerprint
    ok=.true.;message=''
#else
    ok=.false.;message='hybrid RT density update requires MPI'
#endif
  contains
    logical function graph_contains(current,row_position,column)
      type(s_rt_dg_hybrid_state),intent(in)::current
      integer,intent(in)::row_position,column
      integer::q
      graph_contains=.false.
      do q=current%operators%row_offsets(row_position),current%operators%row_offsets(row_position+1)-1
        if(current%operators%column_ids(q)==column)then;graph_contains=.true.;return;endif
      enddo
    end function graph_contains
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
      size(state%local_rows,1)/=nowned.or.size(state%local_rows,2)/=r.or.&
      size(state%sipg_rows,1)/=nowned.or.size(state%sipg_rows,2)/=r)local_bad=1
    if(size(state%basis_values,1)/=r.or.size(state%basis_values,2)/=npoint.or.&
      size(state%density)/=npoint.or.size(state%grid_weights)/=npoint.or.&
      size(state%occupations)/=nocc.or.size(state%eigenvalues)/=r)local_bad=1
    if(.not.state%metric%valid.or.state%metric%global_count/=r.or.&
      size(state%metric%owned_row_ids)/=nowned.or.size(state%metric%row_offsets)/=nowned+1.or.&
      size(state%metric%column_ids)/=nowned*r.or.size(state%metric%values)/=nowned*r.or.&
      size(state%metric%active_rows)/=r.or.size(state%metric%packet_ids)/=r)local_bad=1
    if(.not.state%operators%valid.or.state%operators%global_count/=r.or.&
      state%operators%metric_fingerprint/=state%metric%fingerprint.or.&
      size(state%operators%owned_row_ids)/=nowned.or.size(state%operators%row_offsets)/=nowned+1.or.&
      size(state%operators%column_ids)/=nowned*r.or.size(state%operators%metric_values)/=nowned*r.or.&
      size(state%operators%hamiltonian_values)/=nowned*r.or.&
      size(state%operators%position_values,1)/=3.or.&
      size(state%operators%position_values,2)/=nowned*r)local_bad=1
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
      state%metric%row_offsets(nowned+1)/=nowned*r+1.or.&
      state%operators%row_offsets(nowned+1)/=nowned*r+1)local_bad=1
    if(any(state%metric%column_ids<1).or.any(state%metric%column_ids>r).or.&
      any(state%operators%column_ids<1).or.any(state%operators%column_ids>r))local_bad=1
    do i=1,nowned
      if(state%metric%row_offsets(i)/=(i-1)*r+1.or.state%metric%row_offsets(i+1)/=i*r+1.or.&
        state%operators%row_offsets(i)/=(i-1)*r+1.or.state%operators%row_offsets(i+1)/=i*r+1)local_bad=1
      do j=1,r
        edge=(i-1)*r+j
        if(state%metric%column_ids(edge)/=j.or.state%operators%column_ids(edge)/=j)local_bad=1
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
