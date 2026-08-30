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
    allocate(global_coefficients(state%global_count,state%noccupied),orbital_values(state%noccupied))
    global_coefficients=(0d0,0d0)
    do i=1,size(state%owned_row_ids)
      global_coefficients(int(state%owned_row_ids(i)),:)=state%coefficients(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,global_coefficients,state%global_count*state%noccupied,&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT coefficient redistribution failed';return;endif
    do p=1,size(state%grid_ids)
      orbital_values=matmul(state%basis_values(:,p),global_coefficients)
      state%density(p)=sum(state%occupations*abs(orbital_values)**2)
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(state%density)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0
    if(ok)then;message='';else;message='nonfinite hybrid RT reconstructed density';endif
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
    if(.not.state%valid.or.size(density)/=size(state%grid_ids).or.size(density)/=size(state%grid_weights).or.&
      size(state%basis_values,2)/=size(density).or..not.all(ieee_is_finite(density)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT density update';return;endif
    allocate(new_local(size(state%owned_row_ids),state%global_count),new_h(size(state%owned_row_ids),state%global_count))
    call project_local(state%owned_row_ids,state%grid_ids,density,new_local,callback_ok,callback_message)
    local_bad=merge(0,1,callback_ok.and.all(ieee_is_finite(real(new_local))).and.all(ieee_is_finite(aimag(new_local))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid RT local projection failed: '//trim(callback_message);return;endif
    new_h=state%kinetic_rows+state%nonlocal_rows+new_local+state%sipg_rows
    local_bad=0
    do i=1,size(state%owned_row_ids)
      do j=1,state%global_count
        if(new_h(i,j)==(0d0,0d0))cycle
        if(.not.graph_contains(state,i,j))local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='density update escaped the frozen operator-union envelope';return;endif
    local_hash=0_int64;global_sum=0_int64
    do i=1,size(state%owned_row_ids)
      do j=1,state%global_count
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
end module rt_dg_hybrid_density_update
