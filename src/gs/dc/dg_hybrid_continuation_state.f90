#include "config.h"
module dg_hybrid_continuation_state
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public::s_dg_hybrid_catalog_receipt
    logical::frozen=.false.
    logical::analysis_complete=.false.
    logical::identity_only=.false.
    integer::operation_count=0
    integer::nonidentity_count=0
    integer(int64)::analysis_fingerprint=0_int64
    integer(int64)::catalog_fingerprint=0_int64
    integer(int64)::selection_fingerprint=0_int64
  end type s_dg_hybrid_catalog_receipt

  type,public::s_dg_hybrid_scope_receipt
    logical::valid=.false.
    logical::periodic=.false.
    integer::theory_code=0
    integer::nspin=0
    logical::spinorbit=.false.
    logical::plus_u=.false.
    logical::hse=.false.
    logical::fix_func=.false.
    logical::jm=.false.
    integer,allocatable::xctype(:)
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_scope_receipt

  type,public::s_dg_hybrid_accepted_snapshot
    logical::valid=.false.
    real(real64)::lambda=0d0
    integer::epoch=-1
    integer(int64)::state_fingerprint=0_int64
  end type s_dg_hybrid_accepted_snapshot

  type,public::s_dg_hybrid_continuation_state
    logical::valid=.false.
    logical::lambda_accepted=.false.
    real(real64)::lambda=0d0
    integer::density_epoch=0
    integer::accepted_epoch=-1
    integer(int64)::seed_fingerprint=0_int64
    type(s_dg_hybrid_catalog_receipt)::catalog
    type(s_dg_hybrid_scope_receipt)::scope
    type(s_dg_hybrid_accepted_snapshot)::accepted
    real(real64),allocatable::seed_density(:)
    real(real64),allocatable::mixed_density(:)
  end type s_dg_hybrid_continuation_state

  public::initialize_dg_hybrid_continuation,close_dg_hybrid_selection,build_dg_hybrid_scope_receipt,&
    validate_dg_hybrid_frozen_catalog
contains
  subroutine validate_dg_hybrid_frozen_catalog(comm,state,current,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_continuation_state),intent(in)::state
    type(s_dg_hybrid_catalog_receipt),intent(in)::current
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,state%valid.and.current%frozen.and.current%analysis_complete.and.&
      (current%identity_only.eqv.state%catalog%identity_only).and.&
      current%operation_count==state%catalog%operation_count.and.&
      current%nonidentity_count==state%catalog%nonidentity_count.and.&
      current%analysis_fingerprint==state%catalog%analysis_fingerprint.and.&
      current%catalog_fingerprint==state%catalog%catalog_fingerprint.and.&
      current%selection_fingerprint==state%catalog%selection_fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0
    if(ok)then;message='';else;message='immutable continuation catalog changed';endif
#else
    ok=.false.;message='DG hybrid catalog validation requires MPI'
#endif
  end subroutine validate_dg_hybrid_frozen_catalog

  subroutine build_dg_hybrid_scope_receipt(comm,theory_code,periodic,nspin,spinorbit,plus_u,hse,&
      fix_func,jm,xctype,receipt,ok,message)
    integer,intent(in)::comm,theory_code,nspin,xctype(:)
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    type(s_dg_hybrid_scope_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,local_bad,global_bad
    integer(int64)::minimum_hash,maximum_hash
    ok=.false.;message=''
    local_bad=merge(0,1,theory_code==1.and.periodic.and.nspin==1.and..not.spinorbit.and.&
      .not.plus_u.and..not.hse.and..not.fix_func.and..not.jm.and.size(xctype)>0.and.&
      all(xctype>=1.and.xctype<=4))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='unsupported DG hybrid continuation scope';return
    endif
    receipt%theory_code=theory_code;receipt%periodic=periodic;receipt%nspin=nspin
    receipt%spinorbit=spinorbit;receipt%plus_u=plus_u;receipt%hse=hse
    receipt%fix_func=fix_func;receipt%jm=jm
    allocate(receipt%xctype(size(xctype)));receipt%xctype=xctype
    receipt%fingerprint=int(z'A4093822299F31D0',int64)
    receipt%fingerprint=ieor(ishftc(receipt%fingerprint,7),int(theory_code,int64))
    receipt%fingerprint=ieor(ishftc(receipt%fingerprint,7),int(nspin,int64))
    receipt%fingerprint=ieor(ishftc(receipt%fingerprint,7),merge(1_int64,0_int64,periodic))
    do i=1,size(xctype)
      receipt%fingerprint=ieor(ishftc(receipt%fingerprint,7),int(xctype(i),int64))
    enddo
    if(receipt%fingerprint==0_int64)receipt%fingerprint=1_int64
    call agree_catalog_int64(receipt%fingerprint,comm,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing DG hybrid continuation scope';return
    endif
    receipt%valid=.true.;ok=.true.
#else
    ok=.false.;message='DG hybrid scope receipt requires MPI'
#endif
  end subroutine build_dg_hybrid_scope_receipt

  subroutine close_dg_hybrid_selection(comm,requested_ids,universe_ids,group_action,effective_ids,&
      added_parent,added_operation,fingerprint,ok,message)
    integer,intent(in)::comm,requested_ids(:),universe_ids(:),group_action(:,:)
    integer,allocatable,intent(out)::effective_ids(:),added_parent(:),added_operation(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    logical,allocatable::selected(:)
    integer,allocatable::parent(:),operation(:),seen(:)
    integer::i,j,op,op2,candidate,target,added_count,ierr,local_bad,global_bad
    integer(int64)::minimum_hash,maximum_hash

    ok=.false.;message='';fingerprint=0_int64;local_bad=0
    if(size(universe_ids)<1.or.size(requested_ids)<1.or.size(group_action,1)/=size(universe_ids).or.&
        size(group_action,2)<1)then
      local_bad=1
    elseif(any(universe_ids<=0).or.any(requested_ids<=0))then
      local_bad=1
    endif
    if(local_bad==0)then
      allocate(seen(size(universe_ids)))
      do op=1,size(group_action,2)
        seen=0
        do i=1,size(universe_ids)
          target=group_action(i,op)
          if(target<1.or.target>size(universe_ids))then
            local_bad=1
          else
            seen(target)=seen(target)+1
          endif
        enddo
        if(any(seen/=1))local_bad=1
      enddo
      if(any(group_action(:,1)/=[(i,i=1,size(universe_ids))]))local_bad=1
      do op=1,size(group_action,2)
        do op2=1,size(group_action,2)
          candidate=0
          do j=1,size(group_action,2)
            if(all(group_action(:,j)==group_action(group_action(:,op2),op)))then
              candidate=j;exit
            endif
          enddo
          if(candidate==0)local_bad=1
        enddo
      enddo
      do i=1,size(universe_ids)
        if(count(universe_ids==universe_ids(i))/=1)local_bad=1
      enddo
      do i=1,size(requested_ids)
        if(count(universe_ids==requested_ids(i))/=1)local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid or incomplete symmetry selection action';return
    endif

    fingerprint=int(z'13198A2E03707344',int64)
    do i=1,size(universe_ids)
      fingerprint=ieor(ishftc(fingerprint,7),int(universe_ids(i),int64))
      do op=1,size(group_action,2)
        fingerprint=ieor(ishftc(fingerprint,7),int(group_action(i,op),int64))
      enddo
    enddo
    do i=1,size(requested_ids)
      fingerprint=ieor(ishftc(fingerprint,11),int(requested_ids(i),int64))
    enddo
    call agree_catalog_int64(fingerprint,comm,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing symmetry selection';return
    endif

    allocate(selected(size(universe_ids)),parent(size(universe_ids)),operation(size(universe_ids)))
    selected=.false.;parent=0;operation=0
    do i=1,size(requested_ids)
      j=find_index(requested_ids(i),universe_ids)
      selected(j)=.true.
    enddo
    do
      added_count=count(selected)
      do i=1,size(universe_ids)
        if(.not.selected(i))cycle
        do op=1,size(group_action,2)
          target=group_action(i,op)
          if(selected(target))cycle
          selected(target)=.true.;parent(target)=universe_ids(i);operation(target)=op
        enddo
      enddo
      if(count(selected)==added_count)exit
    enddo
    allocate(effective_ids(count(selected)),added_parent(count(selected)-size(requested_ids)),&
      added_operation(count(selected)-size(requested_ids)))
    j=0;added_count=0
    do i=1,size(universe_ids)
      if(.not.selected(i))cycle
      j=j+1;effective_ids(j)=universe_ids(i)
      if(parent(i)==0)cycle
      added_count=added_count+1
      added_parent(added_count)=parent(i);added_operation(added_count)=operation(i)
    enddo
    do i=1,size(effective_ids)
      fingerprint=ieor(ishftc(fingerprint,13),int(effective_ids(i),int64))
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='DG hybrid selection closure requires MPI';fingerprint=0_int64
#endif
  end subroutine close_dg_hybrid_selection

  subroutine initialize_dg_hybrid_continuation(comm,dc_density,trial_density,core_ids,catalog,scope,&
      state,fingerprint,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::dc_density(:),trial_density(:)
    integer(int64),intent(in)::core_ids(:)
    type(s_dg_hybrid_catalog_receipt),intent(in)::catalog
    type(s_dg_hybrid_scope_receipt),intent(in)::scope
    type(s_dg_hybrid_continuation_state),intent(out)::state
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,local_bad,global_bad,minimum_integer,maximum_integer
    integer,allocatable::ownership(:)
    integer(int64)::minimum_hash,maximum_hash,bits

    ok=.false.;message='';fingerprint=0_int64
    local_bad=merge(0,1,size(dc_density)>0.and.size(trial_density)==size(dc_density).and.&
      all(ieee_is_finite(dc_density)).and.all(ieee_is_finite(trial_density)))
    local_bad=max(local_bad,merge(0,1,catalog%frozen.and.catalog%analysis_complete.and.&
      catalog%operation_count>=1.and.catalog%nonidentity_count>=0.and.&
      catalog%nonidentity_count<catalog%operation_count.and.catalog%analysis_fingerprint/=0_int64.and.&
      catalog%catalog_fingerprint/=0_int64.and.catalog%selection_fingerprint/=0_int64))
    local_bad=max(local_bad,merge(0,1,scope%valid.and.scope%fingerprint/=0_int64))
    if(catalog%identity_only)local_bad=max(local_bad,merge(0,1,catalog%operation_count==1.and.&
      catalog%nonidentity_count==0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid continuation seed or authoritative catalog receipt';return
    endif

    call agree_catalog_integer(catalog%operation_count,comm,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing symmetry operation count';return
    endif
    call agree_catalog_int64(catalog%analysis_fingerprint,comm,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing symmetry provenance';return
    endif
    call agree_catalog_int64(catalog%catalog_fingerprint,comm,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing catalog fingerprint';return
    endif
    call agree_catalog_int64(scope%fingerprint,comm,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing supported-scope receipt';return
    endif

    allocate(ownership(size(dc_density)));ownership=0;local_bad=0
    do i=1,size(core_ids)
      if(core_ids(i)<1_int64.or.core_ids(i)>int(size(dc_density),int64))then
        local_bad=1
      else
        ownership(int(core_ids(i)))=ownership(int(core_ids(i)))+1
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,size(ownership),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='core ownership reduction failed';return;endif
    if(local_bad/=0.or.any(ownership/=1))then;message='invalid distributed core ownership';return;endif

    allocate(state%seed_density(size(dc_density)),state%mixed_density(size(dc_density)))
    state%seed_density=dc_density
    state%mixed_density=dc_density
    state%lambda=0d0
    state%lambda_accepted=.false.
    state%density_epoch=0
    state%accepted_epoch=-1
    state%catalog=catalog
    state%scope=scope
    fingerprint=int(z'243F6A8885A308D3',int64)
    do i=1,size(dc_density)
      bits=transfer(dc_density(i),bits)
      fingerprint=ieor(ishftc(fingerprint,7),bits)
    enddo
    fingerprint=ieor(ishftc(fingerprint,11),catalog%analysis_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,13),catalog%catalog_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,17),catalog%selection_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,19),scope%fingerprint)
    if(fingerprint==0_int64)fingerprint=1_int64
    state%seed_fingerprint=fingerprint
    state%valid=.true.
    ok=.true.
#else
    ok=.false.;message='DG hybrid continuation state requires MPI';fingerprint=0_int64
#endif
  end subroutine initialize_dg_hybrid_continuation

#ifdef USE_MPI
  integer function find_index(value,values) result(index)
    integer,intent(in)::value,values(:)
    integer::i
    index=0
    do i=1,size(values)
      if(values(i)==value)then;index=i;return;endif
    enddo
  end function find_index

  subroutine agree_catalog_integer(value,comm,minimum,maximum,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum,maximum,ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_catalog_integer

  subroutine agree_catalog_int64(value,comm,minimum,maximum,ierr)
    integer(int64),intent(in)::value
    integer,intent(in)::comm
    integer(int64),intent(out)::minimum,maximum
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_catalog_int64
#endif
end module dg_hybrid_continuation_state
