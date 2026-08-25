#include "config.h"
module dg_hybrid_window_distribution
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::prepare_dg_hybrid_window_distribution,redistribute_dg_hybrid_fragment_windows
contains
  subroutine redistribute_dg_hybrid_fragment_windows(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,request_ids,requested_windows,workspace_peak_bytes,&
      fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:)
    integer(int64),intent(in)::box_ids(:),request_ids(:)
    real(real64),intent(in)::box_windows(:,:)
    real(real64),allocatable,intent(out)::requested_windows(:,:)
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,f,target,ierr,local_bad,global_bad,minimum,maximum,allocation_status
    integer,allocatable::fragment_presence(:),local_count(:),global_count(:)
    real(real64),allocatable::local_values(:),global_values(:)
    integer(int64)::bits
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    call agree_integer(global_point_count,minimum,maximum,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='inconsistent requested window point count';return;endif
    call agree_integer(fragment_count,minimum,maximum,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='inconsistent requested window fragment count';return;endif
    local_bad=merge(0,1,global_point_count>0.and.fragment_count>0.and.&
      size(box_windows,1)==size(fragment_ids).and.size(box_windows,2)==size(box_ids).and.&
      all(fragment_ids>=1).and.all(fragment_ids<=fragment_count).and.&
      all(box_ids>=1_int64).and.all(box_ids<=int(global_point_count,int64)).and.&
      all(request_ids>=1_int64).and.all(request_ids<=int(global_point_count,int64)).and.&
      all(ieee_is_finite(box_windows)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid requested window input';return;endif
    allocate(fragment_presence(fragment_count),local_count(fragment_count),global_count(fragment_count),&
      local_values(fragment_count),global_values(fragment_count),requested_windows(fragment_count,size(request_ids)),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate requested window workspace';return;endif
    fragment_presence=0
    do i=1,size(fragment_ids);fragment_presence(fragment_ids(i))=fragment_presence(fragment_ids(i))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,fragment_presence,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(fragment_presence/=1))then
      message='duplicate or missing requested window fragment owner';return
    endif
    requested_windows=0d0;fingerprint=int(z'9B05688C2B3E6C1F',int64)
    do target=1,global_point_count
      local_values=0d0;local_count=0
      do i=1,size(fragment_ids);do j=1,size(box_ids)
        if(box_ids(j)/=int(target,int64))cycle
        f=fragment_ids(i);local_values(f)=local_values(f)+box_windows(i,j);local_count(f)=local_count(f)+1
      enddo;enddo
      call MPI_Allreduce(local_values,global_values,fragment_count,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr==MPI_SUCCESS)&
        call MPI_Allreduce(local_count,global_count,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='requested fragment window reduction failed';return;endif
      do f=1,fragment_count
        if(global_count(f)>0)global_values(f)=global_values(f)/real(global_count(f),real64)
      enddo
      do i=1,size(request_ids)
        if(request_ids(i)==int(target,int64))requested_windows(:,i)=global_values
      enddo
      do f=1,fragment_count
        bits=transfer(global_values(f),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
      enddo
    enddo
    workspace_peak_bytes=24_int64*int(fragment_count,int64)
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(fragment_presence,local_count,global_count,local_values,global_values);ok=.true.
#else
    ok=.false.;message='requested window redistribution requires MPI'
    workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  end subroutine redistribute_dg_hybrid_fragment_windows

  subroutine prepare_dg_hybrid_window_distribution(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,row_action,raw_windows,&
      fragment_action,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:),core_fragment_ids(:)
    integer(int64),intent(in)::box_ids(:),core_ids(:)
    real(real64),intent(in)::box_windows(:,:)
    integer,intent(in)::row_action(:,:)
    real(real64),allocatable,intent(out)::raw_windows(:,:)
    integer,allocatable,intent(out)::fragment_action(:,:)
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,noperation,local_bad,global_bad,i,j,f,op,target,mapped,allocation_status
    integer::minimum_integer,maximum_integer
    integer,allocatable::fragment_presence(:),point_ownership(:),point_fragment(:),local_count(:),global_count(:)
    integer(int64)::bits,minimum_bits,maximum_bits
    real(real64),allocatable::local_values(:),global_values(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    noperation=size(row_action,2)
    call agree_integer(global_point_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent window global point count';return
    endif
    call agree_integer(fragment_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent window fragment count';return
    endif
    call agree_integer(noperation,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent window operation count';return
    endif
    local_bad=merge(0,1,global_point_count>0.and.fragment_count>0.and.noperation>0.and.&
      size(box_windows,1)==size(fragment_ids).and.size(box_windows,2)==size(box_ids).and.&
      size(core_fragment_ids)==size(core_ids).and.&
      all(shape(row_action)==[global_point_count,noperation]).and.&
      all(fragment_ids>=1).and.all(fragment_ids<=fragment_count).and.&
      all(box_ids>=1_int64).and.all(box_ids<=int(global_point_count,int64)).and.&
      all(core_ids>=1_int64).and.all(core_ids<=int(global_point_count,int64)).and.&
      all(core_fragment_ids>=1).and.all(core_fragment_ids<=fragment_count).and.&
      all(ieee_is_finite(box_windows)))
    do op=1,noperation;do i=1,global_point_count
      call agree_integer(row_action(i,op),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid distributed window input';return
    endif
    allocate(fragment_presence(fragment_count),point_ownership(global_point_count),&
      point_fragment(global_point_count),local_count(fragment_count),global_count(fragment_count),&
      local_values(fragment_count),global_values(fragment_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate window distribution workspace';return
    endif
    fragment_presence=0
    do i=1,size(fragment_ids);fragment_presence(fragment_ids(i))=fragment_presence(fragment_ids(i))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,fragment_presence,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(fragment_presence/=1))then
      call cleanup();message='duplicate or missing window fragment ownership';return
    endif
    point_ownership=0;point_fragment=0
    do i=1,size(core_ids)
      point_ownership(int(core_ids(i)))=point_ownership(int(core_ids(i)))+1
      point_fragment(int(core_ids(i)))=core_fragment_ids(i)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,point_ownership,global_point_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr==MPI_SUCCESS)&
      call MPI_Allreduce(MPI_IN_PLACE,point_fragment,global_point_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(point_ownership/=1).or.any(point_fragment<1).or.&
        any(point_fragment>fragment_count))then
      call cleanup();message='duplicate or missing window core ownership';return
    endif
    local_bad=0
    do op=1,noperation
      if(any(row_action(:,op)<1).or.any(row_action(:,op)>global_point_count))local_bad=1
      do target=1,global_point_count
        if(count(row_action(:,op)==target)/=1)local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='window spatial action is not a permutation';return
    endif
    allocate(raw_windows(fragment_count,size(core_ids)),fragment_action(fragment_count,noperation),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate distributed window outputs';return
    endif
    fingerprint=int(z'510E527FADE682D1',int64)
    do target=1,global_point_count
      local_values=0d0;local_count=0
      do i=1,size(fragment_ids);do j=1,size(box_ids)
        if(box_ids(j)/=int(target,int64))cycle
        f=fragment_ids(i);local_values(f)=local_values(f)+box_windows(i,j);local_count(f)=local_count(f)+1
      enddo;enddo
      call MPI_Allreduce(local_values,global_values,fragment_count,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr==MPI_SUCCESS)&
        call MPI_Allreduce(local_count,global_count,fragment_count,MPI_INTEGER,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then
        call cleanup();message='window core-point reduction failed';return
      endif
      do f=1,fragment_count
        if(global_count(f)>0)global_values(f)=global_values(f)/real(global_count(f),real64)
      enddo
      do i=1,size(core_ids)
        if(core_ids(i)==int(target,int64))raw_windows(:,i)=global_values
      enddo
      do f=1,fragment_count
        bits=transfer(global_values(f),bits)
        fingerprint=ieor(ishftc(fingerprint,7),bits)
      enddo
    enddo
    fragment_action=0;local_bad=0
    do op=1,noperation;do target=1,global_point_count
      f=point_fragment(target);mapped=point_fragment(row_action(target,op))
      if(fragment_action(f,op)==0)fragment_action(f,op)=mapped
      if(fragment_action(f,op)/=mapped)local_bad=1
    enddo;enddo
    do op=1,noperation
      do f=1,fragment_count
        if(count(fragment_action(:,op)==f)/=1)local_bad=1
        fingerprint=ieor(ishftc(fingerprint,7),int(fragment_action(f,op),int64))
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='spatial action does not map whole fragments';return
    endif
    workspace_peak_bytes=8_int64*int(global_point_count,int64)+24_int64*int(fragment_count,int64)
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(fragment_presence,point_ownership,point_fragment,local_count,global_count,local_values,global_values)
    ok=.true.
  contains
    subroutine cleanup()
      if(allocated(fragment_presence))deallocate(fragment_presence)
      if(allocated(point_ownership))deallocate(point_ownership)
      if(allocated(point_fragment))deallocate(point_fragment)
      if(allocated(local_count))deallocate(local_count)
      if(allocated(global_count))deallocate(global_count)
      if(allocated(local_values))deallocate(local_values)
      if(allocated(global_values))deallocate(global_values)
      if(allocated(raw_windows))deallocate(raw_windows)
      if(allocated(fragment_action))deallocate(fragment_action)
    end subroutine cleanup
#else
    ok=.false.;message='window distribution requires MPI';workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  end subroutine prepare_dg_hybrid_window_distribution

#ifdef USE_MPI
  subroutine agree_integer(value,minimum,maximum,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum,maximum,ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
#endif
end module dg_hybrid_window_distribution
