#include "config.h"
module dg_hybrid_sparse_metric
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_sparse_metric
    logical::valid=.false.
    integer::global_count=0,numerical_rank=0,max_row_nnz=0
    real(real64)::maximum_value=0d0,condition_estimate=huge(1d0)
    integer(int64)::fingerprint=0_int64,workspace_peak_bytes=0_int64
    integer(int64),allocatable::owned_row_ids(:)
    integer,allocatable::row_offsets(:),column_ids(:),packet_ids(:)
    logical,allocatable::active_rows(:)
    complex(real64),allocatable::values(:)
  end type s_dg_hybrid_sparse_metric
  public::build_dg_hybrid_sparse_metric,apply_dg_hybrid_sparse_metric
contains
  subroutine build_dg_hybrid_sparse_metric(comm,global_count,row_ids,row_offsets,column_ids,values,&
      packet_ids,max_packet_size,selection_fingerprint,packet_fingerprint,complement_fingerprint,&
      tolerance,metric,active_packets,rejected_packets,condition_estimate,workspace_peak_bytes,&
      fingerprint,ok,message)
    integer,intent(in)::comm,global_count,max_packet_size
    integer(int64),intent(in)::row_ids(:)
    integer,intent(in)::row_offsets(:),column_ids(:),packet_ids(:)
    complex(real64),intent(in)::values(:)
    integer(int64),intent(in)::selection_fingerprint,packet_fingerprint,complement_fingerprint
    real(real64),intent(in)::tolerance
    type(s_dg_hybrid_sparse_metric),intent(out)::metric
    integer,allocatable,intent(out)::active_packets(:),rejected_packets(:)
    real(real64),intent(out)::condition_estimate
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,nowned,nnz,npacket,i,j,k,row,target,root,position,local_bad,global_bad
    integer::minimum_integer,maximum_integer,allocation_status,packet,nmember,nactive,nrejected,info,lwork,&
      decision,minimum_decision,maximum_decision
    integer,allocatable::ownership_count(:),owner(:),owner_position(:),members(:),packet_size(:)
    integer(int64)::bits,minimum_bits,maximum_bits,complex_elements,integer_elements,quantized
    logical,allocatable::packet_active(:)
    complex(real64),allocatable::remote_row(:),packet_block(:,:),work(:)
    real(real64),allocatable::eigenvalues(:),rwork(:)
    real(real64)::hermitian_defect,scale,threshold,local_lower,global_lower,local_upper,global_upper,&
      diagonal,radius,quantization_limit
    external::zheev
    ok=.false.;message='';condition_estimate=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64
    local_bad=0;nowned=size(row_ids);nnz=size(values)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call agree_integer(global_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent sparse metric extent';return;endif
    call agree_integer(max_packet_size,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent metric packet bound';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent sparse metric tolerance';return;endif
    call agree_receipt(selection_fingerprint,'selection',local_bad)
    call agree_receipt(packet_fingerprint,'packet',local_bad)
    call agree_receipt(complement_fingerprint,'complement',local_bad)
    if(local_bad/=0)return
    if(global_count<1.or.max_packet_size<1.or.size(packet_ids)/=global_count)local_bad=1
    if(size(row_offsets)/=nowned+1.or.size(column_ids)/=nnz)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance).or..not.finite_vector(values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse metric shape or finite contract';return;endif
    if(tolerance<1d-15.or.tolerance>1d-2.or.any(packet_ids<1))local_bad=1
    do i=1,global_count
      call agree_integer(packet_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent sparse metric packet membership';return
      endif
    enddo
    if(row_offsets(1)/=1.or.row_offsets(nowned+1)/=nnz+1)local_bad=1
    if(any(row_offsets(2:)<row_offsets(:nowned)))local_bad=1
    do i=1,nowned
      do k=row_offsets(i),row_offsets(i+1)-1
        if(column_ids(k)<1.or.column_ids(k)>global_count)local_bad=1
        if(k>row_offsets(i))then
          if(column_ids(k)<=column_ids(k-1))local_bad=1
        endif
      enddo
      if(find_column(int(row_ids(i)),i,row_offsets,column_ids)==0)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse metric CSR contract';return;endif
    npacket=maxval(packet_ids)
    if(npacket>global_count.or.max_packet_size>huge(0)/3)local_bad=1
    if(int(max_packet_size,int64)*int(max_packet_size,int64)>int(huge(0),int64))local_bad=1
    complex_elements=int(nnz,int64)+int(global_count,int64)+&
      int(max_packet_size,int64)*int(max_packet_size,int64)+3_int64*int(max_packet_size,int64)
    integer_elements=3_int64*int(global_count,int64)+int(nowned,int64)+1_int64+int(nnz,int64)+&
      int(global_count,int64)+int(npacket,int64)+int(max_packet_size,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.&
      integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>&
      huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='sparse metric preallocation extent or receipt overflow';return
    endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(ownership_count(global_count),owner(global_count),owner_position(global_count),&
      packet_size(npacket),packet_active(npacket),members(max_packet_size),remote_row(global_count),&
      packet_block(max_packet_size,max_packet_size),eigenvalues(max_packet_size),&
      work(max(1,2*max_packet_size)),rwork(max(1,3*max_packet_size-2)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse metric workspace';return;endif
    packet_size=0
    do i=1,global_count;packet_size(packet_ids(i))=packet_size(packet_ids(i))+1;enddo
    if(any(packet_size<1).or.any(packet_size>max_packet_size))local_bad=1
    ownership_count=0;owner=-1;owner_position=0
    do i=1,nowned
      row=int(row_ids(i));ownership_count(row)=ownership_count(row)+1;owner(row)=rank;owner_position(row)=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='sparse metric ownership count failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='sparse metric owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='sparse metric position reduction failed';return;endif
    if(any(ownership_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='invalid sparse metric ownership or packet size';return;endif
    hermitian_defect=0d0;scale=1d0
    do target=1,global_count
      call broadcast_row(target,remote_row,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='sparse metric row broadcast failed';return;endif
      scale=max(scale,maxval(abs(remote_row)))
      if(abs(aimag(remote_row(target)))>hermitian_defect)hermitian_defect=abs(aimag(remote_row(target)))
      do i=1,nowned
        k=find_column(target,i,row_offsets,column_ids)
        if(k>0)then
          hermitian_defect=max(hermitian_defect,abs(values(k)-conjg(remote_row(int(row_ids(i))))))
        else
          hermitian_defect=max(hermitian_defect,abs(remote_row(int(row_ids(i)))))
        endif
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,hermitian_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.hermitian_defect>100d0*tolerance*scale)then
      call cleanup();message='sparse metric is not Hermitian or has a missing reverse edge';return
    endif
    packet_active=.true.
    do packet=1,npacket
      nmember=0
      do i=1,global_count
        if(packet_ids(i)==packet)then;nmember=nmember+1;members(nmember)=i;endif
      enddo
      packet_block(1:nmember,1:nmember)=(0d0,0d0)
      do i=1,nmember
        call broadcast_row(members(i),remote_row,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric packet row broadcast failed';return;endif
        do j=1,nmember;packet_block(i,j)=remote_row(members(j));enddo
      enddo
      lwork=max(1,2*nmember);call zheev('N','U',nmember,packet_block,max_packet_size,&
        eigenvalues,work,lwork,rwork,info)
      local_bad=merge(0,1,info==0.and.all(ieee_is_finite(eigenvalues(1:nmember))))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='metric packet eigensystem failed';return;endif
      threshold=tolerance*max(1d0,maxval(abs(eigenvalues(1:nmember))))
      decision=merge(1,0,minval(eigenvalues(1:nmember)) < -threshold)
      call agree_integer(decision,minimum_decision,maximum_decision,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_decision/=maximum_decision)then
        call cleanup();message='rank-disagreeing indefinite metric decision';return
      endif
      if(decision==1)then
        call cleanup();message='sparse metric packet is indefinite';return
      endif
      decision=merge(1,0,minval(eigenvalues(1:nmember))<=threshold)
      call agree_integer(decision,minimum_decision,maximum_decision,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_decision/=maximum_decision)then
        call cleanup();message='rank-disagreeing metric rank decision';return
      endif
      if(decision==1)packet_active(packet)=.false.
    enddo
    nactive=count(packet_active);nrejected=npacket-nactive
    if(nactive==0)then;call cleanup();message='all sparse metric packets are rank deficient';return;endif
    allocate(active_packets(nactive),rejected_packets(nrejected),metric%active_rows(global_count),&
      metric%packet_ids(global_count),metric%owned_row_ids(nowned),metric%row_offsets(nowned+1),&
      metric%column_ids(nnz),metric%values(nnz),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse metric outputs';return;endif
    nactive=0;nrejected=0
    do packet=1,npacket
      if(packet_active(packet))then;nactive=nactive+1;active_packets(nactive)=packet
      else;nrejected=nrejected+1;rejected_packets(nrejected)=packet;endif
    enddo
    do i=1,global_count;metric%active_rows(i)=packet_active(packet_ids(i));enddo
    local_lower=huge(1d0);local_upper=0d0
    do i=1,nowned
      row=int(row_ids(i));if(.not.metric%active_rows(row))cycle
      diagonal=0d0;radius=0d0
      do k=row_offsets(i),row_offsets(i+1)-1
        if(.not.metric%active_rows(column_ids(k)))cycle
        if(column_ids(k)==row)then;diagonal=real(values(k));else;radius=radius+abs(values(k));endif
      enddo
      local_lower=min(local_lower,diagonal-radius);local_upper=max(local_upper,diagonal+radius)
    enddo
    call MPI_Allreduce(local_lower,global_lower,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric lower-bound reduction failed';return;endif
    call MPI_Allreduce(local_upper,global_upper,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_lower<=tolerance*max(1d0,global_upper))then
      call cleanup();message='active sparse metric lacks a positive local bound';return
    endif
    condition_estimate=global_upper/global_lower
    metric%valid=.true.;metric%global_count=global_count
    metric%numerical_rank=count(metric%active_rows);metric%condition_estimate=condition_estimate
    metric%maximum_value=scale;metric%max_row_nnz=0
    do i=1,nowned;metric%max_row_nnz=max(metric%max_row_nnz,row_offsets(i+1)-row_offsets(i));enddo
    metric%owned_row_ids=row_ids;metric%row_offsets=row_offsets;metric%column_ids=column_ids
    metric%values=values;metric%packet_ids=packet_ids
    quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*tolerance
    if(scale>quantization_limit)then;call cleanup();message='sparse metric fingerprint range is unsafe';return;endif
    fingerprint=ieor(selection_fingerprint,ishftc(packet_fingerprint,7))
    fingerprint=ieor(fingerprint,ishftc(complement_fingerprint,13))
    do target=1,global_count
      call broadcast_row(target,remote_row,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='metric fingerprint row broadcast failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),int(target,int64))
      fingerprint=ieor(ishftc(fingerprint,9),int(packet_ids(target),int64))
      if(metric%active_rows(target))fingerprint=not(fingerprint)
      do j=1,global_count
        if(abs(remote_row(j))==0d0)cycle
        fingerprint=ieor(ishftc(fingerprint,9),int(j,int64))
        quantized=nint(real(remote_row(j))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
        quantized=nint(aimag(remote_row(j))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    metric%fingerprint=fingerprint;metric%workspace_peak_bytes=workspace_peak_bytes;ok=.true.
#else
    ok=.false.;message='hybrid sparse metric requires MPI';condition_estimate=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;allocate(active_packets(0),rejected_packets(0))
#endif
  contains
#ifdef USE_MPI
    subroutine agree_receipt(value,label,bad)
      integer(int64),intent(in)::value
      character(*),intent(in)::label
      integer,intent(inout)::bad
      call agree_int64(value,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.value==0_int64)then
        message='invalid or inconsistent sparse metric '//trim(label)//' provenance';bad=1
      endif
    end subroutine agree_receipt
    subroutine broadcast_row(global_row,row_values,status)
      integer,intent(in)::global_row
      complex(real64),intent(out)::row_values(:)
      integer,intent(out)::status
      integer::local_position_index,edge
      row_values=(0d0,0d0);root=owner(global_row)
      if(rank==root)then
        local_position_index=owner_position(global_row)
        do edge=row_offsets(local_position_index),row_offsets(local_position_index+1)-1
          row_values(column_ids(edge))=values(edge)
        enddo
      endif
      call MPI_Bcast(row_values,global_count,MPI_DOUBLE_COMPLEX,root,comm,status)
    end subroutine broadcast_row
    subroutine cleanup()
      if(allocated(active_packets))deallocate(active_packets)
      if(allocated(rejected_packets))deallocate(rejected_packets)
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position)
      if(allocated(packet_size))deallocate(packet_size)
      if(allocated(packet_active))deallocate(packet_active)
      if(allocated(members))deallocate(members)
      if(allocated(remote_row))deallocate(remote_row)
      if(allocated(packet_block))deallocate(packet_block)
      if(allocated(eigenvalues))deallocate(eigenvalues)
      if(allocated(work))deallocate(work)
      if(allocated(rwork))deallocate(rwork)
      metric%valid=.false.
    end subroutine cleanup
#endif
  end subroutine build_dg_hybrid_sparse_metric

  subroutine apply_dg_hybrid_sparse_metric(metric,x_global,y_owned,ok,message)
    type(s_dg_hybrid_sparse_metric),intent(in)::metric
    complex(real64),intent(in)::x_global(:)
    complex(real64),allocatable,intent(out)::y_owned(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,k,row,allocation_status
    real(real64)::xscale,safe_scale
    ok=.false.;message=''
    if(.not.metric%valid.or.size(x_global)/=metric%global_count.or..not.finite_vector(x_global))then
      message='invalid sparse metric apply contract';return
    endif
    xscale=maxval(abs(x_global))
    if(metric%maximum_value>0d0.and.xscale>0d0)then
      safe_scale=(huge(1d0)/4d0)/real(max(1,metric%max_row_nnz),real64)/metric%maximum_value
      if(xscale>safe_scale)then;message='sparse metric apply magnitude is unsafe';return;endif
    endif
    allocate(y_owned(size(metric%owned_row_ids)),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate sparse metric apply output';return;endif
    y_owned=(0d0,0d0)
    do i=1,size(metric%owned_row_ids)
      row=int(metric%owned_row_ids(i));if(.not.metric%active_rows(row))cycle
      do k=metric%row_offsets(i),metric%row_offsets(i+1)-1
        if(metric%active_rows(metric%column_ids(k)))&
          y_owned(i)=y_owned(i)+metric%values(k)*x_global(metric%column_ids(k))
      enddo
    enddo
    if(.not.finite_vector(y_owned))then;deallocate(y_owned);message='nonfinite sparse metric result';return;endif
    ok=.true.
  end subroutine apply_dg_hybrid_sparse_metric

  integer function find_column(column,row_position,row_offsets,column_ids)
    integer,intent(in)::column,row_position,row_offsets(:),column_ids(:)
    integer::k
    find_column=0
    do k=row_offsets(row_position),row_offsets(row_position+1)-1
      if(column_ids(k)==column)then;find_column=k;return;endif
    enddo
  end function find_column
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value
    integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module dg_hybrid_sparse_metric
